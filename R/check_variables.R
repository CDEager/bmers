## these are functions for internal use in build_bmer_model

# returns TRUE if a numeric vector is on unit scale, FALSE otherwise
is.scaled <- function(x, tol = 1e-5){
  return(max(abs(x - scale(x)[,1])) < tol)
}


# returns y1 == y2, but with NA == NA
# and comparing factors as characters
eqna <- function(y1,y2){
  if(is.factor(y1)) y1 <- as.character(y1)
  if(is.factor(y2)) y2 <- as.character(y2)
  if(!is.null(dim(y1))&!is.null(dim(y2))){
    if(ncol(y1)==ncol(y2) & min(c(nrow(y1),nrow(y2)))==1){
      if(nrow(y1)==1){
        both.na <- eq <- matrix(NA,nrow(y2),ncol(y2))
      } else {
        both.na <- eq <- matrix(NA,nrow(y1),ncol(y1))
      }
      for(j in 1:ncol(eq)){
        eq[,j] <- y2[,j] == y1[,j]
        both.na[,j] <- is.na(y2[,j]) & is.na(y1[,j])
      }
      eq[is.na(eq)] <- both.na[is.na(eq)]
      eq <- apply(eq,1,all)
      return(eq)
    }
  }
  eq <- y1 == y2
  both.na <- is.na(y1) & is.na(y2)
  eq[is.na(eq)] <- both.na[is.na(eq)]
  return(eq)
}


# remove 'cols' from 'x' and ensure that the result is still a matrix or data.frame
# maintains "terms" attr
safe_rm_cols <- function(x,cols){
  trms <- attr(x,"terms")
  cnms <- colnames(x)[-cols]
  isdf <- is.data.frame(x)
  x <- x[,-cols]
  if(isdf & is.matrix(x)){
    temp <- data.frame(a=rep(NA,nrow(x)))
    temp[,1] <- x
    colnames(temp) <- cnms
    x <- temp
  }
  if(!is.matrix(x) & !is.data.frame(x)){
    if(isdf){
      x <- data.frame(x)
    } else {
      x <- matrix(x,length(x),1)
    }
    colnames(x) <- cnms
  }
  attr(x,"terms") <- trms
  return(x)
}


# create model matrix from frame and formula, set NAs to zero
na_model.matrix <- function(formula,fr){
  naac <- getOption("na.action")
  options(na.action = "na.pass")
  possfail <- tryCatch(mat <- model.matrix(formula,fr), error = function(e) e)
  options(na.action = naac)
  if(inherits(possfail,"error")) stop("Error in forming model matrix")
  mat[is.na(mat)] <- 0
  return(mat)
}


# 'r' is a one-row matrix or data.frame
# 'fr' is a data.frame or matrix with the same number of columns as 'r'
# eqna is called on each column and then the rows are checked
match_row <- function(r,fr){
  if(!is.matrix(r) & !is.data.frame(r)){
    temp <- data.frame(matrix(NA,1,length(r)))
    for(j in 1:length(r)) temp[,j] <- r[j]
    r <- temp
  }
  if((!is.matrix(fr) & !is.data.frame(fr))){
    stop("'r' and 'fr' must be matrices or data.frames")
  }
  if(ncol(r)!=ncol(fr)){
    fr <- fr[,1]
    if(ncol(r)!=ncol(fr)) stop("row and search table have different column numbers")
  }
  return(apply(sapply(1:ncol(r),function(n) eqna(r[,n],fr[,n])),1,all))
}


# returns reference grid for group and effect means and contrasts
make_refgrid <- function(gridlevs,trms){
  g <- expand.grid(gridlevs[trms])
  isnumr <- ispoly <- logical(length(trms))
  names(isnumr) <- names(ispoly) <- trms
  for(i in trms){
    ic <- attr(gridlevs[[i]],"trmclass")
    if(ic=="factor"){
      g[,i] <- factor(g[,i],ordered=FALSE)
      contrasts(g[,i]) <- attr(gridlevs[[i]],"contr")
    } else if(ic=="ordered"){
      g[,i] <- factor(g[,i],ordered=TRUE)
      contrasts(g[,i]) <- attr(gridlevs[[i]],"contr")
    } else if(ic=="poly"){
      isnumr[i] <- ispoly[i] <- TRUE
    } else if(ic=="numeric"){
      isnumr[i] <- TRUE
    }
  }
  if(sum(isnumr)>0){
    if(sum(isnumr)>1){
      all1 <- data.frame(matrix(1,1,sum(isnumr)))
      all0 <- data.frame(matrix(0,1,sum(isnumr)))
      all1 <- which(match_row(all1,g[,isnumr]))
      all0 <- which(match_row(all0,g[,isnumr]))
    } else {
      all1 <- which(eqna(g[,isnumr],1))
      all0 <- which(eqna(g[,isnumr],0))
    }
    g <- g[c(all0,all1),]
    if(!is.data.frame(g)){
      g <- data.frame(g)
      colnames(g) <- trms
    }
  }
  for(i in trms[ispoly]){
    matvals <- attr(gridlevs[[i]],"matvals")
    is0 <- g[,i]==0
    g[,i] <- matrix(0,nrow(g),ncol(matvals))
    g[is0,i] <- matrix(matvals[1,],sum(is0),ncol(matvals),byrow=TRUE)
    g[!is0,i] <- matrix(matvals[2,],sum(!is0),ncol(matvals),byrow=TRUE)
  }
  attr(g,"terms") <- terms(as.formula(paste("~",paste(trms,collapse="*"))))
  return(g)
}


# returns model matrix for highest order interaction of columns in 'g'
# setting NAs to zero.  If drop.lower=TRUE, then all columns except
# those pertaining to the highest order interaction are dropped
make_inter_mat <- function(g,drop.lower=TRUE){
  mat <- na_model.matrix(as.formula(paste("~",paste(colnames(g),collapse="*"))),g)
  if(drop.lower){
    dropcols <- which(attr(mat,"assign")!=max(attr(mat,"assign")))
    mat <- safe_rm_cols(mat,dropcols)
  }
  return(mat)
}


# converts x to an unordered factor and assigns sum contrasts
# with named columns rather than numbered columns for the dummy
# variables.  If the levels are "0" and "1" or "FALSE" and "TRUE",
# the 1/TRUE is set to be 1 for the dummy variable and the 0/FALSE
# is set to be -1 for the dummy variable.  Unused levels are dropped.
named_contr.sum <- function(x){
  x <- factor(x, ordered = FALSE)
  NL <- nlevels(x)
  L <- levels(x)
  if(NL == 2){
    if(all(c("0","1") %in% L)){
      x <- factor(x,levels=c("1","0"))
    } else if(all(c("FALSE","TRUE") %in% L)){
      x <- factor(x,levels=c("TRUE","FALSE"))
    }
    L <- levels(x)
  }
  contrasts(x) <- contr.sum(NL)
  colnames(contrasts(x)) <- L[-NL]
  return(x)
}


# converts x to an ordered factor and assigns orthognal
# polynomnial contrasts where the standard deviation of
# each column in the contrast matrix is 1.
scaled_contr.poly <- function(x){
  x <- factor(x, ordered = TRUE)
  NL <- nlevels(x)
  contr <- scale(contr.poly(NL))[,1:(NL-1)]
  colnames(contr) <- colnames(contr.poly(NL))
  rownames(contr) <- levels(x)
  contrasts(x) <- contr
  return(x)
}


# returns the type of contrasts assigned to a factor x
contr_type <- function(x){
  contr <- contrasts(x)
  NL <- nlevels(x)
  type <- "custom"
  if(all(contr==contr.helmert(NL))){ type <- "helmert"
  } else if(all(contr==contr.poly(NL))){ type <- "poly"
  } else if(all(contr==scale(contr.poly(NL))[,1:(NL-1)])){ type <- "scaled poly"
  } else if(all(contr==contr.sum(NL))){ type <- "sum"
  } else {
    for(i in 1:NL){
      if(all(contr==contr.treatment(NL,base=i))) type <- "treatment"
    }
  }
  return(type)
}


# 'fr' is a model.frame and 'v' is a character vector of columns
# for which a highest-order interaction term is to be evaluated,
# where at least one member of 'v' has NA values.
# First, unordered factors in 'v' are subsetted, and checked for
# NA values in their interaction that do not exist in their main
# effects.  If any levels in the main effect need to be dropped
# for the interaction term, they are dropoped and new sum contrasts
# are set. If 'gridlevs' is supplied, then a reference grid is also
# created, dropping unattested combinations, and creating a copy
# with the newly created contrasts.  returns a list with the following elements:
#    'changed' a logical indicating if unordered factor contrasts were changed
#    'refgrid' a reference grid from make_refgrid() with unattested combinations of
#              unordered factors removed (NULL if 'gridlevs' is NULL)
#    'refgrid.na' the same grid as 'refgrid', but with new contrasts applied
#                 (NULL if 'gridlevs' is NULL)
#    'x' the columns 'v' in 'fr' with new contrasts applied
#    'contr' a list of (possibly altered) contrasts for unordered factors
#            used in the creation of the interaction term.  If no unordered
#            factors in 'v', then NULL
na_interaction <- function(fr,v,gridlevs=NULL){
  nona <- function(y){
    y.nona <- y[apply(is.na(y),1,function(n) !any(n)),]
    newlevs <- lapply(y.nona,unique)
    y[!(y[,1] %in% newlevs[[1]]),1] <- NA
    y[!(y[,2] %in% newlevs[[2]]),2] <- NA
    y[apply(y,1,function(n) any(is.na(n))),1:2] <- NA
    return(y)
  }

  na_interaction_2_unordered <- function(ux){
    new.ux <- nona(ux)
    coll <- xtabs(~.,unique(new.ux))
    drop.f1 <- names(which(rowSums(coll)==1))
    drop.f2 <- names(which(colSums(coll)==1))
    new.ux[new.ux[,1] %in% drop.f1,1] <- NA
    new.ux[new.ux[,2] %in% drop.f2,2] <- NA
    new.ux <- nona(new.ux)
    
    changed <- !is.logical(all.equal(lapply(ux,function(n) sort(unique(n))),
      lapply(new.ux,function(n) sort(unique(n)))))
    
    return(list(changed = changed, ux = new.ux))
  }

  na_interaction_unordered <- function(x,g){
    nc <- ncol(x)
    cnms <- colnames(x)
    ux.temp <- unique(x)
    changed <- FALSE
    
    for(j in cnms){
      x[,j] <- as.character(x[,j])
      ux.temp[,j] <- as.character(ux.temp[,j])
    }
    ux.old <- ux.new <- ux.temp
    
    inter <- na_interaction_2_unordered(ux.temp[,1:2])
    changed <- changed | inter$changed
    ux.new[,1] <- inter$ux[,1]
    ux.new[,2] <- inter$ux[,2]
    ux.temp[,2] <- NA
    temp <- apply(ux.new[,1:2],1,function(n) !any(is.na(n)))
    ux.temp[temp,2] <- apply(ux.new[temp,1:2],1,function(n) paste(n,collapse="."))
    ux.temp <- ux.temp[,-1]

    currcol <- 3
    nc <- nc - 1
    while(nc > 1){
      inter <- na_interaction_2_unordered(ux.temp[,1:2])
      changed <- changed | inter$changed
      ux.new[,currcol] <- inter$ux[,2]
      temp <- is.na(inter$ux[,2])
      ux.new[temp,1:(currcol-1)] <- NA
      ux.temp[temp,1] <- apply(ux.temp[temp,1:2],1,function(n) paste(n,collapse="."))
      ux.temp <- ux.temp[,-2]
      nc <- nc - 1
      currcol <- currcol + 1
    }
    
    for(j in cnms){
      temp <- sort(unique(ux.new[,j]))
      x[!(x[,j] %in% temp),j] <- NA
    }
    
    ux.temp <- unlist(lapply(ux.new,function(n) length(sort(unique(n)))))
    if(any(ux.temp<1)){
      inter <- paste(colnames(x),collapse=":")
      ux.temp <- colnames(x)[which(ux.temp<1)]
      if(length(ux.temp)==2) ux.temp <- paste(ux.temp,collapse= " and ")
      if(length(ux.temp)>2){
        ux.temp <- c(paste(ux.temp[-length(ux.temp)],collapse=", "),ux.temp[length(ux.temp)])
        ux.temp <- paste(ux.temp,collapse=", and ")
      }
      mess <- paste("Due to NAs, the factor(s)",ux.temp,"are singular in the interaction",inter)
      stop(mess)
    }
    
    for(j in cnms){
      x[,j] <- named_contr.sum(x[,j])
      ux.old[,j] <- named_contr.sum(ux.old[,j])
      ux.new[,j] <- named_contr.sum(ux.new[,j])
    }
    
    if(is.null(g)){
      g.na <- NULL
    } else {
      keep <- apply(g[,cnms],1,function(n) any(match_row(n,ux.old)))
      g <- g[keep,]
      g.na <- g
      tona <- apply(g.na[,cnms],1,function(n) !any(match_row(n,ux.new)))
      g.na[tona,] <- NA
      for(j in cnms){
        g.na[,j] <- droplevels(g.na[,j])
        contrasts(g.na[,j]) <- contrasts(ux.new[,j])
      }
    }
    
    contr <- list()
    for(j in cnms){
      contr[[j]] <- contrasts(ux.new[,j])
      dropped <- levels(ux.old[,j])[!(levels(ux.old[,j]) %in% levels(ux.new[,j]))]
      for(i in dropped){
        contr[[j]] <- rbind(contr[[j]],0)
        rownames(contr[[j]])[nrow(contr[[j]])] <- i
      }
      if(any(is.na(ux.old[,j]))){
        contr[[j]] <- rbind(contr[[j]],0)
        rownames(contr[[j]])[nrow(contr[[j]])] <- NA
      }
    }
    
    return(list(changed = changed, x = x, ux.old = ux.old,
      ux.new = ux.new, contr = contr, refgrid = g, refgrid.na = g.na))
  }
  
  if(length(v)<2) stop("only one variable provided")
  if(!is.null(gridlevs)){
    g <- make_refgrid(gridlevs,v)
  } else g <- NULL

  unord <- unlist(lapply(fr,function(n) !is.ordered(n) & is.factor(n)))[v]
  
  if(sum(unord)>1){
    nai <- na_interaction_unordered(fr[,v[unord]],g)
    fr[,v[unord]] <- nai$x
    nai$x <- fr[,v]
  } else {
    nai <- list(changed = FALSE, x = fr[,v], ux.old = NULL, ux.new = NULL,
      contr = NULL, refgrid = g, refgrid.na = g)
  }
  attr(nai$x,"terms") <- terms(as.formula(paste("~",paste(colnames(nai$x),collapse="*"))))
  
  return(nai)
}


# returns fixed effects model matrix with NA values set to zero;
# when variables which have NAs are involved in interactions,
# the interaction terms are evaluated with na_interaction(),
# and any changes to the contrasts are added to the log
na_fixef_mat <- function(fixed_frame,varmat,varlog,gridlevs){
  hasna <- unlist(lapply(fixed_frame,function(n) any(is.na(n))))
  if(hasna[1]) stop("NAs in response variable")
  hasna <- hasna[-1]
  if(!any(hasna)) stop("fixef_na set to TRUE in bmerControl but no fixed effects have NAs")
  
  inter <- varmat.na.inter <- colSums(varmat) > 1
  for(j in 1:ncol(varmat)){
    temp <- any(varmat[hasna,j]>0)
    varmat.na.inter[j] <- varmat.na.inter[j] & temp
  }
  
  modmat <- matrix(1,nrow(fixed_frame),1)
  colnames(modmat)[1] <- "(Intercept)"
  refgrids <- list()
  refmats <- list()
  for(j in colnames(varmat)){
    v <- rownames(varmat)[varmat[,j]>0]
    if(varmat.na.inter[j]){
      ## interaction involving at least one variable with NAs
      nai <- na_interaction(fixed_frame,v,gridlevs)
      refgrids[[j]] <- nai$refgrid
      refmats[[j]] <- make_inter_mat(nai$refgrid.na)
      tempmat <- make_inter_mat(nai$x)
      if(nai$changed){
        cols <- colnames(refmats[[j]])
        separator <- "_"
        new.cols <- paste("nai",cols,sep=separator)
        while(any(new.cols %in% colnames(modmat))){  # shouldn't need to execute
          separator <- paste(separator,"_",sep="")
          new.cols <- paste("nai",cols,sep=separator)
        }
        colnames(tempmat) <- colnames(refmats[[j]]) <- new.cols
        varlog$nainter[[j]] <- nai$contr
      }
      modmat <- cbind(modmat,tempmat)
      
    } else {
      ## interaction involving only variables without NAs or main effect
      refgrids[[j]] <- make_refgrid(gridlevs,v)
      refmats[[j]] <- make_inter_mat(refgrids[[j]])
      modmat <- cbind(modmat,
        make_inter_mat(safe_rm_cols(fixed_frame,which(!(colnames(fixed_frame) %in% v)))))
    }
  }
  
  refmats2 <- list()
  for(j in 1:length(refmats)){
    mumat <- modmat[1:nrow(refmats[[j]]),]
    mumat[,2:ncol(mumat)] <- 0
    mumat[,colnames(refmats[[j]])] <- refmats[[j]]
    if(inter[j]){
      lower <- colnames(attr(attr(refgrids[[j]],"terms"),"factors"))
      lower <- lower[-length(lower)]
      for(k in lower){
        kg <- refgrids[[k]]
        km <- refmats[[k]]
        jk <- safe_rm_cols(refgrids[[j]],which(!(colnames(refgrids[[j]]) %in% colnames(kg))))
        for(i in 1:nrow(kg)){
          jkrows <- which(match_row(kg[i,],jk))
          mumat[jkrows,colnames(km)] <- matrix(km[i,],length(jkrows),ncol(km),byrow=TRUE)
        }
      }
    }
    refmats2[[j]] <- mumat
  }
  refmats <- refmats2
  names(refmats) <- names(refgrids)
  
  for(j in 1:length(refmats)){
    mumat <- refmats[[j]]
    if(any(unlist(lapply(refgrids[[j]],is.numeric)))){
      mumat <- mumat[(nrow(mumat)/2+1):nrow(mumat),] - mumat[1:(nrow(mumat)/2),]
      if(!is.matrix(mumat)){
        mumat <- matrix(mumat,1,length(mumat))
        colnames(mumat) <- colnames(modmat)
      }
      for(k in 1:ncol(refgrids[[j]])){
        if(is.matrix(refgrids[[j]][,k])) refgrids[[j]][,k] <- 1
      }
      temp.cnms <- colnames(refgrids[[j]])
      refgrids[[j]] <- refgrids[[j]][(nrow(refgrids[[j]])/2+1):nrow(refgrids[[j]]),]
      if(!is.data.frame(refgrids[[j]])){
        refgrids[[j]] <- data.frame(refgrids[[j]])
        colnames(refgrids[[j]]) <- temp.cnms
      }
    }
    refmats[[j]] <- mumat
  }
  
  varlog$refgrids <- refgrids
  varlog$refmats <- refmats
  
  return(list(modmat = modmat, varlog = varlog))
}


# returns model matrix to be replicated for each member of random
# effect group 'r' with NA values set to zero.  If within-group factors
# have some levels for which there are no observations when 'r' is not NA,
# then a new factor is created with the name of the group prefixed and the
# new contrasts are added to the log
# TO DO: incorporate the newer utility functions
na_ranef_mat <- function(model_frame,varlog,col.class,r){
  fact <- col.class[-1] %in% c("factor","ordered")
  names(fact) <- names(col.class)[-1]
  
  for(f in names(fact)){
    if(fact[f]){
      tab <- xtabs(~model_frame[,r]+model_frame[,f])
      
      if(any(colSums(tab)==0) & !all(rowSums(tab)==apply(tab,1,max))){
        rfnas <- colnames(tab)[colSums(tab)==0]
        separator <- "_"
        rfname <- paste(r,f,sep="_")
        while(rfname %in% colnames(model_frame)){
          separator <- paste(separator,"_",sep="")
          rfname <- paste(r,f,sep=separator)
        }
        
        model_frame[,rfname] <- model_frame[,f]
        model_frame[model_frame[,rfname] %in% rfnas,rfname] <- NA
        model_frame[,rfname] <- droplevels(model_frame[,rfname])
        model_frame[,rfname] <- named_contr.sum(model_frame[,rfname])
        varlog$random[[r]]$newfacts[[rfname]] <- list(ordered = FALSE, nlevels = nlevels(model_frame[,rfname]),
          levels = levels(model_frame[,rfname]), contr_type = "sum", contr = contrasts(model_frame[,rfname]))

        rownames(varlog$random[[r]]$varmat)[rownames(varlog$random[[r]]$varmat)==f] <- rfname
        temp <- str_split(colnames(varlog$random[[r]]$varmat),":")
        temp <- lapply(temp,function(x) {x[x==f] <- rfname; return(x)})
        temp <- lapply(temp,function(x) paste(x,collapse=":"))
        temp <- unlist(temp)
        colnames(varlog$random[[r]]$varmat) <- temp
      }
    }
  }
  
  fact <- col.class[-1] %in% c("factor","ordered")
  btwn <- btwn_cont <- btwn_fact <- logical(length(fact))
  names(fact) <- names(btwn) <- names(btwn_cont) <- names(btwn_fact) <- rownames(varlog$random[[r]]$varmat)

  for(f in names(fact)){
    if("poly" %in% class(model_frame[,f])){
      tab <- xtabs(~model_frame[,r]+model_frame[,f][,1])
    } else {
      tab <- xtabs(~model_frame[,r]+model_frame[,f])
    }
    if(fact[f]){
      btwn_fact[f] <- all(rowSums(tab)==apply(tab,1,max))
    } else {
      btwn_cont[f] <- all(colSums(apply(tab,1,function(x){x!=0}))==1)
    }
    btwn[f] <- btwn_fact[f] | btwn_cont[f]
  }
  
  varlog$random[[r]]$btwn_fact <- names(btwn_fact)[btwn_fact]
  varlog$random[[r]]$btwn_cont <- names(btwn_cont)[btwn_cont]
  
  slopes <- apply(varlog$random[[r]]$varmat,2,function(x) { any(!btwn & x) & !any(btwn & x)  })
  slopes.inter <- (colSums(varlog$random[[r]]$varmat) > 1) & slopes
  slopes.main <- slopes & !slopes.inter
  
  if(any(slopes.inter)){
    varlog$random[[r]]$slopes <- colnames(varlog$random[[r]]$varmat)[slopes.main]
    mm <- model.matrix(as.formula(paste("~1+",paste(varlog$random[[r]]$slopes,collapse="+"))),model_frame)
    for(i in 1:length(slopes)){
      if(slopes.inter[i]){
        possfail <- tryCatch(nai <- na_interaction(model_frame,
          rownames(varlog$random[[r]]$varmat)[varlog$random[[r]]$varmat[,i]>0]), error = function(e) e)
        if(inherits(possfail,"error")){
          varlog$random[[r]]$nainter$dropped <- c(varlog$random[[r]]$nainter$dropped,names(slopes)[i])
        } else {
          varlog$random[[r]]$slopes <- c(varlog$random[[r]]$slopes,colnames(varlog$random[[r]]$varmat)[i])
          tempmat <- make_inter_mat(nai$x)
          if(nai$changed){
            cols <- colnames(tempmat)
            separator <- "_"
            new.cols <- paste("nai",cols,sep=separator)
            while(any(new.cols %in% colnames(mm))){  # shouldn't ever execute
              separator <- paste(separator,"_",sep="")
              new.cols <- paste("nai",cols,sep=separator)
            }
            colnames(tempmat) <- new.cols
            varlog$random[[r]]$nainter$changed[[colnames(varlog$random[[r]]$varmat)[i]]] <- nai$contr
          }
          mm <- cbind(mm,tempmat)
        }
      }
    }
  } else if(any(slopes.main)){
    varlog$random[[r]]$slopes <- colnames(varlog$random[[r]]$varmat)[slopes]
    mm <- model.matrix(as.formula(paste("~1+",paste(varlog$random[[r]]$slopes,collapse="+"))),model_frame)
  } else {
    varlog$random[[r]]$slopes <- character()
    mm <- model.matrix(~1,model_frame)
  }

  return(list(model_frame = model_frame, varlog = varlog, mm = mm))
}


#### check_variables
## Arguments
# 'formula' a formula passed to bmer() or bmerBuild()
# 'data' a data.frame
# 'family' either gaussian or binomial
# 'control' control parameters, must be class 'bmerControl'
#
## Value: a list with the following elements
# 'family' character string either "gaussian" or "binomial"
# 'fixed_frame' the fixed effects model frame
# 'fixed_varmat' the "factors" matrix of the fixed_frame terms object with the response row removed
# 'x' the fixed effect model matrix
# 'random_frame' the model frame from the right side of the "|" in the formula
# 'random_varmat' the "factors" matrix of the random_frame terms object
# 'model_frame' the fixed_frame and random_frame, with attributes removed
# 'z' a list with the random effects model matrix for each random effect group
# 'messages' a character vector of messages indicating what transformations were applied
# 'warnings' a character vector of warning messages indicating what aspects of the data are not checked to ensure the priors make sense
# 'varlog' a list with information about the variables (see the description for bmerBuild)
# 'stanmod' character string describing the aspects of the model required to write the code
# TO DO: incorporate the newer utility functions
check_variables <- function(formula, data, family, control){
  # check arguments
  if(!is.data.frame(data)) stop("'data' must be a data.frame")
  if(class(control)!="bmerControl") stop("invalid control")
  
  if(is.function(family)) family <- family()$family
  if(!is.character(family)) stop("invalid family")
  if(!(family %in% c("gaussian","binomial"))) stop("only gaussian and binomial families currently supported")
  
  if(class(formula)!="formula") stop("must supply a formula")
  formula <- str_split(paste(formula),"\\|")[[1]]
  if(length(formula)!=2) stop("formula should be specified as 'response ~ fixed_effects | random_effect_groups'")
  fixedform <- as.formula(formula[1])
  ranefgrps <- as.formula(paste("~",formula[2]))
  if(is.null(lhs(fixedform))) stop("formula has no response variable")


  # start logs
  checked <- list(family = family)
  warn <- mess <- character()
  
  if(control@scale_cont){
    mess <- ifelse(family=="gaussian","scaling response and covariates","scaling covariates")
  } else {
    warn <- ifelse(family=="gausian","WARNING: response and covaraites left on original scales","WARNING: covariates left on original scale")
  }
  if(control@set_contr){
    mess <- c(mess,"sum contrasts set for unordered factors","scaled orthogonal polynomial contrasts set for ordered factors")
  } else {
    warn <- c(warn,"WARNING: factor contrasts left at their defaults")
  }
  if(control@fixef_na){
    mess <- c(mess,"any NA values in the fixed effects will be set to zero (i.e. level slashing)")
  }
  if(control@ranef_na){
    mess <- c(mess, "random effect structures will only be applied to observations where the grouping factor is not NA")
  }

  varlog <- list(resp=list(),fact=list(),cont=list(),nainter=list(),random=list(),
    refgrids=list(),refmats=list())

  # get fixed effects model frame
  if(control@fixef_na){
    options(na.action = "na.pass")
  } else {
    options(na.action = "na.omit")
  }
  possfail <- tryCatch(fixed_frame <- model.frame(fixedform,data,drop.unused.levels=TRUE),error=function(e) e)
  if(inherits(possfail,"error")){
    stop("Error in forming fixed effects model frame from data and formula")
  }
  
  if(any(!is.na(str_locate(colnames(fixed_frame),pattern="offset\\(")))){
    stop("offsets not currently supported")
  }
  col.class <- attr(attr(fixed_frame,"terms"),"dataClasses")
  nvals <- unlist(lapply(fixed_frame,function(x){sum(!is.na(unique(x)))}))
  if(any(nvals==1)){
    nvals <- nvals[nvals==1]
    stop(paste("only one value in",names(nvals)))
  }

  
  # check response
  varlog$resp$Name <- colnames(fixed_frame)[1]
  if(family == "gaussian"){
    if(col.class[1] != "numeric") stop("family is gaussian but response is not numeric")
    if(nvals[1]==2) stop("response only has two values but family is gaussian")
    varlog$resp$Class <- "numeric"
    varlog$resp$DataMean <- mean(fixed_frame[,1])
    varlog$resp$DataSD <- sd(fixed_frame[,1])
    if(control@scale_cont & !is.scaled(fixed_frame[,1])){
      fixed_frame[,1] <- scale(fixed_frame[,1])[,1]
    }
    varlog$resp$RegScaled <- is.scaled(fixed_frame[,1])

  } else if(family == "binomial"){
    if(nvals[1]!=2) stop("response has more than two values but family is binomial")
    original <- fixed_frame[,1]
    fixed_frame[,1] <- factor(fixed_frame[,1],ordered=F)
    levs <- levels(fixed_frame[,1])
    varlog$resp$Class <- "binary factor coded as integer 0/1"
    fixed_frame[,1] <- as.integer(as.numeric(fixed_frame[,1])-1)
    col.class[1] <- "integer"
    varlog$resp$Levels <- levs
    if(any(paste(fixed_frame[,1])!=paste(original))){
      varlog$resp$Conversion <- paste(levs[1],"-> 0;",levs[2],"-> 1")
    }
  }


  # check fixed effects
  gridlevs <- list()
  for(f in colnames(fixed_frame)[-1]){
    if(nvals[f]==2 | col.class[f] %in% c("factor","character")){
      if(control@set_contr){
        fixed_frame[,f] <- named_contr.sum(fixed_frame[,f])
      } else if(col.class[f]!="factor"){
        fixed_frame[,f] <- factor(fixed_frame[,f])
      }
      varlog$fact[[f]] <- list(ordered = FALSE, nlevels = nlevels(fixed_frame[,f]), levels = levels(fixed_frame[,f]),
        contr_type = contr_type(fixed_frame[,f]), contr = contrasts(fixed_frame[,f]), hasNAs = any(is.na(fixed_frame[,f])))
      col.class[f] <- "factor"
      gridlevs[[f]] <- varlog$fact[[f]]$levels
      if(any(is.na(fixed_frame[,f]))) gridlevs[[f]] <- c(gridlevs[[f]],NA)
      attr(gridlevs[[f]],"trmclass") <- "factor"
      attr(gridlevs[[f]],"contr") <- varlog$fact[[f]]$contr

    } else if(col.class[f]=="ordered"){
      if(control@set_contr){
        fixed_frame[,f] <- scaled_contr.poly(fixed_frame[,f])
      }
      varlog$fact[[f]] <- list(ordered = TRUE, nlevels = nlevels(fixed_frame[,f]), levels = levels(fixed_frame[,f]),
        contr_type = contr_type(fixed_frame[,f]), contr = contrasts(fixed_frame[,f]), hasNAs = any(is.na(fixed_frame[,f])))
      gridlevs[[f]] <- varlog$fact[[f]]$levels
      if(any(is.na(fixed_frame[,f]))) gridlevs[[f]] <- c(gridlevs[[f]],NA)
      attr(gridlevs[[f]],"trmclass") <- "ordered"
      attr(gridlevs[[f]],"contr") <- varlog$fact[[f]]$contr

    } else if(col.class[f]=="numeric"){
      varlog$cont[[f]] <- list(DataMean = mean(fixed_frame[,f]), DataSD = sd(fixed_frame[,f]), hasNAs = any(is.na(fixed_frame[,f])))
      if(control@scale_cont & !is.scaled(fixed_frame[,f])){
        fixed_frame[,f] <- scale(fixed_frame[,f])[,1]
      }
      varlog$cont[[f]]$RegScaled <- is.scaled(fixed_frame[,f])
      if(nvals[f]<=5) warn <- c(warn,paste("WARNING:",
        f,"is numeric but only has",nvals[f],"unique values; better coded as ordered factor?"))
      gridlevs[[f]] <- c(0,1)
      attr(gridlevs[[f]],"trmclass") <- "numeric"

    } else if("poly" %in% class(fixed_frame[,f])){
      degree <- ncol(fixed_frame[,f])
      coefs <- attr(fixed_frame[,f],"coefs")
      if(is.null(coefs)){
        mu <- mean(fixed_frame[,f][,1])
        sdev <- sd(fixed_frame[,f][,1])
      } else {
        mu <- coefs$alpha[1]
        sdev <- sqrt(coefs$norm2[3] * nrow(fixed_frame) / (nrow(fixed_frame)-1) / coefs$norm2[2])
      }
      gridlevs[[f]] <- c(mu, mu + sdev)
      matvals <- poly(gridlevs[[f]], degree = degree, coefs = coefs, raw = is.null(coefs))
      attr(matvals,"coefs") <- attr(matvals,"degree") <- attr(matvals,"class") <- NULL
      for(d in 1:degree){
        varlog$cont[[paste(f,d,sep="")]] <- list(DataMean = mean(fixed_frame[,f][,d]), DataSD = sd(fixed_frame[,f][,d]),
          coefs = coefs, degree = degree, hasNAs = any(is.na(fixed_frame[,f][,d])))
        if(control@scale_cont & !is.scaled(fixed_frame[,f][,d])){
          fixed_frame[,f][,d] <- scale(fixed_frame[,f][,d])[,1]
          matvals[,d] <- (matvals[,d] - varlog$cont[[paste(f,d,sep="")]]$DataMean) / varlog$cont[[paste(f,d,sep="")]]$DataSD
        }
        varlog$cont[[paste(f,d,sep="")]]$RegScaled <- is.scaled(fixed_frame[,f][,d])
      }
      gridlevs[[f]] <- c(0,1)
      attr(gridlevs[[f]],"matvals") <- matvals
      attr(gridlevs[[f]],"trmclass") <- "poly"

    } else {
      stop(paste(f," has an unsupported class '",class(fixed_frame[,f]),"'",sep=""))
    }
  }
  attr(attr(fixed_frame,"terms"),"dataClasses") <- col.class

  
  # check for collinearity, violations of the interaction hierarchy, and lack of intercept; get X
  varmat <- attr(attr(fixed_frame,"terms"),"factors")[-1,]
  if(!is.matrix(varmat)){ # if only one predictor
    varmat <- matrix(varmat)
    colnames(varmat) <- colnames(attr(attr(fixed_frame,"terms"),"factors"))
    rownames(varmat) <- rownames(attr(attr(fixed_frame,"terms"),"factors"))[-1]
  }
  checked$fixed_frame <- fixed_frame
  checked$fixed_varmat <- varmat
  
  if(any(varmat>1)) stop("model violates interaction hierarchy (not currently supported and not in general recommended)")
  
  if(!control@fixef_na){
    fixedform <- as.formula(paste("~",paste(colnames(varmat),collapse="+")))
    fixed_mat <- model.matrix(fixedform,fixed_frame)
    if(attr(fixed_mat,"assign")[1]!=0) stop("zero-intercept models not currently supported")
    for(j in colnames(varmat)){
      v <- rownames(varmat)[varmat[,j]>0]
      varlog$refgrids[[j]] <- make_refgrid(gridlevs,v)
      tempmat <- make_inter_mat(varlog$refgrids[[j]],FALSE)
      mumat <- fixed_mat[1:nrow(tempmat),]
      mumat[,2:ncol(mumat)] <- 0
      mumat[,colnames(tempmat)] <- tempmat
      if(any(unlist(lapply(varlog$refgrids[[j]],is.numeric)))){
        mumat <- mumat[(nrow(mumat)/2+1):nrow(mumat),] - mumat[1:(nrow(mumat)/2),]
        if(!is.matrix(mumat)){
          mumat <- matrix(mumat,1,length(mumat))
          colnames(mumat) <- colnames(fixed_mat)
        }
        for(k in 1:ncol(varlog$refgrids[[j]])){
          if(is.matrix(varlog$refgrids[[j]][,k])) varlog$refgrids[[j]][,k] <- 1
        }
        temp.cnms <- colnames(varlog$refgrids[[j]])
        varlog$refgrids[[j]] <- varlog$refgrids[[j]][(nrow(varlog$refgrids[[j]])/2+1):nrow(varlog$refgrids[[j]]),]
        if(!is.data.frame(varlog$refgrids[[j]])){
          varlog$refgrids[[j]] <- data.frame(varlog$refgrids[[j]])
          colnames(varlog$refgrids[[j]]) <- temp.cnms
        }
      }
      varlog$refmats[[j]] <- mumat
    }
  } else {
    fixed_mat <- na_fixef_mat(fixed_frame,varmat,varlog,gridlevs)
    varlog <- fixed_mat$varlog
    fixed_mat <- fixed_mat$modmat
  }

  possfail <- tryCatch(temp <- solve(t(fixed_mat) %*% fixed_mat), error = function(e) e)
  if(inherits(possfail, "error")) stop("Collinearity in the fixed effects")
  checked$x <- fixed_mat

  
  # check random effects; get Z
  if(control@ranef_na){
    options(na.action = "na.pass")
  } else {
    options(na.action = "na.omit")
  }
  possfail <- tryCatch(random_frame <- model.frame(ranefgrps,data,drop.unused.levels=TRUE),error=function(e) e)
  if(inherits(possfail,"error")){
    stop("Error in forming random effects model frame from data and formula")
  }
  
  nvals <- unlist(lapply(random_frame,function(x){length(unique(x))}))
  if(any(nvals==1)){
    nvals <- nvals[nvals==1]
    stop(paste("only one value in",names(nvals)))
  }
  if(nrow(random_frame)!=nrow(fixed_frame)){
    stop("There are NA values in either the fixed or random effects but the appropriate bmerControl (fixef_na or ranef_na) parameter is set to FALSE")
  }
  
  if(control@ranef_na){
    hasna <- unlist(lapply(random_frame,function(n) any(is.na(n))))
    if(!any(hasna)) stop("ranef_na set to TRUE in bmerControl but no random effects have NAs")
  }  
  
  for(r in colnames(random_frame)){
    random_frame[,r] <- factor(random_frame[,r],ordered=FALSE)
  }
  checked$random_frame <- random_frame
  checked$random_varmat <- random <- attr(attr(random_frame,"terms"),"factors")

  if(any(rownames(random) %in% rownames(varmat))) stop("factors listed in both fixed and random effects")
  model_frame <- fixed_frame
  model_frame[,colnames(random_frame)] <- random_frame
  attr(model_frame,"terms") <- NULL
  checked$model_frame <- model_frame

  startcol <- ncol(random_frame) + 1
  for(r in 1:ncol(random)){
    random_frame[,ncol(random_frame)+1] <- interaction(random_frame[,rownames(random)[random[,r]>0]],drop=T)
  }
  random_frame <- data.frame(random_frame[,startcol:ncol(random_frame)])
  colnames(random_frame) <- colnames(random)
  random <- colnames(random_frame)
  
  model_frame <- fixed_frame
  for(r in random){
    model_frame[,r] <- random_frame[,r]
    varlog$random[[r]] <- list(nlevels=nlevels(model_frame[,r]), levels = levels(model_frame[,r]),
      hasNAs=any(is.na(model_frame[,r])),newfacts=list(),nainter=list(changed=list(),dropped=character()),
      varmat = varmat)
  }

  checked$z <- list()
  for(r in random){
    if(control@ranef_na){
      ranmat <- na_ranef_mat(model_frame,varlog,col.class,r)
      model_frame <- ranmat$model_frame
      varlog <- ranmat$varlog
      mm <- ranmat$mm
    } else {
      fact <- col.class[-1] %in% c("factor","ordered")
      btwn <- btwn_cont <- btwn_fact <- logical(length(fact))
      names(fact) <- names(btwn) <- names(btwn_cont) <- names(btwn_fact) <- rownames(varlog$random[[r]]$varmat)

      for(f in names(fact)){
        if("poly" %in% class(model_frame[,f])){
          tab <- xtabs(~model_frame[,r]+model_frame[,f][,1])
        } else {
          tab <- xtabs(~model_frame[,r]+model_frame[,f])
        }
        if(fact[f]){
          btwn_fact[f] <- all(rowSums(tab)==apply(tab,1,max))
        } else {
          btwn_cont[f] <- all(colSums(apply(tab,1,function(x){x!=0}))==1)
        }
        btwn[f] <- btwn_fact[f] | btwn_cont[f]
      }
    
      varlog$random[[r]]$btwn_fact <- names(btwn_fact)[btwn_fact]
      varlog$random[[r]]$btwn_cont <- names(btwn_cont)[btwn_cont]
    
      slopes <- apply(varlog$random[[r]]$varmat,2,function(x) { any(!btwn & x) & !any(btwn & x)  })
      varlog$random[[r]]$slopes <- colnames(varlog$random[[r]]$varmat)[slopes]
    
      if(any(slopes)){
        mm <- model.matrix(as.formula(paste("~1+",paste(varlog$random[[r]]$slopes,collapse="+"))),model_frame)
      } else {
        mm <- model.matrix(~1,model_frame)
      }
    }
    
    mm[is.na(mm)] <- 0
    varlog$random[[r]]$cols <- colnames(mm)
    varlog$random[[r]]$Q <- ncol(mm)
    ff <- model.matrix(as.formula(paste("~0+",r)),checked$random_frame)
    ff[is.na(ff)] <- 0
    ff <- ff[,colSums(ff)>0]
    checked$z[[r]] <- Matrix::t(KhatriRao(t(ff),t(mm)))
  }

  Q <- unlist(lapply(varlog$random,function(x) x$Q))
  varlog$random <- varlog$random[c(which(Q>1),which(Q==1))]
  checked$z <- checked$z[c(which(Q>1),which(Q==1))]

  stanmod <- family
  if((control@estimate_scale_beta=="P" & ncol(checked$x)>=10)  |  control@estimate_scale_beta=="yes"){
    stanmod[2] <- "T"
  } else {
    stanmod[2] <- "F"
  }
  if(control@nu_beta==1){
    stanmod[3] <- "T"
  } else {
    stanmod[3] <- "F"
  }
  stanmod[4] <- paste(sum(Q>1))
  stanmod[5] <- paste(sum(Q==1))
  checked$stanmod <- paste(stanmod,collapse="_")
  
  # return a list of the information needed to build the model
  checked$messages <- mess; checked$warnings <- warn; checked$varlog <- varlog
  return(checked)
}
