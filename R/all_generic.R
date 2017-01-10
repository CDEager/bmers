
#' Show a \code{bmerControl}.
#' @param object A \code{bmerControl}.
#' @rdname bmerControl-class
#' @export
setMethod("show",signature(object="bmerControl"),function(object){
  Name <- slotNames(object)
  Value <- NULL
  for(i in Name){
    Value <- c(Value,paste(slot(object,i)))
  }
  summ <- data.frame(Name,Value)
  cat(paste("bmerControl (bmers version ",bmers_version(),")\n",sep=""))
  print(summ,row.names=FALSE,right=FALSE)
})

#' Show a \code{bmerBuildSummary}.
#' @param object A \code{bmerBuildSummary}.
#' @rdname bmerBuildSummary-class
#' @export
setMethod("show",signature(object="bmerBuildSummary"),function(object){
	w <- getOption("width")
	options(width=200)
  intermess <- FALSE
	cat(paste("bmerBuild for model '",object$name,"'",sep=""),"","Call:",sep="\n")
	show(object$call)
	cat("",paste("Observations:",object$N),paste("Parameters (fixed and random):",object$PQ),"",sep="\n")

  y <- data.frame(Family=object$family)
	for(i in 1:length(object$response)){
		if(length(object$response[[i]])>1){
			y[1,names(object$response)[i]] <- paste(object$response[[i]],collapse=", ")
		} else {
			y[1,names(object$response)[i]] <- object$response[[i]]
		}
		if(is.numeric(y[1,ncol(y)])) y[1,ncol(y)] <- format(round(y[1,ncol(y)],3),nsmall=3)
		y[1,ncol(y)] <- paste(y[1,ncol(y)],"   ",sep="")
	}
	colnames(y)[1] <- "  Family"
	colnames(y) <- paste(colnames(y),"   ",sep="")
	cat("Response:",sep="\n")
	print(y,row.names=FALSE,right=FALSE)
	
	cat("","Fixed Effects:","   Variables",sep="\n")
	temp <- data.frame(Factor="",Ordered="",Levels="",Contrasts="",ContainsNAs1="",
		Covariate="",RegScaled="",DataMean="",DataSD="",ContainsNAs2="",stringsAsFactors=FALSE)
  colnames(temp)[c(5,10)] <- "Contains NAs"
	for(i in 1:max(c(length(object$fixed$factors),length(object$fixed$covariates)))){
		if(i <= length(object$fixed$factors)){
			temp[i,1] <- names(object$fixed$factors)[i]
			temp[i,2] <- object$fixed$factors[[i]]$ordered
			temp[i,3] <- object$fixed$factors[[i]]$nlevels
			temp[i,4] <- object$fixed$factors[[i]]$contr_type
      temp[i,5] <- object$fixed$factors[[i]]$hasNAs
		} else {
			temp[i,1:5] <- ""
		}
		if(i <= length(object$fixed$covariates)){
			temp[i,6] <- names(object$fixed$covariates)[i]
			temp[i,7] <- object$fixed$covariates[[i]]$RegScaled
			temp[i,8] <- format(round(object$fixed$covariates[[i]]$DataMean,3),nsmall=3)
			temp[i,9] <- format(round(object$fixed$covariates[[i]]$DataSD,3),nsmall=3)
      temp[i,10] <- object$fixed$covariates[[i]]$hasNAs
		} else {
			temp[i,6:10] <- ""
		}
	}
	if(temp[1,1]==""){
		temp <- temp[,6:10]
    if(!object$control@fixef_na) temp <- temp[,1:4]
	} else if(temp[1,6]==""){
		temp <- temp[,1:5]
    if(!object$control@fixef_na) temp <- temp[,1:4]
	} else if(!object$control@fixef_na){
    temp <- temp[,-c(5,10)]
  }
	for(j in 1:ncol(temp)){
		temp[,j] <- paste(temp[,j],"  ",sep="")
	}
	temp[,1] <- paste("     ",temp[,1],sep="")
	colnames(temp)[1] <- paste("     ",colnames(temp)[1],sep="")
	colnames(temp) <- paste(colnames(temp),"  ",sep="")
	print(temp,row.names=FALSE,right=FALSE)
	
	cat("\n")
	cols <- object$fixed$terms
  cols[cols %in% names(object$fixed$nainter)] <- paste(
    cols[cols %in% names(object$fixed$nainter)],"*",sep="")
	cols <- paste("     ",cols,sep="")
	if(length(cols)%%5 > 0) cols <- c(cols,rep("",5-(length(cols)%%5)))
	colframe <- matrix(cols,length(cols)/5,5)
	colframe <- data.frame(colframe,stringsAsFactors=FALSE)
	colnames(colframe) <- c("  Terms",rep("",4))
	print(colframe,row.names=FALSE,right=FALSE)
  if(length(object$fixed$nainter)>0) intermess <- TRUE
	
	cat("","Random Effects:",sep="\n")
	random <- names(object$random)
	
	for(r in random){
		temp <- data.frame(
      Group = r,
			Levels = object$random[[r]]$nlevels,
			"Between-Group Factors" = paste(object$random[[r]]$btwn_fact,collapse = ", "),
			"Between-Group Covariates" = paste(object$random[[r]]$btwn_cont, collapse = ", "),
      "Contains NAs" = object$random[[r]]$hasNAs,
			stringsAsFactors=FALSE)
		colnames(temp)[1] <- "  Group"
    colnames(temp)[3:5] <- c("Between-Group Factors","Between-Group Covariates","Contains NAs")
    if(temp[1,3]=="") temp[1,3] <- "N/A"
    if(temp[1,4]=="") temp[1,4] <- "N/A"
		colnames(temp) <- paste(colnames(temp),"   ",sep="")
		for(j in 1:5){
			temp[,j] <- paste(temp[,j],"   ",sep="")
		}
		temp[,1] <- paste("  ",temp[,1],sep="")
    if(!object$control@ranef_na) temp <- temp[,1:4]
		print(temp,row.names=FALSE,right=FALSE)
    
    if(length(object$random[[r]]$newfacts)>0){
      cat("\n     Due to the distribution of the NAs in random effect structure, the following unordered factor(s) were created:")
      for(i in 1:length(object$random[[r]]$newfacts)){
        cat(paste("\n      ",names(object$random[[r]]$newfacts)[i],"with contrasts:\n"))
        nfcontr <- object$random[[r]]$newfacts[[i]]$contr
        rownames(nfcontr) <- paste("        ",rownames(nfcontr))
        print(nfcontr)
      }
    }

    if(length(object$random[[r]]$slopes) == 0){
      cat("\n     The maximal random effect structure is random intercepts only\n")
    } else {
      cols <- object$random[[r]]$slopes
      cols[cols %in% names(object$random[[r]]$nainter$changed)] <- paste(
        cols[cols %in% names(object$random[[r]]$nainter$changed)],"*",sep="")
      cols <- paste("     ",cols,sep="")
      if(length(cols)%%5 > 0) cols <- c(cols,rep("",5-(length(cols)%%5)))
      colframe <- matrix(cols,length(cols)/5,5)
      colframe <- data.frame(colframe,stringsAsFactors=FALSE)
      colnames(colframe) <- rep("",5)
      cat("\n     The maximal random effect structure is random intercepts and random slopes for:")
      print(colframe,row.names=FALSE,right=FALSE)
      if(length(object$random[[r]]$nainter$changed)>0) intermess <- TRUE
    }
		cat("\n")
	}
	
  mess <- c(object$messages,object$warnings)
  cat("Messages:",paste("  ",mess),"",sep="\n")
  
  if(intermess) cat("* Interaction terms marked with an asterisk have different contrasts than the main effects involved",
  "\ndue to the distribution of NAs.  In the regression output, these interaction terms are preceded by 'nai_'",
  "\nto indicate the difference in interpretation.  Call get_contrasts to see the contrast changes.\n\n", sep="")

	options(width=w)
})

#' Show a \code{\link{build_summary}} for a \linkS4class{bmerBuild}.
#' @param object A \linkS4class{bmerBuild}.
#' @export
setMethod("show",signature(object="bmerBuild"),function(object){
  show(build_summary(object))
})

#' Show a \code{bmerFitSummary}.
#' @param object A \code{bmerFitSummary}.
#' @rdname bmerFitSummary-class
#' @export
setMethod("show",signature(object="bmerFitSummary"),function(object){
	w <- getOption("width")
	options(width=200)
	cat(paste("bmerFit for model '",object$build$name,"'",sep=""),"","Call:",sep="\n")
	show(object$call)
  
  y <- paste("Response:",object$build$response$Name)
  if("Conversion" %in% names(object$build$response)){
    y <- paste(y," (",object$build$response$Conversion,")",sep="")
  }
  
	cat("",paste("Family:",object$build$family),y,paste("Observations:",object$build$N),
    paste("Parameters (fixed and random):",object$build$PQ),"",sep="\n")

  cat("Random Effects:\n")
  rdf <- list()
  for(r in names(object$random)){
    Q <- length((rsd <- object$random[[r]]$sd))
    rdf[[r]] <- data.frame(Groups = rep("",Q), Effect = names(rsd), SD = format(round(rsd,4),nsmall=4), stringsAsFactors = FALSE)
    rdf[[r]][1,1] <- paste(r," (",nrow(object$random[[r]]$ranef),")",sep="")
    if(Q > 1){
      rcor <- format(round(data.frame(object$random[[r]]$cor),2),nsmall = 2)
      rcor[upper.tri(rcor,diag=T)] <- ""
      for(i in 1:ncol(rcor)) rcor[,i] <- str_pad(rcor[,i],side="left",width=5,pad=" ")
      rdf[[r]][,4:(ncol(rcor)+3)] <- rcor
      rdf[[r]] <- rdf[[r]][,-ncol(rdf[[r]])]
      colnames(rdf[[r]])[4:ncol(rdf[[r]])] <- ""
    }
  }
  cols <- max(unlist(lapply(rdf,ncol)))
  for(r in names(object$random)){
    if(ncol(rdf[[r]])<cols){
      rdf[[r]][,(ncol(rdf[[r]])+1):cols] <- ""
      colnames(rdf[[r]])[4:ncol(rdf[[r]])] <- ""
    }
  }
  rdf <- do.call(rbind,args=rdf)
  if(object$build$family == "gaussian"){
    rdf[nrow(rdf)+1,] <- ""
    rdf[nrow(rdf),1] <- "Residual"
    rdf[nrow(rdf),3] <- format(round(object$sigma_residual,4),nsmall=4)
  }
  if(cols > 3){
    colnames(rdf)[4] <- " Corr"
    if(cols > 4) colnames(rdf)[5:cols] <- ""
  }
  print(rdf,row.names=FALSE,right=FALSE)
  
  fdf <- data.frame(Effect = rownames(object$fixed$coef), object$fixed$coef[,-2])
  colnames(fdf) <- c("Effect","Mean","SD","2.5%","97.5%","Neff","Rhat","P(>0)")
  fdf[,c(2:5,7)] <- round(fdf[,c(2:5,7)],4)
  for(i in c(2,4,5)){
    if(any(fdf[,i]<0)){
      colnames(fdf)[i] <- paste(" ",colnames(fdf)[i],sep="")
    }
  }
  fdf[,c(2:5,7)] <- format(fdf[,c(2:5,7)],nsmall=4)
  fdf[,6] <- format(round(fdf[,6],1),nsmall=1)
  fdf[,8] <- round(fdf[,8],4)
  is0 <- fdf[,8]==0
  is1 <- fdf[,8]==1
  fdf[,8] <- format(fdf[,8],nsmall=4)
  if(any(is0) | any(is1)){
    if(any(!is0 & !is1)){
      fdf[!is0 & !is1,8] <- paste("  ",fdf[!is0 & !is1,8],sep="")
    }
    if(any(is0)){
      fdf[is0,8] <- "< 0.0001"
    }
    if(any(is1)){
      fdf[is1,8] <- "> 0.9999"
    }
  }
  cat("\nFixed Effects:\n")
  print(fdf,row.names=FALSE,right=FALSE)
  cat("\n")

  warn <- c(object$build$warnings,object$warnings)
  if(length(warn)>0) cat(warn,"",sep="\n")
  
  options(width=w)
})

#' Show a \code{\link{fit_summary}} for a \linkS4class{bmerFit}.
#' @param object A \linkS4class{bmerFit}.
#' @export
setMethod("show",signature(object="bmerFit"),function(object){
  show(fit_summary(object))
})
