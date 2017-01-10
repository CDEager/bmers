bmers_version <- function(){
  return("0.1.0")
}

#' Get the \linkS4class{bmerBuild} stored in a \linkS4class{bmerFit}.
#'
#' If the function is called on a \linkS4class{bmerFit} then the \linkS4class{bmerBuild}
#' used to fit the model is returned. If the function is called on a \linkS4class{bmerBuild},
#' then the same build is returned.
#'
#' @param object A \linkS4class{bmerFit} or \linkS4class{bmerBuild}.
#' @export
get_build <- function(object){
  if(class(object)=="bmerFit"){
    return(object@build)
  } else if(class(object)=="bmerBuild"){
    return(object)
  } else stop("Must supply a bmerFit or bmerBuild")
}

#' Extract posterior samples with named dimensions.
#'
#' The function calls \code{\link[rstan]{extract}} and then gives names to the
#' dimensions of the posterior samples.
#'
#' As an example, if there were 4 chains
#' and 1000 post-warmup iterations in each chain, then extracting the random effect
#' estimates for individual members of a random effects group (the \code{gamma}
#' parameter for the group) with \code{M} members and \code{Q} effects will return an array
#' which is \code{4000 x M x Q},
#' where the names of the second dimension are the group member names and the names
#' of the third are the random effect names.
#'
#' @param fit A \linkS4class{bmerFit}.
#' @param pars A character vector of parameter names. See \code{\link{par_names}}.
#' @return A named list where each component contains the posterior samples for a parameter,
#' and each array of posterior samples has named dimensions.
#' @export
named_extract <- function(fit, pars = names(par_names(fit))){
  if(class(fit)!="bmerFit") stop("must supply a bmerFit")
  pars <- extract(object = fit, pars = pars)
  for(p in names(pars)){
    d <- fit@build@pars[[p]]$dims
    if(length(d)>1 | length(d[[1]])>1){
      for(i in 1:length(d)){
        dimnames(pars[[p]])[[i+1]] <- d[[i]]
      }
    }
  }
  return(pars)
}

#' Get factor contrasts.
#'
#' \code{get_contrasts} lists the contrast matrix for each factor in the regression.
#'
#' If there are no factors in the regression, then an empty list is returned.
#' Otherwise, the list will have an element \code{main_effects} which is itself
#' a named list with an entry for the contrast matrix for each factor in the regression.
#' If \code{fixef_na = TRUE} in \code{\link{bmer_control}} and there were interactions
#' involving unordered factors with NA values which required a change in contrasts for
#' the interaction term, then in addition to \code{main_effects} there will be an element
#' \code{na_interactions} with an entry for each interaction where contrasts were changed,
#' with these entries having the same structure as \code{main_effects}.  Similarly,
#' if \code{ranef_na = TRUE} in \code{\link{bmer_control}} and there are fixed effects
#' with levels which only occur when the random effect group is NA, then these factors
#' will have entries structured the same as \code{main_effects} in a sublist named
#' with the random effect group followed by \code{_altered_main_effects}.  Interactions
#' in the slopes for a random effect group which require contrast alterations relative to the way
#' the main effect random slopes are coded are handled in the same way as fixed effect alterations
#' relative the fixed main effects are handled, and stored in a sublist named with the random
#' effect group followed by \code{_na_interactions}.  In the most common use of \code{bmer},
#' where both \code{fixef_na = FALSE} and \code{ranef_na = FALSE}, the result will only contain
#' the \code{main_effects} sublist.  For examples of changes due to NAs, see the vignette.
#'
#' @param object A \linkS4class{bmerBuild} or \linkS4class{bmerFit}.
#' @return A list of contrast matrices organized in named sublists.
#' @export
get_contrasts <- function(object){
  build <- get_build(object)
  v <- build@variables
  contr <- list()
  if(length(v$fact)>0) contr$main_effects <- list()
  if(length(v$nainter)>0) contr$na_interactions <- list()
  for(i in names(v$fact)){
    contr$main_effects[[i]] <- v$fact[[i]]$contr
  }
  for(i in names(v$nainter)){
    contr$na_interactions[[i]] <- v$nainter[[i]]
  }
  for(r in names(v$random)){
    if(length(v$random[[r]]$newfacts)>0){
      n <- paste(r,"altered_main_effects",sep="_")
      contr[[n]] <- list()
      for(i in names(v$random[[r]]$newfacts)){
        contr[[n]][[i]] <- v$random[[r]]$newfacts[[i]]$contr
      }
    }
    if(length(v$random[[r]]$nainter$changed)>0){
      n <- paste(r,"na_interactions",sep="_")
      contr[[n]] <- list()
      for(i in names(v$random[[r]]$nainter$changed)){
        contr[[n]][[i]] <- v$random[[r]]$nainter$changed[[i]]$contr
      }
    }
  }
  return(contr)
}

#' Print or save Stan code.
#'
#' The \code{cat_code} function extracts the Stan code used in a \code{bmer} model and either
#' prints it within R or writes it to a file depending on the value of \code{file}.
#'
#' @param object A \linkS4class{bmerBuild} or \linkS4class{bmerFit}.
#' @param file Name of a file to create and write the code to.  If NULL (the default) then the code is printed within R.
#' @param commented Logical.  If FALSE, then the code used within the \linkS4class{stanmodel} is used.  If TRUE (the recommended default),
#' then a version of the code with additional comments specific to the current model is used (but not different in execution from the raw code).
#' @export
cat_code <- function(object, file = NULL, commented = TRUE){
  build <- get_build(object)
  code <- ifelse(commented,build@code,build@stancode)
  if(!is.null(file)){
    cat(code,file=file,append=FALSE)
  } else cat(code)
}

#' Get parameter names for a \code{bmer} model.
#'
#' @param object A \linkS4class{bmerBuild} or \linkS4class{bmerFit}.
#' @return A named list with an entry for each parameter in the model, where the entries have
#' the names of each dimension in the parameter.
#' @export
par_names <- function(object){
  build <- get_build(object)
  pars <- lapply(build@pars,function(x) x$dims)
  return(pars)
}

#' Retrieve estimates of all effects.
#'
#' The \code{all_comparisons} function computes the posterior distribution of each regression term's effect
#' in the model (i.e. both main effects and interaction terms with the result depending on the term detials).
#' This function can only be used if \code{set_contr = TRUE} and \code{scale_cont = TRUE} in \code{\link{bmer_control}}.
#'
#' In all cases, the estimates returned hold the value of all variables not involved in the term at their average values.
#' For main effects, if the term is a factor, then group means for each factor level are computed
#' and pairwise comparisons of each of these are also computed.  If the main effect term is continuous,
#' then the effect of a one-standard-deviation increase in the effect is computed (which is the same as the raw regression output).
#'
#' Interactions involving only factors are treated in the same way as main effects, with an estimate for each level of
#' the interaction term and pairwise comparisons for all interaction levels.  For interactions involving more than one
#' continuous predictor, the effect of a simultaneous one-standard-deviation increase in all of the continuous predictors
#' involved is computed (which is the same as the raw regression output).  For interactions involving both
#' factors and continuous predictors, the effect of a one-standard-deviation
#' increase in the continuous predictors for each factor level is computed, and pairwise comparisons for group differences
#' in the effect of this one-standard-deviation increase are also computed.
#'
#' For unordered factors involving NAs, the estimate for observations set to NA is also computed.  This may or may not have
#' any real meaning depending on the interpretation of the NA value, and so should be interpreted according to the use of
#' NA in the data.
#'
#' @param fit A \linkS4class{bmerFit}.
#' @param probs A vector of quantiles to compute for all group means and pairwise comparisons.
#' @return A named list with an entry for each term in the regression.  Each element is itself a list with two components.
#' The first, \code{estimates}, contains group means or one-standard-deviation increase effects (see 'Details'). The second,
#' \code{contrasts} contains pairwise comparisons of the distributions in \code{estimates} when applicable.
#' @export
all_comparisons <- function(fit, probs = c(.025,.975)){
  if(class(fit)!="bmerFit") stop("Must supply a bmerFit")
  if(!fit@build@control@set_contr | !fit@build@control@scale_cont){
    stop("all_comparisons can only be used when set_contr = TRUE and scale_cont = TRUE")
  }
  comps <- list()
  g <- fit@build@variables$refgrids
  m <- fit@build@variables$refmats
  bt <- t(extract(fit,pars="beta")$beta)
  if(length(probs)<1) stop("must supply quantiles in 'probs'")
  
  for(n in names(g)){
    comps[[n]] <- list()
    
    gn <- g[[n]]
    nms <- gn
    for(j in 1:ncol(nms)){
      nms[,j] <- paste(colnames(nms)[j],nms[,j],sep="")
    }
    nms <- apply(nms,1,function(x) paste(x,collapse=":"))
    nc <- ncol(gn)
    np <- length(probs)
    gn[,(nc+1):(nc+np+3)] <- 0
    colnames(gn)[(nc+1):ncol(gn)] <- c("Mean","SD","P(>0)",paste(100*probs,"%",sep=""))
    yhat <- m[[n]] %*% bt
    gn[,nc+1] <- apply(yhat,1,mean)
    gn[,nc+2] <- apply(yhat,1,sd)
    gn[,nc+3] <- apply(yhat,1,function(n) mean(n>0))
    gn[,(nc+4):ncol(gn)] <- t(apply(yhat,1,function(x) quantile(x,probs=probs)))
    comps[[n]]$estimates <- gn
    
    if(nrow(gn)>1){
      ncomps <- choose(nrow(gn),2)
      combns <- combn(nrow(gn),2)
      contrhat <- t(apply(combns,2,function(x) yhat[x[1],]-yhat[x[2],]))
      contr <- data.frame(matrix(0,ncomps,np+4))
      contr[,1] <- apply(combns,2,function(x) paste(nms[x],collapse=" - "))
      contr[,2] <- apply(contrhat,1,mean)
      contr[,3] <- apply(contrhat,1,sd)
      contr[,4] <- apply(contrhat,1,function(x) mean(x>0))
      contr[,5:ncol(contr)] <- t(apply(contrhat,1,function(x) quantile(x,probs=probs)))
      colnames(contr) <- c("Contrast","Mean","SD","P(>0)",paste(100*probs,"%",sep=""))
    } else {
      contr <- NULL
    }
    
    comps[[n]]$contrasts <- contr
  }
  
  return(comps)
}

#' Summarize a \linkS4class{bmerBuild}.
#'
#' Returns a list which contains summary information about a \linkS4class{bmerBuild}. While the individual
#' elements of the list can be viewed and used, it is most useful for its default show method,
#' which lays out the important aspects of the model which will be fit when the \linkS4class{bmerBuild}
#' is passed to \code{\link{fit_bmer_build}}.
#'
#' The elements of the returned list are:
#'
#' \strong{call}.  The function call evaluated to create the \linkS4class{bmerBuild}.
#'
#' \strong{name}.  The name given to the model.
#'
#' \strong{family}.  The regression family.
#'
#' \strong{frame}.  The model frame generated from evaluating the formula.  It has attributes "fixed_varmat"
#'     and "random_varmat" indicating how the variables relate to interaction terms.
#'
#' \strong{N}.  The number of observations.
#'
#' \strong{PQ}.  The total number of unconstrained parameters (both fixed and random).
#'
#' \strong{response}.  A list giving information about the dependent variable: name; class; if binomial,
#' the levels of the response factor and which was assigned to 0 and to 1; if Gaussian, the raw mean
#' and standard deviation and a logical indicating whether or not the response is on unit scale in the regression.
#'
#' \strong{fixed}.  A list giving information about the fixed effects, with three elements: \emph{factors}, a named list
#' with an entry for each factor in the regression containing its number of levels, contrasts, etc.; \emph{covariates}, a
#' named list with an entry for each continuous predictor in the regression containing its mean and standard deviation in the
#' raw data and a logical indicating whether or not it was scaled in the regression; and \emph{terms}, a character vector of
#' the main effects and interactions in the regression from the expanded fixed effects formula.  If interactions involving
#' fixed effects with NAs resulted in changing factor contrasts, then additionally there will be a fourth element \emph{nainter},
#' which is a named list with an element for each interaction term for which contrasts were altered, giving the contrasts
#' used in the computation of the interaction term.  See the vignette for more details.
#'
#' \strong{random}.  A named list with an element for each random effects grouping factor, each of which is a list with the
#' number of levels (i.e. group members), their names, whether or not there are NA values, a list of new factor contrasts
#' created (if any), a list of slope interactions where factor contrasts were changed (if any) or dropped due to redundancy
#' (if any), the terms which were evaluated as potential slopes, the terms which were included as slopes, and the names
#' of the random effects model matrix columns.  See the vignette for more details.
#'
#' \strong{code}.  The Stan code used to fit the model, along with model-specific comments.  See \code{\link{cat_code}}.
#'
#' \strong{data}.  The list passed as the \code{data} argument to \code{\link[rstan]{sampling}}.
#'
#' \strong{control}.  The \code{\link{bmer_control}} used to build the model.
#'
#' \strong{messages}.  The messages (not warnings) returned from \code{\link{build_bmer_model}} (if any).
#'
#' \strong{warnings}.  The warnings returned from \code{\link{build_bmer_model}} (if any).
#'
#' @param object A \linkS4class{bmerBuild} or \linkS4class{bmerFit}.
#' @return A list of class \linkS4class{bmerBuildSummary} whose elements are given in 'Details'.
#' @seealso \code{\link{build_summary}},  \code{\link{get_contrasts}},  \code{\link{cat_code}},  \code{\link{par_names}},
#' \code{\link{get_control}}.
#'
#' @export
build_summary <- function(object){
  build <- get_build(object)
	summ <- bmerBuildSummary()
	
	summ$call <- build@call
	summ$name <- build@name
  summ$family <- build@family
  summ$frame <- build@frame
	summ$N <- build@data$N
	summ$PQ <- attr(build@pars,"pq")
	summ$response <- build@variables$resp
	summ$fixed <- list(factors=build@variables$fact,covariates=build@variables$cont,
    terms=colnames(attr(build@frame,"fixed_varmat")))
  if(length(build@variables$nainter)>0) summ$fixed$nainter <- build@variables$nainter
	summ$random <- build@variables$random
	summ$code <- build@code
	summ$data <- build@data
  summ$control <- build@control
	summ$messages <- build@messages
	summ$warnings <- build@warnings
	
	return(summ)
}

#' Summarize a \linkS4class{bmerFit}.
#'
#' Returns a list which contains summary information about a bmerFit. The default printing option
#' gives a summary of the random effects standard deviations and correlations (posterior means), as well
#' as the fixed effects (posterior means, standard deviations, 95\% credible intervals, effect sample sizes,
#' Gelman-Rubin Rhat statistics, and posterior probability that each is positive)
#'
#' The elements of the returned list are:
#'
#' \strong{call}.  The function call evaluated to create the \linkS4class{bmerFit}.
#'
#' \strong{build}.  A list returned by \code{\link{build_summary}} summarizing the \linkS4class{bmerBuild} used to fit the \linkS4class{bmerFit}.
#'
#' \strong{random}.  A named list with an element for each random effect grouping factor, each being a list
#' with the following elements: \emph{ranef}, the intercept and slope estimates for individual member of the grouping factor;
#' \emph{coef}, the coefficient vector for each group member obtained from summing the random intercept and slope estimates
#' with the fixed effects coefficients; \emph{sd}, the standard deviation in the random effects; \emph{cor} the correlation
#' matrix for the random effects; and \emph{cov} the covariance matrix for the random effects.
#'
#' \strong{fixed}.  A list with the following elements: \emph{coef}, the fixed effect estiamtes, including the posterior
#' probability that each is greater than zero (see the vignette for more details); \emph{cor} the correlation in
#' the fixed effects; and \emph{cov}, the covariance matrix for the fixed effects.
#'
#' \strong{fitted.values}.  A vector containing the mean estimate for each observation (predicted response value for Gaussian
#' regression and predicted log-odds of being in category 1 for binomial regression).
#'
#' \strong{warnings}.  The warnings returned by \code{\link{fit_bmer_build}} (if any).
#'
#' @param fit A \linkS4class{bmerFit}.
#' @return A list of class \linkS4class{bmerFitSummary} whose elements are given in 'Details'.
#' @seealso \code{\link{fit_summary}},  \code{\link{get_build}},  \code{\link{named_extract}}, 
#' \code{\link{get_contrasts}},  \code{\link{cat_code}}, \code{\link{par_names}},  \code{\link{all_comparisons}},
#' \code{\link{bmers_summary}},  \code{\link{get_control}}.
#'
#' @export
fit_summary <- function(fit){
  if(class(fit)!="bmerFit") stop("Must supply a bmerFit")
  
  summ <- bmerFitSummary()
  summ$call <- fit@call
  summ$build <- build_summary(fit@build)
  summ$random <- list()
  summ$fixed <- list(coef = rstan::summary(fit,pars="beta",probs=c(.025,.975))$summary)
  summ$fitted.values <- as.numeric(rstan::summary(fit,pars="y_hat")$summary[,1])
  summ$warnings <- fit@warnings
  
  P <- length((Pnames <- fit@build@pars$beta$dims$Effect))
  b <- extract(fit, pars = "beta")[[1]]
  summ$fixed$coef <- cbind(summ$fixed$coef,apply(b,2,function(x) mean(x>0)))
  colnames(summ$fixed$coef)[8] <- "P(>0)"
  rownames(summ$fixed$coef) <- Pnames
  summ$fixed$cor <- cor(b)
  summ$fixed$cov <- cov(b)
  dimnames(summ$fixed$cor) <- dimnames(summ$fixed$cov) <- list(Pnames,Pnames)

  for(r in names(fit@build@variables$random)){
    g <- paste("gamma",r,sep="_"); s <- paste("sigma",r,sep="_"); o <- paste("omega",r,sep="_")
    M <- length((Mnames <- fit@build@pars[[g]]$dims$Member))
    Q <- length((Qnames <- fit@build@pars[[g]]$dims$Effect))
    
    summ$random[[r]]$ranef <- matrix(rstan::summary(fit,pars=g)$summary[,1],M,Q,byrow=TRUE)
    dimnames(summ$random[[r]]$ranef) <- list(Mnames,Qnames)
    summ$random[[r]]$coef <- matrix(colMeans(b),M,P,byrow=TRUE)
    dimnames(summ$random[[r]]$coef) <- list(Mnames,Pnames)
    for(i in Qnames[!(Qnames %in% Pnames)]){
      summ$random[[r]]$coef <- cbind(summ$random[[r]]$coef,0)
      colnames(summ$random[[r]]$coef)[ncol(summ$random[[r]]$coef)] <- i
    }
    for(i in Qnames) summ$random[[r]]$coef[,i] <- summ$random[[r]]$coef[,i] + summ$random[[r]]$ranef[,i]
    rsd <- summ$random[[r]]$sd <- rstan::summary(fit,pars=s)$summary[,1]
    names(summ$random[[r]]$sd) <- Qnames
    if(Q == 1){
      summ$random[[r]]$cor <- NULL
      summ$random[[r]]$cov <- NULL
    } else {
      rcor <- extract(fit, pars = o)[[1]]
      rcor <- apply(rcor,2:3,mean)
      summ$random[[r]]$cor <- rcor
      summ$random[[r]]$cov <- diag(rsd) %*% rcor %*% t(diag(rsd))
      dimnames(summ$random[[r]]$cor) <- dimnames(summ$random[[r]]$cov) <- list(Qnames,Qnames)
    }
  }
  
  if(fit@build@family == "gaussian") summ$sigma_residual <- rstan::summary(fit, pars = "sigma_residual")$summary[,1]
  
  return(summ)
}


#' Get parameter summary from a \linkS4class{bmerFit}.
#'
#' Returns a \linkS4class{stanfit} summary along with the posterior probability that a parameter is on a given interval.
#'
#' @param fit A \linkS4class{bmerFit}.
#' @param pars Names of the parameters to be summarized (defaults to all parameters)
#' @param probs Quantiles to return for each parameter (defaults to a 95\% equal-tailed credible interval)
#' @param interval The interval on which the posterior probability of each parameter is to be computed.  The
#' default is (0, Inf) which returns the posterior probability that each parameter is positive.
#' @return A matrix summarizing the indicated parameters.
#' @export
bmers_summary <- function(fit, pars = names(par_names(fit)), probs = c(.025, .975), interval = c(0, Inf)){
  if(!is.numeric(interval)) stop("'interval' must be numeric")
  if(length(interval)!=2) stop("'interval' must contain an interval where the first element is less than the second")
  if(interval[1] >= interval[2]) stop("'interval' must contain an interval where the first element is less than the second")
  
  summ <- rstan::summary(fit, pars = pars, probs = probs)$summary
  samps <- named_extract(fit, pars = pars)
  
  p <- list()
  for(i in 1:length(samps)){
    d <- length(dim(samps[[i]]))
    if(d){
      p[[i]] <- apply(samps[[i]], 2:d, function(x) mean(x > interval[1] & x < interval[2]))
      if(is.matrix(p[[i]])) p[[i]] <- as.vector(t(p[[i]]))
    } else {
      p[[i]] <- mean(samps[[i]] > interval[1] & samps[[i]] < interval[2])
    }
  }
  p <- do.call(c, args = p)
  
  summ <- cbind(summ, p)
  if(interval[1] == -Inf){
    cols <- paste("P(<",interval[2],")",sep="")
  } else if(interval[2] == Inf){
    cols <- paste("P(>",interval[1],")",sep="")
  } else {
    cols <- paste("P(on (",paste(interval,collapse=","),"))",sep="")
  }
  colnames(summ)[ncol(summ)] <- cols
  
  keep <- unlist(lapply(fit@build@pars[pars], function(x) x$keep))
  cols <- colnames(summ)
  summ <- summ[keep,]
  if(!is.matrix(summ)){
    summ <- matrix(summ, 1, length(summ))
    colnames(summ) <- cols
  }
  
  return(summ)
}


#' See the control parameters applied to a \linkS4class{bmerBuild} or \linkS4class{bmerFit}.
#'
#' @param object A \linkS4class{bmerBuild} or \linkS4class{bmerFit}.
#' @export
get_control <- function(object){
  build <- get_build(object)
  return(build@control)
}
