#' Bayesian mixed effects regression.
#'
#' Fit a Bayesian mixed effects regression with a maximal random effects structure
#' and weakly informative priors in Stan through the rstan interface using NUTS.
#' Additional details not covered here can be found in the vignette.  To get to the
#' vignette, go to the package Index at the bottom of the page and click on
#' "User guides, package vignettes and other documentation".
#'
#' The \code{formula} should have the form \code{response ~ fixed | random_groups}
#' without slopes and intercepts specified. Only the names of the variables which
#' are to be considered random grouping factors should be given to the right of the
#' the pipe '|', separated by '+'. Nesting structure is expressed as in any other
#' case. For example, \code{y ~ x1 + x2 | subject + A/B} would indicate that there
#' are three random effect groups: \code{subject}, \code{'A'}, and \code{'B'} nested
#' within \code{'A'} (that is, \code{A:B}). Random intercepts will be fit for each of
#' these groups.  For each of the fixed effects \code{x1} and \code{x2}, random slopes
#' will be fit for a group if the effect varies within-group (i.e. individual group
#' members have more than one value for the effect).
#'
#' The arguments to \code{bmer} are first passed to \code{build_bmer_model}.  All continuous
#' variables (fixed-effect covariates and, in Gaussian regressions, the response) are scaled
#' (i.e. mean-centered and divided by their standard deviations) so that they have
#' mean 0 and standard deviation 1.  Sum contrasts are set for unordered factors
#' rather than treatment contrasts, and scaled polynomial contrasts (orthogonal 
#' polynomial contrasts where the standard deviation of each column in the contrast
#' matrix is 1) are set for ordered factors rather than polynomial contrasts. For
#' details on the handling of NA values, see the vignette.  The maximal random effect
#' structure for each random effect group is then determined, and the model matrix is
#' generated; see \code{\link{bmer_control}} for details on the priors.
#' Stan code and a list of data to be passed to \code{\link[rstan]{sampling}} are then assembled
#' and stored in an object of class \linkS4class{bmerBuild}.
#'
#' The \linkS4class{bmerBuild} is then passed to \code{fit_bmer_build}, which calls \code{\link[rstan]{sampling}}, first
#' creating a \linkS4class{stanmodel} if required (see the vignette).
#'
#' @param formula A formula describing the model to be fit.  See 'Details'.
#' @param data A data.frame containing variables for the model.
#' @param family Regression family. Currently, gaussian and binomial are supported.
#' @param control A call to \code{\link{bmer_control}}.
#' @param model_name Character. A name given to the model.
#' @param adapt_delta Passed to the control argument for \code{\link[rstan]{sampling}}. See the vignette.
#' @param max_treedepth Passed to the control argument for \code{\link[rstan]{sampling}}. See the vignette.
#' @param ... Additional arguments to be passed to \code{\link[rstan]{sampling}}. See the vignette.
#' @param build An object of class \linkS4class{bmerBuild}.
#' @return The \code{bmer} function begins by calling \code{build_bmer_model}, a modular function which returns a
#' \linkS4class{bmerBuild}, and then passes this onto \code{fit_bmer_build}, which returns a \linkS4class{bmerFit}.
#' If you are new to \code{bmers}, you may want to use the two modular functions manually at first to ensure that
#' \code{build_bmer_model} is interpreting your variables the way you intended before passing
#' the \linkS4class{bmerBuild} to \code{fit_bmer_build}.
#' @export
bmer <- function(formula, data, family = gaussian, control = bmer_control(),
model_name = "anon_model", adapt_delta = 0.8, max_treedepth = 10, ...){

  mc <- match.call()
  build <- build_bmer_model(formula, data, family, control, model_name)
  build@call <- mc
  fit <- fit_bmer_build(build, adapt_delta, max_treedepth, ...)
  fit@call <- mc
  return(fit)
}

#' @rdname bmer
#' @export
build_bmer_model <- function(formula, data, family = gaussian, control = bmer_control(), model_name = "anon_model"){

	mc <- match.call()
  
	if(class(control)!="bmerControl") stop("control must be a call to bmer_control()")
  
	cat("",paste("Building model '",model_name,"' ...",sep=""),"",sep="\n")

  na_action <- getOption("na.action")
	possfail <- tryCatch(checked <- check_variables(formula,data,family,control), error = function(e) e)
  options(na.action = na_action)
  if(inherits(possfail,"error")){
    cat(possfail$message,sep="\n")
    stop("Error in check_variables")
  }
	dat <- bmer_stan_data(checked,control)
  # code to pass to stan
  stancode <- raw_stan_code(checked$stanmod)
  # same code but commented specifically for current build
	code <- commented_stan_code(checked$family,checked$varlog$random,dat$stan_data$P,model_name,control)
  possfail <- tryCatch(temp <- rstan::stanc(model_code = code), error = function(e) e)
  if(inherits(possfail,"error")){
    cat(code)
    cat(possfail$message,sep="\n")
    stop("Error in stan code")
  }
  possfail <- tryCatch(temp <- rstan::stanc(model_code = stancode), error = function(e) e)
  if(inherits(possfail,"error")){
    cat(stancode)
    cat(possfail$message,sep="\n")
    stop("Error in stan code")
  }
  attr(checked$model_frame,"fixed_varmat") <- checked$fixed_varmat
  attr(checked$model_frame,"random_varmat") <- checked$random_varmat
  attr(checked$x,"assign") <- NULL
  
	build <- bmerBuild(call = mc, control = control, family = checked$family, frame = checked$model_frame,
    x = checked$x, z = checked$z, name = model_name, code = code, stancode = stancode,
    data = dat$stan_data, variables = checked$varlog, stanmod = paste("bmers",bmers_version(),"model",checked$stanmod,sep="_"),
		pars = dat$pars, messages = c(checked$messages,dat$messages), warnings = c(checked$warnings,dat$warnings))

	if(control@verbose){
		show(build)
	} else if(length(build@warnings)>0){
		cat(build@warnings,"",sep="\n")
	}

	return(build)
}

#' @rdname bmer
#' @export
fit_bmer_build <- function(build, adapt_delta = 0.8, max_treedepth = 10, ...){

  mc <- match.call()

	if(class(build)!="bmerBuild") stop("must supply a bmerBuild")

	fit <- bmerFit(call = mc, build = build)
	
  possfail <- tryCatch(smdir <- find.package("bmers"), error = function(e) e)  
  
  if(inherits(possfail,"error")){
    cat("Compiling new stanmodel...")
    stanmod <- rstan::stan_model(model_name = build@stanmod, model_code = build@stancode)
  } else {
    currwd <- getwd()
    setwd(smdir)
    smfiles <- list.files()
    if(build@stanmod %in% smfiles){
      cat("Loading previously compiled stanmodel...")
      stanmod <- readRDS(build@stanmod)
    } else {
      cat("Compiling new stanmodel...")
      stanmod <- rstan::stan_model(model_name = build@stanmod, model_code = build@stancode)
      saveRDS(stanmod,file=build@stanmod)
    }
    setwd(currwd)
  }

	cat("",paste("Fitting model '",build@name,"' ...",sep=""),"",sep="\n")
  
	fit.stan <- rstan::sampling(object = stanmod, data = build@data, algorithm = "NUTS",
		pars = unlist(lapply(build@pars,function(x) x$name)), show_messages = FALSE,
		control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth), ...)
  
  p <- build@pars
  for(i in 1:length(p)){
    n2 <- names(p)[i]
    n1 <- p[[i]]$name
    fit.stan@model_pars[fit.stan@model_pars==n1] <- n2
    names(fit.stan@par_dims)[names(fit.stan@par_dims)==n1] <- n2
    fit.stan@sim$pars_oi[fit.stan@sim$pars_oi==n1] <- n2
    names(fit.stan@sim$dims_oi)[names(fit.stan@sim$dims_oi)==n1] <- n2
  }
  rnms <- unlist(lapply(p,function(x) x$rownames))
  names(rnms) <- NULL
  names(fit.stan)[1:length(rnms)] <- rnms
  
  for(i in slotNames(fit.stan)){
    slot(fit,i) <- slot(fit.stan,i)
  }
  
  keep <- c(unlist(lapply(p,function(x) x$keep)),FALSE)
  allpars <- rstan::summary(fit.stan)
  allpars <- allpars$summary[keep,]

	samp <- do.call(rbind,args=get_sampler_params(fit.stan,inc_warmup=FALSE))
	rhat <- allpars[,"Rhat"]
	neff <- allpars[,"n_eff"]
	fit@diagnostics <- list(sampler = samp, rhat = rhat, neff = neff)

	nmtd <- sum(samp[,"treedepth__"] > max_treedepth)
	ndiv <- sum(samp[,"divergent__"])
	neff <- sum(neff < 100)
	rhat <- sum(rhat >= 1.1)

	warn <- character()
	if(nmtd>0) warn <- c(warn,paste("WARNING:",nmtd,"post-warmup transitions exceeded the maximum treedepth; try increasing max_treedepth"))
	if(ndiv>0) warn <- c(warn,paste("WARNING:",ndiv,"divergent transitions post-warmup; try increasing adapt_delta"))
	if(neff>0) warn <- c(warn,paste("WARNING:",neff,"parameters have an effective samples size under 100; try increasing iter"))
	if(rhat>0) warn <- c(warn,paste("WARNING:",rhat,"parameters have an R-hat value over 1.01 (convergence not reached); try increasing warmup and iter"))
	fit@warnings <- warn
	
	if(build@control@verbose){
		show(fit)
	} else if(length(fit@warnings)>0){
		cat("",fit@warnings,"",sep="\n")
	}
	
	return(fit)
}
