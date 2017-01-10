#' An S4 class to store control parameters for \code{\link{bmer}}
#'
#' Stores hyperparameters for the priors, determines the amount of function printing,
#' and gives instructions for how variable scaling, factor contrasts, and NA values are to be handled.
#'
#' It is very important to note that the default weakly-informative priors are set assuming that \code{scale_cont} and
#' \code{set_contr} are both set to TRUE. These parameters should only be set to FALSE if you have a lot of experience
#' with Bayesian modeling and can properly adjust the priors.  See the vignette for more details.
#'
#' @slot verbose Logical. If TRUE, summaries of the \linkS4class{bmerBuild} and \linkS4class{bmerFit} are printed when
#' \code{\link{bmer}} is called. If FALSE (the default), only warning messages (if any) are printed.
#' @slot scale_cont Logical. If TRUE (the default), continuous variables are scaled (i.e. mean-centered and divided by their
#' standard deviation) so that they have mean 0 and standard deviation 1 prior to regression. If FALSE (highly not recommended),
#' continuous variables are left on their original scales.
#' @slot set_contr Logical. If TRUE (the default), sum contrasts are set for unordered factors and scaled orthogonal polynomial
#' contrasts are set for ordered factors prior to regression. If FALSE (highly not recommended), factor contrasts are not checked.
#' @slot estimate_scale_beta Character. One of "P", "yes", or "no". If "yes", the scale of the fixed effect t-distribution prior
#' is estimated with a hyper-prior which is assigned a half-normal distribution with scale \code{sc_beta}. If "no", the prior scale
#' is not estimated with a hyper-prior, and is instead given a value of \code{scale_beta}.  If "P" (the default), then
#' \code{estimate_scale_beta} is interpreted as "yes" if there are 10 or more fixed effect coefficients, and as "no" if there are
#' fewer than 10 fixed effect coefficients.
#' @slot sc_beta Numeric (default 1; must be greater than 0). The scale for the half-normal hyper-prior on the scale of the
#' fixed effect t-distribution prior if a hyper-prior is assigned (depends on \code{estimate_scale_beta}).
#' @slot scale_beta Numeric (default 2; must be greater than 0). The scale for the fixed effect t-distribution prior if a
#' hyper-prior is not assigned (depends on \code{estimate_scale_beta}).
#' @slot nu_beta Numeric (default 5; must be greater than 0). The degrees of freedom for the fixed effect t-distribution prior.
#' Setting to 1 results in a Cauchy prior.
#' @slot sc_q0 Numeric (default is 1.5; must be greater than 0). The scale for the half-normal prior on the standard deviation in
#' the random intercepts.
#' @slot sc_qs Numeric (default is 1; must be greater than 0). The scale for the half-normal prior on the standard deviation in
#' the random slopes.
#' @slot eta_q Numeric (default is 2; must be greater than 0). The shape parameter for the LKJ prior on the random effects
#' correlation matrices. Setting to 1 results in a flat prior.
#' @slot sc_res Numeric (default is 0.5; must be greater than 0). The scale for the half-normal prior on the standard deviation in
#' the residuals.  Only applies to linear (i.e. Gaussian) regressions.
#' @slot fixef_na Logical. If FALSE (the default), observations with NA values in the fixed effects are dropped. If TRUE, then
#' \code{set_contr} must also be TRUE, and NA values will be set to 0 in the model matrix (i.e. level slashing). When an unordered
#' factor with NA values is involved in interactions, the interaction term is checked for higher-order NA values and, if necessary,
#' the contrasts for the factors with NA values are changed for the interaction term, and these new contrasts which apply to the
#' interaction term but not the main effects are logged. In regression output, interaction terms with contrast changes are prefixed
#' with "nai_" so that they can be easily identified. See the vignette for examples.
#' @slot ranef_na Logical. If FALSE (the default), observations with NA values in the random effect grouping factors are dropped.
#' If TRUE, then set_contr must also be TRUE, and NA values will be set to 0 in the model matrix. When a within-group fixed effect
#' factor only takes on a subset of its values when the random grouping factor is not NA, then the contrasts for the factor are
#' changed in the coding of slopes for the random effect grouping factor, and the new contrasts which apply to the slope but not the
#' fixed effect are logged. In the regression output, these slopes are prefixed with the name of the random effect grouping factor
#' so they can be easily identified.
#' @export
bmerControl <- setClass("bmerControl",
  slots = list(
    verbose = "logical",
    scale_cont = "logical",
    set_contr = "logical",
    estimate_scale_beta = "character",
    sc_beta = "numeric",
    scale_beta = "numeric",
    nu_beta = "numeric",
    sc_q0 = "numeric",
    sc_qs = "numeric",
    eta_q = "numeric",
    sc_res = "numeric",
    fixef_na = "logical",
    ranef_na = "logical"),

  prototype = list(
    verbose = FALSE,
    scale_cont = TRUE,
    set_contr = TRUE,
    estimate_scale_beta = "P",
    sc_beta = 1,
    scale_beta = 2,
    nu_beta = 5,
    sc_q0 = 1.5,
    sc_qs = 1,
    eta_q = 2,
    sc_res = 0.5,
    fixef_na = FALSE,
    ranef_na = FALSE),

    validity = function(object){
      messages <- NULL

      for(i in slotNames(object)){
        if(length(slot(object,i))!=1){
          messages <- c(messages,paste(i,"must be of length 1"))
        }
        
        if(i %in% c("sc_beta","scale_beta","nu_beta","sc_q0","sc_qs","eta_q","sc_res")){
          if(slot(object,i)<=0){
            messages <- c(messages,paste(i,"must be positive"))
          }
        }
      }
      
      if(!(object@estimate_scale_beta %in% c("P","yes","no"))){
        messages <- c(messages,"estimate_scale_beta must be one of 'P', 'yes', or 'no'")
      }
      
      if(object@fixef_na & !object@set_contr){
        messages <- c(messages,"If fixef_na is TRUE, then set_contr must also be TRUE")
      }
      
      if(length(messages)>0) return(messages)
      return(TRUE)
    })

#' Control parameters for \code{\link{bmer}}.
#'
#' Provides \code{\link{build_bmer_model}} with hyperparameters for the priors, determines the amount of function printing,
#' and gives instructions for how variable scaling, factor contrasts, and NA values are to be handled.
#'
#' It is very important to note that the default weakly-informative priors are set assuming that \code{scale_cont} and
#' \code{set_contr} are both set to TRUE. These parameters should only be set to FALSE if you have a lot of experience
#' with Bayesian modeling and can properly adjust the priors.
#'
#' @param verbose Logical. If TRUE, summaries of the \linkS4class{bmerBuild} and \linkS4class{bmerFit} are printed when
#' \code{\link{bmer}} is called. If FALSE (the default), only warning messages (if any) are printed.
#' @param scale_cont Logical. If TRUE (the default), continuous variables are scaled (i.e. mean-centered and divided by their
#' standard deviation) so that they have mean 0 and standard deviation 1 prior to regression. If FALSE (highly not recommended),
#' continuous variables are left on their original scales.
#' @param set_contr Logical. If TRUE (the default), sum contrasts are set for unordered factors and scaled orthogonal polynomial
#' contrasts are set for ordered factors prior to regression. If FALSE (highly not recommended), factor contrasts are not checked.
#' @param estimate_scale_beta Character. One of "P", "yes", or "no". If "yes", the scale of the fixed effect t-distribution prior
#' is estimated with a hyper-prior which is assigned a half-normal distribution with scale \code{sc_beta}. If "no", the prior scale
#' is not estimated with a hyper-prior, and is instead given a value of \code{scale_beta}.  If "P" (the default), then
#' \code{estimate_scale_beta} is interpreted as "yes" if there are 10 or more fixed effect coefficients, and as "no" if there are
#' fewer than 10 fixed effect coefficients.
#' @param sc_beta Numeric (default 1; must be greater than 0). The scale for the half-normal hyper-prior on the scale of the
#' fixed effect t-distribution prior if a hyper-prior is assigned (depends on \code{estimate_scale_beta}).
#' @param scale_beta Numeric (default 2; must be greater than 0). The scale for the fixed effect t-distribution prior if a
#' hyper-prior is not assigned (depends on \code{estimate_scale_beta}).
#' @param nu_beta Numeric (default 5; must be greater than 0). The degrees of freedom for the fixed effect t-distribution prior.
#' Setting to 1 results in a Cauchy prior.
#' @param sc_q0 Numeric (default is 1.5; must be greater than 0). The scale for the half-normal prior on the standard deviation in
#' the random intercepts.
#' @param sc_qs Numeric (default is 1; must be greater than 0). The scale for the half-normal prior on the standard deviation in
#' the random slopes.
#' @param eta_q Numeric (default is 2; must be greater than 0). The shape parameter for the LKJ prior on the random effects
#' correlation matrices. Setting to 1 results in a flat prior.
#' @param sc_res Numeric (default is 0.5; must be greater than 0). The scale for the half-normal prior on the standard deviation in
#' the residuals.  Only applies to linear (i.e. Gaussian) regressions.
#' @param fixef_na Logical. If FALSE (the default), observations with NA values in the fixed effects are dropped. If TRUE, then
#' \code{set_contr} must also be TRUE, and NA values will be set to 0 in the model matrix (i.e. level slashing). When an unordered
#' factor with NA values is involved in interactions, the interaction term is checked for higher-order NA values and, if necessary,
#' the contrasts for the factors with NA values are changed for the interaction term, and these new contrasts which apply to the
#' interaction term but not the main effects are logged. In regression output, interaction terms with contrast changes are prefixed
#' with "nai_" so that they can be easily identified. See the vignette for examples.
#' @param ranef_na Logical. If FALSE (the default), observations with NA values in the random effect grouping factors are dropped.
#' If TRUE, then set_contr must also be TRUE, and NA values will be set to 0 in the model matrix. When a within-group fixed effect
#' factor only takes on a subset of its values when the random grouping factor is not NA, then the contrasts for the factor are
#' changed in the coding of slopes for the random effect grouping factor, and the new contrasts which apply to the slope but not the
#' fixed effect are logged. In the regression output, these slopes are prefixed with the name of the random effect grouping factor
#' so they can be easily identified.
#' @return An object of class \linkS4class{bmerControl} which provides control parameters to \code{\link{bmer}}.
#' @export
bmer_control <- function(verbose = FALSE, scale_cont = TRUE, set_contr = TRUE, estimate_scale_beta = "P",
sc_beta = 1, scale_beta = 2, nu_beta = 5, sc_q0 = 1.5, sc_qs = 1, eta_q = 2, sc_res = 0.5, fixef_na = FALSE,
ranef_na = FALSE){
  return(bmerControl(verbose = verbose, scale_cont = scale_cont, set_contr = set_contr,
    estimate_scale_beta = estimate_scale_beta, sc_beta = sc_beta, scale_beta = scale_beta,
    nu_beta = nu_beta, sc_q0 = sc_q0, sc_qs = sc_qs, eta_q = eta_q, sc_res = sc_res,
    fixef_na = fixef_na, ranef_na = ranef_na))
}


#' An S4 class which contains all the information required to fit a \linkS4class{bmerFit}.
#'
#' A \code{bmerBuild} can be summarized using \code{\link{build_summary}}, which is also the default printing method.
#'
#' @seealso \code{\link{build_summary}},  \code{\link{get_contrasts}},  \code{\link{cat_code}},  \code{\link{par_names}}
#'
#' @slot call The function call which created the bmerBuild.
#' @slot control The \linkS4class{bmerControl} used.
#' @slot family The regression family ("gaussian" or "binomial").
#' @slot frame The data.frame containing the fixed effect variables and random grouping factors.
#' @slot x The fixed effect model matrix.
#' @slot z A named list where each element is the random effects sparse model matrix for a single random grouping factor.
#' @slot name The model_name passed to the function call.
#' @slot code The commented Stan code.
#' @slot stancode The uncommented Stan code.
#' @slot stanmod A character string describing the uncommented Stan code. Used to name and identify previously compiled models for future use.
#' @slot data A list to be passed to the \code{data} argument of \code{\link[rstan]{sampling}}.
#' @slot variables A list containing information about the variables and random effects structure.
#' @slot pars A named list containing the names of the paramters and their dimension names.
#' @slot messages A character vector of informational messages (not warnings) returned from \code{\link{build_bmer_model}}.
#' @slot warnings A character vector of warnings returned from \code{\link{build_bmer_model}}.
#' @export
bmerBuild <- setClass("bmerBuild",
	slots = list(
		call = "call",
    control = "bmerControl",
    family = "character",
    frame = "data.frame",
    x = "matrix",
    z = "list",
    name = "character",
    code = "character",
    stancode = "character",
    stanmod = "character",
    data = "list",
    variables = "list",
    pars = "list",
    messages = "character",
    warnings = "character"))

#' An S4 class to hold the summary of a \code{bmerBuild}.
#'
#' For details see \code{\link{build_summary}}.
#'
#' @export
bmerBuildSummary <- setClass("bmerBuildSummary",contains = "list")

#' An S4 class to store a \linkS4class{bmerBuild} and \linkS4class{stanfit}.
#'
#' A \code{bmerFit} can be summarized using \code{\link{fit_summary}}, which is also the default printing method.
#'
#' @seealso \code{\link{fit_summary}},  \code{\link{get_build}},  \code{\link{named_extract}}, 
#' \code{\link{all_comparisons}}, \code{\link{bmers_summary}}
#'
#' @slot call The function call which produced the object.
#' @slot build The \linkS4class{bmerBuild} passed to the \code{\link{fit_bmer_build}} function call.
#' @slot diagnostics A list of diagnostics from the sampler and for Rhat and n_eff.
#' @slot warnings A character vector of warnings returned by \code{\link{fit_bmer_build}} (if any).
#' @export
bmerFit <- setClass("bmerFit",
	slots = list(
    call = "call",
		build = "bmerBuild",
		diagnostics = "list",
		warnings = "character"),
  contains = "stanfit")

#' An S4 class to hold the summary of a \code{bmerFit}.
#'
#' For details see \code{\link{fit_summary}}.
#'
#' @export
bmerFitSummary <- setClass("bmerFitSummary",contains = "list")
