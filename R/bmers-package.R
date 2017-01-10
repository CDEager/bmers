#' bmers: A package for fitting Bayesian mixed effects regressions fit with Stan
#'
#' The \emph{bmers} package fits mixed effects regressions with weakly informative priors and maximal random effects structures
#' in \emph{Stan} with the No U-Turn Sampler (through the \emph{rstan} interface).  The maximal random effects structure
#' is determined automatically, and only the fixed effects structure and the random effects grouping factors need to
#' be specified.  Continuous variables are automatically scaled, and factor contrasts are automatically set; neither
#' of these needs to be done beforehand.  NA values in both the fixed and random effects are supported through the
#' use of sum contrasts with NAs set to zero in the model matrix.  The default weakly informative priors and the scaling
#' and contrasts settings can all be altered through a control variable.
#'
#' @section Vignette:
#' A vignette demonstrating how to use \emph{bmers} can be found
#' by going to the package Index at the bottom of the page and clicking on
#' "User guides, package vignettes and other documentation".
#'
#' @section Fit a Bayesian mixed effects regression:
#' \code{\link{bmer}},  \code{\link{build_bmer_model}},  \code{\link{fit_bmer_build}},  \code{\link{bmer_control}}
#'
#' @section Query a \linkS4class{bmerFit}:
#' \code{\link{get_build}},  \code{\link{named_extract}},  \code{\link{get_contrasts}},  \code{\link{cat_code}},
#' \code{\link{par_names}},  \code{\link{all_comparisons}},  \code{\link{build_summary}},  \code{\link{fit_summary}},
#' \code{\link{bmers_summary}},  \code{\link{get_control}}.
#'
#' @section Query a \linkS4class{bmerBuild}:
#' \code{\link{get_contrasts}},  \code{\link{cat_code}}, 
#' \code{\link{par_names}},  \code{\link{build_summary}},   \code{\link{get_control}}.
#'
#' @docType package
#' @name bmers
#'
#' @import formula.tools
#' @import methods
#' @import stringr
#' @import rstan
#'
#' @importFrom stats as.formula contr.helmert contr.poly contr.sum contr.treatment contrasts
#' @importFrom stats contrasts<- cor cov gaussian model.frame model.matrix poly quantile sd terms xtabs
#' @importFrom Matrix KhatriRao
#' @importFrom utils combn
#'
#' @importMethodsFrom Matrix t
#' @importClassesFrom rstan stanfit stanmodel
#' @importMethodsFrom rstan summary
#'
NULL
