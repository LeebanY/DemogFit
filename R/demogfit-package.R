#' demogfit: Blockwise Composite-Likelihood Demographic Inference
#'
#' Fits closed-form Wilkinson-Herbots (2008, 2012) blockwise
#' composite-likelihood models - Strict Isolation (\code{\link{wh_si}}),
#' Isolation-with-Migration (\code{\link{wh_im}}),
#' Isolation-with-Initial-Migration (\code{\link{wh_iim}}), and Secondary
#' Contact (\code{\link{wh_sc}}) - to the blockwise pairwise-mutation-count
#' distribution derived from a gIMble blockwise site-frequency spectrum
#' (bSFS; see \code{\link{bsfs_to_s_distribution}}).
#'
#' Each model has its own fitting function (\code{\link{fit_si}},
#' \code{\link{fit_im}}, \code{\link{fit_iim}}, \code{\link{fit_sc}}) that
#' returns fitted parameters and the maximized log-likelihood, with
#' user-overridable parameter bounds (see \code{\link{default_bounds}}).
#'
#' \code{\link{fit_demography}} and \code{\link{fit_demography_all}} fit
#' multiple models to one or many species pairs and report which model is
#' best supported - or, at minimum, whether there is evidence for gene
#' flow (Strict Isolation vs. Isolation-with-Migration).
#'
#' @importFrom stats dbinom pgamma ppois optim rmultinom pchisq
#' @importFrom graphics barplot lines legend
#' @importFrom utils modifyList
#' @keywords internal
"_PACKAGE"
