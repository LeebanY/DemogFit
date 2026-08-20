# =============================================================================
# Closed-form blockwise composite-likelihood functions.
#
# All four models below are ported term-for-term from Wilkinson-Herbots
# (2008, 2012), as implemented in the source Mathematica notebooks
# (WilkHeSim2s, WilkHeSimIIMs, pulseADmix). SI is deliberately NOT a
# separate formula: it is the M = 0 special case of wh_im(), exactly as in
# the original notebook.
#
# Two models present in the original notebook - MIG (a stationary,
# no-split migration model) and IIML (the infinite-ancestral-time limit of
# IIM) - have been deliberately removed from this package, since they do
# not correspond to sensible standalone demographic histories for model
# comparison; they were used only as internal diagnostics in the original
# analysis.
# =============================================================================

#' Partial Poisson sum: sum_{l=0}^{k} lambda^l / l!
#' @keywords internal
.partial_exp_sum <- function(lambda, k) {
  exp(lambda) * ppois(k, lambda)
}

#' Eigenvalues of the 2-deme migration matrix used throughout WH's solutions.
#' @keywords internal
.wh_lambdas <- function(M) {
  disc <- -4 * M + (1 + 2 * M)^2  # == 1 + 4*M^2
  disc[disc < 0] <- 0             # guard against floating-point noise at M ~ 0
  list(
    lambda1 = 0.5 * (1 + 2 * M - sqrt(disc)),
    lambda2 = 0.5 * (1 + 2 * M + sqrt(disc))
  )
}

#' Upper incomplete gamma function Gamma(s, x), via R's regularized lower
#' incomplete gamma (pgamma).
#' @keywords internal
.upper_incomplete_gamma <- function(s, x) {
  gamma(s) * pgamma(x, shape = s, lower.tail = FALSE)
}

#' Isolation-with-Migration (IM) model likelihood
#'
#' Blockwise probability of observing \code{k} pairwise differences under a
#' two-population model with symmetric migration at rate \code{M},
#' continuing from divergence time \code{tau} to the present.
#' Wilkinson-Herbots (2008), eq. 24.
#'
#' Setting \code{M = 0} recovers the Strict Isolation model exactly - this
#' is how \code{\link{wh_si}} is defined, and is not an independent formula.
#'
#' @param M scaled symmetric migration rate (M >= 0; M < 0 returns 0 for all k)
#' @param tau scaled divergence time
#' @param theta scaled mutation rate (per block)
#' @param alpha mutation-rate scalar for a second block/data type; set to 1
#'   when analysing a single block type (the default and typical use case)
#' @param k integer vector of blockwise pairwise-difference counts
#' @return numeric vector, the same length as \code{k}, of
#'   P(k | M, tau, theta) under the IM model
#' @examples
#' wh_im(M = 1.5, tau = 0.8, theta = 2, k = 0:20)
#' @export
wh_im <- function(M, tau, theta, alpha = 1, k) {
  if (M < 0) return(rep(0, length(k)))
  L <- .wh_lambdas(M)
  l1 <- L$lambda1; l2 <- L$lambda2

  sum_term <- function(lam, t) .partial_exp_sum((lam + theta) * t, k)

  A1 <- (l1 * theta^k / (l1 + theta)^(k + 1)) *
    (1 - exp(-tau * (l1 + theta)) * sum_term(l1, tau))
  B1 <- (exp(-tau * (l1 + theta)) * (theta * alpha)^k / (1 + theta * alpha)^(k + 1)) *
    .partial_exp_sum((1 / alpha + theta) * tau, k)

  A2 <- (l2 * theta^k / (l2 + theta)^(k + 1)) *
    (1 - exp(-tau * (l2 + theta)) * sum_term(l2, tau))
  B2 <- (exp(-tau * (l2 + theta)) * (theta * alpha)^k / (1 + theta * alpha)^(k + 1)) *
    .partial_exp_sum((1 / alpha + theta) * tau, k)

  p <- (l2 / (l2 - l1)) * (A1 + B1) - (l1 / (l2 - l1)) * (A2 + B2)
  pmax(p, 0)
}

#' Strict Isolation (SI) model likelihood
#'
#' The M = 0 special case of \code{\link{wh_im}}: two populations that
#' diverged at scaled time \code{tau} with no subsequent migration.
#'
#' @param tau scaled divergence time
#' @param theta scaled mutation rate (per block)
#' @param k integer vector of blockwise pairwise-difference counts
#' @return numeric vector, the same length as \code{k}, of
#'   P(k | tau, theta) under the SI model
#' @examples
#' wh_si(tau = 1, theta = 2, k = 0:20)
#' @export
wh_si <- function(tau, theta, k) wh_im(M = 0, tau = tau, theta = theta, alpha = 1, k = k)

#' Isolation-with-Initial-Migration (IIM) model likelihood
#'
#' Symmetric migration during an ancestral period, ceasing \code{tau1} time
#' units before the present (isolation for the most recent \code{tau1},
#' migration for the period stretching back from \code{tau1} to
#' \code{tau0_total} in the past). Wilkinson-Herbots (2012), eq. 29.
#'
#' @param M scaled symmetric migration rate during the ancestral period
#' @param tau1 recent period of isolation (time since migration ceased)
#' @param tau0_total total scaled time back to the ancestral population
#'   merger (must be > \code{tau1}; the migration period itself has
#'   duration \code{tau0_total - tau1})
#' @param theta scaled mutation rate (per block)
#' @param alpha mutation-rate scalar for a second block/data type; set to 1
#'   when analysing a single block type (the default and typical use case)
#' @param k integer vector of blockwise pairwise-difference counts
#' @return numeric vector, the same length as \code{k}, of
#'   P(k | M, tau1, tau0_total, theta) under the IIM model
#' @examples
#' wh_iim(M = 1.2, tau1 = 0.3, tau0_total = 2, theta = 2, k = 0:20)
#' @export
wh_iim <- function(M, tau1, tau0_total, theta, alpha = 1, k) {
  if (M < 0) return(rep(0, length(k)))
  L <- .wh_lambdas(M)
  l1 <- L$lambda1; l2 <- L$lambda2

  sum_term <- function(lam, t) .partial_exp_sum((lam + theta) * t, k)

  term1 <- function(lam) {
    (exp(-theta * tau1) * lam * theta^k / (lam + theta)^(k + 1)) * sum_term(lam, tau1)
  }
  term2 <- function(lam) {
    exp(-lam * (tau0_total - tau1) - theta * tau0_total) *
      ((theta * alpha)^k / (1 + theta * alpha)^(k + 1) * .partial_exp_sum((1 / alpha + theta) * tau0_total, k)
       - (lam * theta^k / (lam + theta)^(k + 1)) * sum_term(lam, tau0_total))
  }

  p <- (l2 / (l2 - l1)) * (term1(l1) + term2(l1)) -
       (l1 / (l2 - l1)) * (term1(l2) + term2(l2))
  pmax(p, 0)
}

#' Secondary Contact (SC) model likelihood
#'
#' A single, recent pulse of admixture of proportion \code{f} at scaled time
#' \code{tau1}, following isolation back to the ancestral divergence time
#' \code{tau0}. Ported from the "pulseADmix" function in the source
#' notebook. \code{tau1} and \code{tau0} are assumed positive (as enforced
#' by the default fitting bounds in \code{\link{fit_sc}}).
#'
#' @param f admixture/migration proportion at the pulse (0 <= f <= 1)
#' @param tau1 scaled time of the (recent) admixture pulse
#' @param tau0 scaled ancestral divergence time (tau0 > tau1)
#' @param theta scaled mutation rate (per block)
#' @param k integer vector of blockwise pairwise-difference counts
#' @return numeric vector, the same length as \code{k}, of
#'   P(k | f, tau1, tau0, theta) under the SC model
#' @examples
#' wh_sc(f = 0.3, tau1 = 0.2, tau0 = 1.5, theta = 2, k = 0:20)
#' @export
wh_sc <- function(f, tau1, tau0, theta, k) {
  s <- 1 + k
  pref <- theta^k / (factorial(k) * (1 + theta)^(1 + k))
  term_tau0 <- exp(tau0) * (f - 1) * (-.upper_incomplete_gamma(s, (1 + theta) * tau0))
  term_tau1 <- exp(tau1) * f       * (-.upper_incomplete_gamma(s, (1 + theta) * tau1))
  p <- pref * (term_tau0 - term_tau1)
  pmax(p, 0)
}
