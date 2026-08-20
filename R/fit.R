# =============================================================================
# Per-model composite-likelihood fitting.
#
# Every fit_*() function below optimizes the composite log-likelihood
#   sum_k count_k * log( model(params, k) )
# via stats::optim(method = "L-BFGS-B"), mirroring the bounded
# FindMaximum[{expr, bounds}, {params}] calls in the source notebook.
#
# Default bounds match those recovered from the source notebook, EXCEPT for
# SC, whose exact fitting bounds could not be located in the notebooks
# supplied (only its closed form and usage string were found) - the
# defaults below are inferred to be consistent with the other models'
# parameter scales and should be checked against any historical SC fits
# you may have, and overridden via the `bounds` argument if they differ.
#
# NOTE ON MODEL COMPARISON: this package deliberately does NOT provide
# AIC()/BIC() methods. SI, IM and IIM form a genuine 1-degree-of-freedom
# nested chain (SI = IM at M = 0; IM = IIM at tau1 = 0) and are compared
# via likelihood-ratio tests in fit_demography(). SC is NOT nested with
# any of the others and is deliberately reported on its own, with no
# statistical comparison to the other three models.
# =============================================================================

#' Default parameter bounds for each model
#'
#' A named list of default lower/upper bounds for \code{\link{fit_si}},
#' \code{\link{fit_im}}, \code{\link{fit_iim}} and \code{\link{fit_sc}}.
#' Pass a (partial) override to the \code{bounds} argument of any fit_*()
#' function to change one or more parameters' bounds; unspecified
#' parameters fall back to these defaults.
#'
#' Tau-type parameters (\code{tau}, \code{tau1}, \code{tau0}) have a lower
#' bound of 0.5 across every model, and \code{M} has an upper bound of 5
#' in IM/IIM - both raised from earlier, much looser defaults after a
#' real fit converged to \code{tau0 = 0.001}, exactly its old lower
#' bound: a boundary-constrained optimum, not a genuine estimate that the
#' ancestral isolation period was near-instantaneous. An unconstrained-
#' looking parameter estimate that happens to sit exactly on its own
#' bound is a red flag, not a result - see \code{fit_*()}'s
#' \code{boundary_hits} field, which now detects and surfaces this
#' automatically regardless of where the bounds are set.
#' These are still just defaults, not hard limits on the underlying
#' model - override via \code{bounds} for a specific pair if you have a
#' real reason to think its true parameter value is more extreme
#' (e.g. a genuinely very recent split, or unusually high migration).
#' \code{f} (SC's admixture proportion) and SC's \code{tau1}/\code{tau0}/
#' \code{theta} are unchanged by the note above - \code{f} isn't on the
#' same scale as \code{M} and doesn't warrant the same ceiling, and SC is
#' not nested with SI/IM/IIM (see \code{\link{fit_demography}}), so its
#' bounds don't need to match theirs for a fair comparison the way
#' SI/IM/IIM's shared parameters do (see next paragraph).
#'
#' Parameters SHARED across SI, IM and IIM are deliberately given
#' IDENTICAL bounds, not independently-chosen ones per model: SI vs IM
#' and IM vs IIM are both tested via likelihood-ratio test in
#' \code{\link{fit_demography}}, and a narrower search range for the
#' simpler (null) model than the richer (alternative) model biases that
#' comparison - the null model could be prevented from reaching a
#' genuinely good fit purely because its own bounds were tighter, making
#' the alternative model look relatively better for a reason that has
#' nothing to do with whether gene flow is actually supported. Each
#' shared parameter uses the UNION of whatever ranges existed for it
#' historically across the models that share it, not an arbitrary pick
#' of one model's prior default: \code{tau} (SI, IM) is
#' \code{c(0.5, 15)}; \code{theta} (SI, IM, IIM) is \code{c(0.001, 10)};
#' \code{M} (IM, IIM) is \code{c(0, 5)}.
#'
#' @format A named list with one element per model (\code{SI}, \code{IM},
#'   \code{IIM}, \code{SC}), each itself a named list of
#'   \code{c(lower, upper)} pairs, one per parameter.
#' @export
default_bounds <- list(
  SI  = list(tau  = c(0.5, 15),  theta = c(1e-3, 10)),
  IM  = list(M    = c(0, 5), tau   = c(0.5, 15), theta = c(1e-3, 10)),
  # IIM bounds widened from the notebook's stated (M<16, tau0-duration<12) after
  # validating against the published 93-pair results: several real pairs' true
  # optima sit beyond those values (M up to ~20, tau0 duration up to ~15) -
  # see package development notes. M's upper bound was subsequently lowered
  # again from 25 to 5 (see the tau/M note above) - if you were relying on
  # the earlier wide M range for real data, override via bounds = list(M = c(0, 25)).
  IIM = list(M    = c(0,    5), tau1  = c(0.5, 3),  tau0 = c(0.5, 16), theta = c(1e-3, 10)),
  # SC bounds are inferred (see file header) - verify against prior fits if available.
  SC  = list(f    = c(1e-3, 1),  tau1  = c(0.5, 3),  tau0 = c(0.5, 15), theta = c(1e-3, 8))
)

#' @keywords internal
.resolve_bounds <- function(model, bounds = NULL) {
  defaults <- default_bounds[[model]]
  if (is.null(bounds) || length(bounds) == 0) return(defaults)
  bad <- setdiff(names(bounds), names(defaults))
  if (length(bad) > 0) {
    stop(sprintf("Unknown bound(s) for model %s: %s. Valid parameters: %s",
                 model, paste(bad, collapse = ", "), paste(names(defaults), collapse = ", ")))
  }
  modifyList(defaults, bounds)
}

#' @keywords internal
.bounds_to_vecs <- function(bounds_list) {
  list(lower = vapply(bounds_list, `[`, numeric(1), 1),
       upper = vapply(bounds_list, `[`, numeric(1), 2))
}

#' Composite negative log-likelihood, generic across models
#' @keywords internal
.nll <- function(par, s_dist, model_fn, fixed_args = list()) {
  args <- c(as.list(par), fixed_args, list(k = s_dist$k))
  p <- do.call(model_fn, args)
  p[p <= 0 | !is.finite(p)] <- .Machine$double.eps
  -sum(log(p) * s_dist$count)
}

#' Latin-hypercube-style stratified starting points
#'
#' Generates \code{n} starting points spread evenly across the bounds: each
#' parameter's range is divided into \code{n} equal-probability strata, and
#' the strata are independently permuted across parameters, so that a
#' small number of starts still covers the full range of every parameter
#' individually (unlike pure independent uniform draws, which can by chance
#' cluster and leave large gaps - see package development notes for a real
#' case where this caused a fit to regress after its bounds were widened).
#' @keywords internal
.lhs_starts <- function(n, lower, upper) {
  d <- length(lower)
  U <- matrix(stats::runif(n * d), nrow = n, ncol = d)
  for (j in seq_len(d)) {
    perm <- sample.int(n)
    U[, j] <- (perm - 1 + U[, j]) / n
  }
  sweep(sweep(U, 2, upper - lower, `*`), 2, lower, `+`)
}

#' @keywords internal
#' Flags any fitted parameter sitting within `tol` (as a fraction of that
#' parameter's own bound RANGE, not an absolute distance - so this is
#' meaningful across parameters on very different scales, e.g. tau vs M)
#' of either its lower or upper bound. This is the signature of a
#' boundary-constrained optimum, not necessarily a genuine converged
#' interior estimate - the optimizer wanted to go further and the bound
#' stopped it, which is a fundamentally different situation from "this is
#' where the likelihood surface actually peaks". A degenerate zero-width
#' bound (lower == upper, i.e. a deliberately fixed parameter) is never
#' flagged - there is no "interior" for such a parameter to be near or
#' far from.
.check_boundary_hits <- function(par, bounds_list, tol = 0.01) {
  hits <- character(0)
  for (nm in names(par)) {
    b <- bounds_list[[nm]]
    if (is.null(b)) next
    lower <- b[1]; upper <- b[2]
    range <- upper - lower
    if (range <= 0) next
    if ((par[[nm]] - lower) <= tol * range) {
      hits[nm] <- "lower"
    } else if ((upper - par[[nm]]) <= tol * range) {
      hits[nm] <- "upper"
    }
  }
  hits
}

#' Generic bounded composite-likelihood fit, with multi-start search
#'
#' Runs \code{stats::optim(method = "L-BFGS-B")} from several starting
#' points spread across the parameter bounds (the bounds midpoint plus
#' \code{n_starts - 1} Latin-hypercube-stratified points) and keeps the
#' result with the highest log-likelihood. This matters in practice:
#' several of these composite-likelihood surfaces (especially IM/IIM/SC,
#' which trade off migration rate against divergence time) have local
#' optima or long flat ridges, and a single start from the bounds midpoint
#' can converge to a substantially worse point than the true maximum on
#' real data.
#' @keywords internal
.fit_model <- function(model_fn, s_dist, bounds_list, fixed_args = list(),
                        start = NULL, model = NULL, call = NULL, n_starts = 20) {
  bv <- .bounds_to_vecs(bounds_list)
  starts <- list((bv$lower + bv$upper) / 2)
  if (!is.null(start)) starts <- c(list(start), starts)
  n_random <- max(0, n_starts - length(starts))
  if (n_random > 0) {
    lhs <- .lhs_starts(n_random, bv$lower, bv$upper)
    for (i in seq_len(n_random)) {
      starts[[length(starts) + 1]] <- stats::setNames(lhs[i, ], names(bv$lower))
    }
  }

  best <- NULL
  for (s0 in starts) {
    fit <- tryCatch(
      stats::optim(
        par = s0, fn = .nll, s_dist = s_dist, model_fn = model_fn,
        fixed_args = fixed_args, method = "L-BFGS-B",
        lower = bv$lower, upper = bv$upper,
        control = list(maxit = 500)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (is.null(best) || fit$value < best$value) best <- fit
  }
  if (is.null(best)) stop("Optimization failed from every starting point tried.")

  boundary_hits <- .check_boundary_hits(best$par, bounds_list)

  structure(
    list(model = model, par = best$par, loglik = -best$value, convergence = best$convergence,
         bounds = bounds_list, boundary_hits = boundary_hits, s_dist = s_dist, model_fn = model_fn,
         fixed_args = fixed_args, call = call, n_starts = length(starts)),
    class = "demogfit_fit"
  )
}

#' @export
print.demogfit_fit <- function(x, ...) {
  model_name <- c(SI = "Strict Isolation", IM = "Isolation-with-Migration",
                  IIM = "Isolation-with-Initial-Migration", SC = "Secondary Contact")[x$model]
  cat(sprintf("%s model\n", model_name %||% x$model %||% "demogfit"))
  cat(sprintf("Log-likelihood: %.4f\n\n", x$loglik))
  cat("Parameter estimates\n")
  for (nm in names(x$par)) cat(sprintf("  %-8s = %.6g\n", nm, x$par[[nm]]))
  cat(sprintf("\nConverged: %s\n", if (x$convergence == 0) "Yes" else sprintf("No (optim code %d)", x$convergence)))
  if (length(x$boundary_hits) > 0) {
    cat("\n[!] WARNING: the following parameter(s) converged AT (or within 1% of) their\n")
    cat("    own bound - this is the signature of a boundary-constrained optimum, not\n")
    cat("    necessarily a genuine converged estimate. Consider whether the bound\n")
    cat("    itself is reasonable, or widen it via `bounds` and refit:\n")
    for (nm in names(x$boundary_hits)) {
      cat(sprintf("      %-8s = %.6g (at %s bound %.6g)\n", nm, x$par[[nm]],
                   x$boundary_hits[[nm]], x$bounds[[nm]][if (x$boundary_hits[[nm]] == "lower") 1 else 2]))
    }
  }
  invisible(x)
}

#' @export
summary.demogfit_fit <- function(object, ...) {
  structure(object, class = c("summary.demogfit_fit", class(object)))
}

#' @export
print.summary.demogfit_fit <- function(x, ...) {
  model_name <- c(SI = "Strict Isolation", IM = "Isolation-with-Migration",
                  IIM = "Isolation-with-Initial-Migration", SC = "Secondary Contact")[x$model]
  cat("Call:\n")
  cat(if (!is.null(x$call)) paste(deparse(x$call), collapse = "\n") else "(not recorded)", "\n\n")
  cat(sprintf("Model: %s\n\n", model_name %||% x$model))
  cat(sprintf("Log-likelihood: %.4f\n", x$loglik))
  cat(sprintf("Total blocks:   %d\n\n", sum(x$s_dist$count)))
  cat("Parameter estimates and bounds\n")
  cat(sprintf("  %-8s %12s %12s %12s\n", "", "estimate", "lower", "upper"))
  for (nm in names(x$par)) {
    b <- x$bounds[[nm]]
    flag <- if (nm %in% names(x$boundary_hits)) "  <-- AT BOUND" else ""
    cat(sprintf("  %-8s %12.6g %12.6g %12.6g%s\n", nm, x$par[[nm]], b[1], b[2], flag))
  }
  cat(sprintf("\nOptimizer:  L-BFGS-B\n"))
  cat(sprintf("Converged:  %s\n", if (x$convergence == 0) "Yes" else sprintf("No (code %d)", x$convergence)))
  if (length(x$boundary_hits) > 0) {
    cat("\n[!] WARNING: parameter(s) marked 'AT BOUND' above converged AT (or within 1%\n")
    cat("    of) their own bound - a boundary-constrained optimum, not necessarily a\n")
    cat("    genuine converged estimate. See ?default_bounds.\n")
  }
  invisible(x)
}

#' Extract fitted parameters from a demogfit_fit object
#' @param object a \code{demogfit_fit} object, as returned by
#'   \code{\link{fit_si}}, \code{\link{fit_im}}, \code{\link{fit_iim}} or
#'   \code{\link{fit_sc}}
#' @param ... unused, present for S3 consistency
#' @return a named numeric vector of fitted parameter estimates
#' @export
coef.demogfit_fit <- function(object, ...) object$par

#' Extract the composite log-likelihood from a demogfit_fit object
#'
#' Returns a plain numeric value. Note: this package does NOT return a
#' classed \code{"logLik"} object and does NOT provide \code{AIC()}/
#' \code{BIC()} methods - SI/IM/IIM are compared via likelihood-ratio tests
#' in \code{\link{fit_demography}} instead (a simpler framework, appropriate
#' since not all four models are mutually nested).
#'
#' @param object a \code{demogfit_fit} object
#' @param ... unused, present for S3 consistency
#' @return a numeric scalar: the maximized composite log-likelihood
#' @export
logLik.demogfit_fit <- function(object, ...) object$loglik

#' Plot observed vs. fitted blockwise mutation-count distributions
#'
#' A simple base-R diagnostic: bars show the observed blockwise
#' mutation-count distribution (\code{s_dist}), the line shows the expected
#' distribution under the fitted model's parameters.
#'
#' @param x a \code{demogfit_fit} object
#' @param ... additional arguments passed to \code{graphics::barplot}
#' @return invisibly, \code{x}
#' @export
plot.demogfit_fit <- function(x, ...) {
  k <- x$s_dist$k
  observed <- x$s_dist$count
  n <- sum(observed)
  args <- c(as.list(x$par), x$fixed_args, list(k = k))
  expected_p <- do.call(x$model_fn, args)
  expected <- expected_p * n

  model_name <- c(SI = "Strict Isolation", IM = "Isolation-with-Migration",
                  IIM = "Isolation-with-Initial-Migration", SC = "Secondary Contact")[x$model]
  bp <- graphics::barplot(observed, names.arg = k, col = "grey85", border = NA,
                           xlab = "pairwise differences per block (k)", ylab = "number of blocks",
                           main = sprintf("%s: observed vs. fitted", model_name %||% x$model), ...)
  graphics::lines(bp, expected, col = "firebrick", lwd = 2, type = "b", pch = 16, cex = 0.6)
  graphics::legend("topright", legend = c("observed", "fitted"),
                    fill = c("grey85", NA), border = c(NA, NA),
                    lty = c(NA, 1), pch = c(NA, 16), col = c(NA, "firebrick"), bty = "n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Fit the Strict Isolation (SI) model
#'
#' Fits \code{\link{wh_si}} to a blockwise mutation-count distribution by
#' maximizing the composite log-likelihood over \code{tau} and \code{theta}.
#'
#' @param s_dist a data.frame with columns \code{k} and \code{count} (see
#'   \code{\link{bsfs_to_s_distribution}})
#' @param bounds optional named list to override one or more of the default
#'   bounds \code{list(tau = c(lower, upper), theta = c(lower, upper))}
#'   (see \code{\link{default_bounds}}); unspecified parameters use defaults
#' @param start optional named numeric vector of starting values
#'   (\code{tau}, \code{theta}); if supplied, used as one of the multi-start
#'   points alongside the bounds midpoint and random draws (see \code{n_starts})
#' @param n_starts number of optimizer starting points to try (the bounds
#'   midpoint, plus \code{n_starts - 1} points from a Latin-hypercube-style
#'   stratified design spread across the bounds); the best (highest-
#'   likelihood) result is kept. Defaults are higher for the
#'   higher-dimensional models (10 for SI's 2 parameters, 20 for IM's 3,
#'   30 for SC's 4, 40 for IIM's 4 - IIM's migration-rate/divergence-time
#'   ridge is the hardest to search reliably). Some of these composite-
#'   likelihood surfaces have local optima or long flat ridges (particularly
#'   trading off migration rate against divergence time), so a fit's
#'   log-likelihood that looks suspiciously low relative to a nested/related
#'   model is a good reason to increase this further.
#' @return an object of class \code{demogfit_fit}: a list with elements
#'   \code{par} (named numeric vector of fitted parameters),
#'   \code{loglik} (the maximized composite log-likelihood),
#'   \code{convergence} (the \code{optim()} convergence code; 0 = success),
#'   and \code{boundary_hits} (named character vector, one entry per
#'   parameter that converged within 1% of its own bound range of either
#'   its lower or upper bound - \code{"lower"} or \code{"upper"} - empty
#'   if none did; a non-empty result means that parameter's estimate is a
#'   boundary-constrained optimum, not necessarily where the likelihood
#'   surface actually peaks - see \code{\link{default_bounds}}). Has
#'   \code{print}, \code{summary}, \code{coef}, \code{logLik} and
#'   \code{plot} methods; both \code{print} and \code{summary} surface any
#'   \code{boundary_hits} as an explicit warning.
#' @examples
#' s_dist <- data.frame(k = 0:5, count = c(120, 40, 15, 6, 2, 1))
#' fit_si(s_dist)
#' fit_si(s_dist, bounds = list(theta = c(0.1, 20)))
#' @export
fit_si <- function(s_dist, bounds = NULL, start = NULL, n_starts = 10) {
  bl <- .resolve_bounds("SI", bounds)
  .fit_model(wh_si, s_dist, bl, start = start, model = "SI", call = match.call(), n_starts = n_starts)
}

#' Fit the Isolation-with-Migration (IM) model
#'
#' Fits \code{\link{wh_im}} to a blockwise mutation-count distribution by
#' maximizing the composite log-likelihood over \code{M}, \code{tau} and
#' \code{theta}.
#'
#' @inheritParams fit_si
#' @param bounds optional named list to override one or more of the default
#'   bounds \code{list(M = c(lower, upper), tau = c(lower, upper),
#'   theta = c(lower, upper))} (see \code{\link{default_bounds}})
#' @return an object of class \code{demogfit_fit} (see \code{\link{fit_si}})
#' @examples
#' s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
#' fit_im(s_dist)
#' fit_im(s_dist, bounds = list(M = c(0, 25)))
#' @export
fit_im <- function(s_dist, bounds = NULL, start = NULL, n_starts = 20) {
  bl <- .resolve_bounds("IM", bounds)
  .fit_model(wh_im, s_dist, bl, fixed_args = list(alpha = 1), start = start,
             model = "IM", call = match.call(), n_starts = n_starts)
}

#' Fit the Isolation-with-Initial-Migration (IIM) model
#'
#' Fits \code{\link{wh_iim}} to a blockwise mutation-count distribution by
#' maximizing the composite log-likelihood over \code{M}, \code{tau1},
#' \code{tau0} (migration-period duration) and \code{theta}. Internally
#' calls \code{\link{wh_iim}} with \code{tau0_total = tau1 + tau0}, so that
#' \code{tau0} here represents the *duration* of the ancestral migration
#' period, not the total time to the ancestral merger - this matches the
#' parametrization used for fitting in the source notebook.
#'
#' @inheritParams fit_si
#' @param bounds optional named list to override one or more of the default
#'   bounds \code{list(M = c(lower, upper), tau1 = c(lower, upper),
#'   tau0 = c(lower, upper), theta = c(lower, upper))}
#'   (see \code{\link{default_bounds}})
#' @return an object of class \code{demogfit_fit} (see \code{\link{fit_si}})
#' @examples
#' s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
#' fit_iim(s_dist)
#' @export
fit_iim <- function(s_dist, bounds = NULL, start = NULL, n_starts = 40) {
  bl <- .resolve_bounds("IIM", bounds)
  wrapped <- function(M, tau1, tau0, theta, k) {
    wh_iim(M = M, tau1 = tau1, tau0_total = tau1 + tau0, theta = theta, alpha = 1, k = k)
  }
  .fit_model(wrapped, s_dist, bl, start = start, model = "IIM", call = match.call(), n_starts = n_starts)
}

#' Fit the Secondary Contact (SC) model
#'
#' Fits \code{\link{wh_sc}} to a blockwise mutation-count distribution by
#' maximizing the composite log-likelihood over \code{f}, \code{tau1},
#' \code{tau0} and \code{theta}.
#'
#' @inheritParams fit_si
#' @param bounds optional named list to override one or more of the default
#'   bounds \code{list(f = c(lower, upper), tau1 = c(lower, upper),
#'   tau0 = c(lower, upper), theta = c(lower, upper))}. NOTE: these
#'   defaults are inferred, not recovered from the source notebook - see
#'   the note in \code{\link{default_bounds}}.
#' @return an object of class \code{demogfit_fit} (see \code{\link{fit_si}}).
#'   SC is not nested with SI/IM/IIM; \code{\link{fit_demography}} reports
#'   its fit for reference but never statistically compares it to the
#'   other three models.
#' @examples
#' s_dist <- data.frame(k = 0:5, count = c(100, 45, 20, 9, 4, 2))
#' fit_sc(s_dist)
#' @export
fit_sc <- function(s_dist, bounds = NULL, start = NULL, n_starts = 30) {
  bl <- .resolve_bounds("SC", bounds)
  .fit_model(wh_sc, s_dist, bl, start = start, model = "SC", call = match.call(), n_starts = n_starts)
}
