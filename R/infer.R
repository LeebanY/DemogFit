# =============================================================================
# Master inference functions: fit multiple models to one (or many) species
# pair(s).
#
# Deliberately simple statistical framework (no AIC/BIC, no multiple
# criteria for the same comparison):
#
#   SI vs IM   - a genuine 1-degree-of-freedom nested comparison (SI is the
#                M = 0 special case of IM). Tested via a likelihood-ratio
#                test: LRT = 2 * (logLik_IM - logLik_SI) ~ chi-sq(1) under
#                H0. This is the headline "is there evidence for gene
#                flow?" result.
#
#   IM vs IIM  - also a genuine 1-degree-of-freedom nested comparison (IM
#                is the tau1 = 0 special case of IIM - i.e. no distinct
#                recent isolation period after ancestral migration ceases).
#                Tested the same way, reported as a secondary result.
#
#   SC         - fit and reported on its own. SC is NOT nested with SI,
#                IM or IIM (it has the same number of free parameters as
#                IIM but a structurally different history - a single
#                admixture pulse rather than a continuous migration rate),
#                so no likelihood-ratio test or other statistical
#                comparison is computed for it. Its log-likelihood and
#                parameter estimates are reported for the user's own
#                reference only.
#
# CAVEAT (documented here and surfaced in print/summary output): these are
# composite-likelihood tests across many blocks. With genome-scale block
# counts (tens of thousands), p-values will become extremely small for
# almost any non-zero effect, including biologically trivial amounts of
# gene flow - the p-value alone does not tell you the effect is large,
# only that it is unlikely to be exactly zero. Always look at the
# log-likelihoods and fitted parameters (e.g. M, or an estimated
# probability of migration) alongside the p-value, not instead of it.
# =============================================================================

#' @keywords internal
.lrt_test <- function(loglik_null, loglik_alt, label_null, label_alt, alpha,
                       supported_msg, not_supported_msg) {
  LRT <- 2 * max(loglik_alt - loglik_null, 0)
  p_value <- stats::pchisq(LRT, df = 1, lower.tail = FALSE)
  supported <- p_value < alpha
  list(
    table = data.frame(model = c(label_null, label_alt),
                        logLik = c(loglik_null, loglik_alt),
                        stringsAsFactors = FALSE),
    LRT = LRT, df = 1L, p_value = p_value,
    supported = supported,
    interpretation = if (supported) supported_msg else not_supported_msg
  )
}

#' @keywords internal
.aic <- function(loglik, k) -2 * loglik + 2 * k

#' @keywords internal
#' Ranks every fitted model by AIC or by raw log-likelihood. Unlike the
#' SI-vs-IM and IM-vs-IIM likelihood-ratio tests above, this ranking makes
#' no nesting assumption, so it is the only place in this package SC is
#' ever compared against the other models numerically.
#'
#' AIC and raw logLik give IDENTICAL rankings for any pair of models with
#' the same number of free parameters (the AIC penalty term is equal on
#' both sides and cancels) - this is exactly the case for SC vs IIM (both
#' 4 parameters: SC has f/tau1/tau0/theta, IIM has M/tau1/tau0/theta), so
#' "aic" and "loglik" only actually disagree once models with different
#' parameter counts (e.g. SI or IM) are included in the same comparison.
.rank_models <- function(fits, method) {
  rows <- lapply(names(fits), function(nm) {
    fit <- fits[[nm]]
    k <- length(fit$par)
    data.frame(model = nm, logLik = fit$loglik, k = k,
               AIC = .aic(fit$loglik, k), stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, rows)
  table <- table[order(if (method == "aic") table$AIC else -table$logLik), ]
  rownames(table) <- NULL

  mixed_k <- length(unique(table$k)) > 1
  caveat <- if (method == "loglik" && mixed_k) {
    paste("NOTE: models being compared have different numbers of free",
          "parameters (raw log-likelihood always favours the more",
          "flexible model, regardless of whether the extra flexibility",
          "is actually supported by the data) - consider comparison_method",
          "= \"aic\" instead, which penalises this.")
  } else {
    NULL
  }

  list(table = table, method = method, best_model = table$model[1], caveat = caveat)
}

#' Fit demographic models to one species pair
#'
#' Fits any combination of the Strict Isolation (SI), Isolation-with-
#' Migration (IM), Isolation-with-Initial-Migration (IIM) and Secondary
#' Contact (SC) models to one pair's blockwise mutation-count distribution.
#'
#' Two likelihood-ratio tests are ALWAYS computed (regardless of
#' \code{comparison_method}), both against a 1-degree-of-freedom
#' chi-squared reference distribution, since both comparisons are
#' genuinely nested:
#' \itemize{
#'   \item \strong{SI vs IM} (if both are fit) - the headline "is there
#'     evidence for gene flow?" test.
#'   \item \strong{IM vs IIM} (if both are fit) - a secondary test of
#'     whether allowing migration to have ceased before the present
#'     improves the fit over continuous migration.
#' }
#' SC is not nested with SI, IM or IIM (same number of free parameters as
#' IIM, but a structurally different history - a single admixture pulse
#' rather than a continuous migration rate), so it is never included in
#' either LRT.
#'
#' Separately, \code{comparison_method} controls an additional, optional
#' ranking across ALL fitted models (including SC) with no nesting
#' assumption:
#' \itemize{
#'   \item \code{"lrt"} (default): only the two LRTs above; \code{$model_comparison}
#'     is \code{NULL}. Identical behaviour to previous versions of this package.
#'   \item \code{"aic"}: ranks every fitted model by AIC (\code{-2*logLik + 2*k}).
#'   \item \code{"loglik"}: ranks every fitted model by raw log-likelihood.
#'     Only meaningful on its own when every fitted model has the same
#'     number of parameters (e.g. comparing just IIM and SC) - a caveat is
#'     added to the result when models of different sizes are mixed, since
#'     raw log-likelihood always favours the more flexible model regardless
#'     of whether that flexibility is actually supported by the data.
#' }
#' AIC and raw log-likelihood give IDENTICAL rankings for any two models
#' with equal parameter counts - in particular, this means \code{"aic"}
#' and \code{"loglik"} always agree on SC vs IIM specifically, and only
#' diverge once SI and/or IM (2 and 3 parameters, vs IIM/SC's 4) are also
#' part of the comparison.
#'
#' @param s_dist a data.frame with columns \code{k} and \code{count} (see
#'   \code{\link{bsfs_to_s_distribution}})
#' @param models character vector of models to fit; must include at least
#'   \code{"SI"} and \code{"IM"} (the minimum needed to test for gene
#'   flow). Defaults to all four: \code{c("SI", "IM", "IIM", "SC")}.
#' @param bounds optional named list, keyed by model name (\code{SI},
#'   \code{IM}, \code{IIM}, \code{SC}), each element itself a bounds list
#'   passed through to the corresponding \code{fit_*()} function - e.g.
#'   \code{list(IM = list(M = c(0, 25)))}. See \code{\link{default_bounds}}.
#' @param no_blocks total number of blocks represented by \code{s_dist};
#'   if omitted, \code{sum(s_dist$count)} is used
#' @param alpha significance threshold for the likelihood-ratio tests
#'   (default 0.05)
#' @param comparison_method one of \code{"lrt"} (default), \code{"aic"} or
#'   \code{"loglik"} - see Details
#' @return an object of class \code{demogfit_result} with elements:
#'   \code{gene_flow} (logical or \code{NA}), \code{gene_flow_test} (list
#'   with \code{table}, \code{LRT}, \code{p_value}, \code{interpretation}),
#'   \code{iim_test} (same structure, \code{NULL} if IIM not fit),
#'   \code{model_comparison} (list with \code{table}, \code{method},
#'   \code{best_model}, \code{caveat}; \code{NULL} unless
#'   \code{comparison_method != "lrt"}), \code{fits} (named list of the
#'   underlying \code{demogfit_fit} objects), \code{no_blocks}, \code{alpha},
#'   \code{comparison_method}, \code{call}
#' @examples
#' set.seed(1)
#' p_true <- wh_im(M = 1.5, tau = 1, theta = 2, k = 0:60)
#' s_dist <- data.frame(k = 0:60, count = as.vector(rmultinom(1, 5000, p_true / sum(p_true))))
#' fit_demography(s_dist)
#' fit_demography(s_dist, models = c("SI", "IM"))
#' fit_demography(s_dist, comparison_method = "aic")$model_comparison
#' @export
fit_demography <- function(s_dist, models = c("SI", "IM", "IIM", "SC"),
                            bounds = list(), no_blocks = NULL, alpha = 0.05,
                            comparison_method = c("lrt", "aic", "loglik")) {
  comparison_method <- match.arg(comparison_method)
  models <- unique(match.arg(models, c("SI", "IM", "IIM", "SC"), several.ok = TRUE))
  if (!all(c("SI", "IM") %in% models)) {
    stop("`models` must include at least 'SI' and 'IM' to test for gene flow.")
  }
  if (is.null(no_blocks)) no_blocks <- sum(s_dist$count)

  fits <- list()
  if ("SI"  %in% models) fits$SI  <- fit_si(s_dist,  bounds = bounds$SI)
  if ("IM"  %in% models) fits$IM  <- fit_im(s_dist,  bounds = bounds$IM)
  if ("IIM" %in% models) fits$IIM <- fit_iim(s_dist, bounds = bounds$IIM)
  if ("SC"  %in% models) fits$SC  <- fit_sc(s_dist,  bounds = bounds$SC)

  gene_flow_test <- .lrt_test(
    fits$SI$loglik, fits$IM$loglik, "SI", "IM", alpha,
    supported_msg     = "Yes - there is statistical support for gene flow (IM fits significantly better than Strict Isolation).",
    not_supported_msg = "No - there is no statistical support for gene flow (Strict Isolation is an adequate description of this pair)."
  )

  iim_test <- NULL
  if ("IIM" %in% names(fits)) {
    iim_test <- .lrt_test(
      fits$IM$loglik, fits$IIM$loglik, "IM", "IIM", alpha,
      supported_msg     = "Yes - IIM (migration that ceased before the present) fits significantly better than continuous IM.",
      not_supported_msg = "No - no additional support for migration having ceased; continuous IM is an adequate description."
    )
  }

  model_comparison <- if (comparison_method != "lrt") .rank_models(fits, comparison_method) else NULL

  structure(
    list(
      gene_flow = gene_flow_test$supported,
      gene_flow_test = gene_flow_test,
      iim_test = iim_test,
      model_comparison = model_comparison,
      fits = fits,
      no_blocks = no_blocks,
      alpha = alpha,
      comparison_method = comparison_method,
      call = match.call()
    ),
    class = "demogfit_result"
  )
}

#' @keywords internal
.large_n_caveat <- function(no_blocks, threshold = 2000) {
  if (no_blocks > threshold) {
    cat(sprintf(
      "\nNote: this pair has %d blocks. With this many blocks, p-values from\ncomposite-likelihood tests become very small for almost any non-zero\neffect, including biologically trivial gene flow. Interpret the p-value\nalongside the log-likelihoods and fitted parameters, not on its own.\n",
      no_blocks))
  }
}

#' @keywords internal
.boundary_hits_across_fits <- function(fits) {
  msgs <- character(0)
  for (model_name in names(fits)) {
    fit <- fits[[model_name]]
    if (length(fit$boundary_hits) == 0) next
    for (nm in names(fit$boundary_hits)) {
      bound_val <- fit$bounds[[nm]][if (fit$boundary_hits[[nm]] == "lower") 1 else 2]
      msgs <- c(msgs, sprintf("%s: %s = %.6g (at %s bound %.6g)",
                               model_name, nm, fit$par[[nm]], fit$boundary_hits[[nm]], bound_val))
    }
  }
  msgs
}

#' @export
print.demogfit_result <- function(x, ...) {
  boundary_msgs <- .boundary_hits_across_fits(x$fits)
  if (length(boundary_msgs) > 0) {
    cat("[!] WARNING: one or more fitted models converged AT (or within 1% of) a\n")
    cat("    parameter bound - a boundary-constrained optimum, not necessarily a\n")
    cat("    genuine converged estimate. This can distort model comparison (a\n")
    cat("    model degenerating to mimic a simpler one at its own boundary can\n")
    cat("    look artificially competitive). Consider widening the relevant bound\n")
    cat("    via `bounds` and refitting:\n")
    for (m in boundary_msgs) cat("     -", m, "\n")
    cat("\n")
  }

  cat("Gene flow (SI vs IM)\n")
  print(x$gene_flow_test$table, row.names = FALSE)
  cat(sprintf("\nLRT = %.4f, df = 1, p-value = %.3g\n", x$gene_flow_test$LRT, x$gene_flow_test$p_value))
  cat(x$gene_flow_test$interpretation, "\n")

  if (!is.null(x$iim_test)) {
    cat("\n--------------------------------------------------\n")
    cat("Isolation with Initial Migration (IM vs IIM)\n")
    print(x$iim_test$table, row.names = FALSE)
    cat(sprintf("\nLRT = %.4f, df = 1, p-value = %.3g\n", x$iim_test$LRT, x$iim_test$p_value))
    cat(x$iim_test$interpretation, "\n")
  }

  if (!is.null(x$fits$SC) && is.null(x$model_comparison)) {
    cat("\n--------------------------------------------------\n")
    cat(sprintf("Secondary Contact: logLik = %.4f\n", x$fits$SC$loglik))
    cat("(not statistically compared to SI/IM/IIM - see ?fit_demography)\n")
  }

  if (!is.null(x$model_comparison)) {
    mc <- x$model_comparison
    cat("\n--------------------------------------------------\n")
    cat(sprintf("Model comparison (ranked by %s)\n", toupper(mc$method)))
    print(mc$table, row.names = FALSE)
    cat(sprintf("\nBest model: %s\n", mc$best_model))
    if (!is.null(mc$caveat)) cat(mc$caveat, "\n")
  }
  .large_n_caveat(x$no_blocks)
  invisible(x)
}

#' @export
summary.demogfit_result <- function(object, ...) {
  structure(object, class = c("summary.demogfit_result", class(object)))
}

#' @export
print.summary.demogfit_result <- function(x, ...) {
  cat("======================================================\n")
  cat("Demographic inference summary\n")
  cat("======================================================\n\n")
  cat(sprintf("Total blocks: %d\n\n", x$no_blocks))

  boundary_msgs <- .boundary_hits_across_fits(x$fits)
  if (length(boundary_msgs) > 0) {
    cat("[!] WARNING: parameter(s) converged AT (or within 1% of) a bound:\n")
    for (m in boundary_msgs) cat("     -", m, "\n")
    cat("    A boundary-constrained optimum is not necessarily a genuine\n")
    cat("    converged estimate - see ?default_bounds.\n\n")
  }

  cat("Gene flow\n")
  for (i in seq_len(nrow(x$gene_flow_test$table))) {
    row <- x$gene_flow_test$table[i, ]
    cat(sprintf("  %-4s logLik = %.4f\n", row$model, row$logLik))
  }
  cat(sprintf("\n  LRT(SI -> IM) = %.4f, p-value = %.3g\n", x$gene_flow_test$LRT, x$gene_flow_test$p_value))
  cat(sprintf("  Evidence for gene flow: %s\n", if (x$gene_flow) "YES" else "NO"))
  cat("  ", x$gene_flow_test$interpretation, "\n", sep = "")

  if (!is.null(x$iim_test)) {
    cat("\n------------------------------------------------------\n\n")
    cat("Isolation with Initial Migration\n")
    for (i in seq_len(nrow(x$iim_test$table))) {
      row <- x$iim_test$table[i, ]
      cat(sprintf("  %-4s logLik = %.4f\n", row$model, row$logLik))
    }
    cat(sprintf("\n  LRT(IM -> IIM) = %.4f, p-value = %.3g\n", x$iim_test$LRT, x$iim_test$p_value))
    cat("  ", x$iim_test$interpretation, "\n", sep = "")
  }

  if (!is.null(x$fits$SC) && is.null(x$model_comparison)) {
    cat("\n------------------------------------------------------\n\n")
    cat("Secondary Contact\n")
    cat(sprintf("  logLik = %.4f\n", x$fits$SC$loglik))
    cat("\n  NOTE: the SC model is not nested within the SI/IM/IIM\n")
    cat("  hierarchy. Its log-likelihood is reported for reference\n")
    cat("  but is not statistically compared by this package.\n")
  }

  if (!is.null(x$model_comparison)) {
    mc <- x$model_comparison
    cat("\n------------------------------------------------------\n\n")
    cat(sprintf("Model comparison (ranked by %s)\n\n", toupper(mc$method)))
    for (i in seq_len(nrow(mc$table))) {
      row <- mc$table[i, ]
      cat(sprintf("  %-4s logLik = %-12.4f k = %d   AIC = %.4f\n",
                   row$model, row$logLik, row$k, row$AIC))
    }
    cat(sprintf("\n  Best model: %s\n", mc$best_model))
    if (!is.null(mc$caveat)) cat("\n  ", mc$caveat, "\n", sep = "")
  }
  .large_n_caveat(x$no_blocks)
  cat("\n======================================================\n")
  invisible(x)
}

#' Fit demographic models across many species pairs
#'
#' Runs \code{\link{fit_demography}} over a named collection of pairs and
#' returns one summary row per pair. The full per-pair \code{demogfit_result}
#' objects are retained as an attribute (\code{"results"}) for anyone who
#' wants to drill into a specific pair's full output (e.g. its \code{plot()}
#' diagnostics).
#'
#' @param pair_data a named list; each element is either a data.frame with
#'   columns \code{k} and \code{count}, or, if \code{convert_bsfs = TRUE},
#'   a raw bSFS table with columns \code{count}, \code{hetA}, \code{hetB},
#'   \code{hetAB}, \code{fixed} (see \code{\link{bsfs_to_s_distribution}})
#' @param no_blocks optional named numeric vector/list of block counts per
#'   pair; if omitted, computed as \code{sum(count)} per pair
#' @param convert_bsfs if \code{TRUE}, run \code{\link{bsfs_to_s_distribution}}
#'   on each element of \code{pair_data} before fitting
#' @param ... additional arguments passed to \code{\link{fit_demography}}
#'   (e.g. \code{models}, \code{bounds}, \code{alpha})
#' @return a data.frame with one row per pair (\code{pair}, \code{no_blocks},
#'   \code{gene_flow}, \code{gene_flow_p_value}, \code{logLik_SI},
#'   \code{logLik_IM}, and if fit, \code{iim_p_value}, \code{logLik_IIM},
#'   \code{logLik_SC}, and if \code{comparison_method != "lrt"} was passed
#'   through \code{...}, \code{best_model}). Pairs that error out are
#'   dropped, with a warning. The full per-pair result objects are
#'   available via \code{attr(result, "results")}.
#' @examples
#' set.seed(1)
#' sim_pair <- function(M, tau, theta, n = 5000) {
#'   p <- wh_im(M = M, tau = tau, theta = theta, k = 0:60)
#'   data.frame(k = 0:60, count = as.vector(rmultinom(1, n, p / sum(p))))
#' }
#' pair_data <- list(
#'   pair_with_gene_flow = sim_pair(M = 1.5, tau = 1,   theta = 2),
#'   pair_no_gene_flow    = sim_pair(M = 0,   tau = 0.8, theta = 2)
#' )
#' fit_demography_all(pair_data, models = c("SI", "IM"))
#' @export
fit_demography_all <- function(pair_data, no_blocks = NULL, convert_bsfs = FALSE, ...) {
  pairs <- names(pair_data)
  full_results <- list()
  rows <- lapply(pairs, function(p) {
    s_dist <- if (convert_bsfs) bsfs_to_s_distribution(pair_data[[p]]) else pair_data[[p]]
    nb <- if (!is.null(no_blocks)) no_blocks[[p]] else NULL
    res <- tryCatch(
      fit_demography(s_dist, no_blocks = nb, ...),
      error = function(e) {
        warning(sprintf("Pair '%s' failed: %s", p, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    full_results[[p]] <<- res
    row <- data.frame(
      pair = p, no_blocks = res$no_blocks,
      gene_flow = res$gene_flow, gene_flow_p_value = res$gene_flow_test$p_value,
      logLik_SI = res$fits$SI$loglik, logLik_IM = res$fits$IM$loglik,
      stringsAsFactors = FALSE
    )
    if (!is.null(res$iim_test)) {
      row$iim_p_value <- res$iim_test$p_value
      row$logLik_IIM <- res$fits$IIM$loglik
    }
    if (!is.null(res$fits$SC)) {
      row$logLik_SC <- res$fits$SC$loglik
    }
    if (!is.null(res$model_comparison)) {
      row$best_model <- res$model_comparison$best_model
    }
    row
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  attr(out, "results") <- full_results
  out
}
