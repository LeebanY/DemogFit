# =============================================================================
# Converting coalescent-scaled parameter estimates (theta, tau, M) into
# real-world units (Ne, divergence time, migration rate per generation).
#
# Convention used (verified against the gIMble methods paper - Laetsch et
# al. 2023, PLoS Genetics - which states the scaled migration rate as
# "M = 4*Ne*me", i.e. the same "4Ne" scalar used for theta = 4*Ne*mu*L
# under the standard diploid-autosomal coalescent; this also matches the
# Ne/divergence-time scaling used in the original analysis this package
# replaces):
#
#   theta = c * Ne * mu * L            =>  Ne = theta / (c * mu * L)
#   M     = c * Ne * m                 =>  m  = M * mu * L / theta
#   1 coalescent time unit = (c / 2) * Ne generations
#                                       =>  T_generations = tau * theta / (2 * mu * L)
#
# where c is `ploidy_scalar` (4 for autosomal diploid loci - the default;
# 2 for haploid/mitochondrial/Y-linked loci; an intermediate value such as
# ~3 for X/Z-linked loci under equal sex ratios) and L is `block_length`.
#
# NOTE: c cancels algebraically in both the T_generations and m formulas
# above - it only affects the standalone Ne estimate. So an imprecise
# ploidy_scalar biases Ne, but not the divergence-time or migration-rate
# conclusions.
# =============================================================================

#' @keywords internal
.select_robust_model <- function(result) {
  if (!isTRUE(result$gene_flow)) return("SI")
  if (is.null(result$iim_test) || !isTRUE(result$iim_test$supported)) return("IM")
  "IIM"
}

#' Convert the most robust model's parameters into real-world units
#'
#' Automatically selects the simplest model in the nested SI/IM/IIM
#' hierarchy not rejected by the likelihood-ratio tests already computed by
#' \code{\link{fit_demography}} - i.e. the same model that workflow already
#' identifies as best-supported - and converts its \code{theta}, \code{tau}
#' (and \code{M}, if applicable) estimates into an effective population
#' size (Ne), a divergence time (in generations, and in years if
#' \code{generation_time_years} is supplied), and a migration rate per
#' generation.
#'
#' The model is chosen automatically and cannot be overridden here: SI is
#' selected if there's no evidence for gene flow, IM if there is evidence
#' for gene flow but no additional evidence of migration having ceased, and
#' IIM if migration having ceased is also supported. SC is never selected,
#' since it is not part of the nested hierarchy these tests apply to (see
#' \code{\link{fit_demography}}) - there is no statistically justified way
#' for this function to decide SC is "more robust" than SI/IM/IIM. If you
#' specifically want SC's parameters in real-world units, you can apply the
#' same formulas documented here to \code{result$fits$SC$par} yourself,
#' but be aware this package deliberately does not automate that, since SC
#' fits are never statistically validated against the other models.
#'
#' @param result a \code{demogfit_result} object, as returned by
#'   \code{\link{fit_demography}} (must include at least SI and IM; IIM is
#'   used if it was also fit)
#' @param mu mutation rate per site per generation (no default - always
#'   specific to your species/locus, so must be supplied explicitly)
#' @param block_length block length in bp: the \code{-l} value used when
#'   running \code{gIMble blocks} (i.e. the length of each block underlying
#'   your bSFS). No default, for the same reason as \code{mu}.
#' @param ploidy_scalar the coalescent scalar relating \code{theta}/\code{M}
#'   to Ne (default 4, for autosomal diploid loci). Use 2 for
#'   haploid/mitochondrial/Y-linked loci, or an intermediate value (e.g.
#'   ~3) for X/Z-linked loci under equal sex ratios. Only affects the
#'   reported Ne - divergence time and migration rate do not depend on it.
#' @param generation_time_years optional: generation time in years. If
#'   supplied, divergence time (and, for IIM, the migration duration) are
#'   also reported in years in addition to generations.
#' @return an object of class \code{demogfit_scaled} with elements:
#'   \code{model} (which of SI/IM/IIM was selected), \code{Ne},
#'   \code{divergence_time_generations}, \code{divergence_time_years}
#'   (\code{NULL} if \code{generation_time_years} not supplied),
#'   \code{migration_rate_per_generation} (0 for SI), and, for IIM only,
#'   \code{migration_duration_generations} / \code{migration_duration_years}.
#'   Has a \code{print} method.
#' @examples
#' set.seed(1)
#' p_true <- wh_im(M = 1.5, tau = 1, theta = 2, k = 0:60)
#' s_dist <- data.frame(k = 0:60, count = as.vector(rmultinom(1, 8000, p_true / sum(p_true))))
#' result <- fit_demography(s_dist)
#' scale_parameters(result, mu = 2.8e-9, block_length = 200)
#' scale_parameters(result, mu = 2.8e-9, block_length = 200, generation_time_years = 0.1)
#' @export
scale_parameters <- function(result, mu, block_length, ploidy_scalar = 4,
                              generation_time_years = NULL) {
  if (!inherits(result, "demogfit_result")) {
    stop("`result` must be a demogfit_result object, as returned by fit_demography().")
  }
  if (missing(mu) || missing(block_length)) {
    stop("`mu` and `block_length` must both be supplied explicitly - there is no sensible universal default across species.")
  }
  if (mu <= 0 || block_length <= 0 || ploidy_scalar <= 0) {
    stop("`mu`, `block_length` and `ploidy_scalar` must all be positive.")
  }

  model <- .select_robust_model(result)
  fit <- result$fits[[model]]
  if (is.null(fit)) {
    stop(sprintf("The selected model (%s) was not fit - make sure `models` in fit_demography() included it.", model))
  }
  theta <- fit$par[["theta"]]
  Ne <- theta / (ploidy_scalar * mu * block_length)

  gens_per_tau <- theta / (2 * mu * block_length)  # ploidy_scalar cancels; see file header

  if (model == "SI") {
    tau_total <- fit$par[["tau"]]
    M <- 0
    migration_duration_gen <- NA_real_
  } else if (model == "IM") {
    tau_total <- fit$par[["tau"]]
    M <- fit$par[["M"]]
    migration_duration_gen <- NA_real_
  } else { # IIM
    tau_total <- fit$par[["tau1"]] + fit$par[["tau0"]]
    M <- fit$par[["M"]]
    migration_duration_gen <- fit$par[["tau0"]] * gens_per_tau
  }

  T_generations <- tau_total * gens_per_tau
  m_per_gen <- if (M > 0) M * mu * block_length / theta else 0

  T_years <- if (!is.null(generation_time_years)) T_generations * generation_time_years else NULL
  migration_duration_years <- if (!is.null(generation_time_years) && !is.na(migration_duration_gen)) {
    migration_duration_gen * generation_time_years
  } else NULL

  structure(
    list(
      model = model, Ne = Ne,
      divergence_time_generations = T_generations,
      divergence_time_years = T_years,
      migration_rate_per_generation = m_per_gen,
      migration_duration_generations = migration_duration_gen,
      migration_duration_years = migration_duration_years,
      mu = mu, block_length = block_length, ploidy_scalar = ploidy_scalar,
      generation_time_years = generation_time_years,
      fit = fit
    ),
    class = "demogfit_scaled"
  )
}

#' @export
print.demogfit_scaled <- function(x, ...) {
  model_name <- c(SI = "Strict Isolation", IM = "Isolation-with-Migration",
                  IIM = "Isolation-with-Initial-Migration")[x$model]
  cat(sprintf("Scaled parameter estimates\n(most statistically robust model selected: %s)\n\n", model_name))
  cat(sprintf("Effective population size (Ne): %.4g\n", x$Ne))
  cat(sprintf("Divergence time: %.4g generations", x$divergence_time_generations))
  if (!is.null(x$divergence_time_years)) cat(sprintf(" (%.4g years)", x$divergence_time_years))
  cat("\n")
  cat(sprintf("Migration rate: %.4g per generation%s\n", x$migration_rate_per_generation,
              if (x$model == "SI") " (Strict Isolation - no migration)" else ""))
  if (x$model == "IIM") {
    cat(sprintf("Duration of ancestral migration: %.4g generations", x$migration_duration_generations))
    if (!is.null(x$migration_duration_years)) cat(sprintf(" (%.4g years)", x$migration_duration_years))
    cat("\n")
  }
  cat(sprintf("\n(mu = %.3g /site/gen, block length = %d bp, ploidy scalar = %g",
              x$mu, x$block_length, x$ploidy_scalar))
  if (!is.null(x$generation_time_years)) cat(sprintf(", generation time = %.3g yr", x$generation_time_years))
  cat(")\n")
  invisible(x)
}
