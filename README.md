# demogfit

An R package for fitting closed-form Wilkinson-Herbots (2008, 2012)
blockwise composite-likelihood demographic models to a pair of diverged
species/populations, from a gIMble blockwise site-frequency spectrum (bSFS).

## Models

| Model | Description | Parameters |
|---|---|---|
| **SI**  | Strict Isolation | `tau`, `theta` |
| **IM**  | Isolation-with-Migration (continuous migration to present) | `M`, `tau`, `theta` |
| **IIM** | Isolation-with-Initial-Migration (migration ceases before present) | `M`, `tau1`, `tau0`, `theta` |
| **SC**  | Secondary Contact (single recent admixture pulse) | `f`, `tau1`, `tau0`, `theta` |

SI, IM and IIM form a genuine 1-degree-of-freedom nested chain (SI is the
`M = 0` special case of IM; IM is the `tau1 = 0` special case of IIM), so
they are compared with likelihood-ratio tests by default. **SC is not
nested with any of the other three models** and is never included in
those tests - its fit is reported on its own.

This package's default model-comparison framework is two honestly-scoped
likelihood-ratio tests (SI vs IM, IM vs IIM), not `AIC()`/`BIC()` - on
purpose, since not all four models are mutually nested. An optional,
explicit AIC/log-likelihood ranking across *all* fitted models (including
SC) is also available - see "Model comparison: LRT vs AIC" below.

## Installation

```r
install.packages("demogfit_0.4.1.tar.gz", repos = NULL, type = "source")
library(demogfit)
```

Or, from the unzipped source directory:

```r
devtools::install("path/to/demogfit")
```

## Quick start

### Fitting one model directly

Every model has its own fitting function, each returning a `demogfit_fit`
object with `print()`, `summary()`, `coef()`, `logLik()` and `plot()`
methods - similar to `lm()` or `glm()`.

```r
library(demogfit)

# s_dist: a data.frame with columns k (pairwise differences per block)
# and count (number of blocks with that many differences)
set.seed(11)
p_true <- wh_im(M = 1.8, tau = 2.5, theta = 1.5, k = 0:40)
s_dist <- data.frame(k = 0:40, count = as.vector(rmultinom(1, 4000, p_true / sum(p_true))))

fit <- fit_im(s_dist)
fit
#> Isolation-with-Migration model
#> Log-likelihood: -8920.8541
#>
#> Parameter estimates
#>   M        = 2.57416
#>   tau      = 2.17991
#>   theta    = 1.62192
#>
#> Converged: Yes

coef(fit)
logLik(fit)
summary(fit)   # adds bounds, call, convergence detail
plot(fit)      # observed vs. fitted blockwise mutation-count distribution
```

### Specifying your own parameter bounds

Every fitting function accepts a `bounds` argument. Only the parameters you
name are overridden; everything else falls back to `default_bounds`.

```r
fit_im(s_dist, bounds = list(M = c(0, 30)))
fit_iim(s_dist, bounds = list(tau1 = c(0, 1), tau0 = c(0, 20)))
```

Inspect the built-in defaults with `default_bounds` (a named list, one
element per model). Parameters shared across SI/IM/IIM (`tau`, `theta`,
`M` where applicable) deliberately have *identical* bounds across those
models - a narrower search range for the simpler model in a
likelihood-ratio test biases the comparison, since it can be prevented
from reaching a genuinely good fit for reasons that have nothing to do
with the actual question being tested.

### Boundary-hit detection

Every `demogfit_fit` object has a `boundary_hits` field: any parameter
that converged within 1% of its own bound range gets flagged by name.
This matters because a parameter estimate sitting exactly on its bound is
the signature of a boundary-constrained optimum, not necessarily a
genuine converged estimate - the optimizer wanted to go further and the
bound stopped it. Both `print()`/`summary()` on an individual fit and
`fit_demography()`'s own top-level output surface this automatically:

```r
result <- fit_demography(s_dist)
summary(result)
#> ======================================================
#> Demographic inference summary
#> ======================================================
#>
#> Total blocks: 4000
#>
#> [!] WARNING: parameter(s) converged AT (or within 1% of) a bound:
#>      - SI: tau = 0.5 (at lower bound 0.5)
#>      - IIM: tau1 = 0.5 (at lower bound 0.5)
#>      - IIM: tau0 = 16 (at upper bound 16)
#>      - SC: tau1 = 0.5 (at lower bound 0.5)
#>     A boundary-constrained optimum is not necessarily a genuine
#>     converged estimate - see ?default_bounds.
#>
#> Gene flow
#>   SI   logLik = -9008.1997
#>   IM   logLik = -8920.8541
#>
#>   LRT(SI -> IM) = 174.6912, p-value = 6.99e-40
#>   Evidence for gene flow: YES
#>   Yes - there is statistical support for gene flow (IM fits significantly
#>   better than Strict Isolation).
#>
#> ------------------------------------------------------
#>
#> Isolation with Initial Migration
#>   IM   logLik = -8920.8541
#>   IIM  logLik = -8944.1615
#>
#>   LRT(IM -> IIM) = 0.0000, p-value = 1
#>   No - no additional support for migration having ceased; continuous IM
#>   is an adequate description.
#>
#> ------------------------------------------------------
#>
#> Secondary Contact
#>   logLik = -8931.0642
#>
#>   NOTE: the SC model is not nested within the SI/IM/IIM
#>   hierarchy. Its log-likelihood is reported for reference
#>   but is not statistically compared by this package.
#> ======================================================
```

This is worth understanding rather than being alarmed by: this `s_dist`
was simulated purely under **IM** (`M = 1.8, tau = 2.5, theta = 1.5`), so
`fit_im()` alone converges cleanly, no warning at all (see above). But
`SI`, `IIM` and `SC` are all being asked to fit data that doesn't
actually match their own structure - `SI` has no migration parameter to
absorb the real gene-flow signal at all, and `IIM`/`SC` have real IM
buried somewhere in their extra flexibility but no natural place to put
it, so their own extra parameters strain toward a boundary trying to
mimic IM as best they can. That's expected, not a bug - and it's exactly
what `boundary_hits` exists to make visible rather than silently returned
as if it were an ordinary result.

### The master function: is there evidence for gene flow?

`fit_demography()` runs the SI-vs-IM and IM-vs-IIM likelihood-ratio tests
automatically (shown in full above). If you only care about the
gene-flow question, fit just SI and IM:

```r
fit_demography(s_dist, models = c("SI", "IM"))
```

### Model comparison: LRT vs AIC

By default (`comparison_method = "lrt"`, unchanged from the examples
above), only the two likelihood-ratio tests run, and SC is never ranked
against the others. Pass `comparison_method = "aic"` or `"loglik"` to
additionally rank *every* fitted model, including SC:

```r
fit_demography(s_dist, comparison_method = "aic")$model_comparison$table
#>   model    logLik k      AIC
#> 1    IM -8920.854 3 17847.71
#> 2    SC -8931.064 4 17870.13
#> 3   IIM -8944.161 4 17896.32
#> 4    SI -9008.200 2 18020.40
```

`"aic"` and `"loglik"` give *identical* rankings for any two models with
equal parameter counts (the AIC penalty term cancels) - in particular,
SC vs IIM always agree between the two methods, since both have 4
parameters; they only diverge once SI/IM (fewer parameters) are also
part of the comparison. `"loglik"` prints an explicit caveat whenever
models of different sizes are being compared this way, since raw
log-likelihood always favours the more flexible model regardless of
whether that flexibility is actually supported by the data.

### A caveat about p-values at genome scale

These are composite-likelihood tests across many (approximately
independent) blocks. With genome-scale block counts - tens of thousands is
typical for a whole-genome pair - the likelihood-ratio test statistic
scales with the number of blocks, so **p-values become extremely small for
almost any non-zero effect, including biologically trivial amounts of gene
flow.** A tiny p-value tells you the effect is unlikely to be exactly
zero - it does not by itself tell you the effect is large or biologically
important. `fit_demography()` prints a reminder of this whenever a pair has
more than 2000 blocks; always look at the log-likelihoods and fitted
parameters (e.g. `M`, or a derived probability of migration) alongside the
p-value, not instead of it.

### From a raw gIMble bSFS

If you have raw bSFS output (one row per distinct `(hetA, hetB, hetAB,
fixed)` combination, with a block `count` for each), convert it to the
`s_dist` format first:

```r
bsfs_table <- data.frame(
  count = c(...), hetA = c(...), hetB = c(...), hetAB = c(...), fixed = c(...)
)
s_dist <- bsfs_to_s_distribution(bsfs_table)
fit_demography(s_dist)
```

### Many species pairs at once

```r
pair_data <- list(
  pair_a = s_dist_a,   # or a raw bSFS table, with convert_bsfs = TRUE
  pair_b = s_dist_b
)
results_table <- fit_demography_all(pair_data)
results_table

# drill into any one pair's full result:
attr(results_table, "results")$pair_a
```

`fit_demography_all()` accepts `comparison_method` too, passed through to
every pair; when set, the returned table gains a `best_model` column.

### Converting to real-world units (Ne, divergence time, migration rate)

`scale_parameters()` takes a `fit_demography()` result and converts the
most statistically robust model's `theta`/`tau`/`M` estimates into an
effective population size, a divergence time, and a migration rate per
generation. The model is chosen automatically - the simplest model in the
SI/IM/IIM chain not rejected by the tests `fit_demography()` already ran -
so you can't accidentally report scaled parameters for a poorly-supported
model. SC is never selected this way, since it isn't part of that nested
hierarchy.

```r
result <- fit_demography(s_dist)
scale_parameters(result, mu = 2.8e-9, block_length = 200)
#> Scaled parameter estimates
#> (most statistically robust model selected: Isolation-with-Migration)
#>
#> Effective population size (Ne): 3.724e+05
#> Divergence time: 3.187e+06 generations
#> Migration rate: 4.02e-07 per generation

# add a generation time (in years) to also get calendar-time estimates:
scale_parameters(result, mu = 2.8e-9, block_length = 200, generation_time_years = 0.1)
```

`mu` (mutation rate per site per generation) and `block_length` (the `-l`
value from your `gIMble blocks` step) must always be supplied explicitly -
there's no sensible universal default across species. `ploidy_scalar`
defaults to 4 (autosomal diploid loci); use 2 for haploid/mitochondrial/
Y-linked data, or an intermediate value (e.g. ~3) for X/Z-linked loci.
Note that `ploidy_scalar` only affects the reported Ne - divergence time
and migration rate are unaffected by it, since it cancels algebraically.

## What changed from the original analysis

This package replaces an original workflow that fit six models (SI, IM,
MIG, IIM, IIML, SC) in Mathematica, using a nested cascade of nine
delta-log-likelihood comparisons re-scaled to a reference of 500 blocks.
Two of those models (MIG, a stationary no-split migration model; IIML, the
infinite-ancestral-time limit of IIM) were diagnostic/boundary models
rather than biologically distinct histories and have been removed
entirely. The remaining four models' closed-form likelihoods are ported
term-for-term from the Wilkinson-Herbots (2008, 2012) solutions used in
the original notebook, but the statistical framework has been simplified
to two honestly-scoped likelihood-ratio tests (SI vs IM, IM vs IIM) plus a
standalone SC fit by default, rather than a nine-way cascade that included
comparing non-nested models against each other - with an optional,
explicit AIC ranking across all four models available when you want it.

## License

MIT - see [LICENSE.md](LICENSE.md).
