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

## Obtaining the input file

The package requires generating a blockwise site-frequency spectrum (bSFS) from the software package gIMble (https://github.com/LohseLab/gIMble). Alongside this I've also provided a Snakemake pipeline called [**DemogFit-Helper**](https://github.com/LeebanY/DemogFit-HELPER). The pipeline is designed to take a genome assembly (fasta format), paired fastq files (forward and reverse reads) from species A and species B (four fastq files total), an annotation file (gff format) and ideally a repeat annotation file (in bed format), all of which must be specified in metadata table files called "pairs.tsv" and "samples.tsv". The pipeline can then be run with default settings (config file can be adjusted) and should generate mapped reads, a variant-specific VCF and perform preprocessing steps for you in gimble. It will spit out a tally (.tsv) file that can then be downloaded and used as input for demog_fit. The pipeline can run as many pairs as you'd like so long as the metadata is properly filled in. 

## Quick start

### Fitting one model directly

Every model has its own fitting function, each returning a `demogfit_fit`
object with `print()`, `summary()`, `coef()`, `logLik()` and `plot()`. 

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
`M` where applicable) deliberately have identical bounds across those
models. A narrower search range for the simpler model in a
likelihood-ratio test biases the comparison, since it can be prevented
from reaching a genuinely good fit for reasons that have nothing to do
with the actual question being tested.

### Boundary-hit detection

Any parameter that converged too close to the bounds gets flagged. 
This matters because a parameter estimate sitting exactly on its bound is
the signature of a boundary-constrained optimum, not necessarily a
genuine converged estimate. The optimizer wanted to go further and the
bound stopped it. The default bounds have been set according to a reasonable
expectation for Drosophila species, but it's important that you make sensible
bound decisions for your system. Both `print()`/`summary()` on an individual fit and
`fit_demography()`'s own top-level output surface this automatically:

A simulated example of what boundary-hitting looks like below: 

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


### The quick-analysis function: is there evidence for gene flow?

If you just want a quick answer about what simple demographic history
best reflects what might be happening in your system, you can use the
'fit_demography()' function:

`fit_demography()` runs the SI-vs-IM and IM-vs-IIM likelihood-ratio tests
automatically (shown in full above). If you only care about the
gene-flow question, fit just SI and IM:

```r
fit_demography(s_dist, models = c("SI", "IM"))
```

### Model comparison: LRT or AIC

By default (`comparison_method = "lrt"`, unchanged from the examples
above), only the two likelihood-ratio tests run, and SC is never ranked
against the others because SC is non-nested with the other models. Pass 
`comparison_method = "aic"` or `"loglik"` to additionally rank *every* 
fitted model, including SC:

```r
fit_demography(s_dist, comparison_method = "aic")$model_comparison$table
#>   model    logLik k      AIC
#> 1    IM -8920.854 3 17847.71
#> 2    SC -8931.064 4 17870.13
#> 3   IIM -8944.161 4 17896.32
#> 4    SI -9008.200 2 18020.40
```

`"aic"` and `"loglik"` give *identical* rankings for any two models with
equal parameter counts (the AIC penalty term cancels). In particular,
SC vs IIM always agree between the two methods, since both have 4
parameters; they only diverge once SI/IM (fewer parameters) are also
part of the comparison. `"loglik"` prints an explicit caveat whenever
models of different sizes are being compared this way, since raw
log-likelihood always favours the more complex model regardless of
whether that flexibility is actually supported by the data.

### From a raw gIMble bSFS

The input for DemogFit is an S distirbution which converted by within-package functions.
If you have raw bSFS output, convert it to the `s_dist` format first:

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
generation. The model is chosen automatically following 'fit_demography()'
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

## Citation(s):
If you have used this package for your work (thanks!), please cite the following:
- Yusuf, L.H., Laetsch, D.R., Lohse, K. and Ritchie, M.G., 2026. Genomic analyses in Drosophila do not support the classic allopatric model of speciation. Evolution Letters, 10(2), pp.186-194.
- Wilkinson-Herbots, H.M., 2012. The distribution of the coalescence time and the number of pairwise nucleotide differences in a model of population divergence or speciation with an initial period of gene flow. Theoretical population biology, 82(2), pp.92-108.
- Wilkinson-Herbots, H.M., 2008. The distribution of the coalescence time and the number of pairwise nucleotide differences in the “isolation with migration” model. Theoretical Population Biology, 73(2), pp.277-288.

## A small note:
Hi all, for transparency sake: much of the package was converted from Mathematica scripts to the R package with the help of Claude (Sonnet/Opus) models. This dramatically sped up the time it took to come up with these packages and I have tried to sanity check the results as much as possible, but LLMs are error-prone and something could've definitely slipped by. Of course, if you do spot something, please do raise an issue and I will try to fix it as quickly as I can. 

## License

MIT - see [LICENSE.md](LICENSE.md).
