# biomAid

[![R-CMD-check](https://github.com/DrJ001/bioMaid/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/DrJ001/bioMaid/actions/workflows/R-CMD-check.yml)
[![Codecov](https://codecov.io/gh/DrJ001/bioMaid/graph/badge.svg)](https://codecov.io/gh/DrJ001/bioMaid)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org/)

Biometric tools to support the analysis of comparative experiment data fitted
with **ASReml-R V4** mixed models. Hold tight theer a lot more tools coming.

---

## Installation

```r
# From GitHub (requires remotes or pak)
pak::pkg_install("DrJ001/bioMaid")

# or
remotes::install_github("DrJ001/bioMaid")
```

> **Note:** The core functions require an [ASReml-R V4](https://vsni.co.uk/software/asreml-r)
> licence. `simTrialData()` is fully standalone.

---

## Functions

### `compare()` — Pairwise comparison criteria

Computes **HSD**, **LSD**, or **Bonferroni**-corrected LSD for predicted values
from an ASReml-R V4 model, optionally within subgroups.

```r
# Tukey HSD for Genotype within each Site
res <- compare(model,
               term  = "Treatment:Site:Genotype",
               by    = "Site",
               type  = "HSD")
```

Two predicted values in the same group are significantly different when their
absolute difference exceeds the returned criterion value.

---

### `waldTest()` — Wald / F-tests on contrasts

Tests linear contrasts of predicted values using predicted information from
`predict.asreml()`. Supports pairwise, custom, and joint zero tests
with optional p-value adjustment.

```r
pred <- predict(model, classify = "Treatment", vcov = TRUE)

waldTest(pred,
         cc     = list(list(coef = c("N0","N1","N2"), type = "con", comp = "pairwise")),
         test   = "F",
         adjust = "fdr",
         df_error = model$nedf)
```

---

### `randomRegressMV()` — Random regression (BLUP-based)

Decomposes Genotype x Environment x (Treatment/Trait) variety BLUP result from
an ASReml-R v4 model to generate **responsiveness (or tolerance) indices** using
a natural genetic regression that forms from Gaussian conditional distribution
theory. Supports four conditioning schemes.

| `type` | Description |
|--------|-------------|
| `"baseline"` | Each non-baseline treatment conditioned on `levs[1]` |
| `"sequential"` | Gram-Schmidt orthogonalisation (diagonal TGmat) |
| `"partial"` | Each treatment conditioned on all others simultaneously |
| `"custom"` | User-specified conditioning sets |

```r
res <- randomRegressMV(model,
                       Env  = "TSite:Variety",
                       levs = c("N0","N1","N2"),
                       type = "sequential")

res$blups    # Site x Variety BLUPs + responsiveness indices
res$TGmat    # Transformed G-matrix
res$beta     # Site-specific regression coefficients
```

---

### `fixedRegressMV()` — Fixed regression (BLUE-based)

The fixed-effects analogue of `randomRegressMV()`. Regresses treatment BLUEs
via OLS within each group and returns **response or tolerance indices** for
every genotype. The same four conditioning schemes are available.

```r
res <- fixedRegressMV(model,
                      term = "Treatment:Site:Genotype",
                      by   = "Site",
                      levs = c("N0","N1","N2"),
                      type = "baseline")

res$blues    # BLUEs + response indices + HSD per group
res$beta     # OLS regression coefficients per group
```

---

### `fast()` — Factor Analytic Selection Tools

Implements the **FAST** (Smith & Cullis 2018) and **iClass** (Smith et al. 2021)
approaches for summarising variety performance from an FA mixed model.

```r
res <- fast(model,
            term   = "fa(Site, 3):Genotype",
            type   = "all",
            ic.num = 2L)

# Columns: Site, Genotype, loads1..k, score1..k, CVE, VE,
#          OP, dev, stab (FAST)
#          iclass, iClassOP, iClassRMSD (iClass)
```

---

### `simTrialData()` — Simulate trial data

Generates a balanced split-plot dataset for testing and demonstration, with
known ground-truth genetic parameters for comparison against model estimates.

```r
out <- simTrialData(nvar       = 20,
                    nsite      = 10,
                    treatments = c("N0","N1","N2"),
                    seed       = 2024,
                    verbose    = TRUE)

out$data          # plot-level data frame
out$params$beta   # true site-specific regression coefficients
```

---

## References

Smith, A.B. & Cullis, B.R. (2018). Plant breeding selection tools built on
factor analytic mixed models for multi-environment trial data.
*Euphytica*, **214**, 143.

Smith, A.B., Norman, A., Kuchel, H. & Cullis, B.R. (2021). Plant variety
selection using interaction classes derived from factor analytic linear mixed
models. *Frontiers in Plant Science*, **12**, 737462.

---

## License

MIT © J Taylor
