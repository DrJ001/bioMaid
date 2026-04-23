# biomAid

[![R-CMD-check](https://github.com/DrJ001/bioMaid/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/DrJ001/bioMaid/actions/workflows/R-CMD-check.yml)
[![Codecov](https://codecov.io/gh/DrJ001/bioMaid/graph/badge.svg)](https://codecov.io/gh/DrJ001/bioMaid)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org/)

Biometric tools to support the analysis of comparative experiment data fitted
with **ASReml-R V4** mixed models. Watch this space, there are a lot more functions coming.

---

## Installation

```r
# From GitHub (requires remotes or pak)
pak::pkg_install("DrJ001/bioMaid")

# or
remotes::install_github("DrJ001/bioMaid")
```

> **Note:** The core modelling functions require an [ASReml-R V4](https://vsni.co.uk/software/asreml-r)
> licence. `simTrialData()` and all `plot_*()` functions are fully standalone.

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
         cc       = list(list(coef = c("N0","N1","N2"), type = "con", comp = "pairwise")),
         test     = "F",
         adjust   = "fdr",
         df_error = model$nedf)
```

---

### `plot_waldTest()` — Forest plot for Wald test contrasts

Produces a publication-ready forest plot from the output of `waldTest()`.
Each contrast is shown as a filled circle with horizontal confidence interval
bars. Points are coloured by **−log₁₀(p)**: non-significant results appear in
grey, with a warm gradient from gold through to dark red as evidence
strengthens. The raw p-value is printed beside each row.

```r
res <- waldTest(pred,
                cc = list(list(coef = c("N0","N1","N2"),
                               type = "con",
                               comp = "pairwise")))

# Simple forest plot
plot_waldTest(res)

# By-group — one facet panel per site (default)
res2 <- waldTest(pred_site,
                 cc = list(list(coef = c("N0","N1","N2"),
                                type = "con", comp = "pairwise")),
                 by = "Site")
plot_waldTest(res2, facet = TRUE)

# All groups on one panel with separator lines
plot_waldTest(res2, facet = FALSE)

# 99% CI, stricter significance threshold
plot_waldTest(res, ci_level = 0.99, alpha = 0.01)
```

| Argument | Description |
|----------|-------------|
| `facet` | `TRUE` (default): one panel per `by`-group; `FALSE`: single panel |
| `ci_level` | Confidence level for CI arms. Default `0.95` |
| `alpha` | Significance threshold for colour-scale break. Default `0.05` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `randomRegress()` — Random regression (BLUP-based)

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
res <- randomRegress(model,
                     Env  = "TSite:Variety",
                     levs = c("N0","N1","N2"),
                     type = "sequential")

res$blups    # Site x Variety BLUPs + responsiveness indices
res$TGmat    # Transformed G-matrix
res$beta     # Site-specific regression coefficients
```

---

### `plot_randomRegress()` — Visualise random regression results

Generates ggplot2 visualisations from `randomRegress()` output. Three plot
types are available, each returned as a ggplot object that can be further
customised with `+`.

```r
# Regression scatter — per-site BLUPs with site-specific β line
plot_randomRegress(res, type = "regress")

# Quadrant plot — efficiency (x) vs responsiveness (y)
plot_randomRegress(res, type = "quadrant")

# G-matrix correlation heatmap
plot_randomRegress(res, type = "gmat")

# Highlight six archetypal varieties (3 top-right, 3 bottom-left)
plot_randomRegress(res, type = "quadrant", highlight = "default")
```

| Argument | Description |
|----------|-------------|
| `type` | `"regress"`, `"quadrant"`, or `"gmat"` |
| `treatments` | Character vector to restrict conditioning pairs plotted |
| `highlight` | `"default"` auto-selects 6 archetypes; `NULL` no highlighting |
| `centre` | `TRUE` adds back within-site means (useful for demo data) |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `fixedRegress()` — Fixed regression (BLUE-based)

The fixed-effects analogue of `randomRegress()`. Regresses treatment BLUEs
via OLS within each group and returns **response or tolerance indices** for
every genotype. The same four conditioning schemes are available.

```r
res <- fixedRegress(model,
                    term = "Treatment:Site:Genotype",
                    by   = "Site",
                    levs = c("N0","N1","N2"),
                    type = "baseline")

res$blues    # BLUEs + response indices + HSD per group
res$beta     # OLS regression coefficients per group
```

---

### `plot_fixedRegress()` — Visualise fixed regression results

Generates ggplot2 visualisations from `fixedRegress()` output. Two plot
types are available, each returned as a ggplot object.

```r
# Regression scatter — BLUEs with OLS regression line per group
plot_fixedRegress(res, type = "regress")

# Quadrant plot — efficiency (x) vs response index (y)
plot_fixedRegress(res, type = "quadrant")

# Custom significance highlighting
plot_fixedRegress(res, type = "quadrant", highlight = "default")

# Raw (uncentred) axes
plot_fixedRegress(res, type = "regress", centre = FALSE)
```

| Argument | Description |
|----------|-------------|
| `type` | `"regress"` or `"quadrant"` |
| `treatments` | Character vector to restrict conditioning pairs plotted |
| `highlight` | `"default"` auto-selects 6 archetypes; `NULL` no highlighting |
| `centre` | `TRUE` (default) subtracts within-group means; `FALSE` raw BLUEs |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `padTrial()` — Extract and pad a field trial layout

Extracts a rectangular sub-trial from a field layout data frame based on a
plot-type classification, trims guard rows outside the bounding box of the
target plot types, and pads any missing interior grid positions with blank rows.
Useful for preparing irregular trial layouts for spatial analysis.

```r
res <- padTrial(data,
                pattern    = "Row:Column",
                match      = "DH",
                split      = "Block",
                pad        = TRUE,
                keep       = "Block",
                fill_value = "Blank",
                type_col   = "Type",
                verbose    = TRUE)
```

| Argument | Description |
|----------|-------------|
| `pattern` | `"Row:Column"` — names of the spatial coordinate columns |
| `match` | Plot type(s) to include (e.g. `"DH"`, `c("DH","Check")`) |
| `split` | Column(s) to split by (e.g. `"Block"`, `c("Site","Block")`); `NULL` = whole dataset |
| `pad` | `TRUE` (default) inserts blank rows for missing grid cells |
| `keep` | Column(s) whose values are carried into padded rows |
| `fill_value` | String written into all character/factor columns of padded rows |
| `verbose` | `TRUE` prints a per-group summary message |

The result is a data frame with an `add` column: `"old"` for original rows,
`"new"` for padded rows. Row and Column are returned as ascending factors.

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
