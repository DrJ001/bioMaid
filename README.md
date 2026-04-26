# biomAid

[![R-CMD-check](https://github.com/DrJ001/biomAid/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/DrJ001/biomAid/actions/workflows/R-CMD-check.yml)
[![Codecov](https://codecov.io/gh/DrJ001/biomAid/graph/badge.svg)](https://codecov.io/gh/DrJ001/biomAid)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org/)

Biometric tools to support the analysis of comparative experiment data fitted
with **ASReml-R V4** mixed models. Watch this space, there are a lot more functions coming.

---

## Vignettes

| Vignette | Description |
|----------|-------------|
| [Wald Tests on Fixed-Effect Contrasts](https://DrJ001.github.io/biomAid/waldTest.html) | Mathematical framework and worked examples for `waldTest()` and `plot_waldTest()` |
| [Multivariate Random Regression](https://DrJ001.github.io/biomAid/randomRegress.html) | Conditioning schemes, efficiency/responsiveness decomposition, and all plot types for `randomRegress()` and `plot_randomRegress()` |

---

## Installation

```r
# From GitHub (requires remotes or pak)
pak::pkg_install("DrJ001/biomAid")

# or
remotes::install_github("DrJ001/biomAid")
```

> **Note:** The core modelling functions require an [ASReml-R V4](https://vsni.co.uk/software/asreml-r)
> licence. `simTrialData()` and all `plot_*()` functions are fully standalone.

---

## Functions

### `compare()` — Pairwise comparison criteria

Computes **HSD**, **LSD**, or **Bonferroni**-corrected LSD for predicted values
from an ASReml-R V4 model, optionally within subgroups. Two predicted values in
the same group are significantly different when their absolute difference exceeds
the returned criterion value.

```r
compare(model, term, by = NULL,
        type  = c("HSD", "LSD", "Bonferroni"),
        pev   = TRUE,
        alpha = 0.05,
        ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | Classify string passed to `predict.asreml()` |
| `by` | Column(s) to split comparisons by. `NULL` = one group |
| `type` | `"HSD"` (default), `"LSD"`, or `"Bonferroni"` |
| `pev` | `TRUE` (default) uses prediction error variance; `FALSE` uses posterior variance |
| `alpha` | Significance level. Default `0.05` |

---

### `waldTest()` — Wald / F-tests on contrasts

Tests linear contrasts of predicted values using prediction error information
from `predict.asreml()`. Supports pairwise, custom contrast matrix, and joint
zero tests with optional p-value adjustment.

```r
waldTest(pred, cc, by = NULL,
         test     = c("Wald", "F"),
         df_error = NULL,
         adjust   = c("none", "bonferroni", "holm", "fdr", "BH", "BY"))
```

| Argument | Description |
|----------|-------------|
| `pred` | List returned by `predict(model, vcov = TRUE)` |
| `cc` | Named list of test specifications (`coef`, `type`, `comp`, `group`) |
| `by` | Column(s) to run tests within. `NULL` = single group |
| `test` | `"Wald"` (default, χ² statistic) or `"F"` (requires `df_error`) |
| `df_error` | Denominator degrees of freedom for F-tests (e.g. `model$nedf`) |
| `adjust` | P-value adjustment method. Default `"none"` |

---

### `plot_waldTest()` — Forest plot for Wald test contrasts

Produces a publication-ready forest plot from the output of `waldTest()`.
Each contrast is shown as a filled circle with horizontal confidence interval
bars. Points are coloured by **−log₁₀(p)**: non-significant results appear in
grey with a warm gradient from gold through to dark red as evidence strengthens.
The raw p-value is printed beside each row.

```r
plot_waldTest(res,
              facet       = TRUE,
              ci_level    = 0.95,
              alpha       = 0.05,
              theme       = ggplot2::theme_bw(),
              return_data = FALSE,
              ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `waldTest()` |
| `facet` | `TRUE` (default): one panel per `by`-group; `FALSE`: single panel |
| `ci_level` | Confidence level for CI arms. Default `0.95` |
| `alpha` | Significance threshold for colour-scale break. Default `0.05` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `randomRegress()` — Random regression (BLUP-based)

Decomposes Genotype × Environment × (Treatment/Trait) variety BLUPs from an
ASReml-R V4 model into **efficiency and responsiveness indices** using a natural
genetic regression derived from Gaussian conditional distribution theory.
Supports four conditioning schemes.

```r
randomRegress(model, Env = "TSite:Variety", levs = NULL,
              type = "baseline", cond = NULL,
              sep = "-", pev = TRUE, ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `Env` | Classify string identifying the G-matrix term. Default `"TSite:Variety"` |
| `levs` | Character vector of treatment levels to decompose |
| `type` | `"baseline"`, `"sequential"`, `"partial"`, or `"custom"` |
| `cond` | User-supplied conditioning list when `type = "custom"` |
| `sep` | Separator between site and treatment in composite labels. Default `"-"` |
| `pev` | `TRUE` (default) uses PEV; `FALSE` uses posterior variance |

---

### `plot_randomRegress()` — Visualise random regression results

Generates ggplot2 visualisations from `randomRegress()` output. Three plot
types are available, each returned as a ggplot object that can be further
customised with `+`.

```r
plot_randomRegress(res,
                   type        = c("regress", "quadrant", "gmat"),
                   treatments  = NULL,
                   highlight   = "default",
                   centre      = FALSE,
                   theme       = ggplot2::theme_bw(),
                   return_data = FALSE,
                   ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `randomRegress()` |
| `type` | `"regress"`, `"quadrant"`, or `"gmat"` |
| `treatments` | Character vector to restrict conditioning pairs plotted. `NULL` = all |
| `highlight` | `"default"` auto-selects 6 archetypes; `NULL` = no highlighting |
| `centre` | `TRUE` adds back within-site means (useful for demo data). Default `FALSE` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `fixedRegress()` — Fixed regression (BLUE-based)

The fixed-effects analogue of `randomRegress()`. Regresses treatment BLUEs via
OLS within each group and returns **efficiency and response indices** for every
genotype. The same four conditioning schemes are available.

```r
fixedRegress(model, term = "Treatment:Genotype",
             by = NULL, levs = NULL,
             type = "baseline", cond = NULL,
             min_obs = NULL, ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | Classify string passed to `predict.asreml()`. Default `"Treatment:Genotype"` |
| `by` | Column(s) defining groups for separate regressions |
| `levs` | Character vector of treatment levels to decompose |
| `type` | `"baseline"`, `"sequential"`, `"partial"`, or `"custom"` |
| `cond` | User-supplied conditioning list when `type = "custom"` |
| `min_obs` | Minimum common genotypes required to fit a regression. `NULL` = auto |

---

### `plot_fixedRegress()` — Visualise fixed regression results

Generates ggplot2 visualisations from `fixedRegress()` output. Two plot types
are available, each returned as a ggplot object that can be further customised
with `+`.

```r
plot_fixedRegress(res,
                  type        = c("regress", "quadrant"),
                  treatments  = NULL,
                  highlight   = "default",
                  centre      = TRUE,
                  theme       = ggplot2::theme_bw(),
                  return_data = FALSE,
                  ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `fixedRegress()` |
| `type` | `"regress"` or `"quadrant"` |
| `treatments` | Character vector to restrict conditioning pairs plotted. `NULL` = all |
| `highlight` | `"default"` auto-selects 6 archetypes; `NULL` = no highlighting |
| `centre` | `TRUE` (default) subtracts within-group means; `FALSE` = raw BLUEs |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `padTrial()` — Extract and pad a field trial layout

Extracts a rectangular sub-trial from a field layout data frame based on a
plot-type classification, trims guard rows outside the bounding box of the
target plot types, and pads any missing interior grid positions with blank rows.
Useful for preparing irregular trial layouts for spatial analysis.

```r
padTrial(data,
         pattern    = "Row:Column",
         match      = "DH",
         split      = "Block",
         pad        = TRUE,
         keep       = split,
         fill_value = "Blank",
         type_col   = "Type",
         verbose    = FALSE)
```

| Argument | Description |
|----------|-------------|
| `data` | Data frame containing the trial layout |
| `pattern` | Colon-separated names of the spatial coordinate columns. Default `"Row:Column"` |
| `match` | Plot type(s) defining the target sub-trial (e.g. `"DH"`, `c("DH","Check")`) |
| `split` | Column(s) to process independently (e.g. `"Block"`); `NULL` = whole dataset |
| `pad` | `TRUE` (default) inserts blank rows for missing grid cells |
| `keep` | Column(s) whose values are carried into padded rows. Defaults to `split` |
| `fill_value` | String written into character/factor columns of padded rows. Default `"Blank"` |
| `type_col` | Name of the plot-type column. Default `"Type"` |
| `verbose` | `TRUE` prints a per-group summary message |

---

### `plot_padTrial()` — Before/after field layout tile map

Visualises the effect of `padTrial()` as a pair of tile maps — **Before** on
top, **After** below. Tiles are coloured by plot type; inserted (padded) cells
appear in light grey. When `data` is supplied the Before panel shows the
complete original layout including guard rows.

```r
plot_padTrial(result,
              data        = NULL,
              type_col    = "Type",
              pattern     = "Row:Column",
              split       = NULL,
              label       = NULL,
              theme       = ggplot2::theme_bw(),
              return_data = FALSE,
              ...)
```

| Argument | Description |
|----------|-------------|
| `result` | Data frame returned by `padTrial()` |
| `data` | Original data passed to `padTrial()`. `NULL` (default) reconstructs Before from result |
| `type_col` | Name of the plot-type column used to colour tiles. Default `"Type"` |
| `pattern` | Colon-separated names of the spatial coordinate columns. Default `"Row:Column"` |
| `split` | Grouping column(s) matching the `split` used in `padTrial()` |
| `label` | Column whose values are printed inside each tile. `NULL` = no labels |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

---

### `fast()` — Factor Analytic Selection Tools

Implements the **FAST** (Smith & Cullis 2018) and **iClass** (Smith et al. 2021)
approaches for summarising variety performance from an FA mixed model.

```r
fast(model, term = "fa(Site, 4):Genotype",
     type   = c("all", "FAST", "iClass"),
     ic.num = 2L,
     ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | FA model term string. Default `"fa(Site, 4):Genotype"` |
| `type` | `"all"` (default), `"FAST"`, or `"iClass"` |
| `ic.num` | Number of interaction class factors (1 to k). Default `2` |

---

### `simTrialData()` — Simulate trial data

Generates a balanced split-plot dataset for testing and demonstration, with
known ground-truth genetic parameters for comparison against model estimates.

```r
simTrialData(nvar           = 20L,
             nsite          = 10L,
             treatments     = c("N0", "N1", "N2"),
             nrep           = 3L,
             rows_per_strip = NULL,
             cols_per_strip = NULL,
             site_mean      = 4500,
             site_sd        = 600,
             treat_effects  = NULL,
             sigma_base     = 250,
             sigr           = NULL,
             beta_min       = NULL,
             beta_max       = NULL,
             rep_sd         = 150,
             row_sd         = 80,
             col_sd         = 60,
             error_sd       = 350,
             sep            = "-",
             variety_prefix = "Var",
             site_prefix    = "Env",
             seed           = NULL,
             outfile        = NULL,
             verbose        = TRUE)
```

| Argument | Description |
|----------|-------------|
| `nvar` | Number of varieties. Default `20` |
| `nsite` | Number of sites. Default `10` |
| `treatments` | Character vector of treatment labels. Default `c("N0","N1","N2")` |
| `nrep` | Number of replicates per site. Default `3` |
| `rows_per_strip` | Rows per treatment strip. `NULL` = auto |
| `cols_per_strip` | Columns per treatment strip. `NULL` = auto |
| `site_mean` | Mean yield across sites. Default `4500` |
| `site_sd` | Site-to-site SD. Default `600` |
| `treat_effects` | Named numeric vector of treatment mean shifts. `NULL` = auto |
| `sigma_base` | Base genetic SD. Default `250` |
| `sigr` | Per-treatment residual SDs. `NULL` = auto |
| `beta_min` | Minimum regression coefficient. `NULL` = auto |
| `beta_max` | Maximum regression coefficient. `NULL` = auto |
| `rep_sd` | Replicate effect SD. Default `150` |
| `row_sd` | Row spatial effect SD. Default `80` |
| `col_sd` | Column spatial effect SD. Default `60` |
| `error_sd` | Residual error SD. Default `350` |
| `sep` | Separator for composite site-treatment labels. Default `"-"` |
| `variety_prefix` | Prefix for variety labels. Default `"Var"` |
| `site_prefix` | Prefix for site labels. Default `"Env"` |
| `seed` | Random seed for reproducibility. `NULL` = no fixed seed |
| `outfile` | Optional path to write the simulated data as a CSV |
| `verbose` | `TRUE` (default) prints simulation progress |

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
