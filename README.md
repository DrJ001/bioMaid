# biomAid

[![R-CMD-check](https://github.com/DrJ001/biomAid/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/DrJ001/biomAid/actions/workflows/R-CMD-check.yml)
[![Codecov](https://app.codecov.io/gh/DrJ001/biomAid/graph/badge.svg)](https://app.codecov.io/gh/DrJ001/biomAid)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://cran.r-project.org/)

---

<img src="man/figures/biomAid_logo.png" align="right" height="150" alt="biomAid"/>

Welcome to bioMaid! This package has been specifically built to provide
biometricians with flexible functions for interpreting and further modelling of results
from complex linear mixed models fitted wth software such as **ASReml-R V4**. Watch this
space, there are a lot more functions coming. 


## Vignettes

| Vignette | Description |
|----------|-------------|
| [Wald Tests on Fixed-Effect Contrasts](https://DrJ001.github.io/biomAid/waldTest.html) | Mathematical framework and worked examples for `waldTest()` and `plot_waldTest()` |
| [Multi-treatment Random Regression](https://DrJ001.github.io/biomAid/randomRegress.html) | Conditioning schemes, efficiency/responsiveness decomposition, and all plot types for `randomRegress()` and `plot_randomRegress()` |
| [Multi-treatment Fixed-Effects Regression](https://DrJ001.github.io/biomAid/fixedRegress.html) | OLS conditioning schemes, efficiency/response index decomposition, and plot types for `fixedRegress()` and `plot_fixedRegress()` |
| [Extracting and Padding Field Trial Layouts](https://DrJ001.github.io/biomAid/padTrial.html) | Step-by-step guide to guard-row removal, missing-plot padding, and Before/After visualisation with `padTrial()` and `plot_padTrial()` |
| [Multiple Comparison Criteria](https://DrJ001.github.io/biomAid/compare.html) | HSD, LSD, and Bonferroni criteria, by-group comparisons, and all three plot types for `compare()` and `plot_compare()` |
| [BLUP Accuracy in Multi-Environment Trials](https://DrJ001.github.io/biomAid/accuracy.html) | Mrode accuracy and Cullis H², supported random structures, and all six plot types for `accuracy()` and `plot_accuracy()` |
| [Simulating Multi-Environment Trials](https://DrJ001.github.io/biomAid/simTrialData.html) | Mathematical framework, balanced/unbalanced/split-plot designs, and all four plot types for `simTrialData()` and `plot_simTrialData()` |
| [Factor Analytic Selection Tools: FAST and iClass](https://DrJ001.github.io/biomAid/fastIC.html) | Mathematical framework, FAST global metrics, iClass interaction classes, and all seven plot types for `fastIC()` and `plot_fastIC()` |

## Function reference

Full documentation for every exported function is available on the package website:

👉 **[https://DrJ001.github.io/biomAid/reference/](https://DrJ001.github.io/biomAid/reference/)**

---

## Installation

```r
# From GitHub (requires remotes or pak)
pak::pkg_install("DrJ001/biomAid")

# or
remotes::install_github("DrJ001/biomAid")
```

> **Note:** The core modelling functions require an [ASReml-R V4](https://vsni.co.uk/software/asreml-r/)
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

### `plot_compare()` — Visualise pairwise comparison results

Four plot types for exploring the output of `compare()`. All types handle
multi-factor `by`-group structure automatically through faceting, and support
an optional interactive [plotly](https://plotly.com/r/) version with hover
tooltips. The result object supports `+` for adding ggplot2 layers before
plotly conversion via `pc_add()`.

```r
plot_compare(res,
             type        = c("dotplot", "errbar", "letters", "heatmap"),
             reference   = NULL,
             interactive = FALSE,
             theme       = ggplot2::theme_bw(),
             return_data = FALSE,
             ...)
```

| Argument | Description |
|----------|-------------|
| `res` | Data frame returned by `compare()` |
| `type` | Plot type. Default `"dotplot"` |
| `reference` | `"dotplot"` only. Character name of a check variety to anchor the criterion band; `NULL` (default) anchors to the top-ranked variety |
| `interactive` | `TRUE` converts to a plotly object with hover tooltips. Requires **plotly**. Default `FALSE` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns the tidy data frame instead of the plot |

| Type | Description |
|------|-------------|
| `"dotplot"` | Varieties sorted by predicted value (highest at top). A shaded band of width equal to the criterion is drawn below the top-ranked variety — points inside the band are not significantly different from the best; points outside are red. `reference` anchors the band to a named check variety instead. |
| `"errbar"` | Each variety shown as a point with a horizontal error bar spanning ±½ criterion. Two varieties are significantly different if and only if their bars do not overlap — an exact visual equivalent of the criterion test. |
| `"letters"` | Compact letter display (CLD) overlaid on the sorted dot plot. Varieties sharing at least one letter are not significantly different. |
| `"heatmap"` | n × n pairwise absolute-difference tile matrix; white × marks significantly different pairs. Useful for large variety sets where a letter display becomes unreadable. |

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
randomRegress(model, term = "us(TSite):Variety", levs = NULL,
              type = "baseline", cond = NULL,
              sep = "-", pev = TRUE, ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | Full random-effect interaction string, e.g. `"fa(TSite, 2):Variety"` or `"corgh(TSite):vm(Variety, giv1)"`. Default `"us(TSite):Variety"` |
| `levs` | Character vector of treatment levels to decompose |
| `type` | `"baseline"`, `"sequential"`, `"partial"`, or `"custom"` |
| `cond` | User-supplied conditioning list when `type = "custom"` |
| `sep` | Separator between treatment and site in composite `TSite` labels. Default `"-"` |
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
                   cond_x      = 1L,
                   theme       = ggplot2::theme_bw(),
                   return_data = FALSE,
                   ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `randomRegress()` |
| `type` | `"regress"`, `"quadrant"`, or `"gmat"` |
| `treatments` | Character vector to restrict conditioning pairs plotted. `NULL` = all |
| `highlight` | `"default"` auto-selects archetypes by distance from origin; character vector of variety names for custom highlights; `NULL` = no highlighting |
| `centre` | `TRUE` adds back within-site means (useful for demo data). Default `FALSE` |
| `cond_x` | `"regress"` only. Positive integer selecting which member of the conditioning set $A_j$ appears on the x-axis (added variable plot). Default `1L` |
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

### `faSummary()` — Factor analytic model summary

Extracts and rotates the factor analytic variance structure from an ASReml-R V4
model, returning one entry for every `fa()` term in the random formula. Provides
the genetic covariance (`Gmat`) and correlation (`Cmat`) matrices, the rotated
loadings, specific variances, variance accounted for, and the genotype BLUPs and
factor score EBLUPs. This is the engine underlying `fastIC()` and the FA path of
`randomRegress()`.

Where a model fits both `fa(Site, k):vm(Variety, G)` and
`fa(Site, k):ide(Variety)`, an additional `"<outer>:<inner>-total"` entry holds
the summed covariance structure and BLUPs.

```r
faSummary(model, term = NULL, blups = TRUE, combine.ide = TRUE)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object with at least one `fa()` random term |
| `term` | FA term(s) to summarise, e.g. `"fa(Site, 3):Variety"`. Default `NULL` = all FA terms found |
| `blups` | Return genotype BLUPs and factor score EBLUPs. Default `TRUE` |
| `combine.ide` | Append the combined `vm()` + `ide()` "total" structure where such a pair exists. Default `TRUE` |

Each `$gammas[[term]]` element contains `Gmat`, `Cmat`, `loads`, `loads_cor`,
`spec_var`, `vaf_env`, `vaf_summary`, `vaf_total` and `k`; each
`$blups[[term]]` element contains `blups` and `scores` data frames.

---

### `plot_faSummary()` — Visualise the FA structure

Five diagnostic plots for the output of `faSummary()`.

```r
plot_faSummary(res,
               type        = c("VAF", "heatmap", "loadings",
                               "regress", "added"),
               term        = NULL,
               order       = c("loading", "asis", "cluster"),
               varieties   = "default",
               n_varieties = 6L,
               tol         = 0.85,
               theme       = ggplot2::theme_bw(),
               return_data = FALSE,
               ...)
```

| Argument | Description |
|----------|-------------|
| `res` | An object of class `faSummary` |
| `type` | Plot type. Default `"VAF"` |
| `term` | Which FA term to plot. Default `NULL` = first term carrying loadings |
| `order` | Environment ordering for `"VAF"` and `"heatmap"`: `"loading"`, `"asis"` or `"cluster"` |
| `varieties` | Varieties shown in `"regress"` / `"added"`; `"default"` picks each factor's extremes |
| `n_varieties` | Maximum number of varieties selected by `"default"`. Default `6L` |
| `tol` | Radius threshold flagging well-explained environments in `"loadings"`. Default `0.85` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns a list with `$plot` and `$data`. Default `FALSE` |

| Type | Description |
|------|-------------|
| `"VAF"` | Stacked 100% bar chart of Variance Accounted For per environment. Each bar is subdivided by FA factor (sequential blues, bottom to top) plus specific variance (grey, top). A dashed line marks the overall mean proportion explained across environments. |
| `"heatmap"` | Genetic correlation between environments on a diverging palette centred at zero. Use `order = "cluster"` to reveal block structure. |
| `"loadings"` | Correlation-scaled loadings as vectors from the origin, one panel per factor pair, with the unit circle for reference. Environments reaching beyond `tol` are drawn solid. Requires k ≥ 2. |
| `"regress"` | Genotype BLUPs against environment loadings, one panel per variety and factor, with a line whose slope is the variety's factor score. |
| `"added"` | As `"regress"`, but BLUPs are adjusted by removing every other factor's contribution — the correct display when k ≥ 2. |

---

### `fastIC()` — Factor Analytic Selection Tools

Implements the **FAST** (Smith & Cullis 2018) and **iClass** (Smith et al. 2021)
approaches for summarising variety performance from an FA mixed model fitted in
ASReml-R V4. Always computes both global FAST metrics (`global_op`,
`global_stab`) and within-class iClass metrics (`iClassOP`, `iClassRMSD`).
The genotype factor may be wrapped in `vm(...)` or `ide(...)` for
relationship-aware models.

The FA decomposition is performed by `faSummary()`. Two of its results are
carried through as attributes of the returned data frame — `vaf_env`
(per-environment variance accounted for by each factor and the specific
variance) and `vaf_summary` (overall proportions across all environments) — for
use by `plot_faSummary(type = "VAF")`.

```r
fastIC(model, term = "fa(Site, 4):Genotype",
       ic.num = 2L,
       ...)
```

| Argument | Description |
|----------|-------------|
| `model` | An ASReml-R V4 model object |
| `term` | FA model term string. Default `"fa(Site, 4):Genotype"` |
| `ic.num` | Number of factors used for iClass sign-pattern classification and iClassOP. Must be < k (i.e. 1 to k − 1) so that the kth factor remains available for iClassRMSD. Default `2` |

---

### `plot_fastIC()` — Visualise FAST and iClass results

Six plot types for exploring the output of `fastIC()`. Covers global
performance and stability, the FA factor structure, within-class metrics, and
cross-class comparisons — together providing a comprehensive view of the GEI
landscape captured by the FA model. The per-environment variance decomposition
lives in `plot_faSummary(type = "VAF")`.

```r
plot_fastIC(res,
            type        = c("fast", "biplot", "CVE",
                            "iclass", "OP.pairs", "OP.variety"),
            highlight   = "default",
            n_highlight = 3L,
            theme       = ggplot2::theme_bw(),
            return_data = FALSE,
            ...)
```

| Argument | Description |
|----------|-------------|
| `res` | Data frame returned by `fastIC()` |
| `type` | Plot type. Default `"fast"` |
| `highlight` | `"default"` auto-selects varieties (by `global_op` / instability for `"fast"`, `"biplot"`, `"CVE"`; by mean iClassOP for `"iclass"`, `"OP.pairs"`, `"OP.variety"`); character vector for explicit names; `NULL` = no annotation |
| `n_highlight` | Maximum number of varieties to highlight automatically. Default `3L` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns a list with `$plot` and `$data`. Default `FALSE` |

| Type | Description |
|------|-------------|
| `"fast"` | Scatter of global Overall Performance (`global_op`) vs global stability (`global_stab`) with quadrant annotations (Broadly adapted / Responsive / Poor & stable / Poor & unstable). |
| `"biplot"` | FA biplot with genotype score points and environment loading arrows for Factors 1 and 2. Arrows coloured by iClass when present. |
| `"CVE"` | Diverging-colour heatmap of the Common Variety Effect (genotype × environment). Environments ordered by iClass then first-factor loading; genotypes by `global_op`. |
| `"iclass"` | Scatter of within-class iClassOP vs iClassRMSD, faceted by iClass with a per-class mean-OP reference line. |
| `"OP.pairs"` | Lower-triangular pairs plot of iClassOP across all iClass levels. Uses `patchwork` when available; falls back to `facet_grid`. Requires ≥ 2 iClass levels. |
| `"OP.variety"` | Line plot of iClassOP across ordered interaction classes. Highlighted varieties drawn in colour over a grey background of all other varieties. |

---

### `simTrialData()` — Simulate trial data

Generates a balanced or unbalanced MET or split-plot dataset with a realistic
genetic covariance structure across environments. The genetic covariance matrix
`G` is either auto-generated from correlation bounds (`G = "auto"`, default)
or supplied directly by the user.

Set `treatments = NULL` for a pure MET (simple RCB, no Treatment column).
Supply `treatments = c("T0", "T1", ...)` for a split-plot design; the
genetic structure then operates over all Treatment x Site (`TSite`) combinations.

```r
simTrialData(nvar        = 20L,
             nsite       = 10L,
             treatments  = NULL,
             nrep        = 2L,
             G           = "auto",
             incidence   = "balanced",
             seed        = NULL,
             verbose     = TRUE,
             sim.options = list())
```

| Argument | Description |
|----------|-------------|
| `nvar` | Number of varieties. Default `20` |
| `nsite` | Number of sites. Default `10` |
| `treatments` | Character vector of treatment labels, or `NULL` for MET-only. Default `NULL` |
| `nrep` | Number of replicates per site. Default `2` |
| `G` | `"auto"` (default) — generate random SPD G from correlation bounds; or a user-supplied ngroup x ngroup SPD covariance matrix |
| `incidence` | `"balanced"` (default) — all varieties at every site; `"unbalanced"` — auto two-tier structure; or a user-supplied `nvar x nsite` matrix of 0/1 |
| `seed` | Random seed. `NULL` = no fixed seed |
| `verbose` | Print design summary and suggested model. Default `TRUE` |
| `sim.options` | Named list of optional controls (see below) |

Key `sim.options` elements (all have built-in defaults):

| Element | Description | Default |
|---------|-------------|---------|
| `site_mean` / `site_sd` | Grand mean and SD of site means | `4500` / `600` |
| `sigma_genetic` | Target mean genetic SD per group | `250` |
| `g_cor_min` / `g_cor_max` | Correlation bounds for auto-generated G | `0.20` / `0.90` |
| `treat_effects` | Fixed treatment effects vector (multi-treatment only) | auto-spaced |
| `error_sd` / `rep_sd` / `row_sd` / `col_sd` | Error SD components | `350`/`150`/`80`/`60` |
| `sep` | Separator for `TSite` labels | `"-"` |
| `variety_prefix` / `site_prefix` | Label prefixes | `"Var"` / `"Env"` |
| `outfile` | Optional CSV output path | `NULL` |

---

### `plot_simTrialData()` — Visualise simulated trial data

Four complementary plot types for exploring the structure of data generated
by `simTrialData()`. Covers the field layout, variety connectivity, true
genetic correlation structure, and the underlying GEI surface — together
providing a complete view of what was simulated before any model is fitted.

```r
plot_simTrialData(res,
                  type        = c("trial", "incidence", "correlation", "blup"),
                  fill        = NULL,
                  label       = NULL,
                  sites       = NULL,
                  sort        = TRUE,
                  ncol        = NULL,
                  theme       = ggplot2::theme_bw(),
                  return_data = FALSE,
                  ...)
```

| Argument | Description |
|----------|-------------|
| `res` | List returned by `simTrialData()` |
| `fill` | `"trial"` only. Column name for tile fill, or `NULL` (default → `Rep`) |
| `label` | `"trial"` only. Column name to overlay as text. `NULL` = no labels |
| `sites` | `"trial"` and `"blup"` only. Character vector of sites to include. `NULL` = all |
| `sort` | Reorder varieties in `"incidence"` (by site count) and `"blup"` (by mean BLUP). Default `TRUE` |
| `ncol` | `"trial"` only. Number of facet columns. `NULL` = ggplot2 default |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns a list with `$plot` and `$data`. Default `FALSE` |


| Type | Description |
|------|-------------|
| `"trial"` | Tiled Row × Column field layout, one panel per site. `fill` maps any data column onto tile colour (`NULL` → `Rep`; numeric → viridis; factor → discrete palette). `label` overlays text from any column. |
| `"incidence"` | Variety × Site presence/absence grid. Each variety label shows its site count `[n]`; each site header shows its variety count `[n]`. `sort = TRUE` orders varieties from most- to least-connected (top to bottom). |
| `"correlation"` | Full heatmap of `cov2cor(params$G)`. Diverging RdBu palette centred at zero; correlation values printed when the matrix is ≤ 15 × 15. |
| `"blup"` | Heatmap of `params$g_arr` (Variety × Group). Diverging RdBu palette centred at zero reveals the true GEI structure. `sort = TRUE` orders varieties by decreasing mean BLUP. |

---

### `accuracy()` — Model-based BLUP accuracy

Computes per-environment mean BLUP accuracy from a fitted ASReml-R V4
mixed model. Two complementary metrics are available:

- **Mrode accuracy** — `r = mean[sqrt(1 - PEV / G_jj)]` per variety,
  averaged over environments.
- **Generalised H²** — `1 - mean(SED²) / (2 G_jj)`, a variety-comparison
  reliability analogue from Cullis et al. (2006)

Supports `fa()`, `diag()`, `corgh()`, `corh()`, `us()`, and single-environment
`id()` random structures. The `only =` argument to `predict.asreml()` is invoked
internally to avoid fixed-effect inflation of the PEV.

```r
accuracy(model,
         term       = NULL,
         metric     = c("accuracy", "gen.H2"),
         pworkspace = "2gb",
         by_variety = FALSE)
```

| Argument | Description |
|----------|-------------|
| `model` | Fitted ASReml-R V4 model object |
| `term` | Full random-effect interaction term, e.g. `"fa(Site, 2):id(Variety)"` or `"corgh(Site):vm(Variety, giv1)"`; auto-detected from the model formula when `NULL` |
| `metric` | One or both of `"accuracy"` and `"gen.H2"` (Cullis H²). Default = both |
| `pworkspace` | Passed to `predict.asreml()`. Default `"2gb"` |
| `by_variety` | `TRUE` returns one row per variety × environment; output includes a `present` column (`TRUE` = observed, `FALSE` = unobserved) |

---

### `plot_accuracy()` — Visualise BLUP accuracy

Six plot types for visualising and comparing accuracy results.
Pass a second accuracy object via `res2` for head-to-head model comparisons.

```r
plot_accuracy(res,
              res2        = NULL,
              type        = c("lollipop", "violin", "heatmap",
                              "dumbbell", "scatter", "diff"),
              metric      = c("accuracy", "gen.H2"),
              label1      = "Model 1",
              label2      = "Model 2",
              theme       = ggplot2::theme_bw(),
              return_data = FALSE,
              ...)
```
| Argument | Description |
|----------|-------------|
| `res` | Data frame returned by `accuracy()` |
| `res2` | Optional second `accuracy()` result for two-model comparison plots (`"dumbbell"`, `"scatter"`, `"diff"`). `NULL` (default) = single-model plots only |
| `type` | Plot type. Default `"lollipop"` |
| `metric` | One or both of `"accuracy"` (Mrode) and `"gen.H2"` (Cullis H²). When both are supplied the plot is split into two facet panels. Default = both |
| `label1` | Label for `res` in two-model plots. Default `"Model 1"` |
| `label2` | Label for `res2` in two-model plots. Default `"Model 2"` |
| `theme` | A ggplot2 theme object. Default `theme_bw()` |
| `return_data` | `TRUE` returns a list with `$plot` and `$data` instead of the plot alone |

| Type | Input | Description |
|------|-------|-------------|
| `"lollipop"` | Group-level | Mean accuracy per environment with ±SD error bars |
| `"violin"` | `by_variety` | Distribution of per-variety accuracies per environment |
| `"heatmap"` | `by_variety` | Variety × environment accuracy grid (viridis fill) |
| `"dumbbell"` | Group-level + `res2` | Paired dots per environment; segment colour = improvement direction |
| `"scatter"` | `by_variety` + `res2` | Per-variety model1 vs model2 accuracy scatter |
| `"diff"` | Group-level + `res2` | Signed accuracy gain (model2 − model1) per environment |

---

## License

MIT © J Taylor
