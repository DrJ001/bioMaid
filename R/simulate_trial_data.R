#' Simulate Multi-Treatment Multi-Environment Plant Breeding Trial Data
#'
#' @description
#' Generates a balanced split-plot field trial dataset suitable for fitting
#' ASReml-R mixed models and testing [randomRegress()].
#'
#' **Experimental design** (per environment):
#' \itemize{
#'   \item Rectangular layout of `nrep * rows_per_strip` rows ×
#'         `length(treatments) * cols_per_strip` columns.
#'   \item `nrep` replicate blocks, each occupying `rows_per_strip` consecutive
#'         rows.
#'   \item Within each replicate: one treatment strip per treatment
#'         (`rows_per_strip × cols_per_strip = nvar` sub-plots), with
#'         treatments randomised to strips independently in every replicate.
#'   \item Varieties randomised to sub-plots within each strip.
#' }
#'
#' **Genetic simulation** uses the same conditional multivariate-normal
#' structure that [randomRegress()] is designed to recover:
#' \deqn{u_1 \sim \mathcal{N}(0,\,\sigma_1^2)}
#' \deqn{u_j \mid u_1 \sim \mathcal{N}(\beta_{j,s}\,u_1,\,\sigma_{r_j}^2)
#'       \quad j = 2,\ldots,k}
#' where \eqn{\beta_{j,s}} varies by site to create genuine
#' Treatment × Genotype × Environment interaction.  The true \eqn{\beta}
#' matrix is returned alongside the data for ground-truth comparison with
#' model estimates.
#'
#' @param nvar        Integer. Number of varieties. Should be composite (not
#'   prime) for a tidy rectangular sub-plot layout; a warning is issued
#'   otherwise.  Default `20`.
#' @param nsite       Integer. Number of environments (sites). Default `10`.
#' @param treatments  Character vector of length ≥ 2 giving treatment labels.
#'   The **first element** is the baseline; all others are regressed onto it
#'   in [randomRegress()].  Default `c("N0", "N1", "N2")`.
#' @param nrep        Integer. Number of complete replicates per environment.
#'   Default `3`.
#' @param rows_per_strip Integer or `NULL`. Rows allocated to each treatment
#'   strip within a replicate block. With `cols_per_strip` must satisfy
#'   `rows_per_strip * cols_per_strip == nvar`. Auto-calculated to the
#'   near-square factor pair when `NULL`.
#' @param cols_per_strip Integer or `NULL`. Columns per treatment strip.
#'   Auto-calculated when `NULL` (see `rows_per_strip`).
#' @param site_mean   Numeric. Grand mean yield across all sites (same units
#'   as the simulated trait). Default `4500`.
#' @param site_sd     Numeric. Standard deviation of site mean yields.
#'   Default `600`.
#' @param treat_effects Numeric vector of length `length(treatments)` giving
#'   the fixed treatment main effects, with `treat_effects[1] = 0` for
#'   identifiability.  Auto-spaced from 0 to `0.25 * site_mean` when `NULL`.
#' @param sigma_base  Numeric. Genetic standard deviation for the baseline
#'   treatment BLUP. Default `250`.
#' @param sigr        Numeric vector of length `length(treatments) - 1`.
#'   Conditional (responsiveness) standard deviations for each non-baseline
#'   treatment given the baseline. Auto-scaled to
#'   `sigma_base * 0.4 * seq_len(ntreat - 1)` when `NULL`.
#' @param beta_min    Numeric vector of length `length(treatments) - 1`.
#'   Lower bound of the uniform distribution from which site-specific
#'   regression coefficients are drawn for each non-baseline treatment.
#'   Auto-calculated when `NULL`.
#' @param beta_max    Numeric vector of length `length(treatments) - 1`.
#'   Upper bound (see `beta_min`). Auto-calculated when `NULL`.
#' @param rep_sd      Numeric. Standard deviation of replicate (block) effects.
#'   Default `150`.
#' @param row_sd      Numeric. Standard deviation of row spatial effects.
#'   Default `80`.
#' @param col_sd      Numeric. Standard deviation of column spatial effects.
#'   Default `60`.
#' @param error_sd    Numeric. Standard deviation of the plot-level residual
#'   (within-sub-plot error). Default `350`.
#' @param sep         Character. Separator used when constructing the composite
#'   `TSite` factor (e.g. `"N0-Env01"` with `sep = "-"`). Must not appear
#'   inside any treatment or site label. Default `"-"`.
#' @param variety_prefix Character. Prefix for variety labels. Default `"Var"`.
#' @param site_prefix    Character. Prefix for site labels. Default `"Env"`.
#' @param seed        Integer or `NULL`. Random seed passed to [set.seed()].
#'   `NULL` uses the current RNG state. Default `NULL`.
#' @param outfile     Character path or `NULL`. If non-`NULL` the data frame
#'   is written to this CSV path via [write.csv()]. Default `NULL`.
#' @param verbose     Logical. Print a design summary and suggested ASReml-R
#'   model call. Default `TRUE`.
#'
#' @return A named list:
#' \describe{
#'   \item{`data`}{Data frame with columns `Site`, `Treatment`, `TSite`,
#'     `Variety`, `Rep`, `Row`, `Column`, `yield`. Ordered by `Site`, `Row`,
#'     `Column` for compatibility with ASReml-R spatial models.}
#'   \item{`params`}{Named list of true simulation parameters for ground-truth
#'     comparison with [randomRegress()] estimates:
#'     \describe{
#'       \item{`beta`}{`nsite × (ntreat-1)` matrix of true site-specific
#'         regression coefficients.}
#'       \item{`sigr`}{Vector of true conditional standard deviations.}
#'       \item{`treat_effects`}{True treatment fixed effects used.}
#'       \item{`site_means`}{True site mean yields used.}
#'       \item{`g_base`}{Standardised variety genetic potentials (length
#'         `nvar`).}
#'       \item{`g_arr`}{`nvar × ntreat × nsite` array of true genetic values.}
#'     }}
#' }
#'
#' @examples
#' \dontrun{
#' # Reproduce the original 3-treatment, 10-site dataset
#' out <- simTrialData(seed = 2024, outfile = "plant_trial_data.csv")
#' head(out$data)
#' out$params$beta   # compare with randomRegress()$beta
#'
#' # Two-treatment, 5-site dataset
#' out2 <- simTrialData(nsite = 5, treatments = c("Dry", "Irr"), seed = 1)
#'
#' # Four treatments, 12 varieties, 4 reps
#' out4 <- simTrialData(nvar = 12, nsite = 8,
#'                      treatments = c("T0","T1","T2","T3"), nrep = 4,
#'                      seed = 99)
#' }
#'
#' @export
simTrialData <- function(nvar           = 20L,
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
                         verbose        = TRUE) {

  # ---- Input validation --------------------------------------------------
  ntreat <- length(treatments)
  if (ntreat < 2L)
    stop("'treatments' must have at least 2 elements.")
  if (any(grepl(sep, treatments, fixed = TRUE)))
    stop("'sep' (\"", sep, "\") must not appear inside any treatment label.")
  if (!is.null(seed))
    set.seed(seed)

  nresp <- ntreat - 1L   # number of non-baseline (responsiveness) treatments

  # ---- Auto layout: near-square factor pair for nvar sub-plots -----------
  .best_dims <- function(n) {
    sq  <- floor(sqrt(n))
    fac <- which(n %% seq_len(sq) == 0L)
    if (!length(fac)) {
      warning("nvar = ", n, " is prime; strip layout will be 1 x ", n,
              ". Consider choosing a composite nvar for a tidier design.")
      return(c(rows = 1L, cols = as.integer(n)))
    }
    r <- max(fac)
    c(rows = as.integer(r), cols = as.integer(n %/% r))
  }

  if (is.null(rows_per_strip) || is.null(cols_per_strip)) {
    dims           <- .best_dims(nvar)
    rows_per_strip <- dims["rows"]
    cols_per_strip <- dims["cols"]
  } else {
    if (rows_per_strip * cols_per_strip != nvar)
      stop("rows_per_strip * cols_per_strip must equal nvar (", nvar, ").")
  }

  nrow_e <- nrep   * rows_per_strip   # total rows per environment
  ncol_e <- ntreat * cols_per_strip   # total columns per environment

  # ---- Default parameter fills -------------------------------------------
  if (is.null(treat_effects))
    treat_effects <- seq(0, site_mean * 0.25, length.out = ntreat)

  if (length(treat_effects) != ntreat)
    stop("'treat_effects' must have length equal to length(treatments).")

  if (is.null(sigr))
    sigr <- sigma_base * 0.4 * seq_len(nresp)

  if (length(sigr) != nresp)
    stop("'sigr' must have length ", nresp, " (one per non-baseline treatment).")

  # Default beta ranges: each subsequent treatment has a higher expected beta
  # beta_j in [0.65 + 0.35*(j-1),  1.30 + 0.65*(j-1)]
  if (is.null(beta_min))
    beta_min <- 0.65 + 0.35 * (seq_len(nresp) - 1L)
  if (is.null(beta_max))
    beta_max <- 1.30 + 0.65 * (seq_len(nresp) - 1L)

  if (length(beta_min) != nresp || length(beta_max) != nresp)
    stop("'beta_min' and 'beta_max' must each have length ", nresp, ".")
  if (any(beta_min >= beta_max))
    stop("All elements of 'beta_min' must be strictly less than 'beta_max'.")

  # ---- Labels ------------------------------------------------------------
  ndig_var  <- nchar(as.character(nvar))
  ndig_site <- nchar(as.character(nsite))
  varieties  <- paste0(variety_prefix, sprintf(paste0("%0", ndig_var,  "d"), seq_len(nvar)))
  sites      <- paste0(site_prefix,    sprintf(paste0("%0", ndig_site, "d"), seq_len(nsite)))

  # ---- Genetic simulation ------------------------------------------------
  # Standardised variety genetic potential (common across all treatments/sites)
  g_base <- rnorm(nvar, 0, 1)

  # Site-specific regression coefficients: beta_js ~ Uniform(min_j, max_j)
  beta_mat <- matrix(NA_real_, nrow = nsite, ncol = nresp,
                     dimnames = list(sites, treatments[-1L]))
  for (j in seq_len(nresp))
    beta_mat[, j] <- runif(nsite, min = beta_min[j], max = beta_max[j])

  # Genetic values: nvar x ntreat x nsite array
  g_arr <- array(NA_real_,
                 dim      = c(nvar, ntreat, nsite),
                 dimnames = list(varieties, treatments, sites))

  for (s in seq_len(nsite)) {
    u_base <- sigma_base * g_base                      # baseline BLUPs
    g_arr[, 1L, s] <- u_base
    for (j in seq_len(nresp))
      g_arr[, j + 1L, s] <- beta_mat[s, j] * u_base + rnorm(nvar, 0, sigr[j])
  }

  # ---- Build plot-level observations -------------------------------------
  site_means <- rnorm(nsite, mean = site_mean, sd = site_sd)
  plots      <- vector("list", nsite * nrep * ntreat * nvar)
  idx        <- 1L

  for (s in seq_len(nsite)) {

    row_eff <- rnorm(nrow_e, 0, row_sd)
    col_eff <- rnorm(ncol_e, 0, col_sd)
    rep_eff <- rnorm(nrep,   0, rep_sd)

    for (r in seq_len(nrep)) {

      treat_order <- sample(ntreat)   # randomise treatments to strips

      for (t_pos in seq_len(ntreat)) {

        t <- treat_order[t_pos]

        rows_in_rep   <- ((r      - 1L) * rows_per_strip + 1L) : (r      * rows_per_strip)
        cols_in_strip <- ((t_pos  - 1L) * cols_per_strip + 1L) : (t_pos  * cols_per_strip)

        sub_grid <- expand.grid(Row = rows_in_rep, Column = cols_in_strip)
        sub_grid <- sub_grid[order(sub_grid$Row, sub_grid$Column), ]

        var_perm <- sample(nvar)   # randomise varieties to sub-plots

        for (p in seq_len(nvar)) {
          v   <- var_perm[p]
          row <- sub_grid$Row[p]
          col <- sub_grid$Column[p]

          yld <- site_means[s]       +   # environment mean
                 treat_effects[t]    +   # treatment fixed effect
                 g_arr[v, t, s]      +   # GxTxE genetic value
                 rep_eff[r]          +   # block effect
                 row_eff[row]        +   # row spatial effect
                 col_eff[col]        +   # column spatial effect
                 rnorm(1L, 0, error_sd)  # plot residual

          plots[[idx]] <- data.frame(
            Site      = sites[s],
            Treatment = treatments[t],
            TSite     = paste(treatments[t], sites[s], sep = sep),
            Variety   = varieties[v],
            Rep       = r,
            Row       = row,
            Column    = col,
            yield     = round(yld, 1L),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
  }

  dat <- do.call(rbind, plots)

  # ---- Factor levels and ordering ----------------------------------------
  tsite_levs    <- sort(unique(dat$TSite))
  dat$Site      <- factor(dat$Site,      levels = sites)
  dat$Treatment <- factor(dat$Treatment, levels = treatments)
  dat$TSite     <- factor(dat$TSite,     levels = tsite_levs)
  dat$Variety   <- factor(dat$Variety,   levels = varieties)
  dat$Rep       <- factor(dat$Rep,       levels = as.character(seq_len(nrep)))

  # Order by Site, Row, Column -- required for ASReml-R spatial residual models
  dat <- dat[order(dat$Site, dat$Row, dat$Column), ]
  rownames(dat) <- NULL

  # ---- Optional CSV output -----------------------------------------------
  if (!is.null(outfile)) {
    write.csv(dat, outfile, row.names = FALSE)
    if (verbose) cat("CSV written to:", outfile, "\n")
  }

  # ---- Verbose summary ---------------------------------------------------
  if (verbose) {
    cat("\n=== simTrialData summary ===\n")
    cat(sprintf("  Environments  : %d\n", nsite))
    cat(sprintf("  Treatments    : %d  ->  %s\n", ntreat, paste(treatments, collapse = ", ")))
    cat(sprintf("  Varieties     : %d\n", nvar))
    cat(sprintf("  Replicates    : %d\n", nrep))
    cat(sprintf("  Layout        : %d rows x %d cols per environment  (%d plots)\n",
                nrow_e, ncol_e, nrow_e * ncol_e))
    cat(sprintf("  Strip size    : %d rows x %d cols = %d sub-plots\n",
                rows_per_strip, cols_per_strip, rows_per_strip * cols_per_strip))
    cat(sprintf("  Total plots   : %d\n", nrow(dat)))
    cat("\n  True beta matrix (site-specific regression coefficients):\n")
    print(round(beta_mat, 3))
    cat("\n  True conditional SDs (sigr):\n")
    print(setNames(round(sigr, 1), treatments[-1L]))
    cat("\n  True treatment fixed effects:\n")
    print(setNames(round(treat_effects, 1), treatments))
    cat("\n  Yield summary by Treatment:\n")
    print(tapply(dat$yield, dat$Treatment, function(x)
      round(c(Min = min(x), Mean = mean(x), Max = max(x), SD = sd(x)), 1)))
    # --- Context-aware variance structure recommendation ------------------
    # us(TSite) has ntsite*(ntsite+1)/2 parameters; this becomes infeasible
    # quickly.  Use corgh() for moderate sizes and fa() for large ones.
    ntsite     <- ntreat * nsite
    n_us_pars  <- ntsite * (ntsite + 1L) / 2L
    var_struct <- if (ntsite <= 8L) {
      "us(TSite)"
    } else if (ntsite <= 16L) {
      "corgh(TSite)"
    } else {
      k <- min(3L, max(1L, as.integer(floor(sqrt(ntsite * 0.15)))))
      paste0("fa(TSite, ", k, ")")
    }
    n_var_pars <- if (ntsite <= 8L) {
      n_us_pars
    } else if (ntsite <= 16L) {
      2L * ntsite - 1L                            # ntsite variances + ntsite-1 correlations
    } else {
      k <- min(3L, max(1L, as.integer(floor(sqrt(ntsite * 0.15)))))
      ntsite * k + ntsite                         # k loadings + ntsite specific variances
    }

    cat("\n  Suggested ASReml-R model:\n")
    cat(sprintf("    ## TSite has %d levels -> %s  (%d variance params vs %d for us())\n",
                ntsite, var_struct, n_var_pars, n_us_pars))
    if (ntsite > 16L) {
      k <- min(3L, max(1L, as.integer(floor(sqrt(ntsite * 0.15)))))
      cat(sprintf("    ## Build FA order up: start fa(TSite,1), then fa(TSite,2) ... fa(TSite,%d)\n", k))
      cat(  "    ## Use AIC/BIC or LRT to choose final order.\n")
    }
    cat("    asreml(\n")
    cat("      fixed    = yield ~ TSite,\n")
    cat(sprintf("      random   = ~ %s:Variety + at(Site):Rep +\n", var_struct))
    cat("                   at(Site):Rep:Treatment,\n")
    cat("      residual = ~ dsum(~ ar1(Row):ar1(Column) | Site),\n")
    cat("      data     = dat\n")
    cat("    )\n")
    cat(sprintf("\n    randomRegress(model, Env = \"TSite:Variety\",\n"))
    cat(sprintf("                    levs = c(%s), sep = \"%s\")\n",
                paste(sprintf('"%s"', treatments), collapse = ", "), sep))
    cat("============================\n\n")
  }

  # ---- Return ------------------------------------------------------------
  list(
    data   = dat,
    params = list(
      beta          = beta_mat,
      sigr          = setNames(sigr, treatments[-1L]),
      treat_effects = setNames(treat_effects, treatments),
      site_means    = setNames(site_means, sites),
      g_base        = setNames(g_base, varieties),
      g_arr         = g_arr
    )
  )
}


# ---- Reproduce the original plant_trial_data.csv -------------------------
if (FALSE) {
  out <- simTrialData(
    nvar       = 20,
    nsite      = 10,
    treatments = c("N0", "N1", "N2"),
    nrep       = 3,
    site_mean  = 4500,
    site_sd    = 600,
    sigma_base = 250,
    error_sd   = 350,
    seed       = 2024,
    outfile    = "plant_trial_data.csv"
  )

  # Quick sanity check: treatment means should increase N0 -> N1 -> N2
  print(tapply(out$data$yield, out$data$Treatment, mean))

  # True betas to compare against randomRegress() estimates
  print(out$params$beta)
}
