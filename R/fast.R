#' Factor Analytic Selection Tools: FAST and iClass Analysis
#'
#' @description
#' A unified implementation of the **FAST** (Factor Analytic Selection Tools;
#' Smith & Cullis 2018) and **iClass** (interaction class; Smith et al. 2021)
#' approaches for summarising variety performance from a Factor Analytic Linear
#' Mixed Model fitted in ASReml-R V4.
#'
#' Both methods build on the rotated FA loadings \eqn{\hat{\bm{\Lambda}}} and
#' score EBLUPs \eqn{\hat{\bm{f}}} returned by [faSummary()].  The
#' **Common Variety Effect** (CVE) for genotype \eqn{g} in environment \eqn{j}
#' is the FA regression prediction:
#'
#' \deqn{
#'   \widehat{\text{CVE}}(g,j) = \sum_{r=1}^{k} \hat{\lambda}_{rj}\,\hat{f}_{rg}
#' }
#'
#' The **Variance Accounted For** (VAF) by each factor in environment \eqn{j}
#' is \eqn{\hat{\lambda}_{rj}^2 / (\sum_r \hat{\lambda}_{rj}^2 + \hat{\psi}_j)},
#' and the specific variance accounts for the remainder; see [plot_faSummary()]
#' with `type = "VAF"` for a per-environment visualisation.
#'
#' @section Unified framework:
#' FAST and iClass are treated as a single framework.  FAST global metrics
#' are computed unconditionally across all environments; iClass metrics
#' partition environments by the sign pattern of the first `ic.num` rotated
#' loadings and compute within-class analogues.  When `ic.num = 1` and all
#' first-factor loadings are positive (the typical outcome after rotation),
#' there is exactly one iClass and `iClassOP` equals `global_op` while
#' `iClassRMSD` equals `global_stab` — the two frameworks agree.  When
#' `ic.num >= 2` the environments split into multiple classes and the iClass
#' metrics diverge from the global ones, revealing crossover GEI that FAST
#' alone cannot detect.
#'
#' @section Global FAST metrics:
#' After rotation the first factor captures the dominant non-crossover
#' genotype-environment interaction (GEI) pattern and typically has
#' all-positive loadings.  The global metrics summarise each genotype across
#' **all** environments:
#'
#' \describe{
#'   \item{**global_op**}{
#'     \eqn{\text{global\_op}(g) = \bar{\lambda}_1 \cdot \hat{f}_{1g}}, where
#'     \eqn{\bar{\lambda}_1 = t^{-1}\sum_j\hat{\lambda}_{1j}}.  On the same
#'     scale as the trait; the genotype's expected performance at the average
#'     environment.}
#'   \item{**global_dev**}{
#'     \eqn{\text{global\_dev}(g,j) = \widehat{\text{CVE}}(g,j) -
#'     \hat{\lambda}_{1j}\hat{f}_{1g}}: residual from the first-factor
#'     regression.  Only present when \eqn{k > 1}.}
#'   \item{**global_stab**}{
#'     \eqn{\text{global\_stab}(g) = \sqrt{t^{-1}\sum_j
#'     \text{global\_dev}(g,j)^2}}: RMSD of the first-factor residual across
#'     all environments.  A small value indicates broad adaptation.
#'     Only present when \eqn{k > 1}.}
#' }
#'
#' @section iClass metrics:
#' Environments are grouped by the sign pattern of their first `ic.num`
#' rotated loadings:
#'
#' \deqn{
#'   \text{iClass}(j) = \bigotimes_{r=1}^{\texttt{ic.num}}
#'   \begin{cases} \text{p} & \hat{\lambda}_{rj} \ge 0 \\
#'                 \text{n} & \hat{\lambda}_{rj} <   0 \end{cases}
#' }
#'
#' Within each iClass \eqn{\omega}:
#'
#' \describe{
#'   \item{**iClassOP**}{
#'     \eqn{\text{iClassOP}(g,\omega) = \sum_{r=1}^{\texttt{ic.num}}
#'     \bar{\lambda}_{r\omega} \cdot \hat{f}_{rg}}, where
#'     \eqn{\bar{\lambda}_{r\omega}} is the mean of factor \eqn{r}'s loadings
#'     across environments in \eqn{\omega}.}
#'   \item{**iClassRMSD**}{
#'     \eqn{\text{iClassRMSD}(g,\omega) = \sqrt{|\omega|^{-1}
#'     \sum_{j\in\omega} \text{dev}_{ic}(g,j)^2}}, where
#'     \eqn{\text{dev}_{ic}(g,j) = \widehat{\text{CVE}}(g,j) -
#'     \sum_{r=1}^{\texttt{ic.num}} \hat{\lambda}_{rj}\hat{f}_{rg}}.
#'     The kth factor is reserved for this residual; therefore
#'     \code{ic.num} must be strictly less than \eqn{k}.}
#' }
#'
#' @param model   An ASReml-R V4 model object containing a Factor Analytic
#'   random term.
#' @param term    Character string identifying the FA random term, written as
#'   `"fa(<EnvFactor>, k):<GenotypeFactor>"`.  The genotype factor may be
#'   wrapped in `vm(...)` (pedigree/genomic models) or `ide(...)` (identity
#'   relationship matrix) for relationship-aware models.
#'   Default `"fa(Site, 4):Genotype"`.
#' @param ic.num  Integer.  Number of factors used to form iClasses and compute
#'   `iClassOP`.  Must be strictly less than \eqn{k} (i.e. between 1 and
#'   \eqn{k - 1}) so that the kth factor remains available for `iClassRMSD`.
#'   Default `2`.
#' @param ...     Additional arguments forwarded to [faSummary()].
#'
#' @return A data frame with one row per environment \eqn{\times} genotype
#'   combination, sorted by `iclass` then environment then genotype, containing:
#'   \describe{
#'     \item{`<EnvFactor>`}{Environment labels.}
#'     \item{`<GenotypeFactor>`}{Genotype labels.}
#'     \item{`loads1`, ..., `loadsK`}{Rotated FA loadings per environment.}
#'     \item{`spec.var`}{Specific (residual) genetic variance per environment.}
#'     \item{`score1`, ..., `scoreK`}{Rotated FA score EBLUPs per genotype.}
#'     \item{`fitted1`, ..., `fittedK`}{Per-factor contributions to CVE:
#'       \eqn{\hat{\lambda}_{rj}\hat{f}_{rg}}.}
#'     \item{`CVE`}{Common Variety Effect (sum of all fitted values).}
#'     \item{`global_op`}{Global Overall Performance (FAST): repeated for
#'       every environment row of a genotype.}
#'     \item{`global_dev`}{Residual from first-factor regression (FAST):
#'       \eqn{\widehat{\text{CVE}} - \hat{\lambda}_{1j}\hat{f}_{1g}}.
#'       Only present when \eqn{k > 1}.}
#'     \item{`global_stab`}{RMSD of `global_dev` across all environments
#'       (FAST stability): repeated for every environment row of a genotype.
#'       Only present when \eqn{k > 1}.}
#'     \item{`iclass`}{Sign-pattern iClass label for the environment.}
#'     \item{`iClassOP`}{Within-iClass Overall Performance for the genotype.}
#'     \item{`iClassRMSD`}{Within-iClass RMSD for the genotype.}
#'   }
#'
#'   Two attributes are attached to the returned data frame for use by
#'   [plot_fastIC()] with `type = "VAF"`:
#'   \describe{
#'     \item{`vaf_env`}{Data frame with one row per environment containing the
#'       proportion of genetic variance accounted for by each factor
#'       (`Factor1`, ..., `FactorK`) and the specific variance (`Specific`),
#'       plus the total genetic variance (`total_var`) for that environment.}
#'     \item{`vaf_summary`}{Data frame with one row per source (`Factor 1`, ...,
#'       `Factor K`, `Specific`) containing the overall proportion of total
#'       genetic variance (`pct_var`) and cumulative proportion (`cum_pct`).}
#'   }
#'
#' @references
#' Smith, A.B. & Cullis, B.R. (2018). Plant breeding selection tools built on
#' factor analytic mixed models for multi-environment trial data.
#' *Euphytica*, 214, 143.
#'
#' Smith, A.B., Norman, A., Kuchel, H. & Cullis, B.R. (2021). Plant variety
#' selection using interaction classes derived from factor analytic linear
#' mixed models: models with independent variety effects.
#' *Frontiers in Plant Science*, 12, 737462.
#'
#' @seealso [faSummary()], [plot_faSummary()], [plot_fastIC()]
#' @export
fastIC <- function(model, term = "fa(Site, 4):Genotype",
                   ic.num = 2L,
                   ...) {

  # ---- Parse the term string -------------------------------------------
  parts  <- strsplit(term, ":")[[1L]]
  fa_idx <- grep("^fa\\s*\\(", parts)
  if (!length(fa_idx))
    stop("'term' must contain a 'fa(...)' component, e.g. \"fa(Site,4):Genotype\".")

  fa_part  <- parts[fa_idx]
  gen_part <- parts[-fa_idx]

  # Environment factor name: extract first argument of fa(...)
  sterm <- trimws(sub(",.*", "", sub("^fa\\s*\\(", "", fa_part)))

  # Genotype factor name: strip vm(...) or ide(...) wrapper if present
  if (any(grepl("^vm\\s*\\(", gen_part))) {
    vm_part <- gen_part[grep("^vm\\s*\\(", gen_part)]
    gterm   <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", vm_part))))
  } else if (any(grepl("^ide\\s*\\(", gen_part))) {
    ide_part <- gen_part[grep("^ide\\s*\\(", gen_part)]
    gterm    <- trimws(sub("[,)].*", "", sub("^ide\\s*\\(", "", ide_part)))
  } else {
    gterm <- trimws(gen_part)
  }

  # ---- FA summary --------------------------------------------------------
  sfa <- faSummary(model, term = term, ...)

  # ---- Scores: pivot long -> wide  (m genotypes × k factors) -----------
  sc_long    <- sfa$blups[[term]]$scores
  score_list <- tapply(sc_long$blupr, sc_long[[sterm]], function(el) el)
  score_mat  <- do.call(cbind, score_list)    # m × k
  k  <- ncol(score_mat)
  m  <- nrow(score_mat)
  colnames(score_mat) <- paste0("score", seq_len(k))

  # Genotype labels, in the order they appear within each Comp block.  Taken
  # from the score EBLUPs themselves so that fastIC() does not depend on the
  # original data object still being in scope.
  genotypes <- if (!is.null(rownames(score_mat))) {
    rownames(score_mat)
  } else if (!is.null(sc_long[[gterm]])) {
    unique(as.character(sc_long[[gterm]]))
  } else {
    levels(eval(model$call$data)[[gterm]])
  }

  # ---- Loadings and specific variances: t × k --------------------------
  loads_mat <- as.matrix(sfa$gammas[[term]]$loads)            # t × k
  spec_var  <- sfa$gammas[[term]]$spec_var                    # length t
  envs      <- rownames(loads_mat)
  t_envs    <- length(envs)
  colnames(loads_mat) <- paste0("loads", seq_len(k))

  # ---- Validate ic.num --------------------------------------------------
  ic.num <- as.integer(ic.num)

  if (k == 1L)
    warning("Only one factor -- iClass with ic.num = 1 places all environments ",
            "in one class.", call. = FALSE)

  if (ic.num < 1L || ic.num >= k)
    stop("'ic.num' must be between 1 and k - 1 = ", k - 1L,
         " (the kth factor must remain available for iClassRMSD).")

  # ---- CVE via matrix multiplication -----------------------------------
  CVE_mat <- score_mat %*% t(loads_mat)        # m × t

  # ---- VAF (Variance Accounted For) ------------------------------------
  # Computed by faSummary(), which owns the FA variance decomposition.
  # Carried through as attributes for plot_faSummary(type = "VAF").
  vaf_env_df  <- sfa$gammas[[term]]$vaf_env
  vaf_summary <- sfa$gammas[[term]]$vaf_summary

  # ---- Build base long-format data frame (environment-major order) -----
  env_rep  <- rep(seq_len(t_envs), each = m)
  geno_rep <- rep(seq_len(m),      times = t_envs)

  out <- data.frame(
    setNames(list(factor(envs[env_rep],       levels = envs)),      sterm),
    setNames(list(factor(genotypes[geno_rep], levels = genotypes)), gterm),
    loads_mat[env_rep,  , drop = FALSE],
    spec.var = spec_var[env_rep],
    score_mat[geno_rep, , drop = FALSE],
    CVE = as.vector(CVE_mat),
    stringsAsFactors = FALSE
  )

  # Per-factor fitted values
  for (r in seq_len(k))
    out[[paste0("fitted", r)]] <-
      out[[paste0("score", r)]] * out[[paste0("loads", r)]]

  # ====================================================================
  # Global FAST metrics  (always computed)
  # ====================================================================

  # global_op = mean(loads1) * score1
  out$global_op <- mean(loads_mat[, 1L]) * out$score1

  if (k > 1L) {
    # global_dev = CVE - fitted1
    fitted1_mat   <- outer(score_mat[, 1L], loads_mat[, 1L])  # m × t
    dev_mat       <- CVE_mat - fitted1_mat
    out$global_dev  <- as.vector(dev_mat)

    # global_stab = RMSD of global_dev across all environments
    out$global_stab <- sqrt(rowMeans(dev_mat^2))[geno_rep]
  }

  # ====================================================================
  # iClass metrics  (always computed)
  # ====================================================================

  # Pre-compute sum of fitted values for factors 1..ic.num
  fitted_ic_mat <- Reduce("+", lapply(seq_len(ic.num), function(r)
    outer(score_mat[, r], loads_mat[, r])))   # m × t

  # Classify environments by sign pattern of first ic.num loadings
  sign_str    <- apply(loads_mat[, seq_len(ic.num), drop = FALSE], 1L,
                       function(el) paste(ifelse(el >= 0, "p", "n"), collapse = ""))
  iclass_levs <- unique(sign_str)
  out$iclass  <- factor(sign_str[as.integer(out[[sterm]])], levels = iclass_levs)

  # Per-(genotype, iClass) OP and RMSD
  n_ic           <- length(iclass_levs)
  iClassOP_mat   <- matrix(NA_real_, m, n_ic,
                            dimnames = list(genotypes, iclass_levs))
  iClassRMSD_mat <- matrix(NA_real_, m, n_ic,
                            dimnames = list(genotypes, iclass_levs))

  for (w in iclass_levs) {
    env_w <- which(sign_str == w)

    mld_w             <- colMeans(loads_mat[env_w, seq_len(ic.num), drop = FALSE])
    iClassOP_mat[, w] <- score_mat[, seq_len(ic.num), drop = FALSE] %*% mld_w

    dev_ic_w              <- CVE_mat[, env_w, drop = FALSE] -
                             fitted_ic_mat[, env_w, drop = FALSE]
    iClassRMSD_mat[, w]   <- sqrt(rowMeans(dev_ic_w^2))
  }

  gi             <- cbind(as.character(out[[gterm]]), as.character(out$iclass))
  out$iClassOP   <- iClassOP_mat[gi]
  out$iClassRMSD <- iClassRMSD_mat[gi]

  # ---- Sort and attach VAF attributes -----------------------------------
  out <- out[do.call(order, out[, c("iclass", sterm, gterm)]), ]
  rownames(out) <- NULL

  attr(out, "vaf_env")     <- vaf_env_df
  attr(out, "vaf_summary") <- vaf_summary

  out
}
