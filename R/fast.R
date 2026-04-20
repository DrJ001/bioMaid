#' Factor Analytic Selection Tools: FAST and iClass Analysis
#'
#' @description
#' A unified implementation of the **FAST** (Factor Analytic Selection Tools;
#' Smith & Cullis 2018) and **iClass** (interaction class; Smith et al. 2021)
#' approaches for summarising variety performance from a Factor Analytic Linear
#' Mixed Model fitted in ASReml-R V4.
#'
#' Both methods build on the rotated FA loadings \eqn{\hat{\bm{\Lambda}}} and
#' score EBLUPs \eqn{\hat{\bm{f}}} returned by [ASExtras4::fa.asreml()].  The
#' **Common Variety Effect** (CVE) for genotype \eqn{g} in environment \eqn{j}
#' is the FA regression prediction:
#'
#' \deqn{
#'   \widehat{\text{CVE}}(g,j) = \sum_{r=1}^{k} \hat{\lambda}_{rj}\,\hat{f}_{rg}
#' }
#'
#' and the **Variety Effect** (VE) adds the environment-specific variance:
#' \eqn{\widehat{\text{VE}}(g,j) = \widehat{\text{CVE}}(g,j) + \hat{\psi}_j}.
#'
#' @section FAST (`type = "FAST"`):
#' After rotation, the first factor captures the dominant non-crossover
#' genotype-environment interaction (GEI) pattern and typically has all-positive
#' loadings.  FAST summarises each genotype by:
#'
#' \describe{
#'   \item{**Overall Performance (OP)**}{
#'     \eqn{\text{OP}(g) = \bar{\lambda}_1 \cdot \hat{f}_{1g}}, where
#'     \eqn{\bar{\lambda}_1 = t^{-1}\sum_j\hat{\lambda}_{1j}} is the mean
#'     first-factor loading.  This is on the same scale as the trait and
#'     represents the genotype's expected performance at the average
#'     environment.}
#'   \item{**Stability (stab)**}{
#'     \eqn{\text{stab}(g) = \sqrt{t^{-1}\sum_j \text{dev}(g,j)^2}}, where
#'     \eqn{\text{dev}(g,j) = \widehat{\text{CVE}}(g,j) - \hat{\lambda}_{1j}
#'     \hat{f}_{1g}} is the residual from the first-factor regression
#'     (i.e.\ the combined contribution of all higher-order bipolar factors).
#'     A small RMSD indicates broad adaptation.}
#' }
#'
#' @section iClass (`type = "iClass"`):
#' When non-trivial crossover GEI is present, a single global OP is
#' misleading.  iClass resolves this by grouping environments into
#' **interaction classes** based on the sign pattern of their first
#' `ic.num` rotated loadings:
#'
#' \deqn{
#'   \text{iClass}(j) = \bigwedge_{r=1}^{k}
#'   \begin{cases} \text{p} & \hat{\lambda}_{rj} \ge 0 \\
#'                 \text{n} & \hat{\lambda}_{rj} <   0 \end{cases}
#' }
#'
#' Within each iClass \eqn{\omega}:
#'
#' \describe{
#'   \item{**iClassOP**}{
#'     \eqn{\text{iClassOP}(g,\omega) = \sum_{r=1}^{k}
#'     \bar{\lambda}_{r\omega} \cdot \hat{f}_{rg}},
#'     where \eqn{\bar{\lambda}_{r\omega}} is the mean of factor \eqn{r}'s
#'     loadings across environments in \eqn{\omega}.  Equals the mean CVE
#'     across environments in \eqn{\omega}.}
#'   \item{**iClassRMSD**}{
#'     \eqn{\text{iClassRMSD}(g,\omega) = \sqrt{|\omega|^{-1}
#'     \sum_{j\in\omega} \text{dev}_{ic}(g,j)^2}}, where
#'     \eqn{\text{dev}_{ic}(g,j) = \widehat{\text{CVE}}(g,j) -
#'     \sum_{r=1}^{k} \hat{\lambda}_{rj}\hat{f}_{rg}} measures
#'     residual crossover GEI within the class.}
#' }
#'
#' @param model   An ASReml-R V4 model object containing a Factor Analytic
#'   random term.
#' @param term    Character string identifying the FA random term, written as
#'   `"fa(<EnvFactor>, k):<GenotypeFactor>"`.  The genotype factor may be
#'   wrapped in `vm(...)` for pedigree/genomic models.
#'   Default `"fa(Site, 4):Genotype"`.
#' @param type    Analysis type.  One of:
#'   \describe{
#'     \item{`"all"` (default)}{Compute both FAST and iClass metrics.}
#'     \item{`"FAST"`}{Compute global Overall Performance and Stability only.}
#'     \item{`"iClass"`}{Compute iClass labels, iClassOP, and iClassRMSD only.}
#'   }
#' @param ic.num  Integer.  Number of factors used to form iClasses and compute
#'   iClassOP.  Must be \eqn{\le k}.  Only used when `type` includes iClass.
#'   Default `2`.
#' @param ...     Additional arguments forwarded to [ASExtras4::fa.asreml()].
#'
#' @return A data frame with one row per environment \eqn{\times} genotype
#'   combination, containing:
#'   \describe{
#'     \item{`<EnvFactor>`}{Environment labels.}
#'     \item{`<GenotypeFactor>`}{Genotype labels.}
#'     \item{`loads1`, ..., `loadsK`}{Rotated FA loadings per environment.}
#'     \item{`spec.var`}{Specific (residual) genetic variance per environment.}
#'     \item{`score1`, ..., `scoreK`}{Rotated FA score EBLUPs per genotype.}
#'     \item{`fitted1`, ..., `fittedK`}{Per-factor contributions to CVE:
#'       \eqn{\hat{\lambda}_{rj}\hat{f}_{rg}}.}
#'     \item{`CVE`}{Common Variety Effect (sum of fitted values).}
#'     \item{`VE`}{Total Variety Effect: `CVE + spec.var`.}
#'     \item{`OP`}{(FAST) Overall Performance  -- same value repeated for each
#'       environment row of a genotype.}
#'     \item{`dev`}{(FAST, \eqn{k > 1}) Residual from first-factor regression.}
#'     \item{`stab`}{(FAST, \eqn{k > 1}) RMSD stability.}
#'     \item{`iclass`}{(iClass) Sign-pattern iClass label.}
#'     \item{`iClassOP`}{(iClass) Within-iClass Overall Performance.}
#'     \item{`iClassRMSD`}{(iClass) Within-iClass RMSD.}
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
#' @seealso [ASExtras4::fa.asreml()]
#' @export
fast <- function(model, term = "fa(Site, 4):Genotype",
                 type   = c("all", "FAST", "iClass"),
                 ic.num = 2L,
                 ...) {

  type <- match.arg(type)

  # ---- Parse the term string -------------------------------------------
  parts   <- strsplit(term, ":")[[1L]]
  fa_idx  <- grep("^fa\\s*\\(", parts)
  if (!length(fa_idx))
    stop("'term' must contain a 'fa(...)' component, e.g. \"fa(Site,4):Genotype\".")

  fa_part  <- parts[fa_idx]
  gen_part <- parts[-fa_idx]

  # Environment factor name: extract first argument of fa(...)
  sterm <- trimws(sub(",.*", "", sub("^fa\\s*\\(", "", fa_part)))

  # Genotype factor name: strip vm(...) wrapper if present
  if (any(grepl("^vm\\s*\\(", gen_part))) {
    vm_part <- gen_part[grep("^vm\\s*\\(", gen_part)]
    gterm   <- trimws(sub(",.*", "", sub("^vm\\s*\\(", "", gsub("\\)", "", vm_part))))
  } else {
    gterm <- trimws(gen_part)
  }

  # ---- Data and fa.asreml -----------------------------------------------
  dat  <- eval(model$call$data)
  sfa  <- fa.asreml(model, ...)

  # ---- Scores: pivot long -> wide (m genotypes * k factors) ---------------
  # In $scores, the sterm column contains factor labels "Comp1","Comp2",...
  sc_long   <- sfa$blups[[term]]$scores
  score_list <- tapply(sc_long$blupr, sc_long[[sterm]], function(el) el)
  score_mat  <- do.call(cbind, score_list)         # m * k
  k  <- ncol(score_mat)
  m  <- nrow(score_mat)
  colnames(score_mat) <- paste0("score", seq_len(k))

  genotypes <- if (!is.null(rownames(score_mat))) rownames(score_mat)
               else levels(dat[[gterm]])

  # ---- Loadings and specific variances: t * k ---------------------------
  loads_mat <- as.matrix(sfa$gammas[[term]]$"rotated loads")  # t * k
  spec_var  <- sfa$gammas[[term]]$"specific var"              # length t
  envs      <- rownames(loads_mat)
  t_envs    <- length(envs)
  colnames(loads_mat) <- paste0("loads", seq_len(k))

  # ---- Validate ic.num --------------------------------------------------
  do_iclass <- type %in% c("iClass", "all")
  do_fast   <- type %in% c("FAST",   "all")

  if (do_iclass) {
    ic.num <- as.integer(ic.num)
    if (ic.num < 1L || ic.num > k)
      stop("'ic.num' must be between 1 and k = ", k, ".")
    if (k == 1L)
      warning("Only one factor  -- iClass with ic.num = 1 places all environments ",
              "in one class.", call. = FALSE)
  }

  # ---- CVE via matrix multiplication ------------------------------------
  # CVE_mat[g,j] = sum_r score_r(g) * loads_r(j) = score_mat %*% t(loads_mat)
  CVE_mat <- score_mat %*% t(loads_mat)       # m * t
  VE_mat  <- CVE_mat + rep(spec_var, each = m) # m * t (broadcast spec_var)

  # ---- Build base long-format data frame (environment-major order) ------
  # Row i encodes: genotype geno_rep[i] in environment env_rep[i]
  # as.vector(m*t matrix) is column-major = environment-major 
  env_rep  <- rep(seq_len(t_envs), each = m)
  geno_rep <- rep(seq_len(m),      times = t_envs)

  out <- data.frame(
    setNames(list(factor(envs[env_rep],      levels = envs)),      sterm),
    setNames(list(factor(genotypes[geno_rep], levels = genotypes)), gterm),
    loads_mat[env_rep,  , drop = FALSE],
    spec.var = spec_var[env_rep],
    score_mat[geno_rep, , drop = FALSE],
    CVE = as.vector(CVE_mat),
    VE  = as.vector(VE_mat),
    stringsAsFactors = FALSE
  )

  # Per-factor fitted values: fitted_r(g,j) = score_r(g) * loads_r(j)
  for (r in seq_len(k))
    out[[paste0("fitted", r)]] <- out[[paste0("score", r)]] * out[[paste0("loads", r)]]

  # ---- Pre-compute sum of first ic.num fitted values (needed for iClass) -
  if (do_iclass && k > 1L) {
    fitted_ic_mat <- Reduce("+", lapply(seq_len(ic.num), function(r)
      outer(score_mat[, r], loads_mat[, r])))   # m * t: sum_r fitted_r for r=1:ic.num
  }

  # ====================================================================
  # FAST metrics
  # ====================================================================
  if (do_fast) {

    # OP = mean(loads1) * score1  [trait-scale non-crossover performance]
    out$OP <- mean(loads_mat[, 1L]) * out$score1

    if (k > 1L) {
      # dev = CVE - fitted1: residual from first-factor regression
      #     = combined contribution of all bipolar (higher-order) factors
      fitted1_mat <- outer(score_mat[, 1L], loads_mat[, 1L])  # m * t
      dev_mat     <- CVE_mat - fitted1_mat
      out$dev  <- as.vector(dev_mat)

      # stab = RMSD per genotype across ALL environments
      # rowMeans(dev^2) is more efficient than tapply on the vector
      out$stab <- sqrt(rowMeans(dev_mat^2))[geno_rep]
    }
  }

  # ====================================================================
  # iClass metrics
  # ====================================================================
  if (do_iclass) {

    # --- Classify environments by sign pattern of first ic.num loadings --
    sign_str  <- apply(loads_mat[, seq_len(ic.num), drop = FALSE], 1L,
                       function(el) paste(ifelse(el >= 0, "p", "n"), collapse = ""))
    iclass_levs <- unique(sign_str)                 # in order of first occurrence
    out$iclass  <- factor(sign_str[as.integer(out[[sterm]])], levels = iclass_levs)

    # --- Per-(genotype, iClass) matrices ---------------------------------
    n_ic <- length(iclass_levs)
    iClassOP_mat   <- matrix(NA_real_, m, n_ic, dimnames = list(genotypes, iclass_levs))
    iClassRMSD_mat <- matrix(NA_real_, m, n_ic, dimnames = list(genotypes, iclass_levs))

    for (w in iclass_levs) {

      env_w <- which(sign_str == w)   # column indices for this iClass

      # iClassOP(g, omega) = sum_{r=1}^{ic.num} mean_{jinomega}(loads_r(j)) * score_r(g)
      #                = score_mat[, 1:ic.num] %*% colMeans(loads_mat[env_w, 1:ic.num])
      mld_w               <- colMeans(loads_mat[env_w, seq_len(ic.num), drop = FALSE])
      iClassOP_mat[, w]   <- score_mat[, seq_len(ic.num), drop = FALSE] %*% mld_w

      # iClassRMSD(g, omega) = sqrt(mean_{jinomega}[dev_ic(g,j)^2])
      # dev_ic = CVE - sum_{r=1}^{ic.num} fitted_r  (residual from ic.num-factor part)
      if (ic.num < k) {
        dev_ic_w <- CVE_mat[, env_w, drop = FALSE] -
                    fitted_ic_mat[, env_w, drop = FALSE]   # m * |env_w|
      } else {
        dev_ic_w <- matrix(0, m, length(env_w))
      }
      iClassRMSD_mat[, w] <- sqrt(rowMeans(dev_ic_w^2))
    }

    # Expand back to row-level using (gterm, iclass) indexing
    gi <- cbind(as.character(out[[gterm]]), as.character(out$iclass))
    out$iClassOP   <- iClassOP_mat[gi]
    out$iClassRMSD <- iClassRMSD_mat[gi]
  }

  # ---- Sort and return ---------------------------------------------------
  sort_by <- if (do_iclass) c("iclass", sterm, gterm) else c(sterm, gterm)
  out <- out[do.call(order, out[, sort_by]), ]
  rownames(out) <- NULL
  out
}
