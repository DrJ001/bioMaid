# ---------------------------------------------------------------------------
# Internal helpers for FA term discovery
# ---------------------------------------------------------------------------

#' Pull the model frame from a fitted ASReml-R object
#' @noRd
.fa_model_frame <- function(model) {
  mf <- model$mf
  if (is.null(mf)) {
    if (is.null(model$RDS))
      stop("'model' is missing its model frame and has no RDS component.",
           call. = FALSE)
    mf <- readRDS(model$RDS)
  }
  mf
}

#' Locate every fa() term in the random formula of a fitted model
#'
#' Returns a list with one element per FA term, each holding the term label
#' plus the outer (environment) and inner (genotype) variable names and the
#' wrapper function applied to the inner variable.
#'
#' @noRd
.fa_find_terms <- function(mf) {
  mt <- attr(mf, "model.terms")$random
  if (is.null(mt))
    stop("'model' has no random terms.", call. = FALSE)

  tobj <- mt$Terms.obj
  fa_rows <- attr(tobj, "specials")$fa
  if (length(fa_rows) == 0L)
    stop("No fa() term found in the random model formula.", call. = FALSE)

  factors <- attr(tobj, "factors")
  vars    <- dimnames(factors)[[1]]
  labels  <- attr(tobj, "term.labels")

  # Every model term (column) that involves an fa() variable (row)
  cols <- unique(unlist(lapply(fa_rows, function(i) which(factors[i, ] > 0))))

  out <- lapply(cols, function(j) {
    parts <- mt$Vars[vars[which(factors[, j] > 0)]]
    fun   <- vapply(parts, function(x) x$Fun, character(1L))
    obj   <- vapply(parts, function(x) x$Obj, character(1L))
    nam   <- vapply(parts, function(x) x$FacNam[1L], character(1L))

    is_fa <- fun == "fa"
    if (sum(!is_fa) != 1L)
      stop("Only a first-order interaction with an fa() term is supported; ",
           "found '", labels[j], "'.", call. = FALSE)

    list(
      label     = labels[j],
      outer     = unname(obj[is_fa]),
      outer_nam = unname(nam[is_fa]),
      inner     = unname(obj[!is_fa]),
      inner_nam = unname(nam[!is_fa]),
      inner_fun = unname(fun[!is_fa])
    )
  })
  names(out) <- labels[cols]

  # Process ide() terms last so their vm() partner is already available
  # when the combined "total" structure is assembled.
  ide <- which(vapply(out, function(x) x$inner_fun == "ide", logical(1L)))
  if (length(ide) > 0L)
    out <- out[c(setdiff(seq_along(out), ide), ide)]

  out
}

# ---------------------------------------------------------------------------
# Internal helpers for the variance-parameter and BLUP calculations
# ---------------------------------------------------------------------------

#' Rotate loadings and build the FA covariance structure for one term
#'
#' Reproduces the ASExtras4 convention exactly: the loading matrix is rotated
#' by the right singular vectors of itself and negated, and the *same*
#' rotation is later applied to the factor scores.  The negation makes the
#' sign of each factor arbitrary but stable, and downstream consumers
#' (notably [fastIC()]'s interaction classes) depend on it.
#'
#' @noRd
.fa_structure <- function(vpar, envs) {
  nenv   <- length(envs)
  gammas <- matrix(vpar, nrow = nenv)
  k      <- ncol(gammas) - 1L

  psi <- gammas[, 1L]
  lam <- gammas[, -1L, drop = FALSE]

  # Rotation is only defined for k > 1; ASExtras4 leaves a single factor as-is
  rot <- NULL
  if (k > 1L) {
    rot <- svd(lam)$v
    lam <- -lam %*% rot
  }

  dimnames(lam) <- list(envs, paste("fac", seq_len(k), sep = "_"))
  names(psi)    <- envs

  gmat <- lam %*% t(lam) + diag(psi, nrow = nenv)
  dimnames(gmat) <- list(envs, envs)
  cmat <- stats::cov2cor(gmat)

  # Correlation-scaled loadings: rows lie inside the unit hypersphere
  lam_cor <- diag(1 / sqrt(diag(gmat)), nrow = nenv) %*% lam
  dimnames(lam_cor) <- dimnames(lam)

  list(k = k, loads = lam, loads_cor = lam_cor, spec_var = psi,
       Gmat = gmat, Cmat = cmat, rotation = rot)
}

#' Variance accounted for, per environment and overall
#'
#' Specific variances are clamped at zero for this calculation only: ASReml-R
#' can return fractionally negative boundary estimates which would otherwise
#' produce negative proportions.  The unclamped values are always reported in
#' `spec_var`.
#'
#' @noRd
.fa_vaf <- function(loads, spec_var, envs, outer_name) {
  k         <- ncol(loads)
  loads_sq  <- loads^2
  spec_clamped <- pmax(0, spec_var)
  total_var <- rowSums(loads_sq) + spec_clamped

  vaf_env <- as.data.frame(loads_sq / total_var)
  names(vaf_env) <- paste0("Factor", seq_len(k))
  vaf_env[[outer_name]] <- envs
  vaf_env$Specific      <- unname(spec_clamped / total_var)
  vaf_env$total_var     <- unname(total_var)
  vaf_env <- vaf_env[, c(outer_name, paste0("Factor", seq_len(k)),
                         "Specific", "total_var"), drop = FALSE]
  rownames(vaf_env) <- NULL

  total_all  <- sum(total_var)
  pct_factor <- colSums(loads_sq) / total_all
  pct_spec   <- sum(spec_clamped) / total_all

  vaf_summary <- data.frame(
    factor  = c(paste0("Factor ", seq_len(k)), "Specific"),
    pct_var = c(unname(pct_factor), pct_spec),
    stringsAsFactors = FALSE
  )
  vaf_summary$cum_pct <- cumsum(vaf_summary$pct_var)

  list(vaf_env = vaf_env, vaf_summary = vaf_summary,
       vaf_total = 100 * sum(pct_factor))
}

#' Split an ASReml-R coefficient label of the form "Site_Env01:Variety_Var01"
#' @noRd
.fa_split_labels <- function(labels) {
  parts <- strsplit(labels, ":", fixed = TRUE)
  strip <- function(x) vapply(
    strsplit(x, "_", fixed = TRUE),
    function(el) paste(el[-1L], collapse = "_"),
    character(1L)
  )
  list(
    outer = strip(vapply(parts, `[`, character(1L), 1L)),
    inner = strip(vapply(parts, `[`, character(1L), 2L))
  )
}

#' Extract BLUPs and rotated factor scores for one FA term
#' @noRd
.fa_blups <- function(model, info, fas, mf, y) {
  k     <- fas$k
  outer_name <- info$outer
  inner_name <- info$inner

  cc <- stats::coef(model, list = TRUE)[[info$label]]
  lab <- .fa_split_labels(dimnames(cc)[[1L]])

  df <- data.frame(blup = as.vector(cc), stringsAsFactors = FALSE)
  df[[outer_name]] <- lab$outer
  df[[inner_name]] <- lab$inner

  # The score EBLUPs are stored against pseudo-levels Comp1 ... Compk
  comps    <- paste0("Comp", seq_len(k))
  is_score <- df[[outer_name]] %in% comps

  scores <- df[is_score, , drop = FALSE]
  score_mat <- matrix(scores$blup, ncol = k)
  if (!is.null(fas$rotation))
    score_mat <- -score_mat %*% fas$rotation
  scores$blupr <- as.vector(score_mat)

  df <- df[!is_score, , drop = FALSE]
  df$regblup <- as.vector(score_mat %*% t(fas$loads))

  # Restrict to genotypes that actually appear in the model frame.  For vm()
  # terms the FA term carries every level of the relationship matrix, which
  # can be a strict superset.
  in_met <- unique(as.character(mf[[inner_name]]))
  df     <- df[df[[inner_name]] %in% in_met, , drop = FALSE]
  scores <- scores[scores[[inner_name]] %in% in_met, , drop = FALSE]

  ord <- function(d) d[order(utils::type.convert(d[[outer_name]], as.is = TRUE),
                             utils::type.convert(d[[inner_name]], as.is = TRUE)), ,
                       drop = FALSE]
  df     <- ord(df)
  scores <- ord(scores)

  # Number of non-missing observations per genotype x environment.  Matched by
  # name rather than position so that it stays correct when level names do not
  # sort lexicographically.
  pres <- tapply(y,
                 list(factor(as.character(mf[[inner_name]])), mf[[outer_name]]),
                 function(x) length(x[!is.na(x)]))
  df$pres <- pres[cbind(match(df[[inner_name]], rownames(pres)),
                        match(df[[outer_name]], colnames(pres)))]

  rownames(df)     <- NULL
  rownames(scores) <- NULL

  df     <- df[, c(outer_name, inner_name, "blup", "regblup", "pres"),
               drop = FALSE]
  scores <- scores[, c(outer_name, inner_name, "blup", "blupr"), drop = FALSE]

  list(blups = df, scores = scores)
}

# ---------------------------------------------------------------------------
# faSummary()
# ---------------------------------------------------------------------------

#' Summarise the factor analytic structure of an ASReml-R model
#'
#' Extracts and rotates the factor analytic (FA) variance structure from a
#' fitted ASReml-R V4 model, returning the genetic covariance and correlation
#' matrices, rotated loadings, specific variances, variance accounted for
#' (VAF), and the genotype BLUPs and factor score EBLUPs.  One entry is
#' produced for every `fa()` term in the random model formula.
#'
#' @details
#' # Rotation
#'
#' An FA(\eqn{k}) structure decomposes the \eqn{t \times t} genetic covariance
#' matrix as \eqn{\bm{G} = \bm{\Lambda}\bm{\Lambda}^\top + \bm{\Psi}}, where
#' \eqn{\bm{\Lambda}} holds the \eqn{t \times k} loadings and \eqn{\bm{\Psi}}
#' is the diagonal matrix of specific variances.  The decomposition is only
#' identified up to an orthogonal rotation, so the raw ASReml-R loadings are
#' rotated by the right singular vectors of \eqn{\bm{\Lambda}} and negated.
#' The identical rotation is applied to the factor scores, so that
#' \eqn{\bm{\Lambda}} and the scores remain conformable and their product --
#' the common variety effect -- is invariant.
#'
#' After rotation the first factor typically captures the dominant
#' non-crossover genotype by environment interaction with all-positive
#' loadings; higher factors are bipolar.  The sign of each factor is arbitrary
#' but reproducible.
#'
#' # Variance accounted for
#'
#' The total genetic variance for environment \eqn{j} is
#' \eqn{\sigma^2_j = \sum_r \lambda_{rj}^2 + \psi_j}.  `vaf_env` reports
#' \eqn{\lambda_{rj}^2 / \sigma^2_j} for each factor together with the
#' specific proportion \eqn{\psi_j / \sigma^2_j}; each row sums to one.
#' `vaf_summary` pools these across environments.
#'
#' # Combined `vm()` and `ide()` terms
#'
#' A common MET model fits both an additive term, `fa(Site, k):vm(Variety, G)`,
#' and a non-additive term, `fa(Site, k):ide(Variety)`.  When `combine.ide` is
#' `TRUE` and such a pair shares the same outer and inner variables, an extra
#' entry named `"<outer>:<inner>-total"` is appended holding the summed
#' covariance matrix and summed BLUPs.  It has no loadings or scores of its
#' own, since the sum of two FA structures is not itself an FA structure.
#'
#' @param model A fitted `asreml` object whose random formula contains at
#'   least one `fa()` term.
#' @param term  Character vector of FA term labels to summarise, e.g.
#'   `"fa(Site, 3):Variety"`.  The default, `NULL`, summarises every FA term
#'   found in the model.
#' @param blups Logical; return the genotype BLUPs and factor score EBLUPs.
#'   Set to `FALSE` to skip the (potentially expensive) call to
#'   [stats::coef()] when only the variance structure is needed.
#' @param combine.ide Logical; when a `vm()` and an `ide()` term share the
#'   same outer and inner variables, append the combined "total" structure.
#'
#' @return An object of class `"faSummary"`: a list with components
#'   \describe{
#'     \item{`gammas`}{Named list, one element per term, each containing
#'       `Gmat`, `Cmat`, `loads`, `loads_cor`, `spec_var`, `vaf_env`,
#'       `vaf_summary`, `vaf_total`, `k`, `env`, `outer`, `inner` and
#'       `inner_fun`.}
#'     \item{`blups`}{Named list, one element per term, each containing a
#'       `blups` data frame (`<outer>`, `<inner>`, `blup`, `regblup`, `pres`)
#'       and a `scores` data frame (`<outer>`, `<inner>`, `blup`, `blupr`),
#'       where `<outer>` holds the pseudo-levels `Comp1` ... `Compk`.  `NULL`
#'       when `blups = FALSE`.}
#'     \item{`terms`}{Character vector of the term labels summarised.}
#'     \item{`call`}{The matched call.}
#'   }
#'
#' @section Note:
#' `pres` counts the non-missing observations for each genotype by environment
#' combination and is `NA`, not `0`, where a combination is entirely absent.
#'
#' @seealso [plot_faSummary()], [fastIC()], [randomRegress()]
#'
#' @examples
#' \dontrun{
#' model <- asreml::asreml(
#'   yield  ~ Site,
#'   random = ~ fa(Site, 3):Variety,
#'   data   = dat
#' )
#' fas <- faSummary(model)
#'
#' fas$gammas[["fa(Site, 3):Variety"]]$loads
#' fas$gammas[["fa(Site, 3):Variety"]]$vaf_summary
#' head(fas$blups[["fa(Site, 3):Variety"]]$scores)
#' }
#'
#' @export
faSummary <- function(model, term = NULL, blups = TRUE, combine.ide = TRUE) {
  if (!inherits(model, "asreml"))
    stop("'model' must be a fitted 'asreml' object.", call. = FALSE)
  if (!requireNamespace("asreml", quietly = TRUE))
    stop("Package 'asreml' is required to summarise a fitted model.",
         call. = FALSE)

  mf <- .fa_model_frame(model)

  lhs <- attr(mf, "traits")$lhs
  if (length(lhs) > 1L)
    stop("No method for multivariate models.", call. = FALSE)
  y <- mf[[lhs]]

  info <- .fa_find_terms(mf)

  if (!is.null(term)) {
    missing_terms <- setdiff(term, names(info))
    if (length(missing_terms) > 0L)
      stop("FA term(s) not found in the model: ",
           paste(missing_terms, collapse = ", "), ".\n  Available: ",
           paste(names(info), collapse = ", "), call. = FALSE)
    info <- info[intersect(names(info), term)]
  }

  vpars <- summary(model, vparameters = TRUE)$vparameters

  gammas_lst <- list()
  blups_lst  <- list()

  for (nm in names(info)) {
    this <- info[[nm]]

    envs <- levels(mf[[this$outer]])
    if (length(envs) == 0L)
      stop("'", this$outer, "' is not a factor.", call. = FALSE)
    if (any(grepl(":", envs, fixed = TRUE)))
      stop("Levels of '", this$outer, "' cannot contain ':'.", call. = FALSE)
    if (any(grepl(":", unique(as.character(mf[[this$inner]])), fixed = TRUE)))
      stop("Levels of '", this$inner, "' cannot contain ':'.", call. = FALSE)

    fas <- .fa_structure(vpars[[nm]], envs)
    vaf <- .fa_vaf(fas$loads, fas$spec_var, envs, this$outer)

    gammas_lst[[nm]] <- list(
      Gmat        = fas$Gmat,
      Cmat        = fas$Cmat,
      loads       = fas$loads,
      loads_cor   = fas$loads_cor,
      spec_var    = fas$spec_var,
      vaf_env     = vaf$vaf_env,
      vaf_summary = vaf$vaf_summary,
      vaf_total   = vaf$vaf_total,
      k           = fas$k,
      env         = envs,
      outer       = this$outer,
      inner       = this$inner,
      inner_fun   = this$inner_fun
    )

    if (blups)
      blups_lst[[nm]] <- .fa_blups(model, this, fas, mf, y)
  }

  if (combine.ide) {
    combined <- .fa_combine_ide(info, gammas_lst, blups_lst, blups)
    gammas_lst <- combined$gammas
    blups_lst  <- combined$blups
  }

  structure(
    list(
      gammas = gammas_lst,
      blups  = if (blups) blups_lst else NULL,
      terms  = names(gammas_lst),
      call   = match.call()
    ),
    class = "faSummary"
  )
}

#' Append the combined vm() + ide() structure where such a pair exists
#' @noRd
.fa_combine_ide <- function(info, gammas_lst, blups_lst, blups) {
  ide_terms <- names(info)[vapply(info, function(x) x$inner_fun == "ide",
                                  logical(1L))]

  for (nm in ide_terms) {
    this <- info[[nm]]
    partner <- names(info)[vapply(info, function(x) {
      x$inner_fun == "vm" && x$outer == this$outer && x$inner == this$inner
    }, logical(1L))]
    if (length(partner) != 1L) next
    if (is.null(gammas_lst[[nm]]) || is.null(gammas_lst[[partner]])) next

    total_nm <- paste0(this$outer, ":", this$inner, "-total")
    gmat <- gammas_lst[[nm]]$Gmat + gammas_lst[[partner]]$Gmat

    gammas_lst[[total_nm]] <- list(
      Gmat      = gmat,
      Cmat      = stats::cov2cor(gmat),
      env       = gammas_lst[[nm]]$env,
      outer     = this$outer,
      inner     = this$inner,
      inner_fun = "total"
    )

    if (blups && !is.null(blups_lst[[nm]]) && !is.null(blups_lst[[partner]])) {
      a <- blups_lst[[nm]]$blups
      b <- blups_lst[[partner]]$blups
      key_a <- paste(a[[this$outer]], a[[this$inner]], sep = "\r")
      key_b <- paste(b[[this$outer]], b[[this$inner]], sep = "\r")
      idx   <- match(key_a, key_b)

      total <- a
      total$blup    <- a$blup    + b$blup[idx]
      total$regblup <- a$regblup + b$regblup[idx]
      blups_lst[[total_nm]] <- list(blups = total, scores = NULL)
    }
  }

  list(gammas = gammas_lst, blups = blups_lst)
}

#' Print a factor analytic summary
#'
#' @param x   An object of class `"faSummary"`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @export
print.faSummary <- function(x, ...) {
  cat("<faSummary>\n")
  for (nm in x$terms) {
    g <- x$gammas[[nm]]
    if (identical(g$inner_fun, "total")) {
      cat("  ", nm, ": combined vm() + ide() structure, ",
          nrow(g$Gmat), " environments\n", sep = "")
    } else {
      cat("  ", nm, ": k = ", g$k, ", ", length(g$env),
          " environments, ", format(g$vaf_total, digits = 4),
          "% variance accounted for\n", sep = "")
    }
  }
  invisible(x)
}
