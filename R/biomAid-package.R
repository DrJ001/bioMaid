#' @keywords internal
"_PACKAGE"

#' @importFrom stats coef contr.helmert lm model.matrix p.adjust pchisq pf
#'   predict qt qtukey residuals rnorm runif sd setNames terms
#' @importFrom utils combn write.csv
NULL

# ---- Internal wrappers for ASReml-R calls --------------------------------
# These thin wrappers live in the biomAid namespace so that
# testthat::local_mocked_bindings() can substitute them in tests without
# requiring an ASReml-R licence in CI.
#
# The FA decomposition itself is implemented natively by faSummary(); the
# ASExtras4 dependency it replaced has been removed.

#' @noRd
.asreml_vparams <- function(model, term) {
  summary(model, vparameters = TRUE)$vparameters[[term]]
}
