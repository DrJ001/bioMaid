#' @keywords internal
"_PACKAGE"

#' @importFrom stats coef contr.helmert lm model.matrix p.adjust pchisq pf
#'   predict qt qtukey residuals rnorm runif sd setNames terms
#' @importFrom utils combn write.csv
NULL

## fa.asreml is from ASExtras4 (in Suggests).
## Suppress R CMD check NOTE about undefined global.
utils::globalVariables("fa.asreml")
