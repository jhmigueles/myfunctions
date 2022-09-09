#' Basic function for cross-sectional linear regression CoDA models
#'
#' @param data data frame or matrix with variables of interest
#' @param compo String with colnames of composition variables or matrix/dataframe with these variables
#' @param outcome String with colname of outcome or numeric vector with the outcome
#' @param covariates String with names of covariates (optional) or matrix/dataframe with these variables
#' @param scale Logical. Whether to standardize results (default = FALSE)
#'
#' @return Linear regression model (lm class)
#' @export
#'
lm_coda = function(data, compo, outcome, covariates = c(), scale = FALSE){
  if (is.character(compo)) compo = data[, compo]
  if (is.character(outcome)) outcome = data[, outcome]
  if (is.character(covariates)) covariates = data[, covariates]
  if (scale == TRUE) {
    compo = scale(compo)
    outcome = scale(outcome)
    covariates = scale(covariates)
  }
  z = LR_ILR(compo)
  model = rep(list(NA), times = length(z))
  for (i in 1:length(z)) {
    if (length(covariates) < 1) model[[i]] = lm(outcome~., data = z[[i]])
    if (length(covariates) > 0) model[[i]] = lm(outcome~., data = cbind(z[[i]], covariates))
  }
  return(model)
}