#' Basic function for cross-sectional logistic regression CoDA models
#'
#' @param data data frame or matrix with variables of interest
#' @param compo String with colnames of composition variables or matrix/dataframe with these variables
#' @param outcome String with colname of outcome or numeric vector with the outcome
#' @param covariates String with names of covariates (optional) or matrix/dataframe with these variables
#' @param scale Logical. Whether to standardize results (default = FALSE)
#'
#' @return Logistic regression model (glm class)
#' @export
#'
logisticReg_coda = function(data, compo, outcome, covariates = c(), scale = FALSE){
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
    if (length(covariates) < 1) model[[i]] = glm(outcome~., data = z[[i]], family = "binomial")
    if (length(covariates) > 0) model[[i]] = glm(outcome~., data = cbind(z[[i]], covariates), family = "binomial")
  }
  return(model)
}