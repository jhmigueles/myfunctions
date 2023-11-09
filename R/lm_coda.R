#' Basic function for cross-sectional linear regression CoDA models
#'
#' @param data data frame or matrix with variables of interest
#' @param compo String with colnames of composition variables or matrix/dataframe with these variables
#' @param outcome String with colname of outcome or numeric vector with the outcome
#' @param covariates String with names of covariates (optional) or matrix/dataframe with these variables
#' @param scale Logical. Whether to standardize results (default = FALSE)
#' @param moderator Optional string with the colname of the interaction term (default = NULL)
#'
#' @return Linear regression model (lm class)
#' @export
#'
lm_coda = function(data, compo, outcome, covariates = c(), 
                   moderator = NULL, scale = FALSE){
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
    if (length(covariates) < 1) dat2fit = z[[i]]
    if (length(covariates) > 0) dat2fit = cbind(z[[i]], covariates)
    if (length(moderator) == 1) {
      ilrMOD = z[[i]][1]
      ilrMODn = "ilr1*"
      dat2fit = cbind(dat2fit, data[, moderator], ilrMOD * data[, moderator])
      colnames(dat2fit)[ncol(dat2fit)] = paste0(ilrMODn, moderator)
    }
    # fit model
    model[[i]] = lm(outcome~., data = dat2fit)
  }
  return(model)
}