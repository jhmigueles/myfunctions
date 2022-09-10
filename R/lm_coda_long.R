#' lm_coda_long
#' @description Run compositional linear regression models using isometric log ratios. Option to perform longitudinal models and moderation analyses implemented.
#'
#' @param covariates Numeric, matrix, or data frame with the covariates to adjust the models for
#' @param scale Logical, whether to standardize or not the output (default = FALSE)
#' @param balance Balance between behaviors to be used in the ILR coordinates
#' @param compo_baseline Variable names for the baseline composition
#' @param compo_followup Variable names for the follow-up composition (optional)
#' @param outcome_baseline Variable name for the baseline outcome
#' @param outcome_followup Variable name for the follow-up outcome (optional)
#' @param moderator Moderator variable, if any (optional)
#' @param longAnalysis Longitudinal analysis to perform, either prospective or changes for now.
#' @param data Data frame with data
#'
#' @return List with models fit for each one of the ILR coordinates as leading coordinate
#' @export
#' @importFrom  compositions gsi.buildilrBase ilr
#'
lm_coda_long = function(compo_baseline,    # Variable names for the baseline composition
                        compo_followup = NULL,    # Variable names for the follow-up composition
                        outcome_baseline,  # Variable name for the baseline outcome
                        outcome_followup = NULL,  # Variable name for the follow-up outcome
                        moderator = NULL,
                        longAnalysis = c("prospective", "change")[1],
                        covariates = NULL,        # Variable names for the covariates
                        scale = FALSE,     # standardize data?
                        balance = NA,      # balance between behaviors to test
                        data){      # dataset

  # Check if it is a cross-sectional analysis 
  if (is.null(compo_followup) & is.null(outcome_followup)) {
    stop("compo_followup and outcome_followup are not specified. 
          Please, use myfunctions::lm_coda for cross-sectional analyses.")
  }
  
  # Define balance between behaviors if not provided ----
  if (!is.matrix(balance) & !is.data.frame(balance)) {
    if (is.na(balance) | is.null(balance)) {
      balance = matrix(data = NA, nrow = length(compo_baseline) - 1, ncol = length(compo_baseline))
      for (i in 1:nrow(balance)) {
        if (i == 1) {
          balance[i,] = c(1, rep(-1, ncol(balance) - 1))
        } else {
          zeroes = sum(balance[i - 1,] == 0) + 1
          balance[i, 1:zeroes] = 0
          balance[i, zeroes + 1] = 1
          balance[i, is.na(balance[i,])] = -1
        }
      }
    } else {
      stop("balance should be a matrix or data frame with the desired time exchange
           between the composition parts. Leave empty for automatic calculation.")
    }
  }
  
  psi = compositions::gsi.buildilrBase(t(balance))
  
  # Select data ----
  data = data[, c(compo_baseline, compo_followup, outcome_baseline, outcome_followup, moderator, covariates)]
  data = data[complete.cases(data),]
  
  # Standardize compositions and get ILR coordinates ----
  if (!is.null(compo_baseline)) {
    old_comps = data[, compo_baseline]
    comp_totals = rowSums(old_comps)
    data[, compo_baseline] = old_comps/matrix(comp_totals, ncol = length(compo_baseline), nrow = nrow(old_comps))
    compo_base = compositions::ilr(data[,compo_baseline], V = psi)
    colnames(compo_base) = paste0("ilr", 1:ncol(compo_base), "_base")
  }
  if (!is.null(compo_followup)) {
    old_comps = data[, compo_followup]
    comp_totals = rowSums(old_comps)
    data[, compo_followup] = old_comps/matrix(comp_totals, ncol = length(compo_followup), nrow = nrow(old_comps))
    compo_fup = compositions::ilr(data[,compo_followup], V = psi)
    colnames(compo_fup) = paste0("ilr", 1:ncol(compo_fup), "_fup")
    # compo change
    compo_change = compo_fup - compo_base
    colnames(compo_change) = paste0("ilr", 1:ncol(compo_change), "_ch")
  }
  
  # Define outcome ----
  if (!is.null(outcome_baseline) & !is.null(outcome_followup)) {
    outcome = data[,outcome_followup] - data[,outcome_baseline]
    outcome = as.numeric(as.matrix(outcome[,1]))
  } else if (is.null(outcome_baseline) | is.null(outcome_followup)) {
    outcome_tmp = c(outcome_baseline, outcome_followup)
    outcome = outcome_tmp[which(!is.null(outcome_tmp))]
    outcome_baseline = NULL # so that it is not included in data2fit later on
    rm(outcome_tmp); gc()
  }
  
  # data to fit in model ----
  if (longAnalysis == "prospective") compo2 = compo_fup
  if (longAnalysis == "change") compo2 = compo_change
  
  if (length(moderator) == 0) {
    dat2fit = cbind(data[, c(outcome, covariates, outcome_baseline)], compo_base, compo2)
    if (colnames(dat2fit)[1] == "") colnames(dat2fit)[1] = outcome
  } else if (length(moderator) == 1) {
    if (!is.numeric(data[, moderator])) stop("Moderator should be a numeric variable in the dataset!")
    if (longAnalysis == "prospective") {
      ilrMOD = compo_base[,1]
      ilrMODn = "ilr1_base*"
    } else if (longAnalysis == "change") {
      ilrMOD = compo_change[,1]
      ilrMODn = "ilr1_ch*"
    }
    dat2fit = cbind(data[, c(outcome, covariates, outcome_baseline, moderator)], compo_base, compo2, ilrMOD * data[, moderator])
    if (colnames(dat2fit)[1] == "") colnames(dat2fit)[1] = outcome
    colnames(dat2fit)[ncol(dat2fit)] = paste0("ilr1_base*", moderator)
  }
  dat2fit = as.data.frame(dat2fit)
  
  # Scaling if required
  if (isTRUE(scale)) dat2fit = scale(dat2fit)
  
  # Fit model ----
  model = lm(dat2fit[, outcome]~., data = dat2fit[, -1])
  
  return(model)
}
