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
  
  # Define balance between behaviors if not provided
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
  }
  
  # Select data
  data = data[, c(compo_baseline, compo_followup, outcome_baseline, outcome_followup, moderator, covariates)]
  data = data[complete.cases(data),]
  
  
  
  # Do we run cross-sectional or longitudinal?
  if (is.null(compo_followup) & is.null(outcome_followup)) {  # then, cross-sectional
    stop("compo_followup and outcome_followup are not specified. 
          Please, use myfunctions::lm_coda for cross-sectional analyses.")
  } else if (!is.null(compo_baseline) & !is.null(compo_followup) & 
             !is.null(outcome_baseline) & !is.null(outcome_followup)) { 
    # longitudinal with repeated measurements of predictor and outcome
    # for now, this is adapted to HealthyMoms 
    # (TO DO: make this piece of code more flexible)

    
    # standardise baseline compositions
    old_comps = data[,compo_baseline]
    comp_totals = rowSums(old_comps)
    data[,compo_baseline] = old_comps/matrix(comp_totals, ncol = length(compo_baseline), nrow = nrow(old_comps))
    
    # standardise followup compositions
    old_comps = data[,compo_followup]
    comp_totals = rowSums(old_comps)
    data[, compo_followup] = old_comps/matrix(comp_totals, ncol = length(compo_followup), nrow = nrow(old_comps))
    
    # Define compositions
    psi = compositions::gsi.buildilrBase(t(balance))
    compo_base = ilr(data[,compo_baseline], V = psi)
    colnames(compo_base) = c("ilr1_base", "ilr2_base", "ilr3_base")
    compo_fup = ilr(data[,compo_followup], V = psi)
    colnames(compo_fup) = c("ilr1_fup", "ilr2_fup", "ilr3_fup")
    compo_change = compo_fup - compo_base
    colnames(compo_change) = c("ilr1_ch", "ilr2_ch", "ilr3_ch")
    
    # Define outcome
    outcome_change = data[,outcome_followup] - data[,outcome_baseline]
    outcome_change = as.numeric(as.matrix(outcome_change[,1]))
    
    # Scaling if required
    if (isTRUE(scale)) {
      data[,outcome_followup] = scale(data[,outcome_followup])
      data[,outcome_baseline] = scale(data[,outcome_baseline])
      data[,outcome_change] = scale(data[,outcome_change])
      if (length(moderator) > 0) data[,moderator] = scale(data[,moderator])
      data[,covariates] = scale(data[,covariates])
    }
    
    # Model
    if (length(covariates) < 1) {
      if (length(moderator) == 0) { # no moderation analysis
        model = lm(outcome_change~., data = cbind(compo_base, compo_change,
                                                  data[,outcome_baseline]))
      } else if (length(moderator) == 1) { # moderation analysis
        if (longAnalysis == "prospective") {
          ind = cbind(compo_base, compo_change,
                      data[,outcome_baseline],
                      data[,moderator],
                      compo_base[,1] * data[,moderator])
          colnames(ind)[ncol(ind)] = paste0("ilr1_base*", moderator)
          model = lm(outcome_change~., data = ind)
        } else if (longAnalysis == "change") {
          ind = cbind(compo_base, compo_change,
                      data[,outcome_baseline],
                      data[,moderator],
                      compo_change[,1] * data[,moderator])
          colnames(ind)[ncol(ind)] = paste0("ilr1_ch*", moderator)
          model = lm(outcome_change~., data = ind)
        }
        
      }
      
    }
    if (length(covariates) > 0) {
      if (length(moderator) == 0) { # no moderation analysis
        model = lm(outcome_change~., data = cbind(compo_base, compo_change,
                                                  data[,outcome_baseline],
                                                  data[,covariates]))
      } else if (length(moderator) == 1) { # moderation analysis
        if (longAnalysis == "prospective") {
          ind = cbind(compo_base, compo_change,
                      data[,outcome_baseline],
                      data[,moderator],
                      data[,covariates],
                      compo_base[,1] * data[,moderator])
          colnames(ind)[ncol(ind)] = paste0("ilr1_base*", moderator)
          model = lm(outcome_change~., data = ind)
        } else if (longAnalysis == "change") {
          ind = cbind(compo_base, compo_change,
                      data[,outcome_baseline],
                      data[,moderator],
                      data[,covariates],
                      compo_change[,1] * data[,moderator])
          colnames(ind)[ncol(ind)] = paste0("ilr1_ch*", moderator)
          model = lm(outcome_change~., data = ind)
        } 
      }
    }
  } else if (!is.null(compo_followup) & 
             (is.null(outcome_baseline) | is.null(outcome_followup))) {
    # Repeated measures of compo_followup, only one measured outcome
    # pick outcome
    outcome_tmp = c(outcome_baseline, outcome_followup)
    outcome = outcome_tmp[which(!is.null(outcome_tmp))]
    rm(outcome_tmp); gc()
    
    # select data
    data = data[, c(compo_baseline, compo_followup, outcome, moderator, covariates)]
    data = data[complete.cases(data),]
    
    # standardise baseline compositions
    old_comps = data[, compo_baseline]
    comp_totals = rowSums(old_comps)
    data[, compo_baseline] = old_comps/matrix(comp_totals, ncol = length(compo_baseline), nrow = nrow(old_comps))
    
    # standardise followup compositions
    old_comps = data[, compo_followup]
    comp_totals = rowSums(old_comps)
    data[, compo_followup] = old_comps/matrix(comp_totals, ncol = length(compo_followup), nrow = nrow(old_comps))
    
    # Define compositions
    psi = compositions::gsi.buildilrBase(t(balance))
    compo_base = compositions::ilr(data[,compo_baseline], V = psi)
    colnames(compo_base) = paste0("ilr", 1:ncol(compo_base), "_base")
    compo_fup = compositions::ilr(data[,compo_followup], V = psi)
    colnames(compo_fup) = paste0("ilr", 1:ncol(compo_fup), "_fup")
    compo_change = compo_fup - compo_base
    colnames(compo_change) = paste0("ilr", 1:ncol(compo_change), "_ch")
    
    # Scaling if required
    if (isTRUE(scale)) {
      data[,outcome] = scale(data[,outcome])
      if (length(moderator) > 0) data[,moderator] = scale(data[,moderator])
      data[,covariates] = scale(data[,covariates])
    }
    
    # data to fit in model
    if (longAnalysis == "prospective") compo2 = compo_fup
    if (longAnalysis == "change") compo2 = compo_change
    
    if (length(moderator) == 0) {
      dat2fit = cbind(data[, c(outcome, covariates)], compo_base, compo2)
      colnames(dat2fit)[1:(1 + length(covariates))] = c(outcome, covariates)
    } else if (length(moderator) == 1) {
      dat2fit = cbind(data[, c(outcome, covariates)], compo_base, compo2, compo_base * data[, moderator])
      colnames(dat2fit)[1:(1 + length(covariates))] = c(outcome, covariates)
      colnames(dat2fit)[ncol(dat2fit)] = paste0("ilr1_base*", moderator)
    }
    
    # model
    model = lm(dat2fit[, outcome]~., data = dat2fit[, -1])
    
  }
  
  return(model)
}
