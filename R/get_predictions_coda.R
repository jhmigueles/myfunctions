#' get_predictions_coda
#' @description Get predictions from CoDA models
#'
#' @param data Data frame with composition (baseline and/or follow up), outcome and covariates
#' @param comp_bl names of composition variables in your dataset (baseline)
#' @param comp_fu names of composition variables in your dataset (follow up, optional)
#' @param comp_names names of composition variables for visualization
#' @param outcome_bl Name of outcome at baseline in your dataset
#' @param outcome_fu Name of outcome at follow up in your dataset (optional)
#' @param covs names of your covariates in your dataset
#' @param comparisons Tow options: "one-v-one" OR "prop-realloc"
#' @param increase name of the behavior to increase (name for visualization)
#' @param decrease name of the behavior to decrease (only in "one-v-one")
#' @param delta minutes to reallocate for prediction (e.g., +30 min)
#' @param model "cross-sectional", "longitudinal", or "change"
#'
#' @return
#' @export
#'
#' @importFrom compositions gsi.buildilrBase ilr
#'
#' @examples
get_predictions_coda = function(data = c(),            # dataframe with COMPOSITION (BASELINE (AND FOLLWOUP), OUTCOME, COVARIATES)
                                comp_bl = c(),    # names of composition variables in your dataset (baseline)
                                comp_fu = c(),    # names of composition variables in your dataset (follow up)
                                comp_names = c(),    # names of composition variables for visualization
                                outcome_bl = c(),    # Name of outcome in your dataset
                                outcome_fu = c(),    # Name of outcome in your dataset
                                covs = c(),    # names of your covariates in your dataset
                                comparisons = "prop-realloc", #  "one-v-one" OR "prop-realloc"
                                increase = c(),    # name of the behavior to increase (name for visualization)
                                decrease = c(),       # name of the behavior to decrease (only in "one-v-one")
                                delta = 30,                  # minutes to reallocate for prediction (e.g., +30 min)
                                model = "change") {          # "cross-sectional", "longitudinal", or "change"


  # define dataset

  if (model == "cross-sectional") {
    df = data[, c(comp_bl, outcome_bl, covs)]
    if (!is.null(comp_fu)) warning("model = 'cross-sectional' BUT you defined a follow-up composition")
    if (!is.null(outcome_fu)) warning("model = 'cross-sectional' BUT you defined a follow-up outcome")

    # rename composition variables
    if (length(comp_names) == length(comp_bl)) {
      colnames(df)[1:length(comp_bl)] = comp_names
      comp_bl = comp_names
    } else {warning("composition names is of different length than the number of
                  elements in the composition. comp_names argument is not used.")}


  } else if (model == "longitudinal" | model == "change") {
    if (is.null(comp_fu) | is.null(outcome_fu)) stop("follow-up composition and/or outcome NOT defined")
    df = data[,c(comp_bl, comp_fu, outcome_bl, outcome_fu, covs)]

    # rename composition variables
    if (length(comp_names) == length(comp_bl) & length(comp_names) == length(comp_fu)) {
      comp_names_BL = paste0(comp_names, "_GW14")
      comp_names_FU = paste0(comp_names, "_GW37")
      colnames(df)[1:length(comp_bl)] = comp_names_BL
      colnames(df)[(length(comp_bl) + 1):(length(comp_bl) + length(comp_fu))] = comp_names_FU
      comp_bl = comp_names
    } else {warning("composition names is of different length than the number of
                  elements in the composition. comp_names argument is not used.")}
  } else {
    stop("model argument is not correctly defined")
  }

  # DEFINE NEW VARIABLES FOR LONGITUDINAL/CHANGE MODELS
  if (model == "longitudinal" | model == "change") {

    ### balance between behaviors
    balance = matrix(c(  1,  -1, -1,  -1,    # MVPA / LPA·SB·Sleep
                         0,  1, -1,  -1,     # LPA / SB·Sleep
                         0,   0, 1,  -1      # SB / Sleep
    ), ncol = 4, byrow = TRUE)

    ### ilr coordinates for baseline, follow up and change
    psi = compositions::gsi.buildilrBase(t(balance))
    ilr_bl = compositions::ilr(df[,comp_names_BL], V = psi)
    colnames(ilr_bl) = c("ilr1_bl", "ilr2_bl", "ilr3_bl")
    ilr_fu = compositions::ilr(df[,comp_names_FU], V = psi)
    colnames(ilr_fu) = c("ilr1_fu", "ilr2_fu", "ilr3_fu")
    ilr_ch = ilr_fu - ilr_bl
    colnames(ilr_ch) = c("ilr1_ch", "ilr2_ch", "ilr3_ch")

    ### change in outcome
    df[,gsub("GW14", "CH", outcome_bl)] = df[,outcome_fu] - df[,outcome_bl]
    outcome_ch = gsub("GW14", "CH", outcome_bl)

    ### redefine dataset and define covariates
    if (model == "longitudinal") {   # then, change composition is adjustment
      df = df[, c(comp_names_BL, outcome_bl, covs, outcome_ch)]
      df = cbind(df, ilr_ch)
      covs = c(covs, colnames(ilr_ch))
    } else if (model == "change") { # then, baseline composition is adjustment
      df = df[, c(comp_names_BL, comp_names_FU, covs, outcome_bl, outcome_ch)]
      df = cbind(df, ilr_bl)
      covs = c(colnames(ilr_bl), covs)
    }
  }

  # RUN PREDICTIONS
  if (model == "cross-sectional") {
    preds = get_plus_minus_changes(dataf = df,
                                   y = outcome_bl,
                                   comps = comp_names,
                                   covars = covs,
                                   deltas = seq(-180, 180, by = 1)/(24*60),
                                   comparisons = comparisons,
                                   alpha = 0.05, verbose = FALSE)
  } else if (model == "longitudinal") {
    preds = get_plus_minus_changes(dataf = df,
                                   y = outcome_ch,
                                   comps = comp_names_BL,
                                   covars = covs,
                                   deltas = seq(-180, 180, by = 1)/(24*60),
                                   comparisons = comparisons,
                                   alpha = 0.05, verbose = FALSE)
  } else if (model == "change") {
    preds = get_plus_minus_changes(dataf = df,
                                   y = outcome_ch,
                                   comps = comp_names_BL,
                                   comps_fup = comp_names_FU,
                                   covars = covs,
                                   deltas = seq(-180, 180, by = 1)/(24*60),
                                   comparisons = comparisons,
                                   alpha = 0.05, verbose = FALSE)
  }

  # GET ESTIMATES
  if (comparisons == "prop-realloc") {
    preds_subset = preds[grep(increase, preds$`comp+`),]
    realloc_time = which(round(preds_subset$delta * 1440) == delta)
    return = preds_subset[realloc_time, c("delta_pred", "ci_lo", "ci_up")]
  } else if (comparisons == "one-v-one") {
    preds_subset = preds[grep(increase, preds$`comp+`),]
    preds_subset = preds_subset[grep(decrease, preds_subset$`comp-`),]
    realloc_time = which(round(preds_subset$delta * 1440) == delta)
    return = preds_subset[realloc_time, c("delta_pred", "ci_lo", "ci_up")]
  }

  # RETURN
  names(return) = c("expected_change", "ci_low", "ci_up")
  cat("\n\n\n")
  cat("Expected change in ",
      ifelse(model == "cross-sectional", yes = outcome_bl, no = outcome_fu),
      "\nfor a ", delta, " min/day reallocation",
      "\nto ", increase, " from ",
      ifelse(comparisons == "prop-realloc", yes = "the rest", no = decrease),
      "\n(", model, " model)\n")
  return(return)
}
