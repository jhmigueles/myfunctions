#' zComposite
#' @description Calculates z composite scores of a set of variables, possibility to stratify by sex implemented
#'
#' @param x Names of the variables to be included in the composite score
#' @param x_fu Names of the variables to be included in the composite score at follow up (it would calculate a z score of the change)
#' @param id Participant identifier variable
#' @param sex_strata Logical. Whether to stratify the composite score by sex
#' @param sex Name of the sex variable in the data frame.
#' @param age Name of the age variable in the dataset. Only used if reference values are provided in \code{refs}
#' @param refs Reference values as, for example, in https://doi.org/10.1016/j.atherosclerosis.2018.10.003 (important the reference values data frame presents the same colnames as the variable names in df)
#' @param df Data frame with variables to calculate the composite score
#'
#' @return Data frame with composite z score calculated
#' @export
#'
#' @examples
zComposite = function(x = c(), x_fu = c(), id = c(), sex_strata = FALSE,
                    sex = "Sex", age = c("Age_I", "Age_II"), refs=NULL, df) {
  # subset data
  df_all = df[, c(id, x, x_fu)]
  if (sex_strata) df_all = df[, c(id, x, x_fu, sex)]
  if (is.matrix(refs) | is.data.frame(refs)) df_all = df[, c(id, x, x_fu, sex, age)]
  # keep only complete cases for z scores
  df = df_all[complete.cases(df_all), ]
  # refs mean reference values for z scores (eg, paper Lars)
  if (is.null(refs)) {
    if (sex_strata) {
      b = which(df[,sex] == 0)
      g = which(df[,sex] == 1)
      if (length(x_fu) > 0) {
        for (i in 1:length(x_fu)) {
          df[b, x_fu[i]] = (df[b, x_fu[i]] - mean(df[b, x[i]])) / sd(df[b, x[i]])
          df[g, x_fu[i]] = (df[g, x_fu[i]] - mean(df[g, x[i]])) / sd(df[g, x[i]])
        }
      }
      df[b, x] = scale(df[b, x])
      df[g, x] = scale(df[g, x])
      cluster_I  = rowMeans(df[, x])
      if (length(x_fu) > 0) cluster_II  = rowMeans(df[, x_fu])
    } else {
      if (length(x_fu) > 0) {
        for (i in 1:length(x_fu)) df[, x_fu[i]] = (df[, x_fu[i]] - mean(df[, x[i]])) / sd(df[, x[i]])
      }
      df[, x] = scale(df[, x])
      cluster_I  = rowMeans(df[, x])
      if (length(x_fu) > 0) cluster_II  = rowMeans(df[, x_fu])
    }
  } else { # if there are reference values...
    df[,age] = floor(df[,age])
    recoding_sex = ifelse(df[,sex] == 0, "b", "g")
    recoding_I.a = floor(df[,age[1]])
    recoding_I = paste0(recoding_sex, recoding_I.a)
    if (length(x_fu) > 0) recoding_II.a = floor(df[,age[2]])
    if (length(x_fu) > 0) recoding_II = paste0(recoding_sex, recoding_II.a)
    recoding_sd = ifelse(recoding_sex == "b", "bsd", "gsd")

    for (i in 1:length(x)) {
      for (child in 1:nrow(df)) {
        mean_ref_I = refs[recoding_I[child], x[i]]
        if (length(x_fu) > 0) mean_ref_II = refs[recoding_II[child], x[i]]
        if (grepl("reverted", x[i])) {
          mean_ref_I = mean_ref_I * -1
          if (length(x_fu) > 0) mean_ref_II = mean_ref_II * -1
        }
        sd_ref = refs[recoding_sd[child], x[i]]
        df[child, x[i]] = (df[child, x[i]] - mean_ref_I) / sd_ref
        if (length(x_fu) > 0) df[child, x_fu[i]] = (df[child, x_fu[i]] - mean_ref_II) / sd_ref
      }
    }
    cluster_I  = rowMeans(df[, x])
    if (length(x_fu) > 0) cluster_II  = rowMeans(df[, x_fu])
  }
  if (length(x_fu) > 0) return(data.frame(ID = df[,id], zCluster_I = scale(cluster_I), zCluster_II = scale(cluster_II)))
  if (length(x_fu) == 0) return(data.frame(ID = df[,id], zCluster_I = scale(cluster_I)))
}
