#' explore
#' @description Calculates descriptives and visualizations of distributions of variables in dataset
#'
#' @param dat Dataset
#' @param idvar Subject identification variable
#' @param outputdir Directory to store the output
#'
#' @return Nothing, it saves a csv and pdf file with descriptives and visualizations
#' @export
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom car qqPlot Boxplot
explore = function(dat = c(), idvar = c(), outputdir = "./"){

  if (is.character(dat)) {
    if(grep(".csv", dat, fixed = TRUE) == 1) {
      dat = read.csv(file=dat)
    } else {print("No .csv file provided")}

  } else {
    dat = dat
  }


  #PACKAGES
  library(car)

  #Cleaning potential extra rows
  dat.tmp = dat[which(is.na(dat[,idvar])==F),]

  #Subsetting only numeric variables
  # numstmp=c()
  # for (i in 1:ncol(dat.tmp)) numstmp[i]=length(unique(dat.tmp[,i]))
  # nums = which(numstmp > 4)
  # 
  # dat = dat.tmp[,nums]

  nums = which(sapply(dat, is.numeric))
  dat = dat[,nums]
  rm(dat.tmp)

  #Output directory
  if(dir.exists(outputdir)==F) dir.create(outputdir)

  # Progess bar
  pb <- txtProgressBar(min = 0, max = ncol(dat)-1, style = 3)


  #QQ-plots
  pdf(file.path(outputdir, "distribution plots.pdf"),
      paper = "a4", width = 7, height = 10)

  for (i in 2:ncol(dat)){
    layout(rbind(c(1, 1), c(2, 3)))
    hist.default(dat[,i], main = colnames(dat[i]))
    qqPlot(dat[,i], envelope = F, ylab = "Observed value")
    Boxplot(dat[,i], ylab = colnames(dat)[i])
    if(!("desc" %in% ls())){
      desc = as.data.frame(matrix(NA, nrow = ncol(dat)-1, ncol = 12))
      colnames(desc) = c("variable","n", "missing", "mean", "sd", "min", "q1", "median", "q3", "max",
                         "outliers", "extreme outliers")
    }

    desc[i-1,] = c(colnames(dat)[i],
                   sum(!is.na(dat[,i])),
                   sum(is.na(dat[,i])),
                   mean(dat[,i], na.rm = T),
                   sd(dat[,i], na.rm = T),
                   quantile(dat[,i], na.rm = T),
                   outliers(dat[,i]),
                   outliers(dat[,i], k = 3))
    setTxtProgressBar(pb, i)
  }
  dev.off()
  write.csv(desc, file.path(outputdir, "descriptives.csv"), row.names = F)
  print(paste(i, "variables explored"))
}
