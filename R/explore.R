#' explore
#' @description Calculates descriptives and visualizations of distributions of variables in dataset
#'
#' @param dat Dataset
#' @param idvar Subject identification variable
#' @param outputdir Directory to store the output
#' @param outputfilename Desired tagname for the excel and the pdf files generated
#'
#' @return Nothing, it saves a csv and pdf file with descriptives and visualizations
#' @export
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom car qqPlot Boxplot
#' @importFrom openxlsx write.xlsx

explore = function(dat = c(), idvar = c(), outputdir = "./", outputfilename = "data"){

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
  dat.tmp = dat[which(is.na(dat[,idvar]) == F),]

  nums = which(sapply(dat, is.numeric))
  dat = dat[,nums]
  rm(dat.tmp)

  #Output directory
  if (dir.exists(outputdir) == F) dir.create(outputdir)

  # Progess bar
  pb <- txtProgressBar(min = 0, max = ncol(dat) - 1, style = 3)


  #QQ-plots
  filename_pdf = ifelse(length(outputfilename > 0), paste0("_", outputfilename, ".pdf"), "")
  pdf(file.path(outputdir, paste0("distribution", filename_pdf)),
      paper = "a4", width = 7, height = 10)

  for (i in 2:ncol(dat)) {
    layout(rbind(c(1, 1), c(2, 3)))
    hist.default(dat[,i], main = colnames(dat[i]))
    car::qqPlot(dat[,i], envelope = F, ylab = "Observed value")
    car::Boxplot(dat[,i], ylab = colnames(dat)[i])
    if (!("desc" %in% ls())) {
      desc = as.data.frame(matrix(NA, nrow = ncol(dat) - 1, ncol = 12))
      colnames(desc) = c("variable","n", "missing", "mean", "sd", "min", "q1", "median", "q3", "max",
                         "outliers", "extreme outliers")
    }

    desc[i - 1,] = c(colnames(dat)[i],
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

  # save to excel file
  filename_xlsx = ifelse(length(outputfilename > 0), paste0("_", outputfilename, ".xlsx"), "")
  openxlsx::write.xlsx(x = desc, file = file.path(outputdir, filename_xlsx),
                       asTable = T, overwrite = TRUE,
                       creator = "Jairo Hidalgo Migueles", sheetName = "Descriptives",
                       firstRow = TRUE, firstCol = TRUE, colWidths = "auto", na.string = " ")
  
  print(paste(i, "variables explored"))
}
