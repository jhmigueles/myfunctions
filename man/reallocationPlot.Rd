% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reallocationPlot.R
\name{reallocationPlot}
\alias{reallocationPlot}
\title{reallocationPlot}
\usage{
reallocationPlot(
  dc_obj = c(),
  data = c(),
  comps = c(),
  comps_fup = NULL,
  comps.names = c(),
  outcome = c(),
  covs = c(),
  comparisons = c("one-v-one", "prop-realloc")[2],
  increase = c(),
  decrease = c(),
  xlab = c(),
  ylab = c(),
  xlim = c(),
  ylim = c(),
  ReallocationLimits = c(-60, 60),
  xaxis.by = 5,
  ribbon = TRUE,
  font = "Times New Roman",
  main = NULL,
  col = "black",
  alpha = 0.2,
  analysis_type = c("cross-sectional", "longitudinal")[1],
  total = NULL
)
}
\arguments{
\item{dc_obj}{Object derived from lm_coda or lm_coda_long}

\item{comps}{String with the name of the composition variables in dataf}

\item{comps_fup}{String with the name of the follow-up composition variables in dataf (optional)}

\item{comps.names}{Names of the composition variables to appear in the plot (it will use \code{comps} if not provided)}

\item{outcome}{String with the name of the outcome in dataf}

\item{covs}{String with the name of the covariates}

\item{comparisons}{Two choices ("one-v-one" or "prop-realloc")}

\item{increase}{What variable in the composition should increase in the visualization (replacement)}

\item{decrease}{What variable/s in the composition should decrease in the visualization (replaced)}

\item{xlab}{X axis label}

\item{ylab}{Y axis label}

\item{xlim}{Limits of the X axis (time reallocation limits)}

\item{ylim}{Limits for the Y axis (autogenerated if not provided)}

\item{ReallocationLimits}{limits for the reallocation time in minutes (vector of 2 numbers)}

\item{xaxis.by}{X axis ticks at every \code{xaxis.by} minutes}

\item{ribbon}{Logical. Whether to plot confidence intervals of the predictions}

\item{font}{Font to use (default = "Times New Roman")}

\item{main}{Title of the plot}

\item{col}{Color for lines and confidence intervals}

\item{alpha}{Alpha for transparency in confidence intervals}

\item{analysis_type}{One of "cross-sectiona", "longitudinal", "change".}

\item{total}{Total time to constrain the composition to.}
}
\value{
reallocation plot with predictions
}
\description{
Visualization of predictions from CoDA model
}
\details{
If follow-up composition is provided, it will run prospective and change models/visualizations
}
