% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logisticReg_coda.R
\name{logisticReg_coda}
\alias{logisticReg_coda}
\title{Basic function for cross-sectional logistic regression CoDA models}
\usage{
logisticReg_coda(data, compo, outcome, covariates = c(), scale = FALSE)
}
\arguments{
\item{data}{data frame or matrix with variables of interest}

\item{compo}{String with colnames of composition variables or matrix/dataframe with these variables}

\item{outcome}{String with colname of outcome or numeric vector with the outcome}

\item{covariates}{String with names of covariates (optional) or matrix/dataframe with these variables}

\item{scale}{Logical. Whether to standardize results (default = FALSE)}
}
\value{
Logistic regression model (glm class)
}
\description{
Basic function for cross-sectional logistic regression CoDA models
}
