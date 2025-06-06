% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_define.R
\name{summary.polymaic}
\alias{summary.polymaic}
\title{Summarises polyMAIC}
\usage{
\method{summary}{polymaic}(object)
}
\arguments{
\item{object}{Object of class polymaic}
}
\value{
A summary of the polymaic object to console
}
\description{
This function takes an object of class 'polymaic' and displays the summary statistics associated within the console.
}
\examples{
## One continuous and one binary variable
TARGS_A    <- c(AGEMEAN = 55, AGESD = 6, SEXFPROP = 0.65)          ## targets
TOLS_AN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005) ## tolerances (narrow)
TOLS_AW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01)         ## tolerances (wide)
VARSTAT_A  <- c(AGE = "MEAN AND SD", SEXF = "MEAN")                ## Variables and measures to match on
TYPE_A     <- c(AGE = "CON", SEXF = "BIN")                         ## Type of variables to match on

polymaic_an <- polymaic(IPD=ScenarioA, VARSTAT=VARSTAT_A, TARGS=TARGS_A, TOLS=TOLS_AN, TYPE=TYPE_A, ID="Rnum")
summary(polymaic_an)

}
