% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_define.R
\name{polyplot}
\alias{polyplot}
\title{Histogram of IPWs}
\usage{
polyplot(object)
}
\arguments{
\item{object}{Object of class polymaic}
}
\value{
A histogram of the re-scaled IPWs
}
\description{
This function uses the ggplot2 package to output a histogram of the re-scaled IPW
}
\examples{
## One continuous and one binary variable
TARGS_A    <- c(AGEMEAN = 55, AGESD = 6, SEXFPROP = 0.65)          ## targets
TOLS_AN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005) ## tolerances (narrow)
TOLS_AW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01)         ## tolerances (wide)
VARSTAT_A  <- c(AGE = "MEAN AND SD", SEXF = "MEAN")                ## Variables and measures to match on
TYPE_A     <- c(AGE = "CON", SEXF = "BIN")                         ## Type of variables to match on

polymaic_an <- polymaic(IPD=ScenarioA, VARSTAT=VARSTAT_A, TARGS=TARGS_A, TOLS=TOLS_AN, TYPE=TYPE_A, ID="Rnum")
polyplot(polymaic_an)

}
