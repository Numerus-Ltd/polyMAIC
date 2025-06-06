#################################################################################################
#  Numerus Ltd
#
# PURPOSE: Defines functions required for the polymaic class and functions to
#          assess performance on polymaic algorithm
#
# Built using R 4.4.3
#################################################################################################

## Output a data frame in the subject number and assigned patient weight

weight_out <- function(beta, id){
  out <- fply_X
  ws <- exp(as.matrix(fply_X[,-1]) %*% as.vector(beta))
  out$WEIGHTS <- fply_SUBJNUM * ws / sum(ws)
  out <- out[,c(id, "WEIGHTS")]
  return(out)
}

##  Output weighted statistic

statistic_out <- function(weights, id){
  DATA <- merge(fply_IPD, weights, by=id)
  out <- vector(mode="numeric", length=fply_STATNUM)
  names(out) <- fply_STATS

  for (i in 1:fply_STATNUM){

    if (str_detect(fply_STATS[i], "MEAN")){

      VAR <- str_remove(fply_STATS[i], "MEAN")
      out[i] <- weighted.mean(x=DATA[,VAR], w=DATA[,"WEIGHTS"])
    }

    else if (str_detect(fply_STATS[i], "PROP")){

      VAR <- str_remove(fply_STATS[i], "PROP")
      out[i] <- weighted.mean(x=DATA[,VAR], w=DATA[,"WEIGHTS"])
    }

    else if (str_detect(fply_STATS[i], "SD")){

      VAR <- str_remove(fply_STATS[i], "SD")
      out[i] <- sqrt(x=wtd.var(x=DATA[,VAR], weights=DATA[,"WEIGHTS"], normwt=T))

    }
  }
  return(out)
}


## Calculate the error between the weighted statistics and the target data

error_out <- function(weights, id){
  WEIGHTED_STATS <- statistic_out(weights=weights, id=id)
  out <- WEIGHTED_STATS - fply_TARGS
  return(out)
}


## Entropy

entropy <- function(weights){
  out <- sum(weights*log(fply_SUBJNUM*weights))
  return(out)
}

## Calculate the proportion of weights greater than 5% of N

ipw_prop <- function(weights, prop){
  out <- sum(weights > (fply_SUBJNUM*prop))/fply_SUBJNUM
  return(out)
}

## Calculate the proportion of weights greater than a fixed number

ipw_fixed <- function(weights, fixed){
  out <- sum(weights > fixed) / fply_SUBJNUM
  return(out)
}


# summary <- function(object){
#   UseMethod("summary")
# }

#' Summarises polyMAIC
#'
#' This function takes an object of class 'polymaic' and displays the summary statistics associated within the console.
#'
#'
#' @param object Object of class polymaic
#' @return A summary of the polymaic object to console
#'
#' @examples
#' ## One continuous and one binary variable
#'TARGS_A    <- c(AGEMEAN = 55, AGESD = 6, SEXFPROP = 0.65)          ## targets
#'TOLS_AN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005) ## tolerances (narrow)
#'TOLS_AW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01)         ## tolerances (wide)
#'VARSTAT_A  <- c(AGE = "MEAN AND SD", SEXF = "MEAN")                ## Variables and measures to match on
#'TYPE_A     <- c(AGE = "CON", SEXF = "BIN")                         ## Type of variables to match on
#'
#'polymaic_an <- polymaic(IPD=ScenarioA, VARSTAT=VARSTAT_A, TARGS=TARGS_A, TOLS=TOLS_AN, TYPE=TYPE_A, ID="Rnum")
#'summary(polymaic_an)
#'
#'@export
#'
summary.polymaic <- function(object){
  cat("############### SUMMARY OF POLYMAIC OBJECT ###############\n")

  ## PRINT EFFECTIVE SAMPLE SIZE
  cat("\nEffective Sample Size (ESS): ", as.character(round(object$ESS,3)), "\n")

  ## PRINT WEIGHTED STATISTICS
  cat("\nWEIGHTED STATISTICS\n")
  print(object$statistic)

  ## PRINT TARGET AGGREGATE STATISTICS
  cat("\nTARGETS \n")
  print(object$target)

  ## PRINT RESIDUAL ERRORS
  cat("\nRESIDUAL ERRORS\n")
  print(round(object$error,4))

  ## PRINT TOLERANCES
  cat("\nUSER-SPECIFIED TOLERANCES\n")
  print(object$tolerance)

  cat("#############################################################\n")
}


#' Histogram of IPWs
#'
#' This function uses the ggplot2 package to output a histogram of the re-scaled IPW
#'
#'
#' @param object Object of class polymaic
#'
#' @return A histogram of the re-scaled IPWs
#'
#' @examples
#' ## One continuous and one binary variable
#'TARGS_A    <- c(AGEMEAN = 55, AGESD = 6, SEXFPROP = 0.65)          ## targets
#'TOLS_AN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005) ## tolerances (narrow)
#'TOLS_AW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01)         ## tolerances (wide)
#'VARSTAT_A  <- c(AGE = "MEAN AND SD", SEXF = "MEAN")                ## Variables and measures to match on
#'TYPE_A     <- c(AGE = "CON", SEXF = "BIN")                         ## Type of variables to match on
#'
#'polymaic_an <- polymaic(IPD=ScenarioA, VARSTAT=VARSTAT_A, TARGS=TARGS_A, TOLS=TOLS_AN, TYPE=TYPE_A, ID="Rnum")
#'polyplot(polymaic_an)
#'
#' @export
#'

polyplot <- function(object){

  df <- object$weights
  colnames(df) <- c("Rnum", "WEIGHTS")


  # Create annotation text
  # txt <- paste0("ESS: ",      round(object$ESS,3),
  #               "\nMean: ",   round(mean(object$weights$WEIGHTS),2),
  #               "\nMin: ",    round(min(object$weights$WEIGHTS),2),
  #               "\nMedian: ", round(median(object$weights$WEIGHTS),2),
  #               "\nMax: ",    round(max(object$weights$WEIGHTS),2))
  #

  # Take min and max to work out binwidth #
  xmin <- floor(min(df[,2], na.rm=T))  # round down minimum observed value
  xmax <- ceiling(max(df[,2], na.rm=T))  # round up maximum observed value



  plot <- ggplot(data=df, aes(WEIGHTS)) +
          geom_histogram(binwidth = (xmax - xmin) / 25 , colour = '#3F5B6E', fill = '#1172C0')  +
          ggtitle("Histogram of Re-Secaled Individual Patient Weights") +
          xlab("Re-Scaled Weights") +
          ylab("Frequency") +
          theme_classic() #+
          #scale_y_continuous(limits = c(0, round(max(ggplot_build(plot)$data[[1]]$count), digits = -1)  )) +
          #scale_x_continuous(limits = c(0, max(ggplot_build(plot)$data[[1]]$x))) +
          #annotate("text",x=xmax - xmax* 0.1, y= max(ggplot_build(plot)$data[[1]]$count) - xmax, label=txt)

  return(plot)
}

## END OF FILE ##
