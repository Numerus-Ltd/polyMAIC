#' Performs polyMAIC.
#'
#' This function takes IPD, targets, tolerances and matching variables and performs polyMAIC.
#'
#'
#' @param IPD data.frame containing IPD from index study containing matching variables
#' @param VARSTAT named character vector of variable (names) to match on with statistics to match on
#' @param TARGS named numeric vector of targets
#' @param TOLS named numeric vector of tolerances
#' @param TYPE  named character vector of data type
#' @param COEFBOUND Parameter space limits for polynomial parameter space (e.g., (-5000, 5000)). Default value 5000.
#' @param ID Unique identifier
#'
#' @return Object of class 'polymaic'
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
#' ## Three continuous and one binary variable
#'TARGS_C    <- c(AGEMEAN = 49.5, AGESD = 7.2, SEXFPROP = 0.683, PAINMEAN = 7.5, PAINSD = 1.5, TIMEDMEAN = 6.1, TIMEDSD=3.8)
#'TOLS_CN    <- c(AGEMEAN = 0.005, AGESD = 0.005, SEXFPROP = 0.0005, PAINMEAN = 0.005, PAINSD = 0.005, TIMEDMEAN = 0.01, TIMEDSD = 0.01)
#'TOLS_CW    <- c(AGEMEAN = 1, AGESD = 0.5, SEXFPROP = 0.01, PAINMEAN = 0.2, PAINSD = 0.1, TIMEDMEAN = 0.2, TIMEDSD = 0.2)
#'VARSTAT_C  <- c(AGE = "MEAN AND SD", SEXF = "MEAN", PAIN = "MEAN AND SD", TIMED = "MEAN AND SD")
#'TYPE_C     <- c(AGE = "CON", SEXF = "BIN", PAIN = "CON", TIMED="CON")
#'
#'polymaic_cn <- polymaic(IPD=ScenarioC, VARSTAT=VARSTAT_C, TARGS=TARGS_C, TOLS=TOLS_CN, TYPE=TYPE_C, ID="Rnum")
#'summary(polymaic_an)
#'
#' @export
#'

polymaic <- function(IPD, VARSTAT, TARGS, TOLS, TYPE, COEFBOUND=5000, ID="USUBJID"){

  ## DATA SET UP -----------------------
  ## chanege to upper case
  VARSTAT <- toupper(VARSTAT)
  TYPE  <- toupper(TYPE)

  fply_IPD     <<- IPD
  fply_TARGS   <<- TARGS
  fply_VARS    <<- names(VARSTAT)    ## extract names of variables to match on
  fply_VARNUM  <<- length(fply_VARS)      ## extract the number of variables
  fply_STATS   <<- names(fply_TARGS)      ## extract names of statistics
  fply_STATNUM <<- length(fply_TARGS)     ## number of statistics to match on
  fply_COEFNUM <<- sum(TYPE=="CON") * 4 + sum(TYPE=="BIN")
  fply_IPD2    <-  fply_IPD[,c(ID,fply_VARS)]  ## extract variables of interest plus USUBJID or id variable
  fply_SUBJNUM <<- nrow(fply_IPD)

  ## simple parameter checks

  for(i in 1:length(TYPE)){
    if(TYPE[i] == "BIN" & VARSTAT[i] == "MEAN AND SD"){
      stop("Matching standard deviation on a binary variable.")
    }
  }

  ## check all matching variables are numeric
  if(!all(sapply(IPD[,fply_VARS], is.numeric))){
    stop("Non-numeric mathing variable detected. Ensure all matching variables are numeric in IPD.")
  }

  ## create standardised data frame for analysis
  ## merge on with USUBJID to ensure records aren't mixed up

  fply_X <<- data.frame(USUBJID = fply_IPD2[,ID])

  ## set up IPD with fourth order variables
  for (i in 1:fply_VARNUM) {
    X_TEMP <- fply_IPD2[,c(ID,fply_VARS[i])]
    X_TEMP[,fply_VARS[i]] <- (fply_IPD[,fply_VARS[i]] - min(fply_IPD[,fply_VARS[i]]))/(max(fply_IPD[,fply_VARS[i]]) - min(fply_IPD[,fply_VARS[i]]))

    ## above binary variables are also standardized in case the use a coding other than 0 1

    if(TYPE[fply_VARS[i]] == "CON"){
      X_TEMP[,paste0(fply_VARS[i],"_2")] <- X_TEMP[,fply_VARS[i]]**2
      X_TEMP[,paste0(fply_VARS[i],"_3")] <- X_TEMP[,fply_VARS[i]]**3
      X_TEMP[,paste0(fply_VARS[i],"_4")] <- X_TEMP[,fply_VARS[i]]**4
    }
    fply_X <<- merge(fply_X, X_TEMP, by=ID)
    ## delete X_TEMP
    rm(X_TEMP)
  }

  ## STANDARDISE TARGS and TOLS -----------------------
  fply_TARGS_ST <<- vector(mode="numeric", length = fply_STATNUM)
  fply_TOLS_ST  <<- vector(mode="numeric", length = fply_STATNUM)

  for (i in 1:fply_STATNUM){

    if (str_detect(fply_STATS[i], "MEAN")){
      VAR <- str_remove(fply_STATS[i], "MEAN")
      ## transformation for mean statistic E((X-a)/(b-a)) = (E(X) - a)/(b-a)
      fply_TARGS_ST[i] <<- (fply_TARGS[i] - min(fply_IPD[,VAR]))/(max(fply_IPD[,VAR]) - min(fply_IPD[,VAR]))
      fply_TOLS_ST[i]  <<- TOLS[i]/(max(fply_IPD[,VAR]) - min(fply_IPD[,VAR]))
    }

    else if (str_detect(fply_STATS[i], "PROP")){
      VAR <- str_remove(fply_STATS[i], "PROP")
      ## transformation for mean statistic E((X-a)/(b-a)) = (E(X) - a)/(b-a)
      fply_TARGS_ST[i] <<- (fply_TARGS[i] - min(fply_IPD[,VAR]))/(max(fply_IPD[,VAR]) - min(fply_IPD[,VAR]))
      fply_TOLS_ST[i]  <<- TOLS[i]/(max(fply_IPD[,VAR]) - min(fply_IPD[,VAR]))
    }

    else if (str_detect(fply_STATS[i], "SD")){
      VAR <- str_remove(fply_STATS[i], "SD")
      ## transformation for SD statistic SD((X -a)/(b-a)) = SD(X)/(b-a)
      fply_TARGS_ST[i] <<- fply_TARGS[i]/(max(fply_IPD[,VAR]) - min(fply_IPD[,VAR]))
      fply_TOLS_ST[i]  <<- TOLS[i]/(max(fply_IPD[,VAR]) - min(fply_IPD[,VAR]))
    }
  }

  names(fply_TARGS_ST) <<- names(fply_TARGS)
  names(fply_TOLS_ST)  <<- names(TOLS)

  ## PERFORM OPTIMISATION ----------
  LOWER <- rep(-COEFBOUND, fply_COEFNUM) ## lower bound for parameter space
  UPPER <- rep( COEFBOUND, fply_COEFNUM) ## upper bound for parameter space
  X0    <- rep(0,    fply_COEFNUM)  ## initial conditions

  tic() ## start timing
  loc_result <- nloptr::slsqp(x0 = X0,
                              fn = ess,
                              gr = grad_ess,
                              lower = LOWER,
                              upper = UPPER,
                              hin = g,
                              hinjac = grad_g,
                              control = list( "xtol_rel"= 1.0e-14,
                                              "ftol_rel"= 1.0e-14,
                                              "maxeval"= 5000),
                              deprecatedBehavior = FALSE) ## call slsqp algorithm for nloptr pacakge

  loc_time    <- toc() ## end timing and save for analytics
  loc_status  <- loc_result$convergence ## convergence status code (produced by nloptr package, see nloptr documentation)
  loc_message <- loc_result$message     ## convergence status message (produced by nloptr package, see nloptr documentation)

  ## Calculate the weights and produce polymaic object.

  WEIGHTS <- weight_out(beta=loc_result$par, id=ID) ## calculate the weights based on the coefficients. Produces data frame with two columns
                                                    ## ID and assigned patient weight

  ESS       <- -1 * loc_result$value ## calculate ESS, note optimiser minimise the negative of the ESS, identical and maximising the positive
  STAT_OUT  <- statistic_out(weights = WEIGHTS, id=ID) ## produce the weighted statistics
  ERROR_OUT <- error_out(weights = WEIGHTS, id=ID) ## produce the errors (residuals) - difference between weighted statistics and targets
  COEF      <- loc_result$par ## save the coefficents

  out <- list(ESS        = ESS,
              weights    = WEIGHTS,
              statistic  = STAT_OUT,
              target     = fply_TARGS,
              error      = ERROR_OUT,
              tolerance  = TOLS,
              coef       = COEF,
              loc_diagnostic  = list(TIME = loc_time$callback_msg, STATUS=loc_status, MESSAGE = loc_message)) ## output desired infomration as a list

  class(out) <- "polymaic" ## set the class of the object

  ## remove global objects
  rm(fply_VARS, fply_VARNUM, fply_STATS, fply_STATNUM, fply_COEFNUM, fply_X, fply_SUBJNUM, fply_TARGS_ST, fply_TOLS_ST, fply_IPD, fply_TARGS, pos=1)

  return(out)
}

## END OF FILE ##
