#################################################################################################
#  Numerus Ltd
#
# PURPOSE: Define set up functions (ESS and Gradients)
#
# Built using R 4.4.3
#################################################################################################



# WEIGHT / ESS -----------------------

## Note: Assumes the index variable is the the first position in the X.frame.

## BETA - coefficients
weight <- function(BETA) {
  ws <- exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA))
  ws_scaled <- fply_SUBJNUM* ws / sum(ws)
  return(ws_scaled)
}


## ESS - Objective function
# Calculates the negative effective sample size
ess <- function(BETA){
  ws <- weight(BETA=BETA)
  ESS <- -sum(ws)^2/sum(ws*ws)
  return(ESS)
}


# INEQUALITY FUNCTIONS -----------------------

## Mean constraint for a given variable
# Calculates the distance beyond the acceptable error

g_mean <- function(VAR, WEIGHTS){
  VARMEAN <- paste0(VAR,"MEAN")
  w.mean <- weighted.mean(fply_X[,VAR], WEIGHTS)
  dist <- abs(w.mean - fply_TARGS_ST[VARMEAN])
  constr <- dist - fply_TOLS_ST[VARMEAN]
  return(constr)
}

## Proportion (mean) constraint for a given variable. Note: function identical to mean function
# Calculates the distance beyond the acceptable error
g_prop <- function(VAR, WEIGHTS){
  VARPROP <- paste0(VAR,"PROP")
  w.mean <- weighted.mean(fply_X[,VAR], WEIGHTS)
  dist <- abs(w.mean - fply_TARGS_ST[VARPROP])
  constr <- dist - fply_TOLS_ST[VARPROP]
  return(constr)
}

## SD constraint for a given variable
# Calculates the distance beyond the acceptable error

g_sd <- function(VAR, WEIGHTS){
  VARSD <- paste0(VAR, "SD")
  w_sd <- sqrt(x=Hmisc::wtd.var(fply_X[,VAR], weights=WEIGHTS, normwt=T))
  dist <- abs(w_sd - fply_TARGS_ST[VARSD])
  constr <- dist - fply_TOLS_ST[VARSD]
  return(constr)
}

#COMBINE MEANS -----------------------

## Takes coefficients and calculate difference to targets for means SDs and
## proportions using lists and conditional programming

g <- function(BETA){
  out <- vector(mode="numeric", length=fply_STATNUM)
  WEIGHTS <- weight(BETA)
  ls1 <- lapply(fply_STATS, create_list, BETA=BETA, WEIGHTS=WEIGHTS) ## create list of statistics and interaction state
  subout <- lapply(ls1, g2) ## calculate weighted statistics
  out <- as.numeric(as.data.frame(subout, row.names=NULL, col.names=fply_STATS))
  return(out)
}

create_list <- function(VAR, BETA, WEIGHTS){
  out <- list(BETA=BETA, VAR=as.character(VAR), WEIGHTS = WEIGHTS)
  return(out)
}

g2 <- function(ls2){
  BETA <- ls2$BETA ## polynomial coefficients
  VAR <- ls2$VAR ## matching variables
  WEIGHTS <- ls2$WEIGHTS

  if (str_detect(VAR, "MEAN")){ ## means
    VAR <- str_remove(VAR, "MEAN")
    out <- g_mean(VAR, WEIGHTS = WEIGHTS)
  }

  else if (str_detect(VAR, "PROP")){ ## proportions
    VAR <- str_remove(VAR, "PROP")
    out <- g_prop(VAR = VAR, WEIGHTS = WEIGHTS)
  }

  else if (str_detect(VAR, "SD")){ ## standard deviations
    VAR <- str_remove(VAR, "SD")
    out <- g_sd(VAR = VAR, WEIGHTS = WEIGHTS)
  }
  return(out)
}

# GRADIENT FUNCTIONS -----------------------

## Caluclate gradient of ESS function
grad_ess<- function(BETA){
  ## set up the vector
  len <- length(BETA)
  grad <- vector(mode="numeric", length = len)

  ## Note: gradient wrt intercept is 0 following derivation
  ## grad[1] <- 0 don't include intercept in X

## For loop to generate the other entries for the gradient
  ## Note: The X matrix is called fply_X.
  ## Note: j indexes coefficients and i indexes patients

  for (j in 1:len) {
    grad[j] <- (sum((exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA)))^2) * 2*sum(as.matrix(fply_X[,-1])[,j]*exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA))) *
                  sum(exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA)))  - sum(exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA)))^2 * 2* sum(as.matrix(fply_X[,-1])[,j]*exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA))^2)) / sum((exp(as.matrix(fply_X[,-1]) %*% as.vector(BETA)))^2)^2

  }
  ### Note: We need a minus sign as we are minimizing the negative of the ESS
  return(-grad)
}


## Mean inequality gradient

grad_mean <- function(BETA, VAR){
  ## set up the vector
  len <-length(BETA)
  grad <- vector(mode="numeric", length = len)

  ## Note:these don't depend on k
  exp_X1Bs  <- exp(as.matrix(fply_X[,-1]) %*% BETA)  ## column vector
  exp_X1Bs_sum <- sum(exp_X1Bs)     ##scalar
  exp_X1Bs_sum_sq <- exp_X1Bs_sum^2 ## scalar

  for (k in 1:len){
    grad[k] <- sum(fply_X[,VAR] * exp_X1Bs*(fply_X[,k+1] * exp_X1Bs_sum - sum(fply_X[,k+1] * exp_X1Bs))) / exp_X1Bs_sum_sq
  }

  # Calculate the weighted mean
  wmean <- weighted.mean(fply_X[,VAR], weight(BETA=BETA))

  if (wmean < fply_TARGS_ST[paste0(VAR, "MEAN")]){
    grad <- -1*grad
  }
  return(grad)
}

## Proportion (mean) inequality gradient Note: function identical to mean function (above)

grad_prop <- function(BETA, VAR){
  ## set up the vector
  len <-length(BETA)
  grad <- vector(mode="numeric", length = len)


  ## Note: these don't depend on k

  exp_X1Bs  <- exp(as.matrix(fply_X[,-1]) %*% BETA) ## column vector
  exp_X1Bs_sum <- sum(exp_X1Bs)     ## scalar
  exp_X1Bs_sum_sq <- exp_X1Bs_sum^2 ## scalar

  for (k in 1:len){
    grad[k] <- sum(fply_X[,VAR] * exp_X1Bs*(fply_X[,k+1] * exp_X1Bs_sum - sum(fply_X[,k+1] * exp_X1Bs))) / exp_X1Bs_sum_sq
  }

  # calculate the weighted mean
  wmean <- weighted.mean(fply_X[,VAR], weight(BETA = BETA))

  if (wmean < fply_TARGS_ST[paste0(VAR, "PROP")]){
    grad <- -1*grad
  }
  return(grad)
}

## SD inequality gradient

grad_sd <- function(BETA, VAR){
  ## set up the vector
  len <-length(BETA)
  pre_grad <- vector(mode="numeric", length = len)

  exp_X1Bs  <- exp(as.matrix(fply_X[,-1]) %*% BETA) ## column vector
  exp_X1Bs_sum <- sum(exp_X1Bs)     ## scalar
  exp_X1Bs_sum_sq <- exp_X1Bs_sum^2 ## scalar


  weights <- weight(BETA)
  samp_sd <- sqrt(wtd.var(fply_X[,VAR], weights, normwt=T)) ## weighted SD

  ## calculate values independent of k
  coef <- fply_SUBJNUM/(fply_SUBJNUM^2 - sum(weights^2))
  coef_term <- sum(weights*(fply_X[,VAR] - weighted.mean(fply_X[,VAR], weights))^2) * 2*fply_SUBJNUM/(fply_SUBJNUM^2 - sum(weights^2))

  for (k in 1:len){
    pre_grad[k] <- sum((exp_X1Bs * (fply_X[,k+1] * exp_X1Bs_sum - sum(fply_X[,k+1]*exp_X1Bs))/exp_X1Bs_sum_sq) * fply_SUBJNUM * (fply_X[,VAR] - weighted.mean(fply_X[,VAR], weights))^2-
                         2 * weights * (fply_X[,VAR] - weighted.mean(fply_X[,VAR], weights)) * sum((exp_X1Bs * (fply_X[,k+1] * exp_X1Bs_sum - sum(fply_X[,k+1]*exp_X1Bs))/exp_X1Bs_sum_sq) * fply_X[,VAR])) +
      coef_term *  sum((exp_X1Bs * (fply_X[,k+1] * exp_X1Bs_sum - sum(fply_X[,k+1]*exp_X1Bs))/exp_X1Bs_sum_sq) * weights)
  }

  grad <- 1/(2*samp_sd) * coef * pre_grad

  if (samp_sd < fply_TARGS_ST[paste0(VAR, "SD")]){
    grad <- -1*grad
  }

  return(grad)
}

# COMBINE GRADIENTS -----------------------

grad_g <- function(BETA){

  out <- matrix(ncol = fply_COEFNUM, nrow=fply_STATNUM)

  for (i in 1:fply_STATNUM){

    if (str_detect(fply_STATS[i], "MEAN")){
      VAR <- str_remove(fply_STATS[i], "MEAN")
      out[i,] <- grad_mean(BETA = BETA, VAR = VAR)
    }

    else if (str_detect(fply_STATS[i], "PROP")){
      VAR <- str_remove(fply_STATS[i], "PROP")
      out[i,] <- grad_prop(BETA = BETA, VAR = VAR)
    }

    else if (str_detect(fply_STATS[i], "SD")){
      VAR <- str_remove(fply_STATS[i], "SD")
      out[i,] <- grad_sd(BETA = BETA, VAR = VAR)
    }
  }

  return(out)
}

## Add the g_minus versions for certain algorithms

g_minus <- function(BETA){
  out <- g(BETA = BETA)
  out2 <- -1 * out
  return(out2)
}

grad_g_minus <- function(BETA){
  out <- grad_g(BETA = BETA)
  out2 <- -1 * out
  return(out2)
}

## END OF FILE ##

