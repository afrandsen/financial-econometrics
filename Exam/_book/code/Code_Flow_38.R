# EXERCISE 1 #
# Author: Andreas Kracht Frandsen
# Date: 2020-01-16
# For documentation please read the Computational Part of the main document.

# Load required packages.
library(quantmod)
library(Rsolnp)

# Obtains the relevant data from the VIX Index. We are only gonna use VIX.Adjusted.
getSymbols(Symbols = "^VIX",
           from    = '2010-01-01',
           to      = '2019-01-01',
           src     = 'yahoo')

# First 5 observation of VIX. 
head(VIX)

# Plot of VIX.Adjusted
plot(y    = VIX$VIX.Adjusted,
     x    = index(VIX),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'VIX',
     main = 'Volatility Index')

## GAS-GAMMA ##

# The following function is the filter for the GAS-GAMMA model.
# Input:  Par (Vector), the parameters in the model, which are omega, alpha, beta and a. 
#         Y (Vector), the returns of the relevant series.
# Output: Output (List), the means (Vector) and log likelihood (Double).
GASGAMMA_Filter <- function(Par, Y) {
  iT = length(Y)
  Mu = numeric(iT)
  
  Omega = Par[1]
  Alpha = Par[2]
  Beta  = Par[3]
  a     = Par[4]
  
  # The first mu is initialized as the unconditional expectation of mu.
  Mu[1] = Omega/(1 - Beta)
  
  # In case mu gets negative or hits zero we adjust it.
  if (Mu[1] < 1e-10) {Mu[1] = 1e-10}
  
  # For the rest (T-1) mu's we use the updating equation.
  for (t in 2:iT) {
    Mu[t] = Omega + Alpha * ((sqrt(a)*(Y[t-1] - Mu[t-1]))/Mu[t-1]) + Beta * Mu[t-1]
    
    # In case mu gets negative or hits zero we adjust it.
    if (Mu[t] < 1e-10) {Mu[t] = 1e-10}
  }
  
  # The log likelihood is computed as shown in the theoretical part.
  LLK = iT * (a * log(a) - lgamma(a)) + (a - 1) * sum(log(Y)) - a * sum(Y / Mu + log(Mu))
  
  Output = list()
  
  Output[["Mu"]]   = Mu
  Output[["LLK"]]  = LLK
  
  return(Output)
  
}

# The following function evaluates the negative log likelihood for further use in the optimization proces.
# Input:  Par (Vector), the parameters in the model, which are omega, alpha, beta and a.
#         Y (Vector), the returns of the relevant series.
# Output: NLL (Double), the negative log likelihood.
NegLogLikelihood <- function(Par, Y) {
  
  Filter = GASGAMMA_Filter(Par, Y)
  
  NLL = -Filter[["LLK"]]
  
  # In case the negative log likelihood isn't finite we adjust it to a large value.
  if (!is.finite(NLL)) {
    NLL = 1e5
  }
  
  return(NLL)
}


# The following function estimates the GAS-GAMMA model by first finding maximum likelihood estimates of our parameters.
# Input:  Y (Vector), the returns of the relevant series.
# Output: Output (List), the optimized parameters (Vector), the BIC (Double) and the filtered values of GASGAMMA_Filter.
Estimate_GASGAMMA <- function(Y) {
  
  # Use the gosolnp from the Rsolnp package to optimize the negative log likelihood. With random initialized starting values.
  optimiser = gosolnp(fun    = NegLogLikelihood,
                      Y      = Y,
                      n.sim  = 20,
                      LB     = c(-0.5, 0.001, 0.01, 0.1),
                      UB     = c(0.5, 1.5, 0.999, 300)
  )
  
  Par = optimiser$pars
  LLK = -tail(optimiser$value, n=1)
  
  # Here we run the filter using the optimal parameter values, to obtain the final estimates of mu.
  FilteredValues = GASGAMMA_Filter(Par, Y)
  
  iT = length(Y)
  
  # Computation of Bayesian Information Criterion, using the fact that we estimate four parameters.
  BIC = (log(iT) * 4 - 2 * LLK)
  
  Output = list()
  
  Output[["Par"]]            = Par
  Output[["BIC"]]            = BIC
  Output[["FilteredValues"]] = FilteredValues
  
  return(Output)
}

## GAS-GAMMA-C ##

GASGAMMA_Filter_c <- function(Par, Y) {
  
  iT = length(Y)
  
  Omega = Par[1]
  Beta  = Par[2]
  a     = Par[3]
  
  # The first mu is initialized as the unconditional expectation of mu.
  Mu = Omega / (1 - Beta)
  
  # In case mu gets negative or hits zero we adjust it.
  if (Mu < 1e-10) {Mu = 1e-10}
  
  # The log likelihood is computed as shown in the theoretical part.
  LLK = iT * (a * log(a) - lgamma(a)) + (a - 1) * sum(log(Y)) - a * sum(Y / Mu + log(Mu))
  
  Output = list()
  
  Output[["Mu"]] = Mu
  Output[["LLK"]] = LLK
  
  return(Output)
  
}

# The following function evaluates the negative log likelihood for further use in the optimization proces.
# Input:  Par (Vector), the parameters in the model, which are omega, alpha, beta and a.
#         Y (Vector), the returns of the relevant series.
# Output: NLL (Double), the negative log likelihood.
NegLogLikelihood_c <- function(Par, Y) {
  
  Filter = GASGAMMA_Filter_c(Par, Y)
  
  NLL = -Filter[["LLK"]]
  
  # In case the negative log likelihood isn't finite we adjust it to a large value.
  if (!is.finite(NLL)) {
    NLL = 1e5
  }
  
  return(NLL)
}


# The following function estimates the GAS-GAMMA-C model by first finding maximum likelihood estimates of our parameters.
# Input:  Y (Vector), the returns of the relevant series.
# Output: Output (List), the optimized parameters (Vector), the BIC (Double) and the filtered values of GASGAMMA_Filter_c.
Estimate_GASGAMMA_c <- function(Y) {
  
  # Use the gosolnp from the Rsolnp package to optimize the negative log likelihood. With random initialized starting values.
  optimiser = gosolnp(fun    = NegLogLikelihood_c,
                      Y      = Y,
                      n.sim = 20,
                      LB  = c(-0.5, 0.01, 0.1),
                      UB  = c(0.5, 0.999, 300)
  )
  
  Par  = optimiser$pars
  LLK = -tail(optimiser$value, n=1)
  
  # Here we run the filter using the optimal parameter values, to obtain the final estimates of mu.
  FilteredValues = GASGAMMA_Filter_c(Par, Y)
  
  iT = length(Y)
  
  # Computation of Bayesian Information Criterion, using the fact that we estimate three parameters.
  BIC = (log(iT) * 3 - 2 * LLK)
  
  Output = list()
  
  Output[["Par"]] = Par
  Output[["BIC"]] = BIC
  Output[["FilteredValues"]] = FilteredValues
  
  return(Output)
}

## MEM-GAMMA ##

# The following function is the filter for the MEM-GAMMA model.
# Input:  Par (Vector), the parameters in the model, which are kappa, eta, phi and a. 
#         Y (Vector), the returns of the relevant series.
# Output: Output (List), the means (Vector) and log likelihood (Double).
MEMGAMMA_Filter <- function(Par, Y) {
  
  iT = length(Y)
  Mu = numeric(iT)
  
  Kappa = Par[1]
  Eta   = Par[2]
  Phi   = Par[3]
  a     = Par[4]
  
  # The first mu is initialized as the unconditional expectation of mu.
  Mu[1] = Kappa/(1 - Eta - Phi)
  
  # In case mu gets negative or hits zero we adjust it.
  if (Mu[1] < 1e-10) {Mu[1] = 1e-10}
  
  # For the rest (T-1) mu's we use the updating equation.
  for (t in 2:iT) {
    Mu[t] = Kappa + Eta * Y[t-1] + Phi * Mu[t-1]
    
    # In case mu gets negative or hits zero we adjust it.
    if (Mu[t] < 1e-10) {Mu[t] = 1e-10}
  }
  
  LLK = iT * (a * log(a) - lgamma(a)) + (a - 1) * sum(log(Y)) - a * sum(Y / Mu + log(Mu))
  
  Output = list()
  
  Output[["Mu"]] = Mu
  Output[["LLK"]] = LLK
  
  return(Output)
  
}

# The following function evaluates the negative log likelihood for further use in the optimization proces.
# Input:  Par (Vector), the parameters in the model, which are omega, alpha, beta and a.
#         Y (Vector), the returns of the relevant series.
# Output: NLL (Double), the negative log likelihood.
NegLogLikelihood_MEM <- function(Par, Y) {
  
  Filter = MEMGAMMA_Filter(Par, Y)
  
  NLL = -Filter[["LLK"]]
  
  # In case the negative log likelihood isn't finite we adjust it to a large value.
  if (!is.finite(NLL)) {
    NLL = 1e5
  }
  
  return(NLL)
}


# The following function estimates the MEM-GAMMA model by first finding maximum likelihood estimates of our parameters.
# Input:  Y (Vector), the returns of the relevant series.
# Output: Output (List), the optimized parameters (Vector), the BIC (Double) and the filtered values of GASGAMMA_Filter.
Estimate_MEMGAMMA <- function(Y) {
  
  # Use the gosolnp from the Rsolnp package to optimize the negative log likelihood. With random initialized starting values.
  optimiser = gosolnp(fun    = NegLogLikelihood_MEM,
                      Y      = Y,
                      n.sim  = 20,
                      LB     = c(0.1, 0.01, 0.01, 0.1),
                      UB     = c(10, 0.99, 0.99, 300)
  )
  
  Par  = optimiser$pars
  LLK = -tail(optimiser$value, n=1)
  
  # Here we run the filter using the optimal parameter values, to obtain the final estimates of mu.
  FilteredValues = MEMGAMMA_Filter(Par, Y)
  
  iT = length(Y)
  
  # Computation of Bayesian Information Criterion, using the fact that we estimate four parameters.
  BIC = (log(iT) * 4 - 2 * LLK)
  
  Output = list()
  
  Output[["Par"]]            = Par
  Output[["BIC"]]            = BIC
  Output[["FilteredValues"]] = FilteredValues
  
  return(Output)
}

# Fit all the models using VIX data.
Fit_GAS_GAMMA = Estimate_GASGAMMA(VIX$VIX.Adjusted)
Fit_GAS_GAMMA_c = Estimate_GASGAMMA_c(VIX$VIX.Adjusted)
Fit_MEM_GAMMA = Estimate_MEMGAMMA(VIX$VIX.Adjusted)

# PLOT OF MEM.
# Control graphics device, to obtain 3 X 1 grid.
par(mfrow = c(3,1))

plot(y    = VIX$VIX.Adjusted,
     x    = index(VIX),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'VIX',
     main = 'i) VIX Index')

plot(y    = Fit_MEM_GAMMA$FilteredValues$Mu,
     x    = index(VIX),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = expression(mu),
     main = 'ii) Mean')

plot(y    = Fit_MEM_GAMMA$FilteredValues$Mu^2/Fit_MEM_GAMMA$Par[4],
     x    = index(VIX),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = expression(sigma^2),
     main = 'iii) Variance')

# EXERCISE 2 #
# Author: Andreas Kracht Frandsen
# Date: 2020-01-16
# For documentation please read the Computational Part of the main document.

# Load required packages.
library(Rsolnp)
library(mvtnorm)
library(quantmod)

# Obtains the relevant updated dataset given in the assignment.
GSPC_DJI <- read.csv2('data/data.csv', sep = ';', dec = ',')

# Just to gather dates related to the dataset. Used for plotting.
getSymbols(Symbols = "^GSPC",
           from    = '2007-01-03',
           to      = '2019-01-01',
           src     = 'yahoo')

# First 5 observations of GSP_DJI.
head(GSPC_DJI)

# Plot of GSPC
plot(y    = GSPC_DJI$GSPC,
     x    = index(head(GSPC$GSPC.Adjusted,-1)),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'Return',
     main = 'SP500 Return Data')

# Plot of DJI
plot(y    = GSPC_DJI$DJI,
     x    = index(head(GSPC$GSPC.Adjusted,-1)),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'Return',
     main = 'DOW Return Data')

## GARCH(1,1) ##

# The following function is the filter for the GARCH(1,1) model.
# Input:  Omega, Alpha and Beta, the parameters in the model (Double). 
#         Y (Vector), the returns of the relevant series.
# Output: Output (List), log likelihood (Double) and the variances (Vector).
GARCHFilter <- function(Y, Omega, Alpha, Beta) {
  
  iT      = length(Y)
  Sigma2 = numeric(iT)
  
  # The first variance is set to the empirical variance of the first 10 % of the observations.
  Sigma2[1] = var(Y[1:round(iT * 0.1)])
  
  # Compute the likelihood of the first observation.
  LLK = dnorm(Y[1], 0, sqrt(Sigma2[1]), log = TRUE)
  
  # For the rest (T-1) observations we use the updating equation.
  for (t in 2:iT) {
    Sigma2[t] = Omega + Alpha * Y[t-1]^2 + Beta * Sigma2[t - 1]
    
    LLK = LLK + dnorm(Y[t], 0, sqrt(Sigma2[t]), log = TRUE)
  }
  
  Output = list()
  
  Output[["LLK"]] = LLK
  Output[["Sigma2"]] = Sigma2
  
  return(Output)
  
}

# The following function evaluates the negative log likelihood for further use in the optimization proces.
# Input:  Par (Vector), the parameters in the model, which are omega, alpha and beta.
#         Y (Vector), the returns of the relevant series.
# Output: -LLK (Double), the negative log likelihood.
ObjectiveFunction <- function(Par, Y) {
  
  Omega = Par[1]
  Alpha = Par[2]
  Beta  = Par[3]
  LLK = GARCHFilter(Y, Omega, Alpha, Beta)$LLK
  
  return(-LLK)
}

# The following function serves as a basis to evaluate the inner part of the inequality constraints that need to be satisfied to impose weak stationarity.
# Input:  Par (Vector), the parameters in the inner part of the inequality constraints, which are alpha and beta.
# Output: Alpha+Beta (Double), the inner part of the inequality constraints.
ineqfun_GARCH_WS <- function(Par, ...) {
  Alpha = Par[2]
  Beta  = Par[3]
  
  return(Alpha + Beta)
}

# The following function estimates the GARCH(1,1) model by first finding maximum likelihood estimates of our parameters.
# Input:  Y (Vector), the returns of the relevant series.
# Output: Output (List), the optimized parameters (Vector), the BIC (Double), the variances (Vector), the log likelihood (Double)
#         and the standardized residuals.

EstimateGARCH <- function(Y, ineqfun_GARCH = ineqfun_GARCH_WS, ineqLB = 0.00, ineqUB = 0.9999) {
  
  # We set starting value for Alpha and Beta and set Omega to target the unconditional variance of the GARCH(1,1) model.
  
  Alpha = 0.125
  Beta  = 0.85
  Omega = var(Y) * (1.0 - Alpha - Beta)
  
  Par = c(Omega, Alpha, Beta)
  
  # Use the solnp from the Rsolnp package to optimize the negative log likelihood.
  # By default we specity ineqLB = 0.00 and ineqUB = 0.9999 in order to match 0 < alpha + beta < 0.9999.
  optimizer = solnp(Par,
                    fun      = ObjectiveFunction,
                    Y        = Y,
                    ineqfun  = ineqfun_GARCH,
                    ineqLB   = ineqLB,
                    ineqUB   = ineqUB,
                    LB       = c(0.00001, 0.0001, 0.0001),
                    UB       = c(10.0, 0.999, 0.999)
  ) 
  
  Par = optimizer$pars
  LLK = -tail(optimizer$values, 1)
  
  # Here we run the filter using the optimal parameter values, to obtain the final estimates of the variance.
  Sigma2 = GARCHFilter(Y, Par[1], Par[2], Par[3])$Sigma2
  
  # Computation of Bayesian Information Criterion.
  iT = length(Y)
  BIC = (-2 * LLK + log(iT) * length(Par))
  
  # Compute standardized residuals.
  st_res <- Y/sqrt(Sigma2)
  
  Output = list()
  
  Output[["Par"]]    = Par
  Output[["LLK"]]    = LLK
  Output[["BIC"]]    = BIC
  Output[["Sigma2"]] = Sigma2
  Output[["st_res"]] = st_res
  
  return(Output)
}

# Fit GARCH(1,1) for GSPC.
Fit_GSPC = EstimateGARCH(GSPC_DJI$GSPC)

# Double check alpha+beta.
sum(Fit_GSPC$Par[-1])

# Fit GARCH(1,1) for DJI.
Fit_DJI = EstimateGARCH(GSPC_DJI$DJI)

# Double check alpha+beta.
sum(Fit_DJI$Par[-1])

# Make one list including both fits.
fit <- list(Fit_GSPC, Fit_DJI)

## DCC AND CCC ##

# The following function is the filter for the DCC (CCC) model.
# Input:  A and B (Double), the parameters in the model, which are a and b.
#         Eta (Matrix), the standardized residuals from GARCH(1,1).
#         Q (Matrix), the unconditional correlation.
# Output: Output (List), the log likelihood (Double) and the correlation matrix R.
DCCFilter <- function(Eta, A, B, Q) {
  
  iN = ncol(Eta)
  iT = nrow(Eta)
  
  Cor = array(0, dim = c(iN, iN, iT))
  aQ  = array(0, dim = c(iN, iN, iT))
  
  ## Initialize to the unconditional correlation.
  Cor[ , , 1] = Q
  aQ[ , , 1]  = Q
  
  # Compute the contribution to the likelihood of the first observation.
  LLK = Eta[1, , drop = FALSE] %*% solve(Cor[,, 1]) %*% t(Eta[1, , drop = FALSE]) - 
    Eta[1, , drop = FALSE]%*% t(Eta[1, , drop = FALSE]) + log(det(Cor[,, 1]))
  
  # For the rest (T-1) observations.
  for (t in 2:iT) {
    # Update the Q matrix.
    aQ[,, t] = Q * (1 - A - B) + A * t(Eta[t - 1, , drop = FALSE]) %*% Eta[t - 1, , drop = FALSE] + 
      B * aQ[,,t - 1]
    
    ## Compute the correlation matrix R.
    Cor[,, t] = diag(sqrt(1/diag(aQ[,, t]))) %*% aQ[,, t] %*% diag(sqrt(1/diag(aQ[,, t]))) 
    
    LLK = LLK + Eta[t, , drop = FALSE] %*% solve(Cor[,, t]) %*% t(Eta[t, , drop = FALSE]) - 
      Eta[t, , drop = FALSE] %*% t(Eta[t, , drop = FALSE]) + log(det(Cor[,, t]))
  }
  
  Output = list()
  
  Output[["LLK"]] = -0.5 * LLK
  Output[["Cor"]] = Cor
  
  return(Output)
}

# The following function estimates the DCC (CCC) model by first finding maximum likelihood estimates of our parameters.
# Input:  Y (Matrix), the returns of the relevant series.
#         fit (List), the fit of the Garch(1,1) models combined.
#         CCC (Boolean), shall the CCC model be computed instead of the DCC.
# Output: Output (List), the optimized parameters (Vector), the BIC (Double), the total log likelihood (Double),
#                        the correlation matrix (Matrix), the standard deviations (Matrix), the parameters of the Garch(1,1) (Vector)
#                        and the standardized residuals (Vector).
Estimate_DCC <- function(Y, fit, CCC = FALSE) {
  
  Eta <- cbind(unlist(fit[[1]]["st_res"]), unlist(fit[[2]]["st_res"]))
  
  
  # Compute unconditional correlation.
  Q = cor(Eta)
  
  
  if(CCC == FALSE){
    
    # Initial parameters of a and b.
    Par = c(0.04, 0.9)
    
    # Use the solnp from the Rsolnp package to optimize the negative log likelihood.
    optimizer = solnp(Par, fun = function(Par, Eta, Q) {
      
      Filter = DCCFilter(Eta, Par[1], Par[2], Q)
      NLLK = -as.numeric(Filter$LLK)
      return(NLLK)
      
    }, ineqfun = function(Par, ...) {
      sum(Par)
    }, ineqLB = 1e-4, ineqUB = 0.999, 
    LB = c(1e-4, 1e-4), UB = c(0.999, 0.999), 
    Eta = Eta, Q = Q)
    
    Par = optimizer$pars
    
    # Likelihood contribution of correlation.
    LLK_C = -tail(optimizer$values, 1)
    
    # Here we run the filter using the optimal parameter values, to obtain the final estimates of the correlation matrix.
    Filter = DCCFilter(Eta, Par[1], Par[2], Q)
  }
  
  else{
    Filter = DCCFilter(Eta, 0, 0, Q)
    
    LLK_C = Filter[["LLK"]]
  }
  
  Sigma = sqrt(cbind(unlist(fit[[1]]["Sigma2"]), unlist(fit[[2]]["Sigma2"])))
  Coef  = cbind(unlist(fit[[1]]["Par"]), unlist(fit[[2]]["Par"]))
  
  # Likelihood contribution of volatility from GARCH(1,1)'s.
  LLK_V = sum(unlist(fit[[1]]["LLK"]), unlist(fit[[2]]["LLK"]))
  
  # Total likelihood.
  LLK = LLK_V + LLK_C
  
  Cor = Filter[["Cor"]]
  
  iT = nrow(Y)
  
  # Computation of Bayesian Information Criterion.
  BIC = log(iT) * 8 - 2 * LLK
  
  Output = list()
  
  Output[["LLK"]]  = LLK
  Output[["Coef"]] = Coef
  
  if(CCC == FALSE){
    Output[["Par"]] = Par
  }
  
  Output[["Sigma"]] = Sigma
  Output[["Cor"]]   = Cor
  Output[["Eta"]]   = Eta
  Output[["BIC"]]   = BIC
  
  return(Output)
  
}

# Fit DCC and CCC for our returns.
Fit_DCC = Estimate_DCC(GSPC_DJI, fit)
Fit_CCC = Estimate_DCC(GSPC_DJI, fit, CCC = TRUE)

## MINIMUM VARIANCE PORTFOLIO ##

# The following function computes the Minimum Variance Portfolio.
# Input:  fit (List), the fit of the DCC or CCC model.
# Output: weight (Array), the optimal portfolio weights for the Minimum Variance Portfolio.
MVP <- function(fit){
  iT = length(fit$Sigma[ , 1])
  iN = ncol(Fit_DCC$Sigma)
  
  D = array(0, dim = c(iN, iN, iT))
  
  SIGMA_INV = array(0, dim = c(iN, iN, iT))
  
  TOP    = array(0, dim = c(1, iN, iT))
  ell    = array(1, dim = c(iN, 1))
  BOTTOM = array(0, dim = c(1,1,iT))
  weight = array(0, dim = c(1, iN, iT))
  
  for (t in 1:iT) {
    D[ , , t] = diag(fit$Sigma[t, ])
    
    SIGMA_INV[ , , t] = solve(D[ , , t]) %*% solve(fit$Cor[ , , t]) %*% solve(D[ , , t])
    
    TOP[ , , t] = SIGMA_INV[ , , t] %*% ell
    BOTTOM[ , , t] = t(ell) %*% SIGMA_INV[ , , t] %*% ell
    
    weight[ , , t] = TOP[ , , t] / BOTTOM[ , , t]
  }
  
  return(weight)
}

# Compute the weights for both DCC and CCC models.
weight_DCC <- MVP(Fit_DCC)
weight_CCC <- MVP(Fit_CCC)

# Plot MVP for DCC
plot(y    = weight_DCC[1,1,],
     x    = index(head(GSPC$GSPC.Adjusted,-1)),
     type ='l',
     ylim =c(-5,5),
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'Weight',
     main = 'Portfolio Weights for S&P500 and DOW using DCC'
)

lines(y    = weight_DCC[1,2,],
      x    = index(head(GSPC$GSPC.Adjusted,-1)),
      type ='l',
      col  = 'red',
      lty = 'dashed')

# Plot MPV for CCC
plot(y    = weight_CCC[1,1,],
     x    = index(head(GSPC$GSPC.Adjusted,-1)),
     type ='l',
     ylim =c(-5,5),
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'Weight',
     main = 'Portfolio Weights for S&P500 and DOW using CCC'
)

lines(y    = weight_CCC[1,2,],
      x    = index(head(GSPC$GSPC.Adjusted,-1)),
      type = 'l',
      col  = 'red',
      lty = 'dashed')

## CoVaR ##

# The following function computes the difference between the Multivariate Gaussian CDF and the squared significance level.
# Input:  CoVar (Double), the CoVaR.
#         VaR (Double), the VaR.
#         sigma (Matrix), the standard deviation matrix.
#         alpha (Double), the significance level.
# Output: target (Double), the value to optimize over.
bi_pnorm_t <- function(CoVaR, VaR, sigma, alpha){
  func <- pmvnorm(upper = c(CoVaR, VaR), sigma = sigma)
  target <- func - alpha^2
}

# The following function computes the CoVaR.
# Input:  fit (List), the fit of either DCC or CCC.
#         alpha (Double), the significance level.
# Output: CoVaR (Vector), the CoVaR over time.
covar <- function(fit, alpha){
  iT <- length(fit$Sigma[,1])
  
  D <- array(0, dim = c(2,2,iT))
  CoVaR <- c()
  
  for (t in 1:iT) {
    D[,,t] = diag(fit$Sigma[t,])
    
    SIGMA = D[,,t] %*% fit$Cor[,,t] %*% D[,,t]
    
    sdY_2 <- sqrt(SIGMA[1, 2])
    
    VaR <- qnorm(alpha, 0, sdY_2)
    
    CoVaR[t] <- uniroot(bi_pnorm_t, interval = c(-10^4, 10), VaR = VaR, sigma = SIGMA, alpha=alpha)[[1]]
  }
  
  return(CoVaR)
}

# Compute the CoVaR at 0.01 and 0.05 significance level for both models.
DCC_CoVaR_1 <- covar(Fit_DCC, 0.01)
DCC_CoVaR_5 <- covar(Fit_DCC, 0.05)

CCC_CoVaR_1 <- covar(Fit_CCC, 0.01)
CCC_CoVaR_5 <- covar(Fit_CCC, 0.05)

# Plot of CoVaR for DCC
plot(y = DCC_CoVaR_1,
     x    = index(head(GSPC$GSPC.Adjusted,-1)),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'CoVaR',
     main = 'CoVaR of DCC Model')
lines(y = DCC_CoVaR_5,
      x    = index(head(GSPC$GSPC.Adjusted,-1)),
      type = 'l',
      xaxs = "i",
      yaxs = "i",
      col = 'red',
      lty = 'dashed')

# Plot of CoVaR for CCC
plot(y= CCC_CoVaR_1,
     x    = index(head(GSPC$GSPC.Adjusted,-1)),
     type = 'l',
     xaxs = "i",
     yaxs = "i",
     xlab = 'Date',
     ylab = 'CoVaR',
     main = 'CoVaR of CCC Model')
lines(y = CCC_CoVaR_5,
      x    = index(head(GSPC$GSPC.Adjusted,-1)),
      type = 'l',
      xaxs = "i",
      yaxs = "i",
      col = 'red',
      lty = 'dashed')
