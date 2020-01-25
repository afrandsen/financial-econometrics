# EXERCISE 1 #
# Author: Flow ID number 38
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
  dLLK = -tail(optimiser$value, n=1)
  
  # Here we run the filter using the optimal parameter values, to obtain the final estimates of mu.
  FilteredValues = MEMGAMMA_Filter(Par, Y)
  
  iT = length(Y)
  
  # Computation of Bayesian Information Criterion, using the fact that we estimate four parameters.
  BIC = (log(iT) * 4 - 2 * dLLK)
  
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
