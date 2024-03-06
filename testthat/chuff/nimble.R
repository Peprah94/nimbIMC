setwd("~/OneDrive - NTNU/GitHub/Thesis/INLA within NIMBLE/Bayesian Lasso")

load("hitters_data.RData")

#Packages
library(nimble)
library(INLA)
library(mvtnorm)
library(MASS)
library(parallel)
library(coda)
library(ggmcmc)


code <- nimbleCode({
 # est_lam ~ dgamma(shape=1, rate=0.0005)
  #Prior for beta1 and beta2
  alpha ~ dnorm(0,1)
  for(j in 1:P){
    beta[j] ~ ddexp(location = 0,rate=est_lam)
  }
  
  prec ~ dgamma(shape=1, rate=0.00005)
  #Fitting the inla with the simulated parameters
  
  # linear model specification
  for(i in 1:N){
    linpred[i] <- alpha + beta[1]*x[i,1] + beta[2]*x[i,2] + beta[3]*x[i,3] + beta[4]*x[i,4] + beta[5]*x[i,5]
    y[i] ~ dnorm(linpred[i],tau=prec ) 
  }

})

## Parameterising the nimble model

#Data
inla_data <- list(y=as.numeric(df$y), 
                  x = df$x)

#Constants
const <- list(N = length(df$y),
              P= ncol(df$x),
              est_lam = 1/0.073
)

# Initial values
idm_inits <- function(){list(beta =rep(0, ncol(df$x),
                                       alpha=0)
)
}

initsList <- idm_inits()

#Putting all together for the creating compilation
modelInfo <- list(
  code = code,
  constants = const,
  data = inla_data,
  inits = initsList
)

#Create the model in nimble
mwtc <- nimbleModel(code,data = inla_data, 
                    constants = const, 
                    inits = initsList)

# Create the model in C
Cmwtc <- compileNimble(mwtc,showCompilerOutput = FALSE) 



mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta", "prec", "alpha"))

Rmcmc <- buildMCMC(mcmcconf, useConjugacy=FALSE)

# Compile 
cmcmc <- compileNimble(Rmcmc, 
                       project = Cmwtc, 
                       resetFunctions = TRUE)

# Run the MCMC
mcmc.out <- runMCMC(cmcmc, 
                    niter = 50000,
                    nchains = 3,
                    nburnin = 20000,
                    #inits = initsList,
                    #thin =10, 
                    setSeed = TRUE, 
                    samples=TRUE, 
                    samplesAsCodaMCMC = TRUE, 
                    summary = TRUE, 
                    WAIC = FALSE)

#Output from the MCMC
output <- mcmc.out$summary
output


mcmc_samples <- ggs(mcmc.out$samples)
ggs_traceplot(mcmc_samples)
ggs_density(mcmc_samples)
ggs_effective(mcmc_samples, proportion = FALSE)
