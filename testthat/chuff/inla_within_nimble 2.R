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

nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE)
nimbleOptions(MCMCsaveHistory = TRUE)

fit.inla <- function(x, y, beta){
  ii  <- get("ii",envir =  parent.frame())
  ii <- assign("ii",ii+1,envir = parent.frame())
  print(ii)
  data <- list(y=y, x=x)
  data$oset = data$x %*% beta
  res = INLA::inla(y ~ -1 + offset(oset), data = data, family = "gaussian",
                   control.predictor = list(compute=TRUE))
  #res = INLA::inla.rerun(res)
  fitted_values = res$summary.fitted.values[,"mean"]
  #intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  ret <- cbind(precision, fitted_values)
  return(ret)
}


nimbleINLA <- nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1) # beta is a vector
  ) {},

  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)
#nimbleINLA(x[1:263,1:5], y[1:263], c(0,0.18, 0,0,0.23))

ii <- 0

code <- nimbleCode({

  alpha ~ dnorm(0,1)


  for(j in 1:P){
    beta[j] ~ ddexp(location = 0,rate=est_lam)
  }
  #Fitting the inla with the simulated parameters
  inla.res[1:N, 1:2] <- nimbleINLA(x[1:N,1:P],y_obs[1:N],beta[1:P])

  #linpred[1:N] <- inla.res[1:100,3]
  sigma <- inla.res[1,1]
  linpred[1:N] <-  inla.res[1:N, 2]+alpha

  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i],tau=sigma )
  }

 # tau <- sigma
  #intercept <- inla.res[1,1]
})

## Parameterising the nimble model

#Data
inla_data <- list(y=as.numeric(df$y),
                  x = df$x,
                  y_obs=as.numeric(df$y))

#Constants
const <- list(N = length(df$y),
              P= ncol(df$x),
              est_lam = 1/0.073
              )

# Initial values
idm_inits <- function(){list(alpha = 0,
                            beta=rep(0,const$P)
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
Cmwtc <- compileNimble(mwtc,showCompilerOutput = FALSE) #Have issues compiling


mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta", "alpha","sigma"))
mcmcconf$removeSamplers(c("beta"))
mcmcconf$addSampler(target=c("beta"),
                    type="RW_block",
                    control=list(adaptInterval=10))
#mcmcconf$removeSamplers("beta")
#mcmcconf$addSampler(c("beta"), "RW_block")
#mcmcconf$printSamplers("beta")

Rmcmc <- buildMCMC(mcmcconf, useConjugacy=FALSE)

# Compile
cmcmc <- compileNimble(Rmcmc,
                       project = Cmwtc,
                       resetFunctions = TRUE)

# Run the MCMC
mcmc.out <- runMCMC(cmcmc,
                    niter = 1000,
                    nchains = 1,
                    nburnin = 200,
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

