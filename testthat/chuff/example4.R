


#Zero Inflated Poisson

library(INLA)
library(pscl)
library(mvtnorm)
library(readr)
library(dplyr)

# Load data
# https://stats.idre.ucla.edu/r/dae/zip/
#zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv")
zinb <- read_csv("testthat/fish.csv")
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})

summary(zinb)

#Maximum likehood estimates
# Full data
d <- zinb

# ZIP
ml.res <- summary(zeroinfl(count ~ child + camper | persons, data = d))
ml.res

#estimates of gamma to use as proposal
muGammaEst <- as.numeric(ml.res$coefficients[[2]][1:2,1])
covGammaEst <- diag(3*as.numeric(ml.res$coefficients[[2]][1:2,2]))

#values for use
x <- cbind(zinb$child,zinb$camper, zinb$persons)

  # zinb%>%
  # dplyr::select(child, camper, persons)

y <- zinb$count

beta = c(1,-0.5)

family = "zeroinflatedpoisson1"
fixedVals <- c("intercept", "beta1", "beta2")
fit.inla <- function(x ,
                     y ,
                     beta,
                     fixedVals,
                     #interInModel,
                     family) {
  data <- data.frame(count = y,
                     child = x[,1],
                     camper = x[,2],
                     persons = x[, 3])

  logit_pi <- beta[1] + beta[2] * data$persons

  # Define pi hyper for likelihood
  hyper_pi <- lapply(logit_pi, function(X) {
    list(hyper = list(prob = list(fixed = TRUE, initial = X)))
  })



  # Define COUNT as diagonal matrix of observed counts
  COUNT <- matrix(NA, nrow = nrow(data), ncol = nrow(data))
  diag(COUNT) <- data$count

  res <- inla(COUNT ~ child + camper,
              data = list(COUNT = COUNT, child = data$child, camper = data$camper),
              family = rep(family, nrow(data)),
              num.threads = "1:1",
              control.fixed = list(prec.intercept = 0.001),
              control.family = hyper_pi
  )
  #precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  beta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
  beta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
  fitted_values = res$mlik[[1]]
  #fitted_values = res$summary.fitted.values[,"mean"]
  ret <- cbind(fitted_values, intercept, beta1, beta2)
  colnames(ret) <- c("mld", fixedVals)
  return(ret)
}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("intercept", "beta1", "beta2")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "zeroinflatedpoisson1")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

# nimbleINLA(x = x,
#            y = y,
#            beta = beta)

### this whole section is available at
## https://r-nimble.org/nimbleExamples/zero_inflated_poisson.html
## with documentation
##


dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), zeroProb = double(),
                 log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if(x != 0) {
      ## return the log probability if log = TRUE
      if(log) return(dpois(x, lambda, log = TRUE) + log(1-zeroProb))
      ## or the probability if log = FALSE
      else return((1-zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1-zeroProb) * dpois(0, lambda, log = FALSE)
    if(log) return(log(totalProbZero))
    return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if(isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
  )))

code <- nimbleCode({

  for( i in 1:3){
    beta[i] ~ dnorm(0, 0.001)
  }

  for( i in 1:2){
    gamma[i] ~ dnorm(0, 0.001)
  }
  #alpha ~ dnorm(0, 0.001)
  for(i in 1:N){
 logit(p[i]) <- gamma[1] + gamma[2]* x[i,3]
  }


  for(i in 1:N){
    log(mu[i]) <- beta[1] + beta[2]*x[i,1] + beta[3]*x[i,2]
  }
  #Fitting the inla with the simulated parameters


  # linear model specification
  for(i in 1:N){
    y[i] ~ dZIP(mu[i], zeroProb = p[i])
  }

})

inla_data <- list(y = as.numeric(y),
                  x = as.matrix(as.data.frame(x))
)

#Constants
const <- list(N = nrow(inla_data$x))

# Initial values
idm_inits <- function(){list(
                             beta = rep(0, 3),
                             gamma = muGammaEst
)
}



ret <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = x,
                               code = code,
                               family = "zeroinflatedpoisson1",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = c("inlamcmc"),
                               inlaMCsampler = "AFSS_INLA_block",
                               samplerControl = list(propCov = covGammaEst,
                                                     mu = muGammaEst),
                               parametersToMonitor = list(mcmc = c("gamma"),
                                                          inla = c("beta")),
                               mcmcConfiguration = list(n.chains = 1,
                                                        n.iterations = 200,
                                                        n.burnin = 50,
                                                        n.thin = 1,
                                                        setSeed = TRUE,
                                                        samples=TRUE,
                                                        samplesAsCodaMCMC = TRUE,
                                                        summary = TRUE,
                                                        WAIC = FALSE)
)








 initsList <- idm_inits()
mwtc <- nimble::nimbleModel(code,
                            data = inla_data,
                            constants =const,
                            inits = initsList)

Cmwtc <- nimble::compileNimble(mwtc,
                               showCompilerOutput = FALSE)

#fixedVals = c("beta[1]", "beta[2]","beta[3]","beta[4]", "sigma")
fixedVals = c("beta")
ret <- buildINLAmodel(mwtc,
                      fam = family,
                      x = x, y= as.numeric(y),
                      control = list(fit.inla = nimbleINLA,
                                     fixedVals = fixedVals))

cwtm <- compileNimble(ret,
                      showCompilerOutput = FALSE)
cwtm$run(beta = beta)
cwtm$mvEWSamples


mwtc <- nimble::nimbleModel(code,
                            data = inla_data,
                            constants = const,
                            inits = initsList)
#cwtm$run(data$x, data$y, c(0,0))
# Create the model in C
Cmwtc <- nimble::compileNimble(mwtc,
                               showCompilerOutput = FALSE) #Have issues compiling

mcmcconf <- nimble::configureMCMC(Cmwtc,
                                  nodes = NULL)

#mcmcconf$removeSampler(c("beta", "a", "sigma"))
fixedVals = "beta"
mcmcconf$addSampler(target = "gamma",
                    type = "AFSS_INLA_block",
                    control = list(fit.inla = nimbleINLA,
                                   x = x,
                                   y = as.numeric(y),
                                   fixedVals = fixedVals,
                                   fam = family,
                                   scale = 1,
                                   sliceMaxSteps = 15,
                                   maxContractions = 20))

mcmcconf$printSamplers()
Rmcmc <- nimble::buildMCMC(mcmcconf)

# Compile
cmcmc <- nimble::compileNimble(Rmcmc,
                               project = Cmwtc,
                               resetFunctions = TRUE)

#cmcmc$run(1000)

cmc.out <- nimble::runMCMC(cmcmc,
                           niter = 5,
                           nchains = 1,
                           nburnin = 1,
                           #inits = initsList,
                           thin = 1,
                           setSeed = TRUE,
                           samples = TRUE,
                           samplesAsCodaMCMC = TRUE,
                           summary = TRUE,
                           WAIC = FALSE)


cmc.out$summary

ggmcmc::ggs(cmc.out$samples)%>%
  ggs_density()
