###MIssing covariates

library(mice)
data(nhanes2)
library(nimble)

# data
d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi)) # finding na's
n.mis <- length(idx.mis) # number of nans
d.mis = data.frame(age = as.numeric(d.mis$age),
              bmi = d.mis$bmi,
              chl = d.mis$chl)%>%
  dplyr::arrange(bmi)

df = list(d.mis = d.mis, idx.mis = idx.mis)

ageModelMat <- model.matrix(~ as.factor(age),
             data = d.mis,
             contrasts.arg = list(age = contrasts(as.factor(d.mis$age), contrasts = TRUE)))

#model matrix covariates for NIMBLE
modelMat <- cbind(ageModelMat[, 2:3], d.mis$bmi)


fitINLAMissingValues <- function(x ,
                                 y ,
                                 beta,
                                 fixedVals,
                                 #interInModel,
                                 family){
  #age = x[,1]
  #bmi = x[,2]
  #chl = y
  missData <- data.frame(age= x[,1],
                         age1 = x[,2],
                         bmi = x[,3],
                         chl = y)
  #subset missing values
  idx.mis <- c(which(is.na(missData[,3])))

  missData[idx.mis,3] = beta



  res = inla(chl ~ 1 + age + age1 + bmi,
             data = missData,
             family = family,
             verbose=FALSE,
             control.predictor = list(compute = TRUE))

  beta0 = res$marginals.fixed[[1]]
  beta1 = res$marginals.fixed[[2]]
  beta2 = res$marginals.fixed[[3]]
  tau = res$marginals.hyperpar[[1]]
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  estbeta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
  estbeta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
  estbeta3 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[4]])
  fitted_values = c(res$mlik[1,1])
  #fitted_values = res$summary.fitted.values[,"mean"]
  ret <- cbind(fitted_values, intercept,estbeta1, estbeta2,
               estbeta3 ,precision)
  colnames(ret) <- c("mld", fixedVals)
  return(ret)
}


nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("beta[1]", "beta[2]","beta[3]","beta[4]", "sigma")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "gaussian")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fitINLAMissingValues'
)

# nimbleINLA(x = modelMat,
#            y = d.mis$chl,
#            beta = rep(10, 9))


code <- nimbleCode({

  for( i in 1:4){
    beta[i] ~ dnorm(0, 0.001)
  }
  #alpha ~ dnorm(0, 0.001)
for(i in 1:n.idx){
  eta[i] ~ dnorm(muMiss, sd = covMiss)
}
  #eta[1: n.idx] ~ dmnorm(muMiss[1:n.idx], cov = covMiss[1:n.idx, 1:n.idx])

  # for(i in 1:17){
  #   xobs[i] <- x[i,2]
  # }
  #
  # for(i in 17:25){
  #   xobs[i] <- x[i,2]
  # }
  #Fitting the inla with the simulated parameters
 # inla.res[1:N, 1:5] <- nimbleINLAMissingValues(x[1:N,1:3], idxMiss[1:n.idx], eta[1: n.idx])

  sigma ~ dgamma(1,0.00005)
  for(i in 1:16){
  linpred[i] <- beta[1] + beta[2]* x[i,1] + beta[3]* x[i,2]  + beta[4]*x[i,3]
  }

  for(i in 17:25){
    linpred[i] <- beta[1] + beta[2]* x[i,1] + beta[3]* eta[i - 16]  + beta[4]*x[i,3]
  }
    #inprod(beta[1:P], x[i, 1:P])


  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

})

inla_data <- list(y = as.numeric(d.mis$chl),
                  x = modelMat
)

#covMiss <- 4 * ()
#Constants
const <- list(N = nrow(inla_data$x),
              n.idx = length(df$idx.mis),
              #muMiss = rep(mean(inla_data$x[,3], na.rm = T), length(df$idx.mis)),
              #covMiss = diag(mean(inla_data$x[,3], na.rm = T), length(df$idx.mis))
              muMiss = mean(inla_data$x[,3], na.rm = T),
              covMiss = 2*sd(inla_data$x[,3], na.rm = T)


)

# Initial values
idm_inits <- function(){list(eta = rep(mean(inla_data$x[,3], na.rm = T), const$n.idx),
                             beta = rep(0, 4)
)
}

x = modelMat; y= as.numeric(d.mis$chl)
nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)
ret <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = x,
                               code = code,
                               family = "gaussian",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = c("inlamcmc"),
                               inlaMCsampler = "RW_INLA_block",
                               samplerControl = list(propCov = diag(4*var(inla_data$x[,3], na.rm = T), 9),
                                                     #interInModel = 0,
                                                     mu = c(rep(mean(inla_data$x[,3], na.rm = T),9)),
                                                     #scale = 4*var(inla_data$x[,3], na.rm = T),
                                                     adaptive = FALSE),
                               parametersToMonitor = list(inla = c("beta", "sigma"),
                                                          mcmc = c("eta")),
                               mcmcConfiguration = list(n.chains = 1,
                                                        n.iterations = 200,
                                                        n.burnin = 20,
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
fixedVals = c("beta", "sigma")
ret <- buildINLAmodel(mwtc,
                      fam = "gaussian",
                      x = modelMat, y = "y", #y= as.numeric(d.mis$chl),
                      control = list(fit.inla = nimbleINLA,
                                     fixedVals = fixedVals))

cwtm <- compileNimble(ret,
                      showCompilerOutput = FALSE)
cwtm$run(beta = rep(10, 9))
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
obs <- c("y")
mcmcconf$addSampler(target = "eta",
                    type = "RW_INLA_block",
                    control = list(fit.inla = nimbleINLA,
                                   x = modelMat,
                                   y = obs,
                                   fixedVals = fixedVals,
                                   fam = "gaussian",
                                   scale = 1))

mcmcconf$printSamplers()
Rmcmc <- nimble::buildMCMC(mcmcconf)

# Compile
cmcmc <- nimble::compileNimble(Rmcmc,
                               project = Cmwtc,
                               resetFunctions = TRUE)

#cmcmc$run(1000)

cmc.out <- nimble::runMCMC(cmcmc,
                           niter = 2000,
                           nchains = 2,
                           nburnin = 500,
                           #inits = initsList,
                           thin = 1,
                           setSeed = TRUE,
                           samples = TRUE,
                           samplesAsCodaMCMC = TRUE,
                           summary = TRUE,
                           WAIC = FALSE)

cmc.out$summary$all.chains

ggmcmc::ggs(cmc.out$samples)%>%
  ggs_density()
