library('nimble')
library(nimbleSMC)
library(myphdthesis)
set.seed(1)


mvS <- modelValues(stateSpaceModel)

#devtools::install_github("Peprah94/myphdthesis")


load("~/Documents/GitHub/myphdthesis/example1BaselineSMC2.RData")
load("~/Documents/GitHub/myphdthesis/example1SimData2.RData")
load("~/Documents/GitHub/myphdthesis/example1ReducedIn45A2.RData")

nIterations = 5
nBurnin = 2
nChains = 2
nThin = 1

# ## define the model
message("Running baseline model")
stateSpaceCode <- nimbleCode({
  x[1] ~ dnorm(mu0, 1)
  y[1] ~ dnorm(x[1], 1)
  for(i in 2:t){
    x[i] ~ dnorm(x[i-1] * a + b, 1)
    y[i] ~ dnorm(x[i] * c, 1)
  }
  a ~ dunif(0, 1)
  b ~ dnorm(0, 1)
  c ~ dnorm(1,1)
  mu0 ~ dnorm(0, 1)
})
#
# ## define data, constants, and initial values
data <- list(
  #   #y = c(0.213, 1.025, 0.314, 0.521, 0.895, 1.74, 0.078, 0.474, 0.656, 0.802)
  y = simData$y
)
constants <- list(
  t = 50
)
inits <- list(
  a = reducedCase2$mcmcSamplesAndSummary$summary$all.chains[1,1],
  b = reducedCase2$mcmcSamplesAndSummary$summary$all.chains[1,2],
  mu0= reducedCase2$mcmcSamplesAndSummary$summary$all.chains[1,4],
  c = reducedCase2$mcmcSamplesAndSummary$summary$all.chains[1,3]
)
#
#
# ## build the model
stateSpaceModel <- nimbleModel(stateSpaceCode,
                               data = data,
                               constants = constants,
                               inits = inits,
                               check = FALSE)

model = stateSpaceModel
latent = "x"
MCMCconfiguration = list(target = c('a', 'b', 'c', 'mu0'),
                         additionalPars = "x",
                         n.iter = nIterations,
                         n.chains = nChains,
                         n.burnin = nBurnin,
                         n.thin = nThin)
weights = example1ReducedModel$weights
unweightedLatentSamples = example1ReducedModel$mvWS
weightedLatentSamples = example1ReducedModel$mvEWS
loglike = example1ReducedModel$logLike
pfControl = list(saveAll = TRUE,
                 #lookahead = "mean",
                 smoothing = FALSE,
                 mcmc = TRUE,
                 M = 5,
                 iNodePrev = 45)
pfType = NULL#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
 #list of controls for particle filter
nParFiltRun = NULL
postReducedMCMC = example1ReducedModel$mcmcSamplesAndSummary
mvSamplesEst =  example1ReducedModel$mvSamplesEst
control = pfControl

mcmccontrol = list(latents = latent,
               target = target,
               pfControl = list(saveAll = TRUE, M = M, iNodePrev = iNodePrev),
               pfNparticles = nParFiltRun,
               pfType = "bootstrapUpdate",
               weights = weights,
               unweightedSamples= unweightedLatentSamples,
               weightedSamples = weightedLatentSamples,
               mvSamplesEst = mvSamplesEst,
               logLikeVals = loglike)












mvConf = modelValuesConf(vars = c('a', 'b', 'c'),
                         type = c('double', 'int', 'double'),
                         size = list(a = 2, b =c(2,2), c = 1) )
customMV = modelValues(mvConf, m = 2)


m = getsize(reducedCase2$compiledParticleFilterEst$particleFilterEst$mvWSamples)
names <- reducedCase2$compiledParticleFilterEst$particleFilterEst$mvWSamples$varNames
size <- reducedCase2$compiledParticleFilterEst$particleFilterEst$mvWSamples$sizes
type <- c("double", "double", "double")

namesEWS <- names[!names %in% c("bootLL", "wts")]
index <- which(names == namesEWS)
mvEWS <- modelValues(modelValuesConf(vars = namesEWS,
                                           types = type[index],
                                           sizes = size[index]))
if(smoothing){
  size$wts <- 1
  size$bootLL <- 1
}
mvWS  <- modelValues(modelValuesConf(vars = names,
                                           types = type,
                                           sizes = size))

resize(mvWS, m)
resize(mvEW, m)

for(i in 1: m){
mvEW[[namesEWS]][[i]] <- weightedLatentSamples[i,]
mvWS[[namesEWS]][[i]] <- unweightedLatentSamples[i,]
}





