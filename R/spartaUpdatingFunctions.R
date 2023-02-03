#library(usethis)
#usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).


#This function runs MCMC using particle filter MCMC and returns weights

spartaNimWeights <- function(model, #nimbleModel
                             MCMCconfiguration = NULL, #configuration for MCMC
                             pfType = NULL,#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                             pfControl = list(saveAll = TRUE,
                                              lookahead = "mean",
                                              smoothing = FALSE), #list of controls for particle filter
                             nParFiltRun = NULL, #Number of PF runs
                             latent #the latent variable
                             #newData,
                             #weights,
                             #weightedLatentSamples,
                             #unweightedLatentSamples
){

  timeStart1 <- Sys.time()
  target = MCMCconfiguration[["target"]]
  additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
  n.iter = MCMCconfiguration[["n.iter"]]
  n.chains = MCMCconfiguration[["n.chains"]]
  n.burnin = MCMCconfiguration[["n.burnin"]]
  if(is.null(nParFiltRun)) nParFiltRun = 10000

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")
  if(!is.null(pfType)){
    if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
#particleFilter is used to run the MCMC
#particleFilterEst is used for returning the weights at the posterior
    #values of the top-level nodes
    if(pfType == "auxiliary"){
      particleFilter <- myphdthesis::buildAuxiliaryFilterNew(model,
                                                             latent,
                                                             control = pfControl)

        particleFilterEst <- myphdthesis::buildAuxiliaryFilterNew(estimationModel,
                                                             latent,
                                                             control = pfControl)
    }

    if(pfType == "bootstrap"){
      particleFilter <- nimbleSMC::buildBootstrapFilter(model,
                                                        latent,
                                                        control = pfControl)

        particleFilterEst <-  nimbleSMC::buildBootstrapFilter(estimationModel,
                                                              latent,
                                                              control = pfControl)
    }
  }else{
    particleFilter <- myphdthesis::buildAuxiliaryFilterNew(model,
                                                           latent,
                                                           control = pfControl)

      particleFilterEst <- myphdthesis::buildAuxiliaryFilterNew(estimationModel,
                                                           latent,
                                                           control = pfControl)
  }

  message("Compiling the particle filter")
  #compiling the model
  compiledParticleFilter <- compileNimble(model,  particleFilter)

  #Loglikelihood of last run and the Effective sample sizes
  message("Running the particle filter")
  logLik <-   compiledParticleFilter$particleFilter$run(m = nParFiltRun)
  ESS <-   compiledParticleFilter$particleFilter$returnESS()


  message("Setting up the MCMC Configuration")
  #model <- model$newModel(replicate = TRUE)
  modelMCMCconf <- nimble::configureMCMC(model, nodes = NULL)

  if(is.null(pfType)) pfType = "auxiliary" #set auxiliary as default pfType

  modelMCMCconf$addSampler(target = target,
                           type = 'RW_PF_block',
                           control = list(latents = latent,
                                          pfControl = list(saveAll = TRUE),
                                          pfNparticles = nParFiltRun,
                                          pfType = pfType))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- nimble::buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = TRUE)

  message("Running the PF MCMC")
  #run MCMC
  timeStart2 <- Sys.time()
  mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                              niter = n.iter,
                              nchains = n.chains,
                              nburnin = n.burnin,
                              #inits = initsList,
                              #thin = 5,
                              setSeed = TRUE,
                              samples=TRUE,
                              samplesAsCodaMCMC = TRUE,
                              summary = TRUE,
                              WAIC = FALSE)
timeEnd <- Sys.time()

timetaken1 <- as.numeric(round(timeEnd - timeStart1, 2))
timetaken2 <- as.numeric(round(timeEnd - timeStart2, 2))

message(paste("Estimating weights at posteror values of ", target))
#posteriorEstimates <- mcmc.out$summary$all.chains[target, 'Mean']

for(i in 1:length(target)){
estimationModel[[target[i]]] <- mcmc.out$summary$all.chains[target[i], 'Mean']
}

message("Compiling the estimation particle filter")
#compiling the model
compiledParticleFilterEst <- compileNimble(estimationModel,  particleFilterEst)

#Loglikelihood of last run and the Effective sample sizes
message("Running the estimation particle filter")
logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun)
ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()


#save weights and samples
message("Extracting the weights and samples from particle fiter")
weights <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "wts")
unweightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, latent)
weightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvEWSamples, latent)
logLike <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "auxlog")


  message("Returning the results")
  #list to return
  returnList <- list(weights = weights,
                     logLike = logLike,
                     unweightedSamples = unweightedSamples,
                     weightedSamples = weightedSamples,
                     particleFilter = particleFilter,
                     mcmcSamplesAndSummary = mcmc.out,
                     timeTakenAll = timetaken1,
                     timeTakenRun = timetaken2)

  return(returnList)
}


###########################################

spartaNimUpdates <- function(model, #nimbleModel
                             MCMCconfiguration = NULL, #configuration for MCMC
                             pfType = NULL,#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                             pfControl = NULL, #list of controls for particle filter
                             nParFiltRun = NULL, #Number of PF runs
                             #updatePFControl = list(iNodePrev = NULL, M = NULL),
                             latent, #the latent variable
                             #newData,
                             weights,
                             weightedLatentSamples,
                             unweightedLatentSamples,
                             loglike
){
  timeStart1 <- Sys.time()
  target = MCMCconfiguration[["target"]]
  additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
  n.iter = MCMCconfiguration[["n.iter"]]
  n.chains = MCMCconfiguration[["n.chains"]]
  n.burnin = MCMCconfiguration[["n.burnin"]]
  iNodePrev = pfControl[["iNodePrev"]]
  M = pfControl[["M"]]
  #iNodePrev = updatePFControl[["iNodePrev"]]
  #M = updatePFControl[["M"]]
  if(is.null(nParFiltRun)) nParFiltRun = 10000

  #create new model for weights
  estimationModel <- model$newModel(replicate = TRUE)

  message("Building particle filter for model")
  if(!is.null(pfType)){
    # if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
    particleFilter <- myphdthesis::buildAuxiliaryFilterUpdate(model,
                                                              latent,
                                                             mvWSamplesWTSaved = weights,
                                                              mvWSamplesXSaved = unweightedLatentSamples,
                                                              mvEWSamplesXSaved = weightedLatentSamples,
                                                              logLikeVals = loglike,
                                                              control = pfControl)

    particleFilterEst <- myphdthesis::buildAuxiliaryFilterUpdate(estimationModel,
                                                                    latent,
                                                                    mvWSamplesWTSaved = weights,
                                                                    mvWSamplesXSaved = unweightedLatentSamples,
                                                                    mvEWSamplesXSaved = weightedLatentSamples,
                                                                    logLikeVals = loglike,
                                                                    control = pfControl)
  }else{
    particleFilter <- myphdthesis::buildAuxiliaryFilterUpdate(model,
                                                              latent,
                                                              mvWSamplesWTSaved = weights,
                                                              mvWSamplesXSaved = unweightedLatentSamples,
                                                              mvEWSamplesXSaved = weightedLatentSamples,
                                                              logLikeVals = loglike,
                                                              control = pfControl)

    particleFilterEst <- myphdthesis::buildAuxiliaryFilterUpdate(estimationModel,
                                                                    latent,
                                                                    mvWSamplesWTSaved = weights,
                                                                    mvWSamplesXSaved = unweightedLatentSamples,
                                                                    mvEWSamplesXSaved = weightedLatentSamples,
                                                                    logLikeVals = loglike,
                                                                    control = pfControl)
  }

  message("Compiling the particle filter")
  #compiling the model
  compiledParticleFilter <- nimble::compileNimble(model,  particleFilter)

  #Loflikelihood of last run and the Effective sample sizes
  message("Running the particle filter")
  logLik <-   compiledParticleFilter$particleFilter$run(m = nParFiltRun)
  ESS <-   compiledParticleFilter$particleFilter$returnESS()

  #save weights and samples
  message("Extracting the weights and samples from particle fiter")
  weightsNew <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, "wts")
  unweightedSamples <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, latent)
  weightedSamples <- as.matrix( compiledParticleFilter$particleFilter$mvEWSamples, latent)
  loglike <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, "auxlog")
  message("Setting up the MCMC Configuration")
  #newModel <- model$newModel(replicate = TRUE)
  modelMCMCconf <- nimble::configureMCMC(model, nodes = NULL)

  pfType = 'auxiliaryUpdate'


  mvWSamplesWTSaved = weights
  mvWSamplesXSaved = unweightedLatentSamples
  mvEWSamplesXSaved = weightedLatentSamples

  modelMCMCconf$addSampler(target = target,
                           type = 'RW_PF_blockUpdate',
                           control = list(latents = latent,
                                          pfControl = list(saveAll = TRUE, M = M, iNodePrev = iNodePrev),
                                          pfNparticles = nParFiltRun,
                                          pfType = pfType,
                                          weights = weights,
                                          unweightedSamples= unweightedLatentSamples,
                                          weightedSamples = weightedLatentSamples,
                                          logLikeVals = loglike))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- nimble::buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = TRUE)
  timeStart2 <- Sys.time()
  message("Running the PF MCMC")
  #run MCMC
  mcmc.out <- nimble::runMCMC(compiledList$modelMCMC,
                              niter = n.iter,
                              nchains = n.chains,
                              nburnin = n.burnin,
                              #inits = initsList,
                              #thin = 5,
                              setSeed = TRUE,
                              samples=TRUE,
                              samplesAsCodaMCMC = TRUE,
                              summary = TRUE,
                              WAIC = FALSE)
  timeEnd <- Sys.time()

  timetaken1 <- as.numeric(round(timeEnd - timeStart1, 2))
  timetaken2 <- as.numeric(round(timeEnd - timeStart2, 2))

  message(paste("Estimating weights at posteror values of ", target))
  #posteriorEstimates <- mcmc.out$summary$all.chains[target, 'Mean']

  for(i in 1:length(target)){
    estimationModel[[target[i]]] <- mcmc.out$summary$all.chains[target[i], 'Mean']
  }

  message("Compiling the estimation particle filter")
  #compiling the model
  compiledParticleFilterEst <- compileNimble(estimationModel,  particleFilterEst)

  #Loglikelihood of last run and the Effective sample sizes
  message("Running the estimation particle filter")
  logLik <-   compiledParticleFilterEst$particleFilterEst$run(m = nParFiltRun)
  ESS <-   compiledParticleFilterEst$particleFilterEst$returnESS()


  #save weights and samples
  message("Extracting the weights and samples from particle fiter")
  weights <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "wts")
  unweightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, latent)
  weightedSamples <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvEWSamples, latent)
  logLike <- as.matrix( compiledParticleFilterEst$particleFilterEst$mvWSamples, "auxlog")



  message("Returning the results")
  #list to return
  returnList <- list(weights = weightsNew,
                     unweightedSamples = unweightedSamples,
                     weightedSamples = weightedSamples,
                     logLike = loglike,
                     particleFilter = particleFilter,
                     mcmcSamplesAndSummary = mcmc.out,
                     timeTakenAll = timetaken1,
                     timeTakenRun = timetaken2)

  return(returnList)
}


