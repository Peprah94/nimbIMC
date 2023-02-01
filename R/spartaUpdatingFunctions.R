#library(usethis)
#usethis::edit_r_environ()
#when the tab opens up in R studio, add this to the 1st line: R_MAX_VSIZE=100Gb (or whatever memory you wish to allocate).


#This function runs MCMC using particle filter MCMC and returns weights

spartaNimWeights <- function(model, #nimbleModel
                            MCMCconfiguration = NULL, #configuration for MCMC
                            pfType = NULL,#Either 'auxiliary' or 'bootstrap'. Defaults to auxiliary
                            pfControl = list(saveAll = TRUE, lookahead = "mean", smoothing = FALSE), #list of controls for particle filter
                            nParFiltRun = NULL, #Number of PF runs
                            latent #the latent variable
                            #newData,
                            #weights,
                            #weightedLatentSamples,
                            #unweightedLatentSamples
                            ){

target = MCMCconfiguration[["target"]]
additionalPars = MCMCconfiguration[["additionalPars"]] #other dependent variables you seek to monitor
n.iter = MCMCconfiguration[["n.iter"]]
n.chains = MCMCconfiguration[["n.chains"]]
n.burnin = MCMCconfiguration[["n.burnin"]]
if(is.null(nParFiltRun)) nParFiltRun = 10000


  message("Building particle filter for model")
  if(!is.null(pfType)){
    if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")

     if(pfType == "auxiliary"){
  particleFilter <- nimbleSMC::buildAuxiliaryFilter(model,
                                                       latent,
                                                       control = pfControl)
  }

  if(pfType == "bootstrap"){
    particleFilter <- nimbleSMC::buildBootstrapFilter(model,
                                                      latent,
                                                      control = pfControl)
  }
  }else{
    particleFilter <- nimbleSMC::buildAuxiliaryFilter(model,
                                                      latent,
                                                      control = pfControl)
  }

  message("Compiling the particle filter")
  #compiling the model
  compiledParticleFilter <- compileNimble(model,  particleFilter)

  #Loflikelihood of last run and the Effective sample sizes
  message("Running the particle filter")
  logLik <-   compiledParticleFilter$particleFilter$run(m = nParFiltRun)
  ESS <-   compiledParticleFilter$particleFilter$returnESS()

  #save weights and samples
  message("Extracting the weights and samples from particle fiter")
  weights <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, "wts")
  unweightedSamples <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, latent)
  weightedSamples <- as.matrix( compiledParticleFilter$particleFilter$mvEWSamples, latent)

  message("Setting up the MCMC Configuration")
  model <- model$newModel(replicate = TRUE)
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

  message("Returning the results")
  #list to return
  returnList <- list(weights = weights,
                     unweightedSamples = unweightedSamples,
                     weightedSamples = weightedSamples,
                     particleFilter = particleFilter,
                     mcmcSamplesAndSummary = mcmc.out)

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
                             unweightedLatentSamples
){

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


  message("Building particle filter for model")
  if(!is.null(pfType)){
   # if(!pfType %in% c("auxiliary", "bootstrap")) stop("Function currently works for auxiliary and bootstap Particle filters")
particleFilter <- myphdthesis::buildAuxiliaryFilterUpdate(model,
                                                          latent,
                                                          mvWSamplesWTSaved = weights,
                                                          mvWSamplesXSaved = unweightedLatentSamples,
                                                          mvEWSamplesXSaved = weightedLatentSamples,
                                                          control = pfControl)
  }else{
    particleFilter <- myphdthesis::buildAuxiliaryFilterUpdate(model,
                                                              latent,
                                                              mvWSamplesWTSaved = weights,
                                                              mvWSamplesXSaved = unweightedLatentSamples,
                                                              mvEWSamplesXSaved = weightedSamples,
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
  weights <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, "wts")
  unweightedSamples <- as.matrix( compiledParticleFilter$particleFilter$mvWSamples, latent)
  weightedSamples <- as.matrix( compiledParticleFilter$particleFilter$mvEWSamples, latent)

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
                                         weightedSamples = weightedLatentSamples))

  modelMCMCconf$addMonitors(additionalPars)

  message("Building and compiling the PF MCMC")
  ## build and compile pMCMC sampler
  modelMCMC <- nimble::buildMCMC(modelMCMCconf)
  compiledList <- nimble::compileNimble(model,
                                        modelMCMC,
                                        resetFunctions = TRUE)

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

  message("Returning the results")
  #list to return
  returnList <- list(weights = weights,
                     unweightedSamples = unweightedSamples,
                     weightedSamples = weightedSamples,
                     particleFilter = particleFilter,
                     mcmcSamplesAndSummary = mcmc.out)

  return(returnList)
}


