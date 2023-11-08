#' Fitting INLA within NIMBLE for Data Generating process
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
INLAWiNimDataGenerating <- function(data,
                                    covariate,
                                    code,
                                    family,
                                    modelData,
                                    modelConstants,
                                    modelInits,
                                    nimbleINLA,
                                    inlaMCMC = c("inla", "mcmc", "inlamcmc"),
                                    inlaMCsampler = "RW_INLA_block",
                                    samplerControl = list(proposal = "normal",
                                                          initMean = NULL,
                                                          initCov = NULL,
                                                          nCores = 2),
                                    parametersToMonitor = list(mcmc = c("mcmc"),
                                                               inla = c("inla"),
                                                               additionalPars = NULL),
                                    mcmcConfiguration = list(n.chains = 1,
                                                             n.iterations = 10,
                                                             n.burnin = 0,
                                                             n.thin = 1,
                                                             setSeed = TRUE,
                                                             samples=TRUE,
                                                             samplesAsCodaMCMC = TRUE,
                                                             summary = TRUE,
                                                             WAIC = FALSE)){



  #extract necessary variables
  x <- covariate # must be a matrix
  y <- data # must be a vector
  family <- family
  fixedVals <- parametersToMonitor$inla
  target <- parametersToMonitor$mcmc

  #set up initial value
  initsList <- modelInits()

  #Create the model in nimble
  mwtc <- nimble::nimbleModel(code,
                              data = modelData,
                              constants = modelConstants,
                              inits = initsList)

  # Create the model in C
  Cmwtc <- nimble::compileNimble(mwtc,
                                 showCompilerOutput = FALSE) #Have issues compiling

  #create list to return
  retList <- list()
  # Fit INLAMCMC

  if(inlaMCMC %in% "inlamcmc"){
    mcmcconf <- nimble::configureMCMC(Cmwtc,
                                      nodes = NULL)

    # setting sampler controls
    samplerControl$fit.inla = nimbleINLA
    samplerControl$x = x
    samplerControl$y = y
    samplerControl$fixedVals = fixedVals
    samplerControl$fam = family

    # mcmc configuration
    mcmcconf$addSampler(target = target,
                        type = inlaMCsampler,
                        control = samplerControl)

    mcmcconf$printSamplers(executionOrder = TRUE)
    mcmcconf$addMonitors(target)

    if(!is.null(parametersToMonitor$additionalPars)){
      mcmcconf$addMonitors(parametersToMonitor$additionalPars)
    }

    #build model
    Rmcmc <- nimble::buildMCMC(mcmcconf)

    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = Cmwtc,
                                   resetFunctions = TRUE)
    startTime <- Sys.time()
    mcmc.out <- nimble::runMCMC(cmcmc,
                                niter = mcmcConfiguration[["n.iterations"]],
                                nchains = mcmcConfiguration[["n.chains"]],
                                nburnin = mcmcConfiguration[["n.burnin"]],
                                #inits = initsList,
                                thin = mcmcConfiguration[["n.thin"]],
                                setSeed = mcmcConfiguration[["setSeed"]],
                                samples = mcmcConfiguration[["samples"]],
                                samplesAsCodaMCMC = mcmcConfiguration[["samplesAsCodaMCMC"]],
                                summary = mcmcConfiguration[["summary"]],
                                WAIC = mcmcConfiguration[["WAIC"]])
    endTime <- Sys.time()
    timeTaken <- difftime(endTime, startTime, units = "secs")
    #as.numeric(endTime - startTime)

    ret <- list(mcmc.out = mcmc.out,
                timeTaken = timeTaken)

    #save inlamcmc results
    retList$inlamcmc <- ret
  }

  if(inlaMCMC %in% "mcmc"){
    mcmcconf <- nimble::configureMCMC(Cmwtc,
                                      monitors = c(target, fixedVals))


    mcmcconf$printSamplers()

    # Add new samplers
    mcmcconf$removeSampler(target)
    # mcmcconf$removeSampler("beta")
    #mcmcconf$addSampler("beta", type = inlaMCsampler)
    mcmcconf$addSampler(target, type = inlaMCsampler,
                        control = samplerControl)

    mcmcconf$printSamplers(executionOrder = TRUE)

    if(!is.null(parametersToMonitor$additionalPars)){
      mcmcconf$addMonitors(parametersToMonitor$additionalPars)
    }
    #build model
    Rmcmc <- nimble::buildMCMC(mcmcconf, useConjugacy = FALSE)

    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = Cmwtc,
                                   resetFunctions = TRUE)
    startTime <- Sys.time()
    mcmc.out <- nimble::runMCMC(cmcmc,
                                niter = mcmcConfiguration[["n.iterations"]],
                                nchains = mcmcConfiguration[["n.chains"]],
                                nburnin = mcmcConfiguration[["n.burnin"]],
                                #inits = initsList,
                                thin = mcmcConfiguration[["n.thin"]],
                                setSeed = mcmcConfiguration[["setSeed"]],
                                samples = mcmcConfiguration[["samples"]],
                                samplesAsCodaMCMC = mcmcConfiguration[["samplesAsCodaMCMC"]],
                                summary = mcmcConfiguration[["summary"]],
                                WAIC = mcmcConfiguration[["WAIC"]])
    endTime <- Sys.time()
    timeTaken <- difftime(endTime, startTime, units = "secs")

    ret <- list(mcmc.out = mcmc.out,
                timeTaken = timeTaken)

    #save results for MCMC
    retList$mcmc <-  ret
  }

  if(inlaMCMC %in%"importanceSampling"){
    samplerControl$nimbleINLA <- nimbleINLA
    samplerControl$family <- family
    samplerControl$fixedVals <- fixedVals
    samplerControl$timeIndex <- ceiling(mcmcConfiguration[["n.iterations"]]/(mcmcConfiguration[["n.chains"]]))
    samplerControl$nSteps <- mcmcConfiguration[["n.chains"]]
    rr <- inlaIS(mwtc,
                 family,
                 x,
                 y,
                 target,
                 control = samplerControl)

    #compile Model
    compileModel <- compileNimble(mwtc, rr)

    startTime <- Sys.time()
    out <- compileModel$rr$run()
    endTime <- Sys.time()
    timeTaken <- difftime(endTime, startTime, units = "secs")

    mcmc.matrix <- as.matrix(compileModel$rr$mvEWSamples)

    indices <- lapply(as.list(c("gamma", "wts", "logLike")), function(x){
      ret <- startsWith(colnames(mcmc.matrix), x)
      ret <- which(ret == TRUE)
      return(ret)
    })%>%
      do.call("c", .)

    nburnin = mcmcConfiguration[["n.burnin"]]
    res <- mcmc.matrix[-(1:nburnin), -indices]

    #as samplesAsCodaMCMC
    samplesList1 <- coda::as.mcmc.list(coda::as.mcmc(res))

    #save samples and summary
    mcmc.out <- list()
    #Estimating summary
    summaryObject <- lapply(samplesList1, samplesSummary)
    names(summaryObject) <- paste0('chain', 1)
    summaryObject$all.chains <- samplesSummary(do.call('rbind', samplesList1))
    mcmc.out$samples <- samplesList1
    mcmc.out$summary <- summaryObject

    ret <- list(mcmc.out = mcmc.out,
                timeTaken = timeTaken,
                ess = out)

    #save inlamcmc results
    retList$isINLA <- ret

    #colnames(mcmc.matrix)%in%grepl(target, colnames(mcmc.matrix))
    #[, mwtc$expandNodeNames(nodes = c(target, fixedVals))]
    #run function
    # out <-
  }
  return(retList)
}


#' Fitting INLA within NIMBLE
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
INLAWiNim <- function(data,
                      code,
                      fam,
                      modelData,
                      modelConstants,
                      modelInits,
                      parametersToMonitor = c("beta", "sigma", "intercept"),
                      mcmcControl = NULL,
                      mcmcSamplerChange = FALSE,
                      parametersForSamplerChange = NULL,
                      newSampler = NULL,
                      newSamplerControl = NULL,
                      mcmcConfiguration = list(n.chains = 1,
                                               n.iterations = 10,
                                               n.burnin = 0,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = TRUE)){




  initsList <- modelInits()

  #Create the model in nimble
  if(!is.null(newSampler) && newSampler == "hmc"){
    mwtc <- nimble::nimbleModel(code,
                                data = modelData,
                                constants = modelConstants,
                                inits = initsList,
                                buildDerivs = TRUE)

    Cmwtc <- nimble::compileNimble(mwtc,
                                   showCompilerOutput = FALSE)

    mcmcconf <- nimbleHMC::configureHMC(mwtc,
                                        monitors = parametersToMonitor,
                                        enableWAIC = FALSE)

    Rmcmc <- nimble::buildMCMC(mcmcconf)

    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = mwtc)
  }else{
    mwtc <- nimble::nimbleModel(code,
                                data = modelData,
                                constants = modelConstants,
                                inits = initsList)

    # Create the model in C
    Cmwtc <- nimble::compileNimble(mwtc,
                                   showCompilerOutput = FALSE) #Have issues compiling

    if(!is.null(mcmcControl)){
      mcmcconf <- nimble::configureMCMC(Cmwtc,
                                        monitors = parametersToMonitor,
                                        control = mcmcControl,
                                        enableWAIC = FALSE)
    }else{
      mcmcconf <- nimble::configureMCMC(Cmwtc,
                                        monitors = parametersToMonitor,
                                        enableWAIC = FALSE
      )
    }

    if(mcmcSamplerChange == TRUE){
      mcmcconf$removeSamplers(parametersForSamplerChange)
      mcmcconf$addSampler(target = parametersForSamplerChange,
                          type = newSampler,
                          control = newSamplerControl)
      mcmcconf$printSamplers()
      Rmcmc <- nimble::buildMCMC(mcmcconf)
    }
    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = Cmwtc,
                                   resetFunctions = TRUE)
  }

  # Run the MCMC
  startTime <- Sys.time()
  mcmc.out <- nimble::runMCMC(cmcmc,
                              niter = mcmcConfiguration[["n.iterations"]],
                              nchains = mcmcConfiguration[["n.chains"]],
                              nburnin = mcmcConfiguration[["n.burnin"]],
                              #inits = initsList,
                              thin = mcmcConfiguration[["n.thin"]],
                              setSeed = mcmcConfiguration[["setSeed"]],
                              samples = mcmcConfiguration[["samples"]],
                              samplesAsCodaMCMC = mcmcConfiguration[["samplesAsCodaMCMC"]],
                              summary = mcmcConfiguration[["summary"]],
                              WAIC = mcmcConfiguration[["WAIC"]])

  endTime <- Sys.time()
  timeTaken <- difftime(endTime, startTime, units = "secs")
  #Output from the MCMC
  output <- mcmc.out$summary
  output

  #cmcmc$

  # MCMCtrace(object = mcmc.out,
  #           pdata = FALSE,
  #           ind = TRUE,
  #           Rhat = TRUE, # add Rhat
  #           n.eff = TRUE, # add eff sample size
  #           params = "theta")



  #set.seed(1)
  #cmcmc$run(niter = 50)
  #cmcmc$Rmcmc$run(niter=2400, nburnin=0)
  #mcmcconf$printSamplers()
  if(fam == "gaussian"){
    #scales <- cmcmc$samplerFunctions[[1]]$getScaleHistory()
    #accept <- cmcmc$samplerFunctions[[1]]$getAcceptanceHistory()
    scales <- NA
    accept <- NA
    prop_history <- NA
    #prop_history <- cmcmc$samplerFunctions[[1]]$getPropCovHistory()
  }else{
    scales <- NA
    accept <- NA
    prop_history <- NA
  }


  returnList = list(output=output,
                    scales=scales,
                    accept=accept,
                    prop_history=prop_history,
                    mcmc.out=mcmc.out,
                    timeTaken = timeTaken)


  return(returnList)
}
