
#' Fitting INLA within NIMBLE for Data Generating process
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
INLAWiNimDataGeneratingTargetDivide <- function(data,
                                    covariate,
                                    code,
                                    family,
                                    modelData,
                                    modelConstants,
                                    modelInits,
                                    nimbleINLA,
                                    inlaMCMC = c("inla", "mcmc", "inlamcmc"),
                                    inlaMCsampler = "RW_INLA_block",
                                    samplerControl = list(),
                                    parametersToMonitor = list(mcmc = c("mcmc"),
                                                               mcmc2 = c("mcmc"),
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
  targetMCMC <- parametersToMonitor$mcmcINLA
  additionalPars <- parametersToMonitor$additionalPars
  #betaWts <- parametersToMonitor$betaWts
  #latentWts <- parametersToMonitor$latentWts

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
                                      nodes = c(target))

    # setting sampler controls
    samplerControl$fit.inla = nimbleINLA
    samplerControl$x = x
    samplerControl$y = y
    samplerControl$fixedVals = fixedVals
    samplerControl$fam = family
    #samplerControl$targetMCMC = targetMCMC

    # mcmc configuration
    mcmcconf$addSampler(target = targetMCMC,
                        type = inlaMCsampler,
                        control = samplerControl)
    #mcmcconf$addSampler(target = target[1],
     #                   type = "RW_block")

    mcmcconf$printSamplers(executionOrder = TRUE)
    mcmcconf$addMonitors(c(target, targetMCMC))

    if(!is.null(parametersToMonitor$additionalPars)){
      mcmcconf$addMonitors(parametersToMonitor$additionalPars)
    }

    #build model
    Rmcmc <- nimble::buildMCMC(mcmcconf)

    # Compile
    cmcmc <- nimble::compileNimble(Rmcmc,
                                   project = Cmwtc)
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
                                      monitors = c(target, fixedVals, targetMCMC))

    mcmcconf$printSamplers()
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
    if(!is.null(additionalPars)){samplerControl$additionalPars <- additionalPars}
    if(is.null(samplerControl$latentIsDependent)) samplerControl$latentIsDependent <- FALSE

    # there are instances, such as spatio-temporal effects where we want to return samples of linear effect
    # we then need to extract such variables
    # if(!is.null(samplerControl$spatioTemporal)){
    #   samplerControl$linearPred <- Cmwtc$getDependencies(nodes = fixedVals, includeData = FALSE, self = FALSE, determOnly = TRUE)
    #   samplerControl$returnLinearPred <- TRUE
    # } else {
    #   samplerControl$returnLinearPred <- FALSE
    #   samplerControl$linearPred <- NULL
    #   }


      # if(!is.null(latentWts)) samplerControl$latentWts <- latentWts
   # if(!is.null(betaWts)) samplerControl$betaWts <- betaWts
    rr <- inlaISmultiple(mwtc,
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


  # Resample latent vars with their weights
    additionalPars <- colnames(res)[startsWith(colnames(res), additionalPars)]
    latentNames <- mwtc$expandNodeNames(target)[mwtc$isDiscrete(target)]
    latentRes <- res[ , colnames(res) %in% c(latentNames, additionalPars, "latentWeights[1]")]

    #resample latent variable with their weights
    #Extract weights and obtain sampling index
    latentWtsIndex <- ncol(latentRes)
    maxWts <- max(latentRes[, latentWtsIndex ])
    nWeights <- exp(latentRes[, latentWtsIndex ] - maxWts)
    nWeights <- nWeights/sum(nWeights)
    latentSampIndex <- sample(1:nrow(res), prob = nWeights, replace = TRUE)
    latentRes <- latentRes[latentSampIndex, -latentWtsIndex]

    otherRes <- res[ , !colnames(res) %in% c(latentNames, additionalPars, "latentWeights[1]", "betaWts[1]")]

    #put all together
    res <- cbind(otherRes, latentRes)

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


#create matrix

