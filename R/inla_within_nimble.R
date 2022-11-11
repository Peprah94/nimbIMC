#' Fitting INLA within NIMBLE
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
fit.inla <- function(x ,
                     y ,
                     beta,
                     interInModel,
                     family
                     ){
  #ii  <- get("ii",envir =  parent.frame())
  #ii <- assign("ii",ii+1,envir = parent.frame())
  #print(ii)
  data <- list(y=y, x=x)
  data$oset = data$x %*% (beta)

  if(interInModel == 1){
formula =  1 + offset(data$oset)
  }else{
formula = - 1 + offset(data$oset)
  }

  res = INLA::inla(y ~ formula,
                   data = data,
                   family = family,
                   verbose=FALSE,
                   control.predictor = list(compute = TRUE))
  fitted_values = res$summary.fitted.values[,"mean"]
  #intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
 if(interInModel == 1){
   intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
 }else{
   intercept = 1
 }
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  ret <- cbind(intercept, precision, fitted_values)
  return(ret)
}



#' Compiling fit.inla in C+
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    interInModel = double(0, default = 1),
    family = character(0, default = "gaussian")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

#' Conditional INLA for missing covariate
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return dataframe
#' @export
fitINLAMissingValues <- function(missData, indexMiss,  eta){
  age = missData[,1]
  bmi = missData[,2]
  chl = missData[,3]
missData <- data.frame(age= age, bmi = bmi, chl = chl)
  missData[indexMiss,2] = eta

  res = inla(chl ~ 1 + bmi + as.factor(age), data = missData)

  beta0 = res$marginals.fixed[[1]]
  beta1 = res$marginals.fixed[[2]]
  beta2 = res$marginals.fixed[[3]]
  tau = res$marginals.hyperpar[[1]]
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  estbeta1 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[2]])
  estbeta2 = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[3]])
  fitted_values = res$summary.fitted.values[,"mean"]
  ret <- cbind(intercept, precision, estbeta1, estbeta2,
               fitted_values)
  return(ret)
}


#' Compiling fitINLAMissingValues in C+
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export
nimbleINLAMissingValues <- nimble::nimbleRcall(
  prototype = function(
     missData = double(2), #x is a matrix
    indexMiss = double(1),
    eta=double(1)
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fitINLAMissingValues'
)

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
                      parametersToMonitor = c("beta", "tau", "intercept"),
                      mcmcConfiguration = list(n.chains = 1,
                                               n.iterations = 10,
                                               n.burnin = 0,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = FALSE)){



initsList <- modelInits()
#initsList <- idm_inits()


#Create the model in nimble
mwtc <- nimble::nimbleModel(code,
                            data = modelData,
                    constants = modelConstants,
                    inits = initsList)

# Create the model in C
Cmwtc <- nimble::compileNimble(mwtc,
                               showCompilerOutput = FALSE) #Have issues compiling


mcmcconf <- nimble::configureMCMC(Cmwtc, monitors = parametersToMonitor,
                                  control = list(scale = 0.75,
                                                 adaptive = TRUE,
                                                 propCov='identity',
                                                 adaptInterval=50))

Rmcmc <- nimble::buildMCMC(mcmcconf)

# Compile
cmcmc <- nimble::compileNimble(Rmcmc,
                       project = Cmwtc, resetFunctions = TRUE)

#MCMC Configurations


# Run the MCMC
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
                  mcmc.out=mcmc.out)
return(returnList)
}


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
                      code,
                      fam,
                      modelData,
                      modelConstants,
                      modelInits,
                      parametersToMonitor = c("omega"),
                      mcmcConfiguration = list(n.chains = 1,
                                               n.iterations = 10,
                                               n.burnin = 0,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = FALSE)){



  initsList <- modelInits()
  #initsList <- idm_inits()


  #Create the model in nimble
  mwtc <- nimble::nimbleModel(code,
                              data = modelData,
                              constants = modelConstants,
                              inits = initsList)

  # Create the model in C
  Cmwtc <- nimble::compileNimble(mwtc,
                                 showCompilerOutput = FALSE) #Have issues compiling

  mcmcconf <- nimble::configureMCMC(Cmwtc,
                                    monitors = parametersToMonitor,
                                    print = TRUE,
                                    useConjugacy=FALSE)

  mcmcconf$removeSamplers(parametersToMonitor)
  #mcmcconf$addSampler(c("omega[1,1:2]"), "myRW_dirichlet")
  #mcmcconf$addSampler(c("omega[2,1:2]"), "myRW_dirichlet")
  mcmcconf$addSampler(c("omega"), "myRW_dirichlet")

  Rmcmc <- nimble::buildMCMC(mcmcconf)

  # Compile
  cmcmc <- nimble::compileNimble(Rmcmc,
                                 project = Cmwtc,
                                 resetFunctions = TRUE)

  #MCMC Configurations


  # Run the MCMC
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
                    mcmc.out=mcmc.out)
  return(returnList)
}




