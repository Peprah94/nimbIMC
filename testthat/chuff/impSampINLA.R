# Importance sampling

##  Contains code to run bootstrap particle filters.
##  We have a build function (buildBootstrapFilter),
##  and step function.
bootStepVirtualUpdate <- nimbleFunctionVirtual(
  run = function(m = integer(),
                 iterRun = integer(0),
                 storeModelValues = double(1),
                 threshNum=double(),
                 prevSamp = logical()) {
    returnType(double(1))
  },
  methods = list(
    returnESS = function() {
      returnType(double())
    }
  )
)

# Bootstrap filter as specified in Doucet & Johnasen '08,
# uses weights from previous time point to calculate likelihood estimate.
bootFStepUpdate <- nimbleFunction(
  name = 'bootFStepUpdate',
  contains = bootStepVirtualUpdate,
  setup = function(model,
                   mvEWSamples,
                   mvWSamples,
                   fixedVals,
                   # inlaModel,
                   x,
                   y,
                   # interInModel,
                   fam
                   ) {

    # setting up parameters
    N <- length(fixedVals)
    ess <- 0

    # calculate dependencies
    calc_thisNode_deps <- model$get

    resamplerFunctionList <- nimbleFunctionList(resamplerVirtual)
    #defaultResamplerFlag <- FALSE
    if(resamplingMethod == 'default'){
      resamplerFunctionList[[1]] <- residualResampleFunction()
      defaultResamplerFlag <- TRUE
    }
    if(resamplingMethod == 'residual')
      resamplerFunctionList[[1]] <- residualResampleFunction()
    if(resamplingMethod == 'multinomial')
      resamplerFunctionList[[1]] <- multinomialResampleFunction()
    if(resamplingMethod == 'stratified')
      resamplerFunctionList[[1]] <- stratifiedResampleFunction()
    if(resamplingMethod == 'systematic')
      resamplerFunctionList[[1]] <- systematicResampleFunction()
  },
  run = function(m = integer(),
                 beta = double(1),
                 prevSamp = logical()) {
    returnType(double(1))

    res <- nimbleINLA(x, y, beta= beta, fixedVals,  family = fam)

    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    llEst <- numeric(m, init=FALSE)

    llEst <- res[1,1]

    for(i in 1:m){
    # For now, simulate beta's from priors
    model$simulate(beta)

    ## The logProbs of calc_thisNode_self are, correctly, not calculated.
    nimCopy(model, mvWSamples, nodes = beta, nodesTo = beta, row = i)

    #calculate
    wts[i]  <- model$calculate(calc_thisNode_deps)
    }
    if(lll == -Inf){
      copy(mvEWSamples, model, nodes = fixedVals, row = 1)
    }else{
      saveResults(fixedVals, res)
      copy( model, mvEWSamples, nodes = fixedVals, row = 1)
    }



  },
  methods = list(
    returnESS = function(){
      returnType(double(0))
      return(ess)
    }
  )
)

#' Create an updated bootstrap particle filter algorithm to estimate log-likelihood.
#'
#'@description Create an updated bootstrap particle filter algorithm for a given NIMBLE state space model.
#'
#' @param model A nimble model object, typically representing a state
#'  space model or a hidden Markov model.
#' @param nodes  A character vector specifying the latent model nodes
#'  over which the particle filter will stochastically integrate to
#'  estimate the log-likelihood function.  All provided nodes must be stochastic.
#'  Can be one of three forms: a variable name, in which case all elements in the variable
#'  are taken to be latent (e.g., 'x'); an indexed variable, in which case all indexed elements are taken
#'  to be latent (e.g., 'x[1:100]' or 'x[1:100, 1:2]'); or a vector of multiple nodes, one per time point,
#'  in increasing time order (e.g., c("x[1:2, 1]", "x[1:2, 2]", "x[1:2, 3]", "x[1:2, 4]")).
#' @param mvSamplesEst  A modelValue object contained posterior samples from the reduced model using MCMC.
#' @param control  A list specifying different control options for the particle filter.  Options are described in the details section below.
#' @author Kwaku Peprah Adjei
#' @details
#'
#' Each of the \code{control()} list options are described in detail here:
#' \describe{
#'  \item{thresh}{ A number between 0 and 1 specifying when to resample: the resampling step will occur when the
#'   effective sample size is less than \code{thresh} times the number of particles. Defaults to 0.8. Note that at the last time step, resampling will always occur so that the \code{mvEWsamples} \code{modelValues} contains equally-weighted samples.}
#'  \item{resamplingMethod}{The type of resampling algorithm to be used within the particle filter. Can choose between \code{'default'} (which uses NIMBLE's \code{rankSample()} function),  \code{'systematic'}, \code{'stratified'}, \code{'residual'}, and \code{'multinomial'}.  Defaults to \code{'default'}.  Resampling methods other than \code{'default'} are currently experimental.}
#'  \item{saveAll}{Indicates whether to save state samples for all time points (TRUE), or only for the most recent time point (FALSE)}
#'  \item{smoothing}{Decides whether to save smoothed estimates of latent states, i.e., samples from f(x[1:t]|y[1:t]) if \code{smoothing = TRUE}, or instead to save filtered samples from f(x[t]|y[1:t]) if \code{smoothing = FALSE}.  \code{smoothing = TRUE} only works if \code{saveAll = TRUE}.}
#'  \item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.
#'  Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
#'  \item{initModel}{A logical value indicating whether to initialize the model before running the filtering algorithm.  Defaults to TRUE.}
#' }
#'  \item{iNodePrev}{An integer specifying the number of years used to fit the reduced model.}
#' }
#'
#'  The updated bootstrap filter starts by copying saved MCMC results into a
#'  modelValue object, and calculate weights for times t<= iNodePrev. For time iNode>iNodePrev,
#'  the updated bootstrap filter generates a sample of estimates from the
#'  prior distribution of the latent states of a state space model.  At each of these time point, these particles are propagated forward
#'  by the model's transition equation.  Each particle is then given a weight
#'  proportional to the value of the observation equation given that particle.
#'  The weights are used to draw an equally-weighted sample of the particles at this time point.
#'  The algorithm then proceeds
#'  to the next time point.  Neither the transition nor the observation equations are required to
#'  be normal for the bootstrap filter to work.
#'
#'  The resulting specialized particle filter algorithm will accept a
#'  single integer argument (\code{m}, default 10,000), which specifies the number
#'  of random \'particles\' to use for estimating the log-likelihood.  The algorithm
#'  returns the estimated log-likelihood value, and saves
#'  unequally weighted samples from the posterior distribution of the latent
#'  states in the \code{mvWSamples} modelValues object, with corresponding logged weights in \code{mvWSamples['wts',]}.
#'  An equally weighted sample from the posterior can be found in the \code{mvEWSamples} \code{modelValues} object.
#'
#'  Note that if the \code{thresh} argument is set to a value less than 1, resampling may not take place at every time point.
#'  At time points for which resampling did not take place, \code{mvEWSamples} will not contain equally weighted samples.
#'  To ensure equally weighted samples in the case that \code{thresh < 1}, we recommend resampling from \code{mvWSamples} at each time point
#'  after the filter has been run, rather than using \code{mvEWSamples}.
#'
#' @section \code{returnESS()} Method:
#'  Calling the \code{returnESS()} method of a bootstrap filter after that filter has been \code{run()} for a given model will return a vector of ESS (effective
#'  sample size) values, one value for each time point.
#'
#' @export
#'
#' @family smc update methods
#' @references Michaud, N., de Valpine, P., Turek, D., Paciorek, C. J., & Nguyen, D. (2021). Sequential Monte Carlo methods in the nimble and nimbleSMC R packages. \emph{Journal of Statistical Software}. 100, 1-39.
#' @examples
#' ## For illustration only.
#' stateSpaceCode <- nimbleCode({
#'   x0 ~ dnorm(0, var = 1)
#'   x[1] ~ dnorm(a* x0, var = 1)
#'   y[1] ~ dnorm(x[1]*c, var = .5)
#'   for(t in 2:10){
#'     x[t] ~ dnorm(a * x[t-1], var = 1)
#'     y[t] ~ dnorm(x[t]*c, var = .5)
#'   }
#' })
#'
#' model <- nimbleModel(code = exampleCode, data = list(y = rnorm(10)),
#'                      inits = list(x0 = 0, x = rnorm(10)))
#' my_BootF <- buildBootstrapFilterUpdate(estimationModel,
#' latent,
#' mvSamplesEst = mvSamplesEst,
#' target = target,
#' control = pfControl)
#' ## Now compile and run, e.g.,
#' ## targetAsScalar <- estimationModel$expandNodeNames(target, returnScalarComponents = TRUE)
#' ## compiledParticleFilter <- compileNimble(estimationModel,  particleFilterEst)
#' ## logLik <- compiledParticleFilter$particleFilterEst$run(m = 2000, iterRun = 1, storeModelValues = values(estimationModel, targetAsScalar))
#' ## ESS <- compiledParticleFilter$particleFilterEst$returnESS()
#' ## boot_X <- as.matrix(compiledParticleFilter$particleFilterEst$mvEWSamples, 'x')

inlaIS <- nimbleFunction(
  name = 'inlaIS',
  setup = function(model,
                   nodes,
                   mvSamplesEst = list(),
                   target,
                   control = list()) {

    #control list extraction
    thresh <- control[['thresh']]
    saveAll <- control[['saveAll']]
    smoothing <- control[['smoothing']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    initModel <- control[['initModel']]
    iNodePrev <- control[['iNodePrev']] #how many time points were used to fit the reduced model
    #initModel <- control[['initModel']]
    #M <- control[['M']]
    resamplingMethod <- control[['resamplingMethod']]
    if(is.null(thresh)) thresh <- 0.8
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(initModel)) initModel <- TRUE
    if(is.null(resamplingMethod)) resamplingMethod <- 'default'
    if(!(resamplingMethod %in% c('default', 'multinomial', 'systematic', 'stratified',
                                 'residual')))
      stop('buildBootstrapFilter: "resamplingMethod" must be one of: "default", "multinomial", "systematic", "stratified", or "residual". ')
    ## latent state info
    targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    #storeModelValues <- values(model, targetNodesAsScalar)
    nodes <- findLatentNodes(model, nodes, timeIndex)
    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1)
      stop('buildBootstrapFilter: sizes or dimensions of latent states varies.')
    vars <- model$getVarNames(nodes =  nodes)

    my_initializeModel <- initializeModel(model, silent = silent)

    if(0 > thresh || 1 < thresh || !is.numeric(thresh))
      stop('buildBootstrapFilter: "thresh" must be between 0 and 1.')
    if(!saveAll & smoothing)
      stop("buildBootstrapFilter: must have 'saveAll = TRUE' for smoothing to work.")

    if(!is.numeric(iNodePrev) || iNodePrev < 0)
      stop("buildBootstrapFilter: must have 'iNodePrev' numeric and greater than 0")
    ## Create mv variables for x state and sampled x states.  If saveAll=TRUE,
    ## the sampled x states will be recorded at each time point.
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]

    if(saveAll){
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))

      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      latent <- names
      names <- c(names, "wts", "bootLL")
      type <- c(type, "double", "double")
      size$wts <- length(dims)
      size$bootLL <- length(dims)
      ##  Only need one weight per particle (at time T) if smoothing == TRUE.
      if(smoothing){
        size$wts <- 1
        size$bootLL <- 1
      }
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))

    }else{
      names <- sapply(modelSymbolObjects, function(x)return(x$name))
      type <- sapply(modelSymbolObjects, function(x)return(x$type))
      size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])

      mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      latent <- names
      names <- c(names, "wts", "bootLL")
      type <- c(type, "double", "double")
      size$wts <- 1
      size$bootLL <- 1
      mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
      names <- names[1]
    }


    bootStepFunctions <- nimbleFunctionList(bootStepVirtualUpdate)
    for(iNode in seq_along(nodes)){
      bootStepFunctions[[iNode]] <- bootFStepUpdate(model, mvEWSamples, mvWSamples,
                                                    nodes, iNode, names, saveAll,
                                                    smoothing, resamplingMethod,
                                                    silent,
                                                    iNodePrev, latent,
                                                    target, mvSamplesEst)
    }

    essVals <- rep(0, length(nodes))
    lastLogLik <- -Inf
    runTime <- 1
  },
  run = function(m = integer(default = 10000),
                 beta = double(1)
  ) {
    returnType(double())

    if(initModel) my_initializeModel$run()
    resize(mvWSamples, m)
    resize(mvEWSamples, m)
    threshNum <- ceiling(thresh*m)
    logL <- 0
    ## prevSamp indicates whether resampling took place at the previous time point.
    prevSamp <- 1
    for(iNode in seq_along(bootStepFunctions)) {
      if(iNode == length(bootStepFunctions)){
        threshNum <- m  ## always resample at last time step so mvEWsamples is equally-weighted
      }
      out <- bootStepFunctions[[iNode]]$run(m, iterRun, storeModelValues, threshNum, prevSamp)
      logL <- logL + out[1]
      prevSamp <- out[2]
      #print(iNode)
      essVals[iNode] <<- bootStepFunctions[[iNode]]$returnESS()
      if(logL == -Inf) {lastLogLik <<- logL; return(logL)}
      if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
      if(logL == Inf)  {lastLogLik <<- -Inf; return(-Inf)}
    }
    lastLogLik <<- logL
    runTime <<- iterRun
    return(logL)
  },
  methods = list(
    getLastLogLik = function() {
      return(lastLogLik)
      returnType(double())
    },
    setLastLogLik = function(lll = double()) {
      lastLogLik <<- lll
    },
    getLastTimeRan = function() {
      return(runTime)
      returnType(integer())
    },
    setLastTimeRan = function(lll = integer()) {
      runTime <<- lll + 1
    },
    returnESS = function(){
      returnType(double(1))
      return(essVals)
    }
  )
)
