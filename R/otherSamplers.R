#Define binary sampler for occupancy models

####################################################################
### multivariate binary Gibbs sampler ###########################################
####################################################################
# The intention is to define the INLA function and use it
# when z is accepted
#' @rdname samplers
#' @export
sampler_binary_blockWithINLA <- nimbleFunction(
  name = 'sampler_binary',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    extraVars <- extractControlElement(control, 'extraVars',  character())
    inlaModel <- extractControlElement(control, 'inlaModel',  NULL)
    x <- extractControlElement(control, 'x',  double(2))
    y <- extractControlElement(control, 'y',  double(1))
    N <- nrow(x)
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    d <- length(targetAsScalar)
    #number of extra vars to save from INLA posterior samples
    nVars <- length(extraVars)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
    ## checks
    if(length(targetAsScalar) < 1)  stop('cannot use binary sampler on less than one target node')
    #if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
  },
  run = function() {
    currentLogProb <- model$getLogProb(calcNodes)
    #model[[target]] <<- 1 - model[[target]]
    #This is the difference
    model$simulate(target)
    otherLogProbPrior <- model$calculate(target)
    if(otherLogProbPrior == -Inf) {
      otherLogProb <- otherLogProbPrior
    } else {
      otherLogProb <- otherLogProbPrior + model$calculate(calcNodesNoSelf)
    }
    acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
    jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
    if(jump) {
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      saveINLA(x,y)
       nimCopy(from = model, to = mvSaved, row = 1, nodes = extraVars, logProb = FALSE)
      } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = extraVars, logProb = FALSE)
    }
  },
  methods = list(
    reset = function() { },
    saveINLA = function(x = double(2), y = double(1)){
      inlaRes[1:N, 1:(nVars+1)] <- inlaModel(x, y)
      values(model, extraVars) <<- inlaRes[1, 2:(nVars+1)]
    }
  )
)

####################################################################
### multivariate binary Gibbs sampler ###########################################
####################################################################

#' @rdname samplers
#' @export
sampler_binary_block <- nimbleFunction(
  name = 'sampler_binary',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    d <- length(targetAsScalar)
    #number of extra vars to save from INLA posterior samples
    #nVars <- length(extraVars)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; calcNodesNoSelf <- ccList$calcNodesNoSelf; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch
    ## checks
    if(length(targetAsScalar) < 1)  stop('cannot use binary sampler on less than one target node')
    #if(!model$isBinary(target))     stop('can only use binary sampler on discrete 0/1 (binary) nodes')
  },
  run = function() {
    currentLogProb <- model$getLogProb(calcNodes)
    #model[[target]] <<- 1 - model[[target]]
    #This is the difference
    model$simulate(target)
    otherLogProbPrior <- model$calculate(target)
    if(otherLogProbPrior == -Inf) {
      otherLogProb <- otherLogProbPrior
    } else {
      otherLogProb <- otherLogProbPrior + model$calculate(calcNodesNoSelf)
    }
    acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
    jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
    if(jump) {
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    }
  },
  methods = list(
    reset = function() { }
  )
)
