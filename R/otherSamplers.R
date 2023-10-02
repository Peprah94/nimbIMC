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
