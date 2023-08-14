
####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################

mcmc_determineCalcAndCopyNodes <- function(model, target) {
  targetExpanded <- model$expandNodeNames(target)
  modelPredictiveNodes <- model$getNodeNames(predictiveOnly = TRUE)
  targetExpandedPPbool <- targetExpanded %in% modelPredictiveNodes
  targetAllPP <- all(targetExpandedPPbool)
  targetAnyPP <- any(targetExpandedPPbool)
  ## if a particular sampler is assigned *jointly to PP and non-PP* nodes, then we're going to bail
  ## out and quit, if the option MCMCusePredictiveDependenciesInCalculations == FALSE.
  ## this is an extreme corner-case, which I think will lead to problems.
  if(targetAnyPP && !targetAllPP && !getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
    stop('cannot assign samplers jointly to posterior predictive (PP) nodes and non-PP nodes, when MCMCusePredictiveDependenciesInCalculations option is FALSE', call. = FALSE)
  ## if the sampler calling this, itself, is operating exclusively on posterior predictive nodes,
  ## then regardless of how the rest of the model is being sampled (w.r.t. inclusion of posterior predictive nodes),
  ## we'll include 'self' and all stochastic dependencies (the full markov blanket) in the calculations,
  ## which necessarily are taking place entirely within a posterior predictive network of nodes.
  ## this should lead to correct behaviour (consistent samples and joint posteriors) in all cases.
  if(targetAllPP) {
    ## when sampler is operating only on posterior predictive nodes,
    ## then always include all predictive dependencies:
    calcNodes <- model$getDependencies(target, includePredictive = TRUE)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE, includePredictive = TRUE)
    ##calcNodesPPomitted <- character()
  } else {
    ## usual case:
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ##calcNodesPPomitted <- setdiff(model$getDependencies(target, includePredictive = TRUE), calcNodes)
  }
  ## copyNodes:
  copyNodes <- model$getDependencies(target, self = FALSE)
  isStochCopyNodes <- model$isStoch(copyNodes)
  copyNodesDeterm <- copyNodes[!isStochCopyNodes]
  copyNodesStoch <- copyNodes[isStochCopyNodes]
  ##
  ccList <- list(
    calcNodes = calcNodes,
    calcNodesNoSelf = calcNodesNoSelf,
    ##calcNodesPPomitted = calcNodesPPomitted,
    copyNodesDeterm = copyNodesDeterm,
    copyNodesStoch = copyNodesStoch
  )
  return(ccList)
}


#' @rdname samplers
#' @export
#######################################################################################
### RW_PF_block, does a block RW, but using a particle filter likelihood function #####
#######################################################################################
# used to sample with AFSS using INLA embedded code that returns marginal likehood to write the sampler
# for continuous and discrete random variables
# reference: Li, Y., Linero, A., & Walker, S. G. (2023). Latent uniform samplers on multivariate binary spaces. Statistics and Computing, 33(5), 102.
#' @rdname samplers
#' @export
sampler_AFSS_INLA_block <- nimbleFunction(
  name = 'sampler_AFSS_INLA_block',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    x <- extractControlElement(control, 'x',  double())
    y <- extractControlElement(control, 'y',  character())
    fixedVals <- extractControlElement(control, 'fixedVals',  double())
    fam <- extractControlElement(control, 'fam',  "gaussian")
    interVal <- extractControlElement(control, 'interInModel',  1)
    existingINLA          <- extractControlElement(control, 'fit.inla',                   NULL)
    ## control list extraction
    widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
    targetMCMC               <- extractControlElement(control, 'targetMCMC',              NULL)
    maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            100)
    adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
    adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
    adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
    adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
    maxContractions        <- extractControlElement(control, 'maxContractions',          1000)
    maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
    eps <- 1e-15

    #Extract y values
    #nodesY <- model$expandNodeNames(y)
    #y <- nimble::values(model, nodesY)
    #y <- c(model[[y]])
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    #calcNodesDepStage1 <- model$getVarNames(nodes = calcNodesDepStage)
    #calcNodesDepStage <- calcNodesDepStage[!model$isData(calcNodesDepStage)]
    #Note that model likelihood will be estimated from INLA
    # so we take the data out
    ## numeric value generation
    d                  <- length(targetAsScalar)
    discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
    anyDiscrete        <- any(discrete)
    gammaMatrix        <- diag(d)         # matrix of orthogonal bases
    if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
    widthVecOriginal   <- widthVec
    nExpansions        <- rep(0, d)       # number of expansions
    nContracts         <- rep(0, d)       # number of contractions
    adaptFactorMaxIterOriginal <- adaptFactorMaxIter
    factorCounter      <- 0               # number of iterations since last factor adaptation
    factorTimesAdapted <- 0               # number of times factors have adapted
    empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
    empirCov           <- diag(d)
    allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
    widthCounter       <- 0               # number of iterations since last width adaptation
    adaptWidthMaxIterOriginal <- adaptWidthMaxIter
    adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
    widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
    ## checks
    if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
    if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
      stop('sliceWidths must be a numeric vector')
    if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)

    #set up for INLA estimated vals
    if(length(model$expandNodeNames(fixedVals)) > 0){
      latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)


    my_particleFilter <- buildINLAmodel(model, fam, x,y = y,control = list(fit.inla = existingINLA,fixedVals = fixedVals))

  #Target values for inla
  if(is.null(targetMCMC)){targetINLA <- targetAsScalar
  }else{
    targetMCMCasScalar <- model$expandNodeNames(targetMCMC, returnScalarComponents = TRUE)
      targetINLA <- c(targetMCMCasScalar, targetAsScalar)}

    #targetVal <- values(model, targetAsScalar)
    particleMV <- my_particleFilter$mvEWSamples
out <- -Inf
    print(latentSamp)
  },
  run = function() {
    #out <<- my_particleFilter$run(beta = values(model,  target))
    maxContractionsReached <- FALSE
    for(i in 1:d) {
      eigenVec <- gammaMatrix[, i]
      width <- widthVec[i]
      u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
      x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
      Lbound <- -1.0 * runif(1, 0, 1) * width
      Rbound <- Lbound + width
      L <- x0 + Lbound * eigenVec
      R <- x0 + Rbound * eigenVec
      maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
      maxStepsR <- maxSteps - 1 - maxStepsL
      lp <- setAndCalculateTarget(L) #+ out
      while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
        Lbound <- Lbound - width
        L <- x0 + Lbound * eigenVec
        lp <- setAndCalculateTarget(L) #+ out
        maxStepsL <- maxStepsL - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      lp <- setAndCalculateTarget(R) + out
      while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
        Rbound <- Rbound + width
        R <- x0 + Rbound * eigenVec
        lp <- setAndCalculateTarget(R) #+ out
        maxStepsR <- maxStepsR - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
      x1 <- x0 + prop * eigenVec
      lp <- setAndCalculateTarget(x1)
      numContractions <- 0
      while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
        ## The checks for Rbound - Lbound small and max number of contractions are
        ## for cases where model is in invalid state and lp calculations are NA/NaN or where
        ## interval contracts to zero
        if(prop < 0) { Lbound <- prop }
        else         { Rbound <- prop }
        nContracts[i] <<- nContracts[i] + 1
        prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
        x1 <- x0 + prop * eigenVec
        lp <- setAndCalculateTarget(x1) #+ out
        numContractions <- numContractions + 1
      }
      if(Rbound - Lbound <= eps | numContractions == maxContractions)
        maxContractionsReached <- TRUE
    }
    #print(out)
    if(maxContractionsReached) {
      if(maxContractionsWarning)
        cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = fixedValsDep, logProb = TRUE)
    } else {
      out <<- my_particleFilter$run(beta = values(model,  targetINLA))
      copy(particleMV, model, fixedVals, fixedVals, row = 1)
      #calculate(model, fixedValsDep)
      model$calculate()
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    if(allWidthsAdapted == 0)   adaptWidths()
    if(adaptFactorMaxIter > 0)  adaptFactors()

  },
  methods = list(
    setAndCalculateTarget = function(targetValues = double(1)) {
      if(anyDiscrete == 1)
        for(i in 1:d)
          if(discrete[i] == 1)   targetValues[i] <- floor(targetValues[i])
      values(model, target) <<- targetValues
      lp <- model$calculate(calcNodesProposalStage)
      if(lp == -Inf) return(lp)
      lp <- lp + model$calculate(calcNodesDepStage) #+ out
      returnType(double())
      return(lp)
    },
    adaptFactors = function() {
      adaptFactorMaxIter <<- adaptFactorMaxIter - 1
      factorCounter <<- factorCounter + 1
      empirSamp[factorCounter, 1:d] <<- values(model, target)
      if(factorCounter == adaptFactorInterval) {
        for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
        gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
        factorTimesAdapted <<- factorTimesAdapted + 1
        factorCounter      <<- 0
        nExpansions        <<- rep(0, d)
        nContracts         <<- rep(0, d)
        allWidthsAdapted   <<- 0
        widthCounter       <<- 0
        adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
        adaptWidthInterval <<- 1
        widthIndicatorVec  <<- rep(1, d)
      }
    },
    adaptWidths = function() {
      adaptWidthMaxIter <<- adaptWidthMaxIter - 1
      widthCounter <<- widthCounter + 1
      if(widthCounter == adaptWidthInterval) {
        for(i in 1:d) {
          if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
            if(nExpansions[i] == 0)   nExpansions[i] <<- 1
            widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
            widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
            adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
            nExpansions[i] <<- 0
            nContracts[i] <<- 0
            if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
              widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
          }
        }
        allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
        widthCounter     <<- 0
      }
      if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
        allWidthsAdapted <<- 1
    },
    reset = function() {
      gammaMatrix        <<- diag(d)
      empirCov           <<- diag(d)
      widthVec           <<- widthVecOriginal
      nExpansions        <<- rep(0, d)
      nContracts         <<- rep(0, d)
      adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
      factorCounter      <<- 0
      factorTimesAdapted <<- 0
      allWidthsAdapted   <<- 0
      widthCounter       <<- 0
      adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
      adaptWidthInterval <<- 1
      widthIndicatorVec  <<- rep(1, d)
    }
  )
)


#' @rdname samplers
#' @export
#######################################################################################
### RW_PF_block, does a block RW, but using a particle filter likelihood function #####
#######################################################################################
# used to sample with AFSS using INLA embedded code that returns marginal likehood to write the sampler
# for continuous and discrete random variables
# reference: Li, Y., Linero, A., & Walker, S. G. (2023). Latent uniform samplers on multivariate binary spaces. Statistics and Computing, 33(5), 102.
#' @rdname samplers
#' @export
sampler_AFSS_INLA_blockV2 <- nimbleFunction(
  name = 'sampler_AFSS_INLA_blockV2',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    x <- extractControlElement(control, 'x',  double())
    y <- extractControlElement(control, 'y',  character())
    fixedVals <- extractControlElement(control, 'fixedVals',  double())
    fam <- extractControlElement(control, 'fam',  "gaussian")
    interVal <- extractControlElement(control, 'interInModel',  1)
    existingINLA          <- extractControlElement(control, 'fit.inla',                   NULL)
    ## control list extraction
    widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
    targetMCMC               <- extractControlElement(control, 'targetMCMC',              NULL)
    extraVars               <- extractControlElement(control, 'extraVars',              NULL)
    maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            100)
    adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
    adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
    adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
    adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
    maxContractions        <- extractControlElement(control, 'maxContractions',          1000)
    maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
    eps <- 1e-15

    #Extract y values
    #nodesY <- model$expandNodeNames(y)
    #y <- nimble::values(model, nodesY)
    #y <- c(model[[y]])
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    #calcNodesDepStage1 <- model$getVarNames(nodes = calcNodesDepStage)
    #calcNodesDepStage <- calcNodesDepStage[!model$isData(calcNodesDepStage)]
    #Note that model likelihood will be estimated from INLA
    # so we take the data out
    ## numeric value generation
    d                  <- length(targetAsScalar)
    discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
    anyDiscrete        <- any(discrete)
    gammaMatrix        <- diag(d)         # matrix of orthogonal bases
    if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
    widthVecOriginal   <- widthVec
    nExpansions        <- rep(0, d)       # number of expansions
    nContracts         <- rep(0, d)       # number of contractions
    adaptFactorMaxIterOriginal <- adaptFactorMaxIter
    factorCounter      <- 0               # number of iterations since last factor adaptation
    factorTimesAdapted <- 0               # number of times factors have adapted
    empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
    empirCov           <- diag(d)
    allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
    widthCounter       <- 0               # number of iterations since last width adaptation
    adaptWidthMaxIterOriginal <- adaptWidthMaxIter
    adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
    widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
    ## checks
    if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
    if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
      stop('sliceWidths must be a numeric vector')
    if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)

    #set up for INLA estimated vals
    if(length(model$expandNodeNames(fixedVals)) > 0){
      latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)


    my_particleFilter <- buildINLAmodelV2(model, fam, x,y = y,control = list(fit.inla = existingINLA,fixedVals = fixedVals))

    extraVal <- model$expandNodeNames(extraVars, returnScalarComponents = TRUE)
    #Target values for inla
    if(is.null(targetMCMC)){targetINLA <- targetAsScalar
    }else{
      targetMCMCasScalar <- model$expandNodeNames(targetMCMC, returnScalarComponents = TRUE)
      targetINLA <- c(targetMCMCasScalar, targetAsScalar)}

    #targetVal <- values(model, targetAsScalar)
    particleMV <- my_particleFilter$mvEWSamples
    out <- -Inf
    print(latentSamp)
  },
  run = function() {
    #out <<- my_particleFilter$run(beta = values(model,  target))
    maxContractionsReached <- FALSE
    for(i in 1:d) {
      eigenVec <- gammaMatrix[, i]
      width <- widthVec[i]
      u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
      x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
      Lbound <- -1.0 * runif(1, 0, 1) * width
      Rbound <- Lbound + width
      L <- x0 + Lbound * eigenVec
      R <- x0 + Rbound * eigenVec
      maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
      maxStepsR <- maxSteps - 1 - maxStepsL
      lp <- setAndCalculateTarget(L) #+ out
      while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
        Lbound <- Lbound - width
        L <- x0 + Lbound * eigenVec
        lp <- setAndCalculateTarget(L) #+ out
        maxStepsL <- maxStepsL - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      lp <- setAndCalculateTarget(R) + out
      while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
        Rbound <- Rbound + width
        R <- x0 + Rbound * eigenVec
        lp <- setAndCalculateTarget(R) #+ out
        maxStepsR <- maxStepsR - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
      x1 <- x0 + prop * eigenVec
      lp <- setAndCalculateTarget(x1)
      numContractions <- 0
      while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
        ## The checks for Rbound - Lbound small and max number of contractions are
        ## for cases where model is in invalid state and lp calculations are NA/NaN or where
        ## interval contracts to zero
        if(prop < 0) { Lbound <- prop }
        else         { Rbound <- prop }
        nContracts[i] <<- nContracts[i] + 1
        prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
        x1 <- x0 + prop * eigenVec
        lp <- setAndCalculateTarget(x1) #+ out
        numContractions <- numContractions + 1
      }
      if(Rbound - Lbound <= eps | numContractions == maxContractions)
        maxContractionsReached <- TRUE
    }
    #print(out)
    if(maxContractionsReached) {
      if(maxContractionsWarning)
        cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = fixedValsDep, logProb = TRUE)
    } else {
      out <<- my_particleFilter$run(beta = values(model,  targetINLA), extraVars = values(model, extraVal))
      copy(particleMV, model, fixedVals, fixedVals, row = 1)
      #calculate(model, fixedValsDep)
      model$calculate()
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    if(allWidthsAdapted == 0)   adaptWidths()
    if(adaptFactorMaxIter > 0)  adaptFactors()

  },
  methods = list(
    setAndCalculateTarget = function(targetValues = double(1)) {
      if(anyDiscrete == 1)
        for(i in 1:d)
          if(discrete[i] == 1)   targetValues[i] <- floor(targetValues[i])
      values(model, target) <<- targetValues
      lp <- model$calculate(calcNodesProposalStage)
      if(lp == -Inf) return(lp)
      lp <- lp + model$calculate(calcNodesDepStage) #+ out
      returnType(double())
      return(lp)
    },
    adaptFactors = function() {
      adaptFactorMaxIter <<- adaptFactorMaxIter - 1
      factorCounter <<- factorCounter + 1
      empirSamp[factorCounter, 1:d] <<- values(model, target)
      if(factorCounter == adaptFactorInterval) {
        for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
        gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
        factorTimesAdapted <<- factorTimesAdapted + 1
        factorCounter      <<- 0
        nExpansions        <<- rep(0, d)
        nContracts         <<- rep(0, d)
        allWidthsAdapted   <<- 0
        widthCounter       <<- 0
        adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
        adaptWidthInterval <<- 1
        widthIndicatorVec  <<- rep(1, d)
      }
    },
    adaptWidths = function() {
      adaptWidthMaxIter <<- adaptWidthMaxIter - 1
      widthCounter <<- widthCounter + 1
      if(widthCounter == adaptWidthInterval) {
        for(i in 1:d) {
          if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
            if(nExpansions[i] == 0)   nExpansions[i] <<- 1
            widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
            widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
            adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
            nExpansions[i] <<- 0
            nContracts[i] <<- 0
            if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
              widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
          }
        }
        allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
        widthCounter     <<- 0
      }
      if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
        allWidthsAdapted <<- 1
    },
    reset = function() {
      gammaMatrix        <<- diag(d)
      empirCov           <<- diag(d)
      widthVec           <<- widthVecOriginal
      nExpansions        <<- rep(0, d)
      nContracts         <<- rep(0, d)
      adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
      factorCounter      <<- 0
      factorTimesAdapted <<- 0
      allWidthsAdapted   <<- 0
      widthCounter       <<- 0
      adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
      adaptWidthInterval <<- 1
      widthIndicatorVec  <<- rep(1, d)
    }
  )
)

####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################

mcmc_determineCalcAndCopyNodes <- function(model, target) {
  targetExpanded <- model$expandNodeNames(target)
  modelPredictiveNodes <- model$getNodeNames(predictiveOnly = TRUE)
  targetExpandedPPbool <- targetExpanded %in% modelPredictiveNodes
  targetAllPP <- all(targetExpandedPPbool)
  targetAnyPP <- any(targetExpandedPPbool)
  ## if a particular sampler is assigned *jointly to PP and non-PP* nodes, then we're going to bail
  ## out and quit, if the option MCMCusePredictiveDependenciesInCalculations == FALSE.
  ## this is an extreme corner-case, which I think will lead to problems.
  if(targetAnyPP && !targetAllPP && !getNimbleOption('MCMCusePredictiveDependenciesInCalculations'))
    stop('cannot assign samplers jointly to posterior predictive (PP) nodes and non-PP nodes, when MCMCusePredictiveDependenciesInCalculations option is FALSE', call. = FALSE)
  ## if the sampler calling this, itself, is operating exclusively on posterior predictive nodes,
  ## then regardless of how the rest of the model is being sampled (w.r.t. inclusion of posterior predictive nodes),
  ## we'll include 'self' and all stochastic dependencies (the full markov blanket) in the calculations,
  ## which necessarily are taking place entirely within a posterior predictive network of nodes.
  ## this should lead to correct behaviour (consistent samples and joint posteriors) in all cases.
  if(targetAllPP) {
    ## when sampler is operating only on posterior predictive nodes,
    ## then always include all predictive dependencies:
    calcNodes <- model$getDependencies(target, includePredictive = TRUE)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE, includePredictive = TRUE)
    ##calcNodesPPomitted <- character()
  } else {
    ## usual case:
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ##calcNodesPPomitted <- setdiff(model$getDependencies(target, includePredictive = TRUE), calcNodes)
  }
  ## copyNodes:
  copyNodes <- model$getDependencies(target, self = FALSE)
  isStochCopyNodes <- model$isStoch(copyNodes)
  copyNodesDeterm <- copyNodes[!isStochCopyNodes]
  copyNodesStoch <- copyNodes[isStochCopyNodes]
  ##
  ccList <- list(
    calcNodes = calcNodes,
    calcNodesNoSelf = calcNodesNoSelf,
    ##calcNodesPPomitted = calcNodesPPomitted,
    copyNodesDeterm = copyNodesDeterm,
    copyNodesStoch = copyNodesStoch
  )
  return(ccList)
}


#' @rdname samplers
#' @export
#######################################################################################
### RW_PF_block, does a block RW, but using a particle filter likelihood function #####
#######################################################################################

#' @rdname samplers
#' @export
sampler_AFSS_INLAlatent_block <- nimbleFunction(
  name = 'sampler_AFSS_INLAlatent_block',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    x <- extractControlElement(control, 'x',  double())
    y <- extractControlElement(control, 'y',  character())
    fixedVals <- extractControlElement(control, 'fixedVals',  double())
    fam <- extractControlElement(control, 'fam',  "gaussian")
    interVal <- extractControlElement(control, 'interInModel',  1)
    existingINLA          <- extractControlElement(control, 'fit.inla',                   NULL)
    ## control list extraction
    widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
    maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            100)
    adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
    adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
    adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
    adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
    maxContractions        <- extractControlElement(control, 'maxContractions',          1000)
    maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
    eps <- 1e-15

    #Extract y values
    #nodesY <- model$expandNodeNames(y)
    #y <- nimble::values(model, nodesY)
    #y <- c(model[[y]])
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    #calcNodesDepStage1 <- model$getVarNames(nodes = calcNodesDepStage)
    #calcNodesDepStage <- calcNodesDepStage[!model$isData(calcNodesDepStage)]
    #Note that model likelihood will be estimated from INLA
    # so we take the data out
    ## numeric value generation
    d                  <- length(targetAsScalar)
    discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
    anyDiscrete        <- any(discrete)
    gammaMatrix        <- diag(d)         # matrix of orthogonal bases
    if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
    widthVecOriginal   <- widthVec
    nExpansions        <- rep(0, d)       # number of expansions
    nContracts         <- rep(0, d)       # number of contractions
    adaptFactorMaxIterOriginal <- adaptFactorMaxIter
    factorCounter      <- 0               # number of iterations since last factor adaptation
    factorTimesAdapted <- 0               # number of times factors have adapted
    empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
    empirCov           <- diag(d)
    allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
    widthCounter       <- 0               # number of iterations since last width adaptation
    adaptWidthMaxIterOriginal <- adaptWidthMaxIter
    adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
    widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
    ## checks
    if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
    if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
      stop('sliceWidths must be a numeric vector')
    if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)

    #set up for INLA estimated vals
    if(length(model$expandNodeNames(fixedVals)) > 0){
      latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)


    my_particleFilter <- buildINLAmodel(model, fam, x,y = y,control = list(fit.inla = existingINLA,
                                                                           fixedVals = fixedVals))
    #targetVal <- values(model, targetAsScalar)
    particleMV <- my_particleFilter$mvEWSamples
    out <- -Inf
    print(latentSamp)
  },
  run = function() {
    #out <<- my_particleFilter$run(beta = values(model,  target))
    maxContractionsReached <- FALSE
    for(i in 1:d) {
      eigenVec <- gammaMatrix[, i]
      width <- widthVec[i]
      u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
      x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
      Lbound <- -1.0 * runif(1, 0, 1) * width
      Rbound <- Lbound + width
      L <- x0 + Lbound * eigenVec
      R <- x0 + Rbound * eigenVec
      maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
      maxStepsR <- maxSteps - 1 - maxStepsL
      lp <- setAndCalculateTarget(L) #+ out
      while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
        Lbound <- Lbound - width
        L <- x0 + Lbound * eigenVec
        lp <- setAndCalculateTarget(L) #+ out
        maxStepsL <- maxStepsL - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      lp <- setAndCalculateTarget(R) + out
      while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
        Rbound <- Rbound + width
        R <- x0 + Rbound * eigenVec
        lp <- setAndCalculateTarget(R) #+ out
        maxStepsR <- maxStepsR - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
      x1 <- x0 + prop * eigenVec
      lp <- setAndCalculateTarget(x1)
      numContractions <- 0
      while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
        ## The checks for Rbound - Lbound small and max number of contractions are
        ## for cases where model is in invalid state and lp calculations are NA/NaN or where
        ## interval contracts to zero
        if(prop < 0) { Lbound <- prop }
        else         { Rbound <- prop }
        nContracts[i] <<- nContracts[i] + 1
        prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
        x1 <- x0 + prop * eigenVec
        lp <- setAndCalculateTarget(x1) #+ out
        numContractions <- numContractions + 1
      }
      if(Rbound - Lbound <= eps | numContractions == maxContractions)
        maxContractionsReached <- TRUE
    }
    #print(out)
    if(maxContractionsReached) {
      if(maxContractionsWarning)
        cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = fixedValsDep, logProb = TRUE)
    } else {
      out <<- my_particleFilter$run(beta = values(model,  targetAsScalar))
      copy(particleMV, model, fixedVals, fixedVals, row = 1)
      #calculate(model, fixedValsDep)
      model$calculate()
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    if(allWidthsAdapted == 0)   adaptWidths()
    if(adaptFactorMaxIter > 0)  adaptFactors()

  },
  methods = list(
    setAndCalculateTarget = function(targetValues = double(1)) {
      if(anyDiscrete == 1)
        for(i in 1:d)
          if(discrete[i] == 1)   targetValues[i] <- floor(targetValues[i])
      values(model, target) <<- targetValues
      lp <- model$calculate(calcNodesProposalStage)
      if(lp == -Inf) return(lp)
      lp <- lp + model$calculate(calcNodesDepStage) #+ out
      returnType(double())
      return(lp)
    },
    adaptFactors = function() {
      adaptFactorMaxIter <<- adaptFactorMaxIter - 1
      factorCounter <<- factorCounter + 1
      empirSamp[factorCounter, 1:d] <<- values(model, target)
      if(factorCounter == adaptFactorInterval) {
        for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
        gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
        factorTimesAdapted <<- factorTimesAdapted + 1
        factorCounter      <<- 0
        nExpansions        <<- rep(0, d)
        nContracts         <<- rep(0, d)
        allWidthsAdapted   <<- 0
        widthCounter       <<- 0
        adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
        adaptWidthInterval <<- 1
        widthIndicatorVec  <<- rep(1, d)
      }
    },
    adaptWidths = function() {
      adaptWidthMaxIter <<- adaptWidthMaxIter - 1
      widthCounter <<- widthCounter + 1
      if(widthCounter == adaptWidthInterval) {
        for(i in 1:d) {
          if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
            if(nExpansions[i] == 0)   nExpansions[i] <<- 1
            widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
            widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
            adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
            nExpansions[i] <<- 0
            nContracts[i] <<- 0
            if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
              widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
          }
        }
        allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
        widthCounter     <<- 0
      }
      if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
        allWidthsAdapted <<- 1
    },
    reset = function() {
      gammaMatrix        <<- diag(d)
      empirCov           <<- diag(d)
      widthVec           <<- widthVecOriginal
      nExpansions        <<- rep(0, d)
      nContracts         <<- rep(0, d)
      adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
      factorCounter      <<- 0
      factorTimesAdapted <<- 0
      allWidthsAdapted   <<- 0
      widthCounter       <<- 0
      adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
      adaptWidthInterval <<- 1
      widthIndicatorVec  <<- rep(1, d)
    }
  )
)



# used to sample with AFSS using INLA embedded in BUGS code
#  for binary samples
# reference: Li, Y., Linero, A., & Walker, S. G. (2023). Latent uniform samplers on multivariate binary spaces. Statistics and Computing, 33(5), 102.
#' @rdname samplers
#' @export
sampler_AF_slice_binary <- nimbleFunction(
  name = 'sampler_AF_slice_binary',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
    maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            10)
    adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
    adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
    adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
    adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
    maxContractions        <- extractControlElement(control, 'maxContractions',          500)
    maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
    eps <- 1e-15
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    out <- 0#model$calculate(calcNodesProposalStage)
    ## numeric value generation
    d                  <- length(targetAsScalar)
    discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
    anyDiscrete        <- any(discrete)
    gammaMatrix        <- diag(d)         # matrix of orthogonal bases
    if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
    widthVecOriginal   <- widthVec
    nExpansions        <- rep(0, d)       # number of expansions
    nContracts         <- rep(0, d)       # number of contractions
    adaptFactorMaxIterOriginal <- adaptFactorMaxIter
    factorCounter      <- 0               # number of iterations since last factor adaptation
    factorTimesAdapted <- 0               # number of times factors have adapted
    empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
    empirCov           <- diag(d)
    allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
    widthCounter       <- 0               # number of iterations since last width adaptation
    adaptWidthMaxIterOriginal <- adaptWidthMaxIter

    adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
    widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
    ## checks
    if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
    if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
      stop('sliceWidths must be a numeric vector')
    if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)
  },
  run = function() {
    maxContractionsReached <- FALSE
    for(i in 1:d) {
      eigenVec <- gammaMatrix[, i]
      width <- widthVec[i]
      u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
      x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
      Lbound <- -1.0 * runif(1, 0, 1) * width
      Rbound <- Lbound + width
      L <- x0 + Lbound * eigenVec
      R <- x0 + Rbound * eigenVec
      maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
      maxStepsR <- maxSteps - 1 - maxStepsL
      lp <- setAndCalculateTarget(L, out)
      while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
        Lbound <- Lbound - width
        L <- x0 + Lbound * eigenVec
        lp <- setAndCalculateTarget(L, out)
        maxStepsL <- maxStepsL - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      lp <- setAndCalculateTarget(R, out)
      while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
        Rbound <- Rbound + width
        R <- x0 + Rbound * eigenVec
        lp <- setAndCalculateTarget(R, out)
        maxStepsR <- maxStepsR - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
      x1 <- x0 + prop * eigenVec
      lp <- setAndCalculateTarget(x1, out)
      numContractions <- 0
      while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
        ## The checks for Rbound - Lbound small and max number of contractions are
        ## for cases where model is in invalid state and lp calculations are NA/NaN or where
        ## interval contracts to zero
        if(prop < 0) { Lbound <- prop }
        else         { Rbound <- prop }
        nContracts[i] <<- nContracts[i] + 1
        prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
        x1 <- x0 + prop * eigenVec
        lp <- setAndCalculateTarget(x1, out)
        numContractions <- numContractions + 1
      }
      if(Rbound - Lbound <= eps | numContractions == maxContractions)
        maxContractionsReached <- TRUE
    }
    if(maxContractionsReached) {
      if(maxContractionsWarning)
        cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    } else {
      out <<- 0#model$calculate(calcNodesProposalStage)
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    }
    if(allWidthsAdapted == 0)   adaptWidths()
    if(adaptFactorMaxIter > 0)  adaptFactors()
  },
  methods = list(
    setAndCalculateTarget = function(targetValues = double(1), out = double(0)) {
      if(anyDiscrete == 1)
        for(i in 1:d){
          if(discrete[i] == 1){
            if(targetValues[i] > 0){
              targetValues[i] = 1
            }else{
              targetValues[i] = 0
              }
            }
          }#floor(targetValues[i])
      values(model, target) <<- targetValues
      lp <- model$calculate(calcNodesProposalStage)
      if(lp == -Inf) return(lp)
      lp <- lp + model$calculate(calcNodesDepStage)
      returnType(double())
      return(lp)
    },
    adaptFactors = function() {
      adaptFactorMaxIter <<- adaptFactorMaxIter - 1
      factorCounter <<- factorCounter + 1
      empirSamp[factorCounter, 1:d] <<- values(model, target)
      if(factorCounter == adaptFactorInterval) {
        for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
        gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
        factorTimesAdapted <<- factorTimesAdapted + 1
        factorCounter      <<- 0
        nExpansions        <<- rep(0, d)
        nContracts         <<- rep(0, d)
        allWidthsAdapted   <<- 0
        widthCounter       <<- 0
        adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
        adaptWidthInterval <<- 1
        widthIndicatorVec  <<- rep(1, d)
      }
    },
    adaptWidths = function() {
      adaptWidthMaxIter <<- adaptWidthMaxIter - 1
      widthCounter <<- widthCounter + 1
      if(widthCounter == adaptWidthInterval) {
        for(i in 1:d) {
          if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
            if(nExpansions[i] == 0)   nExpansions[i] <<- 1
            widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
            widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
            adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
            nExpansions[i] <<- 0
            nContracts[i] <<- 0
            if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
              widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
          }
        }
        allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
        widthCounter     <<- 0
      }
      if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
        allWidthsAdapted <<- 1
    },
    reset = function() {
      gammaMatrix        <<- diag(d)
      empirCov           <<- diag(d)
      widthVec           <<- widthVecOriginal
      nExpansions        <<- rep(0, d)
      nContracts         <<- rep(0, d)
      adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
      factorCounter      <<- 0
      factorTimesAdapted <<- 0
      allWidthsAdapted   <<- 0
      widthCounter       <<- 0
      adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
      adaptWidthInterval <<- 1
      widthIndicatorVec  <<- rep(1, d)
    }
  )
)


# used to sample with AFSS using INLA embedded code that returns marginal likehood to write the sampler
# and for binary samples
# reference: Li, Y., Linero, A., & Walker, S. G. (2023). Latent uniform samplers on multivariate binary spaces. Statistics and Computing, 33(5), 102.
#' @rdname samplers
#' @export
sampler_AFSS_INLA_block_binary <- nimbleFunction(
  name = 'sampler_AFSS_INLA_block_binary',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    x <- extractControlElement(control, 'x',  double())
    y <- extractControlElement(control, 'y',  character())
    fixedVals <- extractControlElement(control, 'fixedVals',  double())
    fam <- extractControlElement(control, 'fam',  "gaussian")
    interVal <- extractControlElement(control, 'interInModel',  1)
    existingINLA          <- extractControlElement(control, 'fit.inla',                   NULL)
    ## control list extraction
    widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
    maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            10)
    adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
    adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
    adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
    adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
    maxContractions        <- extractControlElement(control, 'maxContractions',          500)
    maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
    eps <- 1e-15

    #Extract y values
    #nodesY <- model$expandNodeNames(y)
    #y <- nimble::values(model, nodesY)
    #y <- c(model[[y]])
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    #calcNodesDepStage1 <- model$getVarNames(nodes = calcNodesDepStage)
    #calcNodesDepStage <- calcNodesDepStage[!model$isData(calcNodesDepStage)]
    #Note that model likelihood will be estimated from INLA
    # so we take the data out
    ## numeric value generation
    d                  <- length(targetAsScalar)
    discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
    anyDiscrete        <- any(discrete)
    gammaMatrix        <- diag(d)         # matrix of orthogonal bases
    if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
    widthVecOriginal   <- widthVec
    nExpansions        <- rep(0, d)       # number of expansions
    nContracts         <- rep(0, d)       # number of contractions
    adaptFactorMaxIterOriginal <- adaptFactorMaxIter
    factorCounter      <- 0               # number of iterations since last factor adaptation
    factorTimesAdapted <- 0               # number of times factors have adapted
    empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
    empirCov           <- diag(d)
    allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
    widthCounter       <- 0               # number of iterations since last width adaptation
    adaptWidthMaxIterOriginal <- adaptWidthMaxIter
    adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
    widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
    ## checks
    if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
    if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
      stop('sliceWidths must be a numeric vector')
    if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)

    #set up for INLA estimated vals
    if(length(model$expandNodeNames(fixedVals)) > 0){
      latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)


    my_particleFilter <- buildINLAmodel(model, fam, x,y = y,control = list(fit.inla = existingINLA,
                                                                           fixedVals = fixedVals))
    #targetVal <- values(model, targetAsScalar)
    particleMV <- my_particleFilter$mvEWSamples
    out <- -Inf
    print(latentSamp)
  },
  run = function() {
    #out <<- my_particleFilter$run(beta = values(model,  target))
    maxContractionsReached <- FALSE
    for(i in 1:d) {
      eigenVec <- gammaMatrix[, i]
      width <- widthVec[i]
      u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
      x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
      Lbound <- -1.0 * runif(1, 0, 1) * width
      Rbound <- Lbound + width
      L <- x0 + Lbound * eigenVec
      R <- x0 + Rbound * eigenVec
      maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
      maxStepsR <- maxSteps - 1 - maxStepsL
      lp <- setAndCalculateTarget(L) #+ out
      while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
        Lbound <- Lbound - width
        L <- x0 + Lbound * eigenVec
        lp <- setAndCalculateTarget(L) #+ out
        maxStepsL <- maxStepsL - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      lp <- setAndCalculateTarget(R) + out
      while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
        Rbound <- Rbound + width
        R <- x0 + Rbound * eigenVec
        lp <- setAndCalculateTarget(R) #+ out
        maxStepsR <- maxStepsR - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
      x1 <- x0 + prop * eigenVec
      lp <- setAndCalculateTarget(x1)
      numContractions <- 0
      while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
        ## The checks for Rbound - Lbound small and max number of contractions are
        ## for cases where model is in invalid state and lp calculations are NA/NaN or where
        ## interval contracts to zero
        if(prop < 0) { Lbound <- prop }
        else         { Rbound <- prop }
        nContracts[i] <<- nContracts[i] + 1
        prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
        x1 <- x0 + prop * eigenVec
        lp <- setAndCalculateTarget(x1) #+ out
        numContractions <- numContractions + 1
      }
      if(Rbound - Lbound <= eps | numContractions == maxContractions)
        maxContractionsReached <- TRUE
    }
    #print(out)
    if(maxContractionsReached) {
      if(maxContractionsWarning)
        cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = fixedValsDep, logProb = TRUE)
    } else {
      out <<- my_particleFilter$run(beta = values(model,  targetAsScalar))
      copy(particleMV, model, fixedVals, fixedVals, row = 1)
      #calculate(model, fixedValsDep)
      model$calculate()
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
      copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    if(allWidthsAdapted == 0)   adaptWidths()
    if(adaptFactorMaxIter > 0)  adaptFactors()

  },
  methods = list(
    setAndCalculateTarget = function(targetValues = double(1)) {
      if(anyDiscrete == 1)
        for(i in 1:d){
          if(discrete[i] == 1){
            if(targetValues[i] > 0){
              targetValues[i] = 1
            }else{
              targetValues[i] = 0
            }
          }
        }
      values(model, target) <<- targetValues
      lp <- model$calculate(calcNodesProposalStage)
      if(lp == -Inf) return(lp)
      lp <- lp + model$calculate(calcNodesDepStage) #+ out
      returnType(double())
      return(lp)
    },
    adaptFactors = function() {
      adaptFactorMaxIter <<- adaptFactorMaxIter - 1
      factorCounter <<- factorCounter + 1
      empirSamp[factorCounter, 1:d] <<- values(model, target)
      if(factorCounter == adaptFactorInterval) {
        for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
        gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
        factorTimesAdapted <<- factorTimesAdapted + 1
        factorCounter      <<- 0
        nExpansions        <<- rep(0, d)
        nContracts         <<- rep(0, d)
        allWidthsAdapted   <<- 0
        widthCounter       <<- 0
        adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
        adaptWidthInterval <<- 1
        widthIndicatorVec  <<- rep(1, d)
      }
    },
    adaptWidths = function() {
      adaptWidthMaxIter <<- adaptWidthMaxIter - 1
      widthCounter <<- widthCounter + 1
      if(widthCounter == adaptWidthInterval) {
        for(i in 1:d) {
          if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
            if(nExpansions[i] == 0)   nExpansions[i] <<- 1
            widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
            widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
            adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
            nExpansions[i] <<- 0
            nContracts[i] <<- 0
            if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
              widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
          }
        }
        allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
        widthCounter     <<- 0
      }
      if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
        allWidthsAdapted <<- 1
    },
    reset = function() {
      gammaMatrix        <<- diag(d)
      empirCov           <<- diag(d)
      widthVec           <<- widthVecOriginal
      nExpansions        <<- rep(0, d)
      nContracts         <<- rep(0, d)
      adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
      factorCounter      <<- 0
      factorTimesAdapted <<- 0
      allWidthsAdapted   <<- 0
      widthCounter       <<- 0
      adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
      adaptWidthInterval <<- 1
      widthIndicatorVec  <<- rep(1, d)
    }
  )
)

# Sample discrete and continuous variable using INLA model but with distibutions that take only positive values like Poisson and binomial
#' @rdname samplers
#' @export
sampler_AF_slice_positive <- nimbleFunction(
  name = 'sampler_AF_slice_positive',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    widthVec               <- extractControlElement(control, 'sliceWidths',              'oneVec')
    maxSteps               <- extractControlElement(control, 'sliceMaxSteps',            10)
    adaptFactorMaxIter     <- extractControlElement(control, 'sliceAdaptFactorMaxIter',  15000)
    adaptFactorInterval    <- extractControlElement(control, 'sliceAdaptFactorInterval', 200)
    adaptWidthMaxIter      <- extractControlElement(control, 'sliceAdaptWidthMaxIter',   512)
    adaptWidthTolerance    <- extractControlElement(control, 'sliceAdaptWidthTolerance', 0.1)
    maxContractions        <- extractControlElement(control, 'maxContractions',          500)
    maxContractionsWarning <- extractControlElement(control, 'maxContractionsWarning',   TRUE)
    eps <- 1e-15
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    ccList <- mcmc_determineCalcAndCopyNodes(model, target)
    calcNodes <- ccList$calcNodes; copyNodesDeterm <- ccList$copyNodesDeterm; copyNodesStoch <- ccList$copyNodesStoch   # not used: calcNodesNoSelf
    finalTargetIndex <- max(match(model$expandNodeNames(target), calcNodes))
    if(!is.integer(finalTargetIndex) | length(finalTargetIndex) != 1 | is.na(finalTargetIndex[1]))   stop('problem with target node in AF_slice sampler')
    calcNodesProposalStage <- calcNodes[1:finalTargetIndex]
    calcNodesDepStage <- calcNodes[-(1:finalTargetIndex)]
    out <- 0#model$calculate(calcNodesProposalStage)
    ## numeric value generation
    d                  <- length(targetAsScalar)
    discrete           <- sapply(targetAsScalar, function(x) model$isDiscrete(x))
    anyDiscrete        <- any(discrete)
    gammaMatrix        <- diag(d)         # matrix of orthogonal bases
    if(is.character(widthVec) && widthVec == 'oneVec')   widthVec <- rep(1,d)
    widthVecOriginal   <- widthVec
    nExpansions        <- rep(0, d)       # number of expansions
    nContracts         <- rep(0, d)       # number of contractions
    adaptFactorMaxIterOriginal <- adaptFactorMaxIter
    factorCounter      <- 0               # number of iterations since last factor adaptation
    factorTimesAdapted <- 0               # number of times factors have adapted
    empirSamp          <- matrix(0, nrow=adaptFactorInterval, ncol=d)   # matrix of posterior samples
    empirCov           <- diag(d)
    allWidthsAdapted   <- 0               # indicates whether all widths have finished adapting
    widthCounter       <- 0               # number of iterations since last width adaptation
    adaptWidthMaxIterOriginal <- adaptWidthMaxIter

    adaptWidthInterval <- 1               # interval to adapt widths; doubles each time widths are adaptated
    widthIndicatorVec  <- rep(1, d)       # indicator of which widths are still adapting
    ## checks
    if(d <= 1)                         stop('AF_slice sampler must be used on at least two target nodes')
    if(!inherits(widthVec, 'numeric') && !inherits(widthVec, 'integer'))
      stop('sliceWidths must be a numeric vector')
    if(length(widthVec) != d)          stop('sliceWidths must have length = ', d)
  },
  run = function() {
    maxContractionsReached <- FALSE
    for(i in 1:d) {
      eigenVec <- gammaMatrix[, i]
      width <- widthVec[i]
      u <- model$getLogProb(calcNodes) - rexp(1, 1)   # generate (log)-auxiliary variable: exp(u) ~ uniform(0, exp(lp))
      x0 <- values(model, target)                      # create random interval (L,R), of width 'width', around current value of target
      Lbound <- -1.0 * runif(1, 0, 1) * width
      Rbound <- Lbound + width
      L <- x0 + Lbound * eigenVec
      R <- x0 + Rbound * eigenVec
      maxStepsL <- floor(runif(1, 0, 1) * maxSteps)    # randomly allot (maxSteps-1) into maxStepsL and maxStepsR
      maxStepsR <- maxSteps - 1 - maxStepsL
      lp <- setAndCalculateTarget(L, out)
      while(maxStepsL > 0 & !is.nan(lp) & lp >= u) {   # step L left until outside of slice (max maxStepsL steps)
        Lbound <- Lbound - width
        L <- x0 + Lbound * eigenVec
        lp <- setAndCalculateTarget(L, out)
        maxStepsL <- maxStepsL - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      lp <- setAndCalculateTarget(R, out)
      while(maxStepsR > 0 & !is.nan(lp) & lp >= u) {   # step R right until outside of slice (max maxStepsR steps)
        Rbound <- Rbound + width
        R <- x0 + Rbound * eigenVec
        lp <- setAndCalculateTarget(R, out)
        maxStepsR <- maxStepsR - 1
        nExpansions[i] <<- nExpansions[i] + 1
      }
      prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
      x1 <- x0 + prop * eigenVec
      lp <- setAndCalculateTarget(x1, out)
      numContractions <- 0
      while((is.nan(lp) | lp < u) & Rbound - Lbound > eps & numContractions < maxContractions) {   # must be is.nan()
        ## The checks for Rbound - Lbound small and max number of contractions are
        ## for cases where model is in invalid state and lp calculations are NA/NaN or where
        ## interval contracts to zero
        if(prop < 0) { Lbound <- prop }
        else         { Rbound <- prop }
        nContracts[i] <<- nContracts[i] + 1
        prop <- Lbound + runif(1, 0, 1) * (Rbound - Lbound)
        x1 <- x0 + prop * eigenVec
        lp <- setAndCalculateTarget(x1, out)
        numContractions <- numContractions + 1
      }
      if(Rbound - Lbound <= eps | numContractions == maxContractions)
        maxContractionsReached <- TRUE
    }
    if(maxContractionsReached) {
      if(maxContractionsWarning)
        cat("Warning: AF slice sampler reached maximum number of contractions in at least one dimension.\n")
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    } else {
      out <<- 0#model$calculate(calcNodesProposalStage)
      ##model$calculate(calcNodesPPomitted)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesProposalStage, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodesStoch, logProbOnly = TRUE)
    }
    if(allWidthsAdapted == 0)   adaptWidths()
    if(adaptFactorMaxIter > 0)  adaptFactors()
  },
  methods = list(
    setAndCalculateTarget = function(targetValues = double(1), out = double(0)) {
      if(anyDiscrete == 1)
        for(i in 1:d){
          if(discrete[i] == 1){
            if(targetValues[i] < 0){
              targetValues[i] = 0
            }else{
              targetValues[i] = floor(targetValues[i])
            }
          }
        }#floor(targetValues[i])
      values(model, target) <<- targetValues
      lp <- model$calculate(calcNodesProposalStage)
      if(lp == -Inf) return(lp)
      lp <- lp + model$calculate(calcNodesDepStage)
      returnType(double())
      return(lp)
    },
    adaptFactors = function() {
      adaptFactorMaxIter <<- adaptFactorMaxIter - 1
      factorCounter <<- factorCounter + 1
      empirSamp[factorCounter, 1:d] <<- values(model, target)
      if(factorCounter == adaptFactorInterval) {
        for(i in 1:d)   empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
        empirCov <<- (t(empirSamp) %*% empirSamp) / (adaptFactorInterval - 1)
        gammaMatrix <<- eigen(empirCov)$vectors  # replace old factors with new factors
        factorTimesAdapted <<- factorTimesAdapted + 1
        factorCounter      <<- 0
        nExpansions        <<- rep(0, d)
        nContracts         <<- rep(0, d)
        allWidthsAdapted   <<- 0
        widthCounter       <<- 0
        adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
        adaptWidthInterval <<- 1
        widthIndicatorVec  <<- rep(1, d)
      }
    },
    adaptWidths = function() {
      adaptWidthMaxIter <<- adaptWidthMaxIter - 1
      widthCounter <<- widthCounter + 1
      if(widthCounter == adaptWidthInterval) {
        for(i in 1:d) {
          if(widthIndicatorVec[i] == 1) {   # widths that are still adapting
            if(nExpansions[i] == 0)   nExpansions[i] <<- 1
            widthAdaptRatio <- nExpansions[i] / (nExpansions[i] + nContracts[i])
            widthVec[i] <<- widthVec[i] * 2 * widthAdaptRatio
            adaptWidthInterval <<- 2 * adaptWidthInterval   # double width adapt interval
            nExpansions[i] <<- 0
            nContracts[i] <<- 0
            if(adaptWidthInterval > 16)  # once adapt interval is large enough, determine whether adaptation is finished
              widthIndicatorVec[i] <<- (abs(widthAdaptRatio - .5) > adaptWidthTolerance)  # equals 1 if adaptation isn't finished
          }
        }
        allWidthsAdapted <<- 1 - ceiling(mean(widthIndicatorVec))  # equals 1 only if all slice adapt indicators are 0
        widthCounter     <<- 0
      }
      if(adaptWidthMaxIter <= 0)  # alternatively, if max iters have been reached, stop adapting
        allWidthsAdapted <<- 1
    },
    reset = function() {
      gammaMatrix        <<- diag(d)
      empirCov           <<- diag(d)
      widthVec           <<- widthVecOriginal
      nExpansions        <<- rep(0, d)
      nContracts         <<- rep(0, d)
      adaptFactorMaxIter <<- adaptFactorMaxIterOriginal
      factorCounter      <<- 0
      factorTimesAdapted <<- 0
      allWidthsAdapted   <<- 0
      widthCounter       <<- 0
      adaptWidthMaxIter  <<- adaptWidthMaxIterOriginal
      adaptWidthInterval <<- 1
      widthIndicatorVec  <<- rep(1, d)
    }
  )
)
