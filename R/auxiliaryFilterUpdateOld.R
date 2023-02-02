setSavedWts <-  function(mv, currIndx, indx, mvWSamplesXSaved){
      mv['x',indx][currIndx] <- mvWSamplesXSaved[indx, currIndx]
      return(mv)
}






auxFStepUpdateOld <- nimbleFunction(
  name = 'auxFStepUpdate',
  contains = auxStepVirtual,
  setup = function(model, mvEWSamples, mvWSamples, nodes, iNode, names,
                   saveAll, smoothing, lookahead, resamplingMethod,
                   silent = TRUE, iNodePrev, mvWSamplesWTSaved,
                   mvWSamplesXSaved, mvEWSamplesXSaved
                   ) {
    notFirst <- iNode != 1
    last <- iNode == length(nodes)
    prevNode <- nodes[if(notFirst) iNode-1 else iNode]

    notFirst <- iNode != 1
    modelSteps <- particleFilter_splitModelSteps(model, nodes, iNode, notFirst)
    prevDeterm <- modelSteps$prevDeterm
    calc_thisNode_self <- modelSteps$calc_thisNode_self
    calc_thisNode_deps <- modelSteps$calc_thisNode_deps

    ## prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
    thisNode <- nodes[iNode]
    ## thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
    ## thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
    ## t is the current time point.
    t <- iNode
    ## Get names of x and xs node for current and previous time point,
    ## will be different depending on whether we are saving all time points
    ## or only the most recent.
    if(saveAll == 1){
      allPrevNodes <- model$expandNodeNames(nodes[1:(iNode-1)])
      prevXName <- prevNode
      thisXName <- thisNode
      currInd <- t
      prevInd <- t-1
      if(smoothing == TRUE){
        currInd <- 1
        prevInd <- 1
      }
    }
    else{
      allPrevNodes <- names # dummy value -- not used.
      prevXName <- names
      thisXName <- names
      currInd <- 1
      prevInd <- 1
    }

    auxFuncList <- nimbleFunctionList(auxFuncVirtual)
    allLatentNodes <- model$expandNodeNames(calc_thisNode_self, sort = TRUE) ## They should already be sorted, but sorting here is a failsafe.
    numLatentNodes <- length(allLatentNodes)

    if(lookahead == "mean"){
      for(i in 1:numLatentNodes)
        auxFuncList[[i]] <- auxLookFunc(model, allLatentNodes[i])
    }
    else{
      for(i in 1:numLatentNodes)
        auxFuncList[[i]] <- auxSimFunc(model,  allLatentNodes)
    }
    ess <- 0
    resamplerFunctionList <- nimbleFunctionList(resamplerVirtual)
    defaultResamplerFlag <- FALSE
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
  run = function(m = integer()) {
    returnType(double())
    auxll <- numeric(m, init=FALSE)
    auxWts <- numeric(m, init=FALSE)
    wts <- numeric(m, init=FALSE)
    ids <- integer(m, 0)
    ll <- numeric(m, init=FALSE)

    ## This is the look-ahead step, not conducted for first time-point.
    # if( t < iNodePrev +1){
    #
    #
    #   ## Get p(y_t+1 | x_t+1).
    # ## Use log-sum-exp trick to avoid underflow.
    # #maxWt <- max(wts)
    # #normWts <- exp(wts - maxWt)/sum(exp(wts - maxWt))
    # #ess <<- 1/sum(normWts^2)
    # for(i in 1:m){
    #   ## Save weights for use in next timepoint's look-ahead step.
    #   nimCopy(mvWSamples, model, prevXName, prevNode, row=i)
    #   if(lookahead == "mean"){
    #     for(j in 1:numLatentNodes)
    #       auxFuncList[[j]]$lookahead()
    #   } else auxFuncList[[1]]$lookahead()
    #   }
    #   for(i in 1: m){
    #   mvWSamples['wts', i][currInd] <<- mvWSamplesWTSaved[i, currInd]
    #   mvWSamples['x', i][currInd] <<- mvWSamplesXSaved[i, currInd]
    #   mvEWSamples['x', i][currInd] <<- mvEWSamplesXSaved[i, currInd]
    #   wts[i] <-  mvWSamplesWTSaved[i, currInd]
    #
    #   }
    #   if(lookahead == "mean"){
    #     for(j in 1:numLatentNodes)
    #       auxFuncList[[j]]$lookahead()
    #   } else auxFuncList[[1]]$lookahead()
    #   copy(model, mvEWSamples, nodes = thisNode, nodesTo = thisXName, row=i)
    #   copy(mvEWSamples, mvWSamples, thisXName, thisXName, row = i,  rowTo = i)
    # }
    #}
    ##########
    if(notFirst){
      for(i in 1:m) {
        if(smoothing == 1){
          ## smoothing is only allowed if saveAll is TRUE, so this should be ok.
          ## i.e., mvEWSamples have been resampled.
          copy(mvEWSamples, mvWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = i, rowTo=i)
        }
        nimCopy(mvWSamples, model, prevXName, prevNode, row=i)
        model$calculate(prevDeterm)
        ## The lookahead steps may include determ and stoch steps.
        if(lookahead == "mean"){
          for(j in 1:numLatentNodes)
            auxFuncList[[j]]$lookahead()
        } else auxFuncList[[1]]$lookahead()

        ## Get p(y_t+1 | x_t+1).
        auxll[i] <- model$calculate(calc_thisNode_deps)
        if(is.nan(auxll[i])){
          return(-Inf)
        }
        ## Multiply (on log scale) by weight from time t.
        if(t > iNodePrev){
          auxWts[i] <- auxll[i] + mvWSamples['wts',i][prevInd]
        }else{
          auxWts[i] <- mvWSamplesWTSaved[i, currInd]
        }
      }
      ## Normalize weights and resample, using log-sum-exp trick to avoid underflow.
      maxWt <- max(auxWts)
      normAuxWts <- exp(auxWts - maxWt)/sum(exp(auxWts - maxWt))
      if(defaultResamplerFlag == TRUE){
        rankSample(normAuxWts, m, ids, silent)
      }
      else{
        ids <- resamplerFunctionList[[1]]$run(normAuxWts)
      }
    }

    for(i in 1:m) {
      if(notFirst) {
        nimCopy(mvWSamples, model, nodes = prevXName, nodesTo = prevNode,
                row = ids[i])
        model$calculate(prevDeterm)
      }
      # Simulate from x_t+1 | x_t.
      #if(t > iNodePrev){
      model$simulate(calc_thisNode_self)
        nimCopy(model, mvEWSamples, nodes = thisNode, nodesTo = thisXName, row=i)
     # }else{

       # setSavedWts(mvWSamples, currInd, i, mvWSamplesXSaved)
       # nimCopy(mvWSamples, mvEWSamples, nodes = thisNode, nodesTo = thisXName, row=i)
      #}

      ## Get p(y_t+1 | x_t+1).
      ll[i] <- model$calculate(calc_thisNode_deps)
      if(is.nan(ll[i])){
        return(-Inf)
      }

      if(t > iNodePrev){
      if(notFirst){
        ## Construct weight following step 4 of paper.
        wts[i] <- ll[i]-auxll[ids[i]]
      }
      else{
        ## First step has no auxiliary weights.
        wts[i] <- ll[i]
      }
      }else{
        wts[i] <- mvWSamplesWTSaved[i, currInd]
      }
      }
    #}
    ## Use log-sum-exp trick to avoid underflow.
    maxWt <- max(wts)
    normWts <- exp(wts - maxWt)/sum(exp(wts - maxWt))
    ess <<- 1/sum(normWts^2)
    for(i in 1:m){
      if(t > iNodePrev){
      ## Save weights for use in next timepoint's look-ahead step.
      mvWSamples['wts', i][currInd] <<- log(normWts[i])
      }else{
        mvWSamples['wts', i][currInd] <<- mvWSamplesWTSaved[i, currInd]
      }
    }
    if(defaultResamplerFlag == TRUE){
      rankSample(normWts, m, ids, silent)
    }
    else{
      ids <- resamplerFunctionList[[1]]$run(normWts)
    }
    for(i in 1:m){
      nimCopy(mvEWSamples, mvWSamples, thisXName, thisXName, row = i,  rowTo = i)
    }

    if(saveAll | last) {
      for(i in 1:m) {
        if(smoothing == 1){
          copy(mvWSamples, mvEWSamples, nodes = allPrevNodes,
               nodesTo = allPrevNodes, row = ids[i], rowTo=i)
        }

        copy(mvWSamples, mvEWSamples, thisXName, thisXName, ids[i], i)
      }
    }
#}
    ## Calculate likelihood p(y_t+1 | y_1:t) as in equation (3) of paper.
    ## Use log-sum-exp trick to avoid underflow.
    if(notFirst){
      maxWt <- max(wts)
      maxAuxWt <- max(auxWts)
      outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m) + log(sum(exp(auxWts - maxAuxWt))) + maxAuxWt
    } else {
      maxWt <- max(wts)
      outLL <- log(sum(exp(wts - maxWt))) + maxWt - log(m)
    }
    return(outLL)

  },
  methods = list(
    returnESS = function() {
      returnType(double(0))
      return(ess)
    }
  )
)


buildAuxiliaryFilterUpdateOld <- nimbleFunction(
  name = 'buildAuxiliaryFilterUpdate',
  setup = function(model, nodes, mvWSamplesWTSaved,
                   mvWSamplesXSaved, mvEWSamplesXSaved,control = list()) {

    ## Control list extraction.
    saveAll <- control[['saveAll']]
    smoothing <- control[['smoothing']]
    silent <- control[['silent']]
    timeIndex <- control[['timeIndex']]
    lookahead <- control[['lookahead']]
    iNodePrev <- control[['iNodePrev']]
    initModel <- control[['initModel']]
   # mvWSamplesWTSaved <- control[['mvWSamplesWTSaved1']]
    #mvWSamplesXSaved <- control[['mvWSamplesXSaved1']]
   # mvEWSamplesXSaved <- control[['mvEWSamplesXSaved1']]
    M <- control[['M']]
    resamplingMethod <- control[['resamplingMethod']]
    if(is.null(silent)) silent <- TRUE
    if(is.null(saveAll)) saveAll <- FALSE
    if(is.null(smoothing)) smoothing <- FALSE
    if(is.null(lookahead)) lookahead <- 'simulate'
    if(is.null(initModel)) initModel <- TRUE
    if(!saveAll & smoothing) stop("must have saveAll = TRUE for smoothing to
                                  work")
    if(lookahead == "mean"){
      errors <- sapply(model$expandNodeNames(nodes), function(node){
        tryCatch(getParam(model, node, 'mean'), error=function(a){
          return("error")})})
      if(any(errors == "error", na.rm=TRUE))
        stop("cannot use 'mean' lookahead for this model, try 'simulate'")
    }
    else if(lookahead != "simulate"){
      stop("lookahead argument must be either 'simulate' or 'mean'")
    }
    if(is.null(resamplingMethod)) resamplingMethod <- 'default'
    if(!(resamplingMethod %in% c('default', 'multinomial', 'systematic', 'stratified',
                                 'residual')))
      stop('resamplingMethod must be one of: "default", "multinomial", "systematic",
           "stratified", or "residual". ')

    ## Latent state info.
    nodes <- findLatentNodes(model, nodes, timeIndex)

    dims <- lapply(nodes, function(n) nimDim(model[[n]]))
    if(length(unique(dims)) > 1)
      stop('sizes or dimensions of latent states varies')
    vars <- model$getVarNames(nodes =  nodes)

    my_initializeModel <- initializeModel(model, silent = silent)


    # Create mv variables for x state and sampled x states.  If saveAll=TRUE,
    # the sampled x states will be recorded at each time point.
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]


     if(saveAll){
    #
       names <- sapply(modelSymbolObjects, function(x)return(x$name))
       type <- sapply(modelSymbolObjects, function(x)return(x$type))
       size <- lapply(modelSymbolObjects, function(x)return(x$size))
       mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                  sizes = size))

       names <- c(names, "wts")
       type <- c(type, "double")
       size$wts <- length(dims)
      ## Only need one weight per particle (at time T) if smoothing == TRUE.
       if(smoothing){
         size$wts <- 1
       }
       mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                  types = type,
                                                  sizes = size))

     }
    else{
       names <- sapply(modelSymbolObjects, function(x)return(x$name))
       type <- sapply(modelSymbolObjects, function(x)return(x$type))
       size <- lapply(modelSymbolObjects, function(x)return(x$size))
      size[[1]] <- as.numeric(dims[[1]])

       mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))

      names <- c(names, "wts")
       type <- c(type, "double")
       size$wts <- 1
       mvWSamples  <- modelValues(modelValuesConf(vars = names,
                                                 types = type,
                                                 sizes = size))
     }

    #for(currInd in 1:iNodePrev){
    # for(i in 1:1000){
    #  mvWSamples['wts', i][currInd] <-  cModel$mvWSamples['wts', i][currInd]
   #   mvWSamples['x',i ][currInd] <-  cModel$mvWSamples['x', i][currInd]
    #  mvEWSamples['x',i ][currInd] <-  cModel$mvEWSamples['x', i][currInd]
     #}
    #}

   # mvWSamples <-  modelValues(cModel$mvWSamples)
   # mvEWSamples <- cModel$mvEWSamples

   #for(iNode in 1:iNodePrev){
      # mvWSamples[['x']] <-  mvWSamplesUpdate[['x']]
       #mvWSamples[['wts']] <-  mvWSamplesUpdate[['wts']]
       #mvEWSamples[['x']] <- mvEWSamplesUpdate[['x']]
       #mvEWSamples[['wts']] <- mvWSamplesUpdate[['wts']]
#}
    #for(iNode in 1:iNodePrev){

      #   prevXName <-  nodes[iNode]
      #   thisXName <-  nodes[iNode]
      #
      #
      #   #allPrevNodes <- names # dummy value -- not used.
      #   #prevXName <- names
      #   #thisXName <- names
      #
      #   for(i in 1:M){
      # nimCopy(mvWSamplesUpdate, mvWSamples, nodes = prevXName, nodesTo = thisXName,
      #      row = i)
      #
      # nimCopy(mvEWSamplesUpdate,mvEWSamples, nodes = thisXName, nodesTo = thisXName, row=i)
      #   }

    #}


    names <- names[1]
    auxStepFunctions1 <- nimbleFunctionList(auxStepVirtual)
    #for(iNode in 1:iNodePrev){
    #  auxStepFunctions1[[iNode]] <- cModel$auxStepFunctions[[iNode]]
    #}

    # for(iNode in seq_along(nodes))
    #   auxStepFunctions1[[iNode]] <- auxFStepUpdate(model, mvEWSamples, mvWSamples,
    #                                         nodes, iNode, names, saveAll,
    #                                         smoothing, lookahead,
    #                                         resamplingMethod, silent, iNodePrev,mvWSamplesWTSaved,
    #                                         mvWSamplesXSaved, mvEWSamplesXSaved)

    for(iNode in seq_along(nodes))
      auxStepFunctions1[[iNode]] <- auxFStepUpdate(model, mvEWSamples, mvWSamples,
                                            nodes, iNode, names, saveAll,
                                            smoothing, lookahead,
                                            resamplingMethod, silent,
                                            iNodePrev, mvWSamplesWTSaved,
                                            mvWSamplesXSaved, mvEWSamplesXSaved)

    essVals <- rep(0, length(nodes))
    lastLogLik <- -Inf
  },
  run = function(m = integer(default = 10000)) {
    returnType(double())
    declare(logL, double())
    if(initModel == TRUE) my_initializeModel$run()
    resize(mvEWSamples, m)
    resize(mvWSamples, m)
    logL <- 0
    for(iNode in seq_along(auxStepFunctions1)) {
      logL <- logL + auxStepFunctions1[[iNode]]$run(m)
    essVals[iNode] <<- auxStepFunctions1[[iNode]]$returnESS()
    #essVals <<- 0

      ## When all particles have 0 weight, likelihood becomes NAN
      ## this happens if top-level params have bad values - possible
      ## during pmcmc for example.
      if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
      if(logL == -Inf) {lastLogLik <<- logL; return(logL)}
      if(logL == Inf) {lastLogLik <<- -Inf; return(-Inf)}
    }
    lastLogLik <<- logL
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
    returnESS = function(){
      returnType(double(1))
      return(essVals)
    }
  )
)



# ##################
#
# buildAuxiliaryFilterUpdate1 <- nimbleFunction(
#   name = 'buildAuxiliaryFilterUpdate1',
#   setup = function(model, nodes,mvWSamplesWTSaved,
#                    mvWSamplesXSaved, mvEWSamplesXSaved, control = list()) {
#
#     ## Control list extraction.
#     saveAll <- control[['saveAll']]
#     smoothing <- control[['smoothing']]
#     silent <- control[['silent']]
#     timeIndex <- control[['timeIndex']]
#     lookahead <- control[['lookahead']]
#     initModel <- control[['initModel']]
#     iNodePrev <- control[['iNodePrev']]
#     resamplingMethod <- control[['resamplingMethod']]
#     if(is.null(silent)) silent <- TRUE
#     if(is.null(saveAll)) saveAll <- FALSE
#     if(is.null(smoothing)) smoothing <- FALSE
#     if(is.null(lookahead)) lookahead <- 'simulate'
#     if(is.null(initModel)) initModel <- TRUE
#     if(!saveAll & smoothing) stop("must have saveAll = TRUE for smoothing to
#                                   work")
#     if(lookahead == "mean"){
#       errors <- sapply(model$expandNodeNames(nodes), function(node){
#         tryCatch(getParam(model, node, 'mean'), error=function(a){
#           return("error")})})
#       if(any(errors == "error", na.rm=TRUE))
#         stop("cannot use 'mean' lookahead for this model, try 'simulate'")
#     }
#     else if(lookahead != "simulate"){
#       stop("lookahead argument must be either 'simulate' or 'mean'")
#     }
#     if(is.null(resamplingMethod)) resamplingMethod <- 'default'
#     if(!(resamplingMethod %in% c('default', 'multinomial', 'systematic', 'stratified',
#                                  'residual')))
#       stop('resamplingMethod must be one of: "default", "multinomial", "systematic",
#            "stratified", or "residual". ')
#
#     ## Latent state info.
#     nodes <- findLatentNodes(model, nodes, timeIndex)
#
#     dims <- lapply(nodes, function(n) nimDim(model[[n]]))
#     if(length(unique(dims)) > 1)
#       stop('sizes or dimensions of latent states varies')
#     vars <- model$getVarNames(nodes =  nodes)
#
#     my_initializeModel <- initializeModel(model, silent = silent)
#
#
#     # Create mv variables for x state and sampled x states.  If saveAll=TRUE,
#     # the sampled x states will be recorded at each time point.
#     modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]
#     if(saveAll){
#
#       names <- sapply(modelSymbolObjects, function(x)return(x$name))
#       type <- sapply(modelSymbolObjects, function(x)return(x$type))
#       size <- lapply(modelSymbolObjects, function(x)return(x$size))
#       mvEWSamples <- modelValues(modelValuesConf(vars = names,
#                                                  types = type,
#                                                  sizes = size))
#
#       names <- c(names, "wts")
#       type <- c(type, "double")
#       size$wts <- length(dims)
#       ## Only need one weight per particle (at time T) if smoothing == TRUE.
#       if(smoothing){
#         size$wts <- 1
#       }
#       mvWSamples  <- modelValues(modelValuesConf(vars = names,
#                                                  types = type,
#                                                  sizes = size))
#
#     }
#     else{
#       names <- sapply(modelSymbolObjects, function(x)return(x$name))
#       type <- sapply(modelSymbolObjects, function(x)return(x$type))
#       size <- lapply(modelSymbolObjects, function(x)return(x$size))
#       size[[1]] <- as.numeric(dims[[1]])
#
#       mvEWSamples <- modelValues(modelValuesConf(vars = names,
#                                                  types = type,
#                                                  sizes = size))
#
#       names <- c(names, "wts")
#       type <- c(type, "double")
#       size$wts <- 1
#       mvWSamples  <- modelValues(modelValuesConf(vars = names,
#                                                  types = type,
#                                                  sizes = size))
#     }
#
#     names <- names[1]
#     auxStepFunctions <- nimbleFunctionList(auxStepVirtual)
#     for(iNode in seq_along(nodes))
#       auxStepFunctions[[iNode]] <- auxFStep(model, mvEWSamples, mvWSamples,
#                                             nodes, iNode, names, saveAll,
#                                             smoothing, lookahead,
#                                             resamplingMethod, silent)
#
#     essVals <- rep(0, length(nodes))
#     lastLogLik <- -Inf
#   },
#   run = function(m = integer(default = 10000)) {
#     returnType(double())
#     declare(logL, double())
#     if(initModel == TRUE) my_initializeModel$run()
#     resize(mvEWSamples, m)
#     resize(mvWSamples, m)
#     logL <- 0
#     for(iNode in seq_along(auxStepFunctions)) {
#       logL <- logL + auxStepFunctions[[iNode]]$run(m)
#       essVals[iNode] <<- auxStepFunctions[[iNode]]$returnESS()
#
#       ## When all particles have 0 weight, likelihood becomes NAN
#       ## this happens if top-level params have bad values - possible
#       ## during pmcmc for example.
#       if(is.nan(logL)) {lastLogLik <<- -Inf; return(-Inf)}
#       if(logL == -Inf) {lastLogLik <<- logL; return(logL)}
#       if(logL == Inf) {lastLogLik <<- -Inf; return(-Inf)}
#     }
#     lastLogLik <<- logL
#     return(logL)
#   },
#   methods = list(
#     getLastLogLik = function() {
#       return(lastLogLik)
#       returnType(double())
#     },
#     setLastLogLik = function(lll = double()) {
#       lastLogLik <<- lll
#     },
#     returnESS = function(){
#       returnType(double(1))
#       return(essVals)
#     }
#   )
# )
