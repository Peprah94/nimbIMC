sampler_BASE <- nimble::nimbleFunctionVirtual(
  name = 'sampler_BASE',
  methods = list(
    reset = function() { }
  )
)

myRW_dirichlet <- nimbleFunction(
  name = 'myRW_dirichlet',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scaleOriginal       <- extractControlElement(control, 'scale',               1)
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    # d <- length(targetAsScalar)

    d <- length(model$expandNodeNames(target))
    dr <- length(targetAsScalar)/d

    #targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    expandTarget <- model$expandNodeNames(target)
    target1 <- expandTarget[1]
    target2 <- expandTarget[2]
    thetaVec         <- rep(0, d)
    scaleVec         <- rep(scaleOriginal, d)
    thetaVec2         <- rep(0, d)
    scaleVec2       <- rep(scaleOriginal, d)
    timesRan         <- 0
    timesAcceptedVec <- rep(0, d)
    timesAdapted     <- 0
    optimalAR        <- 0.44
    gamma1           <- 0
    logMHR <- rep(0,d)
    propValue <- rep(0, d)
    logMHR2 <- rep(0,d)
    propValue2 <- rep(0, d)
    ## checks
    #if(length(model$expandNodeNames(target)) > 1)    stop('RW_dirichlet sampler only applies to one target node')
    #if(model$getDistribution(target) != 'ddirch')    stop('can only use RW_dirichlet sampler for dirichlet distributions')
  },
  run = function() {

    for(j in 1:dr){
      if(j == 1){
        if(thetaVec[1] == 0)   thetaVec <<- values(model, target1)   ## initialization
        alphaVec <- model$getParam(target1, 'alpha')
        for(i in 1:d) {
          currentValue <- thetaVec[i]
          propLogScale <- rnorm(1, mean = 0, sd = scaleVec[i])
          propValue[i] <<- currentValue * exp(propLogScale)
          if(propValue[i] != 0) {
            thetaVecProp <- thetaVec
            thetaVecProp[i] <- propValue[i]
            values(model, target1) <<- thetaVecProp / sum(thetaVecProp)
            logMHR[i] <<- alphaVec[i]*propLogScale + currentValue - propValue[i]}
        }
      }

      if(j == 2){
        if(thetaVec2[1] == 0)   thetaVec2 <<- values(model, target2)   ## initialization
        alphaVec <- model$getParam(target2, 'alpha')
        for(i in 1:d) {
          currentValue <- thetaVec2[i]
          propLogScale <- rnorm(1, mean = 0, sd = scaleVec2[i])
          propValue2[i] <<- currentValue * exp(propLogScale)
          if(propValue2[i] != 0) {
            thetaVecProp2 <- thetaVec2
            thetaVecProp2[i] <- propValue2[i]
            values(model, target2) <<- thetaVecProp2 / sum(thetaVecProp2)
            logMHR2[i] <<- alphaVec[i]*propLogScale + currentValue - propValue2[i]}
        }
      }
    }

    eldf <-  model$calculateDiff(calcNodesNoSelf)

    for(j in 1:dr){
      if(j == 1){
        for(i in 1:d) {
          newlogMHR <- logMHR[i] + eldf
          if(propValue[1] != 0){
            jump <- decide(newlogMHR)
          }else{
            jump <- FALSE
          }
          #if(adaptive & jump)   timesAcceptedVec[i] <<- timesAcceptedVec[i] + 1
          if(jump) {
            thetaVec <<- thetaVecProp
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target1, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
          } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target1, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
          }
          model$calculate(target1)                                                             ## update target logProb
          nimCopy(from = model, to = mvSaved, row = 1, nodes = target1, logProbOnly = TRUE)    ##
        }
      }

      if(j == 2){
        for(i in 1:d) {
          newlogMHR <- logMHR2[i] + eldf
          if(propValue2[i] != 0){
            jump <- decide(newlogMHR)
          }else{
            jump <- FALSE
          }
          #if(adaptive & jump)   timesAcceptedVec[i] <<- timesAcceptedVec[i] + 1
          if(jump) {
            thetaVec2 <<- thetaVecProp2
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target2, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
          } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target2, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
          }
          model$calculate(target2)                                                             ## update target logProb
          nimCopy(from = model, to = mvSaved, row = 1, nodes = target2, logProbOnly = TRUE)    ##
        }
      }
    }
    if(adaptive) {
      timesRan <<- timesRan + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRateVec <- timesAcceptedVec / timesRan
        timesAdapted <<- timesAdapted + 1
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        adaptFactorVec <- exp(10 * gamma1 * (acceptanceRateVec - optimalAR))
        scaleVec <<- scaleVec * adaptFactorVec
        timesRan <<- 0
        timesAcceptedVec <<- numeric(d, 0)
      }
    }

  },
  methods = list(
    reset = function() {
      thetaVec         <<- numeric(d, 0)
      scaleVec         <<- numeric(d, scaleOriginal)
      thetaVec2         <<- numeric(d, 0)
      scaleVec2         <<- numeric(d, scaleOriginal)
      timesRan         <<- 0
      timesAcceptedVec <<- numeric(d, 0)
      timesAdapted     <<- 0
      gamma1           <<- 0
    }
  )
)








dmydirch <- nimble::nimbleFunction(
  run = function(x = double(1), alpha = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    logProb <- sum(lgamma(alpha)) - lgamma(sum(alpha)) +
      sum((alpha -1) * log(x))
    if(log) return(logProb)
    else return(exp(logProb))
  })

rmydirch <- nimble::nimbleFunction(
  run = function(n = integer(0), alpha = double(1)) {
    returnType(double(1))
    if(n != 1) print("rdirch only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(p)
  })
