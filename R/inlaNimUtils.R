
####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################



#' @rdname samplers
#' @export
#######################################################################################
### RW_PF_block, does a block RW, but using a INLA built function #####
#######################################################################################

#' @rdname samplers
#' @export
sampler_RW_INLA_block <- nimbleFunction(
  name = 'sampler_RW_INLA_block',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',             FALSE)
    adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',       FALSE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',        200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent',  0.8)
    x                   <- extractControlElement(control, 'x',  double())
    y                   <- extractControlElement(control, 'y',  character())
    targetMCMC               <- extractControlElement(control, 'targetMCMC',              NULL)
    mu                  <- extractControlElement(control, 'mu',  NULL)
    #obsVars <- extractControlElement(control, 'obsVar',  character())
    fixedVals           <- extractControlElement(control, 'fixedVals',  double())
    fam                 <- extractControlElement(control, 'fam',  "gaussian")
    interVal            <- extractControlElement(control, 'interInModel',  1)
    scale               <- extractControlElement(control, 'scale',                1)
    propCov             <- extractControlElement(control, 'propCov',              'identity')
    existingINLA        <- extractControlElement(control, 'fit.inla',                   NULL)
    m                   <- extractControlElement(control, 'pfNparticles',         1000)
    filterType          <- extractControlElement(control, 'pfType',               'bootstrap')
    filterControl       <- extractControlElement(control, 'pfControl',            list())
    optimizeM           <- extractControlElement(control, 'pfOptimizeNparticles', FALSE)
    #Extract y values
    #obsData <- model$expandNodeNames(y)
    #y <- c(model[["y"]])
    ## node list generation
    targetAsScalar      <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes           <- model$getDependencies(target)
    if(length(fixedVals) > 0){
    latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)
    MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)

    topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)
    target <- model$expandNodeNames(target)
    ## numeric value generation
    optimizeM     <- as.integer(optimizeM)
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    prevLL        <- 0
    nVarEsts      <- 0
    itCount       <- 0
    d <- length(targetAsScalar)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    chol_propCov_scale <- scale * chol_propCov
    empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)

    # if(is.null(mu)){
    #   muTarget <- nimble::values(model,target)
    # }else{
    #     muTarget <- mu
    #   }
    storeParticleLP <- -Inf
    storeLLVar  <- 0
    nVarReps <- 7    ## number of LL estimates to compute to get each LL variance estimate for m optimization
    mBurnIn  <- 15   ## number of LL variance estimates to compute before deciding optimal m
    if(optimizeM)   m <- 3000
    ## nested function and function list definitions
    my_setAndCalculate <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, target, calcNodes)
    my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)

    #yVals <- values(model, obsData)
    my_particleFilter <- buildINLAmodel(model,
                                        fam,
                                        x,
                                        y = y,
                                        control = list(fit.inla = existingINLA,
                                                       fixedVals = fixedVals))

                                                       #Target values for inla
if(is.null(targetMCMC)){
  targetVal <- nimble::values(model, targetAsScalar)
}else{
targetMCMCasScalar <- model$expandNodeNames(targetMCMC, returnScalarComponents = TRUE)

targetVal <- nimble::values(model, c(targetMCMCasScalar, targetAsScalar))
   }

    particleMV <- my_particleFilter$mvEWSamples

    print(latentSamp)

    ## checks
    if(!inherits(propCov, 'matrix'))                    stop('propCov must be a matrix\n')
    if(!inherits(propCov[1,1], 'numeric'))              stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))                         stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))                           stop('propCov matrix must be symmetric')
    if(length(targetAsScalar) < 2)                      stop('less than two top-level targets; cannot use RW_PF_block sampler, try RW_PF sampler')
   # if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
  },
  run = function() {
    storeParticleLP <<- my_particleFilter$getLastLogLik()
    modelLP0 <- storeParticleLP + getLogProb(model, target)
    propValueVector <- generateProposalVector()
    my_setAndCalculate$run(propValueVector)
    targetVal <<- values(model, targetAsScalar)
    particleLP <- my_particleFilter$run(beta = targetVal)#,interInModel = interVal)
    modelLP1 <- particleLP + getLogProb(model, target)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
    if(!jump) {
      my_particleFilter$setLastLogLik(storeParticleLP)
    }
     if(jump & latentSamp) {
       ## if we jump, randomly sample latent nodes from pf output and put
       ## into model so that they can be monitored
    #   index <- ceiling(runif(1, 0, m))
       copy(particleMV, model, fixedVals, fixedVals, row = 1)
       calculate(model, fixedValsDep)
       copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
     else if(!jump & latentSamp) {
      ## if we don't jump, replace model latent nodes with saved latent nodes
       copy(from = mvSaved, to = model, nodes = fixedValsDep, row = 1, logProb = TRUE)
     }
    ##if(jump & !resample)  storeParticleLP <<- particleLP
    if(jump & optimizeM) optimM()
    if(adaptive)     adaptiveProcedure(jump)
  },
  methods = list(
    optimM = function() {
      tempM <- 15000
      declare(LLEst, double(1, nVarReps))
      if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
        for(i in 1:nVarReps)
          LLEst[i] <- my_particleFilter$run(beta = targetVal)#, interInModel = interVal)
        ## next, store average of var estimates
        if(nVarEsts == 1)
          storeLLVar <<- var(LLEst)/mBurnIn
        else {
          LLVar <- storeLLVar
          LLVar <- LLVar + var(LLEst)/mBurnIn
          storeLLVar <<- LLVar
        }
        nVarEsts <<- nVarEsts + 1
      }
      else {  # once enough var estimates have been taken, use their average to compute m
        m <<- m*storeLLVar/(0.92^2)
        m <<- ceiling(m)
        storeParticleLP <<- my_particleFilter$run(targetVal)
        optimizeM <<- 0
      }
    },
    generateProposalVector = function() {
      propValueVector <- rmnorm_chol(1, values(model, targetAsScalar), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$getGamma1()
          for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
          empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
        }
        chol_propCov_scale <<- chol_propCov * scale
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      chol_propCov_scale <<- chol_propCov * scale
      storeParticleLP <<- -Inf
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      my_calcAdaptationFactor$reset()
    }
  )
)



####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################



#' @rdname samplers
#' @export
#######################################################################################
### RW_PF_block, does a block RW, but using a INLA built function #####
#######################################################################################

#' @rdname samplers
#' @export
sampler_RW_INLA_blockV2 <- nimbleFunction(
  name = 'sampler_RW_INLA_blockV2',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',             FALSE)
    adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',       FALSE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',        200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent',  0.8)
    x                   <- extractControlElement(control, 'x',  double())
    y                   <- extractControlElement(control, 'y',  character())
    extraVars          <- extractControlElement(control, 'extraVars',              NULL)
    mu                  <- extractControlElement(control, 'mu',  NULL)
    #obsVars <- extractControlElement(control, 'obsVar',  character())
    fixedVals           <- extractControlElement(control, 'fixedVals',  double())
    fam                 <- extractControlElement(control, 'fam',  "gaussian")
    interVal            <- extractControlElement(control, 'interInModel',  1)
    scale               <- extractControlElement(control, 'scale',                1)
    propCov             <- extractControlElement(control, 'propCov',              'identity')
    existingINLA        <- extractControlElement(control, 'fit.inla',                   NULL)
    m                   <- extractControlElement(control, 'pfNparticles',         1000)
    filterType          <- extractControlElement(control, 'pfType',               'bootstrap')
    filterControl       <- extractControlElement(control, 'pfControl',            list())
    optimizeM           <- extractControlElement(control, 'pfOptimizeNparticles', FALSE)

    #Retuen scalar component of target variable
    targetAsScalar      <- model$expandNodeNames(target, returnScalarComponents = TRUE)

    # Return extraVars simulated from MCMC that will be needed with y for INLA
    extraVarsAsScalar <- model$expandNodeNames(extraVars, returnScalarComponents = TRUE)
    # Get dependencies of target vars
    calcNodes           <- model$getDependencies(target)
    if(length(fixedVals) > 0){
      latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)
    MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)

    topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)
    target <- model$expandNodeNames(target)
    ## numeric value generation
    optimizeM     <- as.integer(optimizeM)
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    prevLL        <- 0
    nVarEsts      <- 0
    itCount       <- 0
    d <- length(targetAsScalar)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    chol_propCov_scale <- scale * chol_propCov
    empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)

    storeParticleLP <- -Inf
    storeLLVar  <- 0
    nVarReps <- 7    ## number of LL estimates to compute to get each LL variance estimate for m optimization
    mBurnIn  <- 15   ## number of LL variance estimates to compute before deciding optimal m
    if(optimizeM)   m <- 3000
    ## nested function and function list definitions
    my_setAndCalculate <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, target, calcNodes)
    my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)

    #yVals <- values(model, obsData)
    my_particleFilter <- buildINLAmodelV2(model,
                                        fam,
                                        x,
                                        y = y,
                                        control = list(fit.inla = existingINLA,
                                                       fixedVals = fixedVals))

    #Target values for inla

    #get values of extra Vars
      targetVal <- nimble::values(model, targetAsScalar)
      extraVal <- nimble::values(model, extraVarsAsScalar)

      #return saved values of INLA
    particleMV <- my_particleFilter$mvEWSamples

    print(latentSamp)

    ## checks
    if(!inherits(propCov, 'matrix'))                    stop('propCov must be a matrix\n')
    if(!inherits(propCov[1,1], 'numeric'))              stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))                         stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))                           stop('propCov matrix must be symmetric')
    if(length(targetAsScalar) < 2)                      stop('less than two top-level targets; cannot use RW_PF_block sampler, try RW_PF sampler')
    # if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
  },
  run = function() {
    storeParticleLP <<- my_particleFilter$getLastLogLik()
    modelLP0 <- storeParticleLP + getLogProb(model, target)
    propValueVector <- generateProposalVector()
    my_setAndCalculate$run(propValueVector)
    targetVal <<- values(model, targetAsScalar)
    extraVal <<- values(model, extraVarsAsScalar)
    particleLP <- my_particleFilter$run(beta = targetVal, extraVars = extraVal)#,interInModel = interVal)
    modelLP1 <- particleLP + getLogProb(model, target)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
    if(!jump) {
      my_particleFilter$setLastLogLik(storeParticleLP)
    }
    if(jump & latentSamp) {
      ## if we jump, randomly sample latent nodes from pf output and put
      ## into model so that they can be monitored
      #   index <- ceiling(runif(1, 0, m))
      copy(particleMV, model, fixedVals, fixedVals, row = 1)
      calculate(model, fixedValsDep)
      copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    else if(!jump & latentSamp) {
      ## if we don't jump, replace model latent nodes with saved latent nodes
      copy(from = mvSaved, to = model, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    ##if(jump & !resample)  storeParticleLP <<- particleLP
    if(jump & optimizeM) optimM()
    if(adaptive)     adaptiveProcedure(jump)
  },
  methods = list(
    optimM = function() {
      tempM <- 15000
      declare(LLEst, double(1, nVarReps))
      if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
        for(i in 1:nVarReps)
          LLEst[i] <- my_particleFilter$run(beta = values(model, targetAsScalar), extraVars = values(model, extraVarsAsScalar))#, interInModel = interVal)
        ## next, store average of var estimates
        if(nVarEsts == 1)
          storeLLVar <<- var(LLEst)/mBurnIn
        else {
          LLVar <- storeLLVar
          LLVar <- LLVar + var(LLEst)/mBurnIn
          storeLLVar <<- LLVar
        }
        nVarEsts <<- nVarEsts + 1
      }
      else {  # once enough var estimates have been taken, use their average to compute m
        m <<- m*storeLLVar/(0.92^2)
        m <<- ceiling(m)
        storeParticleLP <<- my_particleFilter$run(targetVal, extraVal)
        optimizeM <<- 0
      }
    },
    generateProposalVector = function() {
      propValueVector <- rmnorm_chol(1, values(model, targetAsScalar), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$getGamma1()
          for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
          empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
        }
        chol_propCov_scale <<- chol_propCov * scale
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      chol_propCov_scale <<- chol_propCov * scale
      storeParticleLP <<- -Inf
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      my_calcAdaptationFactor$reset()
    }
  )
)

####################################################################
### virtual nimbleFunction template, included for ALL samplers #####
####################################################################



#' @rdname samplers
#' @export
#######################################################################################
### RW_PF_block, does a block RW, but using a INLA built function #####
#######################################################################################

#' @rdname samplers
#' @export
sampler_RW_INLAlatent_block <- nimbleFunction(
  name = 'sampler_RW_INLAlatent_block',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target,  control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',             FALSE)
    adaptScaleOnly      <- extractControlElement(control, 'adaptScaleOnly',       FALSE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',        200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent',  0.8)
    x                   <- extractControlElement(control, 'x',  double())
    y                   <- extractControlElement(control, 'y',  character())
    mu                  <- extractControlElement(control, 'mu',  NULL)
    #obsVars <- extractControlElement(control, 'obsVar',  character())
    fixedVals           <- extractControlElement(control, 'fixedVals',  double())
    fam                 <- extractControlElement(control, 'fam',  "gaussian")
    interVal            <- extractControlElement(control, 'interInModel',  1)
    scale               <- extractControlElement(control, 'scale',                1)
    propCov             <- extractControlElement(control, 'propCov',              'identity')
    existingINLA        <- extractControlElement(control, 'fit.inla',                   NULL)
    m                   <- extractControlElement(control, 'pfNparticles',         1000)
    filterType          <- extractControlElement(control, 'pfType',               'bootstrap')
    filterControl       <- extractControlElement(control, 'pfControl',            list())
    optimizeM           <- extractControlElement(control, 'pfOptimizeNparticles', FALSE)
    #Extract y values
    #obsData <- model$expandNodeNames(y)
    #y <- c(model[["y"]])
    ## node list generation
    targetAsScalar      <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes           <- model$getDependencies(target)
    if(length(fixedVals) > 0){
      latentSamp <- TRUE
    }else{
      latentSamp <- FALSE
    }
    fixedValsDep <- model$getDependencies(fixedVals)
    MCMCmonitors <- tryCatch(parent.frame(2)$conf$monitors, error = function(e) e)

    topParams <- model$getNodeNames(stochOnly=TRUE, includeData=FALSE, topOnly=TRUE)
    target <- model$expandNodeNames(target)
    ## numeric value generation
    optimizeM     <- as.integer(optimizeM)
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    prevLL        <- 0
    nVarEsts      <- 0
    itCount       <- 0
    d <- length(targetAsScalar)
    if(is.character(propCov) && propCov == 'identity')     propCov <- diag(d)
    propCovOriginal <- propCov
    chol_propCov <- chol(propCov)
    chol_propCov_scale <- scale * chol_propCov
    empirSamp <- matrix(0, nrow=adaptInterval, ncol=d)

    # if(is.null(mu)){
    #   muTarget <- nimble::values(model,target)
    # }else{
    #   muTarget <- mu
    # }
    storeParticleLP <- -Inf
    storeLLVar  <- 0
    nVarReps <- 7    ## number of LL estimates to compute to get each LL variance estimate for m optimization
    mBurnIn  <- 15   ## number of LL variance estimates to compute before deciding optimal m
    if(optimizeM)   m <- 3000
    ## nested function and function list definitions
    my_setAndCalculate <- setAndCalculate(model, target)
    my_decideAndJump <- decideAndJump(model, mvSaved, target, calcNodes)
    my_calcAdaptationFactor <- calcAdaptationFactor(d, adaptFactorExponent)

    #yVals <- values(model, obsData)
    my_particleFilter <- buildINLAmodel(model,
                                        fam,
                                        x,
                                        y = y,
                                        control = list(fit.inla = existingINLA,
                                                       fixedVals = fixedVals))
    targetVal <- nimble::values(model, targetAsScalar)
    particleMV <- my_particleFilter$mvEWSamples

    print(latentSamp)

    ## checks
    if(!inherits(propCov, 'matrix'))                    stop('propCov must be a matrix\n')
    if(!inherits(propCov[1,1], 'numeric'))              stop('propCov matrix must be numeric\n')
    if(!all(dim(propCov) == d))                         stop('propCov matrix must have dimension ', d, 'x', d, '\n')
    if(!isSymmetric(propCov))                           stop('propCov matrix must be symmetric')
    if(length(targetAsScalar) < 2)                      stop('less than two top-level targets; cannot use RW_PF_block sampler, try RW_PF sampler')
    # if(any(target%in%model$expandNodeNames(latents)))   stop('PMCMC \'target\' argument cannot include latent states')
  },
  run = function() {
    storeParticleLP <<- my_particleFilter$getLastLogLik()
    modelLP0 <- storeParticleLP + getLogProb(model, target)
    propValueVector <- generateProposalVector()
    my_setAndCalculate$run(propValueVector)
    targetVal <<- values(model, targetAsScalar)
    particleLP <- my_particleFilter$run(beta = targetVal)#,interInModel = interVal)
    modelLP1 <- particleLP + getLogProb(model, target)
    jump <- my_decideAndJump$run(modelLP1, modelLP0, 0, 0)
    if(!jump) {
      my_particleFilter$setLastLogLik(storeParticleLP)
    }
    if(jump & latentSamp) {
      ## if we jump, randomly sample latent nodes from pf output and put
      ## into model so that they can be monitored
      #   index <- ceiling(runif(1, 0, m))
      copy(particleMV, model, fixedVals, fixedVals, row = 1)
      calculate(model, fixedValsDep)
      copy(from = model, to = mvSaved, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    else if(!jump & latentSamp) {
      ## if we don't jump, replace model latent nodes with saved latent nodes
      copy(from = mvSaved, to = model, nodes = fixedValsDep, row = 1, logProb = TRUE)
    }
    ##if(jump & !resample)  storeParticleLP <<- particleLP
    if(jump & optimizeM) optimM()
    if(adaptive)     adaptiveProcedure(jump)
  },
  methods = list(
    optimM = function() {
      tempM <- 15000
      declare(LLEst, double(1, nVarReps))
      if(nVarEsts < mBurnIn) {  # checks whether we have enough var estimates to get good approximation
        for(i in 1:nVarReps)
          LLEst[i] <- my_particleFilter$run(beta = targetVal)#, interInModel = interVal)
        ## next, store average of var estimates
        if(nVarEsts == 1)
          storeLLVar <<- var(LLEst)/mBurnIn
        else {
          LLVar <- storeLLVar
          LLVar <- LLVar + var(LLEst)/mBurnIn
          storeLLVar <<- LLVar
        }
        nVarEsts <<- nVarEsts + 1
      }
      else {  # once enough var estimates have been taken, use their average to compute m
        m <<- m*storeLLVar/(0.92^2)
        m <<- ceiling(m)
        storeParticleLP <<- my_particleFilter$run(targetVal)
        optimizeM <<- 0
      }
    },
    generateProposalVector = function() {
      propValueVector <- rmnorm_chol(1, values(model, target), chol_propCov_scale, 0)  ## last argument specifies prec_param = FALSE
      returnType(double(1))
      return(propValueVector)
    },
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(!adaptScaleOnly)     empirSamp[timesRan, 1:d] <<- values(model, target)
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        adaptFactor <- my_calcAdaptationFactor$run(acceptanceRate)
        scale <<- scale * adaptFactor
        ## calculate empirical covariance, and adapt proposal covariance
        if(!adaptScaleOnly) {
          gamma1 <- my_calcAdaptationFactor$getGamma1()
          for(i in 1:d)     empirSamp[, i] <<- empirSamp[, i] - mean(empirSamp[, i])
          empirCov <- (t(empirSamp) %*% empirSamp) / (timesRan-1)
          propCov <<- propCov + gamma1 * (empirCov - propCov)
          chol_propCov <<- chol(propCov)
        }
        chol_propCov_scale <<- chol_propCov * scale
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    reset = function() {
      scale   <<- scaleOriginal
      propCov <<- propCovOriginal
      chol_propCov <<- chol(propCov)
      chol_propCov_scale <<- chol_propCov * scale
      storeParticleLP <<- -Inf
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      my_calcAdaptationFactor$reset()
    }
  )
)

