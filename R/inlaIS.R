# Adaptive multiple Importance Sampling with INLA

##  Contains code to IS
##  We have a build function (buildBootstrapFilter),
##  and step function.
importanceSamplingStepVirtual <- nimbleFunctionVirtual(
  run = function(meanBeta = double(1),
                 sigmaBeta = double(2),
                 #t = integer(0),
                 prevSamp = integer(0)) {
    returnType(double(0))
  },
  methods = list(
    updateBetaMean=function(){
      returnType(double(1))
      return(mu)
    },
    updateBetaSigma=function(){
      returnType(double(2))
      return(sigma)
    }
  )
)



#FUnction to convert mvEWSamples to matrix
convertToMatrix <- function(mvEWSamples, target){
  ret <- as.matrix(mvEWSamples, target)
  return(ret)
}

# nimbleconvertToMatrix <- nimble::nimbleRcall(
#   prototype = function(
# mvEWSamples = character(),
# target = character()
#   ) {},
#   returnType = double(2), # outcome is a vector
#   Rfun = 'convertToMatrix'
# )

# Sample values and run
impSampINLAstep <- nimbleFunction(
  name = 'impSampINLAstep',
  contains = importanceSamplingStepVirtual,
  setup = function(model,
                   mvEWSamples,
                   mvWSamples,
                   fixedVals,
                   iNode,
                   x, #covariates
                   y, #response variable
                   # interInModel,
                   fam,
                   proposal, #proposal distribution
                   beta, #indicates whether t = 0 or not
                  timeIndex,
                  vars,
                  nCores,
                  dfTdist
  ) {

    # setting up parameters
    N <- length(fixedVals) #length of INLA parameters
    ess <- 0

    #store simulated values of beta
    betaNames <- model$expandNodeNames(nodes = beta)
    nBetaSims <- length(betaNames)
    betaVals <- rep(0, length = length(betaNames))
    mu <- rep(0, length = nBetaSims)
    sigma <- diag(nBetaSims) #initialise sigma for updating
    nugget <- seq(1, iNode, 1) #Cummulative index for N's
    nNugget <- length(nugget) #length of N's


    indInc <- (iNode - 1)*timeIndex
    m <- timeIndex
    betaEsts <- matrix(0, nrow = m, ncol = nBetaSims)
    betaEstsUpd <- matrix(0, nrow = m+indInc, ncol = (nBetaSims+1))
    rts <- matrix(0, nrow = m+indInc, ncol = nBetaSims)
    rtsUpd <- matrix(0, nrow = m+indInc, ncol = nBetaSims)
    muEsts <- matrix(0, nrow = m+indInc, ncol = nBetaSims)
    #cummN <- timeIndex
    #T <- length(cummN)
    #No resampling will be done here
    #but we will save weights
  },
  run = function(meanBeta = double(1),
                 sigmaBeta = double(2),
                 #t = integer(0),
                 prevSamp= integer(0)) {
    returnType(double(0))
#print(1)
    wts <- numeric(m, init=FALSE)
    gamma <- numeric(m, init=FALSE)
    nugs <- numeric(nNugget, init = FALSE)

#copy results from mvEWSamples to mvWSamples
    # for(i in 1:m){
    #   nimCopy(mvEWSamples, mvWSamples, nodes = "gamma", row = i)
    #   nimCopy(mvEWSamples, mvWSamples, nodes = "logLike", row = i)
    #   nimCopy(mvEWSamples, mvWSamples, nodes = "wts", row = i)
    # }
   #mvWSamples <<-  mvEWSamples
#print(2)
Nt <- m # Number of samples at Nt
    for(i in 1:m){
      # For now, simulate beta's from proposal distribution
      if(proposal == "normal"){
        betaVals <<- rmnorm_chol(1, meanBeta, chol(sigmaBeta), prec_param = FALSE)
      }
      if(proposal == "studentT"){
        betaVals <<- rmvt_chol(1, mu = meanBeta, chol(sigmaBeta), df= dfTdist,prec_param = FALSE)
      }
      if(proposal == "prior"){
        model$simulate(nodes = betaNames)
        betaVals <<- values(model, betaNames)
      }
        #}
     # for(i in 1:m){
        betaEsts[i, 1:nBetaSims] <<- betaVals#mvEWSamples[beta, (i+indInc)] #nimbleconvertToMatrix(mvEWSamples, beta)

    }
#print(3)
      res <- nimbleINLA(x, y, beta= betaEsts, fixedVals,  family = fam, nCores = nCores)

      #Calulate weights
      #nugget[1:T] <- timeIndex * dmnorm_chol(betaVals, meanBeta, chol(sigmaBeta), prec_param = FALSE)
      #T <- length(nugget)
      #if(T ==1){
      #  gamma[i] <- iNode*timeIndex * dmnorm_chol(betaVals, meanBeta, chol(sigmaBeta), prec_param = FALSE)
      #}else{
        #Calulate weights
      for(i in 1:m){
        betaVals <<- betaEsts[i, 1:nBetaSims]
      for(j in 1:iNode){
        nugs[j] <- timeIndex * dmnorm_chol(betaVals, meanBeta, chol(sigmaBeta), prec_param = FALSE)
      }
        gamma[i] <- sum(nugs[1:iNode])
#print(4)
     # }
      #print(gamma)
      values(model, beta) <<- betaVals
      loglike <- res[i,1]*model$calculate(beta)
      #print(loglike)
      wts[i] <- (loglike)/(gamma[i]/(iNode*timeIndex))
      #print(99)
      mvWSamples["logLike",i][iNode] <<- loglike
      #print(999)
      mvEWSamples["logLike",i][iNode] <<- loglike
#print(5)
      # save results for current time
      saveResults(fixedVals, res, ind = i)
      k <- i + indInc
      nimCopy(model, mvEWSamples, nodes = beta, nodesTo = beta, row = k)
      nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = k)
#print(6)
      mvEWSamples["gamma",i][iNode] <<- gamma[i]
      mvEWSamples["wts", i][iNode] <<- wts[i]
#print(7)

#Saving output for updated mean and sd
betaEstsUpd[k,1:nBetaSims] <<- betaVals
betaEstsUpd[k,(nBetaSims+1)] <<- wts[i]

    if(iNode > 1){
        for(t in 1:(iNode-1)){
          indx <- i + (t-1)*m
          #print(i)
         # print(t)
          nimCopy(mvEWSamples, model, nodes = beta, nodesTo = beta, row = indx)
          betaValsNew <- values(model, beta)
          if(prevSamp == 1){
            if(proposal == "normal"){
              priorDist <- dmnorm_chol(betaValsNew, meanBeta, chol(sigmaBeta), prec_param = FALSE)
            }
            if(proposal == "studentT"){
              priorDist <- dmvt_chol(betaValsNew, mu = meanBeta, chol(sigmaBeta), df= dfTdist,prec_param = FALSE)
            }
            if(proposal == "prior"){
              priorDist <- model$calculate(beta)
              #model$simulate(nodes = betaNames)
              #betaVals <<- values(model, betaNames)
            }
        #priorDist <- dmnorm_chol(betaValsNew, meanBeta, chol(sigmaBeta), prec_param = FALSE)
        ret <-  mvEWSamples["gamma",i][t] #trying to trick nimble
        gammaUpd <- ret + (m * priorDist)
        mvEWSamples["gamma",i][t] <<-  gammaUpd
        mvEWSamples["wts",i][t] <<- (mvEWSamples["logLike",i][t])/(gammaUpd/(m * iNode))
        #print(indx)
         }
          betaEstsUpd[indx,1:nBetaSims] <<- betaValsNew
          betaEstsUpd[indx,(nBetaSims+1)] <<- mvEWSamples["wts",i][t]
        }


      #save results

      #calculate
      #wts[i]  <- model$calculate(calc_thisNode_deps)
    }#else{
     # nimCopy(mvEWSamples, model, nodes = beta, nodesTo = beta, row = i)
     # betaValsNew <- values(model, beta)
     # betaEstsUpd[i,1:nBetaSims] <<- betaValsNew
     # betaEstsUpd[i,(nBetaSims+1)] <<- wts[i]
    #}
}
    #Normalised weights for Effective Sample size
    nWeights <- wts/sum(wts)
    ess <<- 1/sum(nWeights^2)


    # mu and sigma updates
    K <- m + indInc
    #print(K)
    nWeightsUpd <- betaEstsUpd[1:K,(nBetaSims+1)]/sum(betaEstsUpd[1:K,(nBetaSims+1)])
    #print(nWeightsUpd)
    for(j in 1:nBetaSims){
    for(i in 1:K){
      muEsts[i, j] <<- nWeightsUpd[i] * betaEstsUpd[i, j]
    }
      mu[j] <<- sum(muEsts[1:K, j])
    }
    #for(i in 1:nBetaSims){
   ## mu[i] <<- #(t(nWeightsUpd[1:K]) %*% betaEstsUpd[1:K,1:nBetaSims])[1, 1:nBetaSims]
#}
#print(10000)
for(j in 1:nBetaSims){
     for(i in 1:K){
       #rr <- nWeightsUpd[i] * (betaEstsUpd[i, 1:nBetaSims] - mu[1:nBetaSims])
    rts[i, j] <<-  nWeightsUpd[i] * (betaEstsUpd[i, j] - mu[j])
    rtsUpd[i,j] <<- (betaEstsUpd[i, j] - mu[j])
     }
}
    sigma <<-  t(rts)%*%rtsUpd
      #(t(rts[1:K, 1:nBetaSims]))%*%((betaEstsUpd[1:K, 1:nBetaSims] - mu[1:nBetaSims]))

    print(mu)
    print(sigma)
    #save the values

lll <- res[1,1]

    # if(lll == -Inf){
    #   copy(mvEWSamples, model, nodes = fixedVals, row = 1)
    # }else{
    #   saveResults(fixedVals, res)
    #   copy( model, mvEWSamples, nodes = fixedVals, row = 1)
    # }

    return(ess)

  },
  methods = list(
    saveResults = function(fixedVals = character(1),
                           res = double(2),
                           ind = integer(0)){ #which index to subset
      n <- length(fixedVals)
      vals <- numeric(n, init = FALSE)
      #r <- character(0)
      #if(n > 1){
      for(i in seq_along(fixedVals)){
        # r <- fixedVals[i]
        vals[i] <- res[ind, i + 1]
        #model[[r]] <<- vals[i]
      }
      # }
      #for(i in seq_along(fixedVals)){
      #model[[fixedVals[i]]] <<- vals[i]
      #values(model, fixedVals[i]) <<- vals[i]
      #}
      values(model, fixedVals) <<- c(vals)
      #else{
      # vals <<- res[1, 2]
      #}
      return(vals)
      returnType(double(1))
    },
    updateBetaMean=function(){
      returnType(double(1))
      return(mu)
    },
    updateBetaSigma=function(){
      returnType(double(2))
      return(sigma)
    }
  )
)



#

inlaIS <- nimbleFunction(
  name = 'inlaIS',
  setup = function(model,
                   fam, x, y,
                   target, #beta
                   control = list()) {

    #control list extraction
    #proposal <- control[['proposal']]
    #initMean <- control[['initMean']]
    #initCov <- control[['initCov']]
    #silent <- control[['silent']]
    #timeIndex <- control[['timeIndex']]
    #initModel <- control[['initModel']]
   # x <- control[['x']] #how many time points were used to fit the reduced model
    #initModel <- control[['initModel']]
    #M <- control[['M']]
    nimbleINLA <- extractControlElement(control, 'nimbleINLA',  NULL)
    fixedVals <- extractControlElement(control, 'fixedVals',  double())
    proposal <- extractControlElement(control, 'proposal',  character())
    initMean <- extractControlElement(control, 'initMean',  NULL)
    initCov <- extractControlElement(control, 'initCov',  NULL)
    initModel <- extractControlElement(control, 'initModel',  TRUE)
    timeIndex <- extractControlElement(control, 'timeIndex',  double()) #Nt
    nSteps <- extractControlElement(control, 'nSteps',  integer()) #Number of steps at each direction
    nCores <- extractControlElement(control, 'nCores',  NULL)
    dfTdist <- extractControlElement(control, 'dfTdist',  NULL)

    nTarget <- length(model$expandNodeNames(nodes = target))
    if(is.null(initMean)) initMean <- rep(0, nTarget)
    if(is.null(initCov)) initCov <- diag(nTarget)
    if(is.null(nCores)) nCores <- 1
    if(is.null(dfTdist)) dfTdist <- 1

    if(!proposal %in% c("normal", "studentT", "prior")) stop("Proposal distribution must be either 'normal', 'student T' and 'prior' distributions.")

    yExpand <- model$expandNodeNames(y, returnScalarComponents = TRUE)
    y <- model[[y]]
    #y <- c(nimble::values(model, yExpand))
    my_initializeModel <- initializeModel(model, silent = TRUE)
    #save posterior samples
    #modelVals = modelValues(model, m = 1)
    vars <- model$getVarNames(nodes = c(fixedVals, target))
    modelSymbolObjects <- model$getSymbolTable()$getSymbolObjects()[vars]

    names <- sapply(modelSymbolObjects, function(x)return(x$name))
    type <- sapply(modelSymbolObjects, function(x)return(x$type))
    size <- lapply(modelSymbolObjects, function(x)return(x$size))
    size <- lapply(size, function(x){
      if(length(x) == 0){
        #length(model$vars)
        ret <- 1 #c(1, timeIndex)
       # return(ret)
      }else{
        ret <- x #c(x, timeIndex)
      }
      return(ret)
    } )

    #size$beta <- 4
    # Add names and dimensions for wts and gamma
    names <- c(names, "wts", "gamma","logLike")
    type <- c(type, "double", "double", "double")
    size$wts <- nSteps
    size$gamma <- nSteps
    size$logLike <- nSteps
#print(1)
    #model values to save results
    mvEWSamples <- modelValues(modelValuesConf(vars = names,
                                               types = type,
                                               sizes = size))

    mvWSamples <- modelValues(modelValuesConf(vars = names,
                                               types = type,
                                               sizes = size))

    fixedVals <- model$expandNodeNames(fixedVals)

    multiple <- TRUE
    #print(multiple)
    if(length(model$expandNodeNames(fixedVals)) == 1) multiple = FALSE
    # vals <- rep(0, length(fixedVals))
    #for(i in seq_along(fixedVals)){
    # values(model, fixedVals[i]) <- c(res[1,fixedVals[i]])
    #}
    #vals <- c(res[1,fixedVals])
    impSampINLAstepFnx <- nimbleFunctionList(importanceSamplingStepVirtual)
    #for(iNode in seq_along(nodes)){
    beta <- target
    for(iNode in 1:nSteps){
    impSampINLAstepFnx[[iNode]] <- impSampINLAstep(model,
                                                mvEWSamples,
                                                mvWSamples,
                                                fixedVals,
                                                iNode,
                                                x, #covariates
                                                y, #response variable
                                                # interInModel,
                                                fam,
                                                proposal, #proposal distribution
                                                beta, #variable names for beta)
                                               timeIndex,
                                               vars,
                                               nCores,
                                               dfTdist
)
    }
    #}
    #essVals <- rep(0, length(nodes))
    meanBeta <- initMean; sigmaBeta <- initCov
    cummN <- timeIndex
    m <- timeIndex
    #print(m)
   # T <- length(cummN)
    #M <- sum(cummN[1:T])
#essVals <- rep(0, nSteps)
pp <- 0
    lastLogLik <- -Inf
  },
  run = function( ) {
    returnType(integer(0))
    declare(essVals, double())
    if(initModel) my_initializeModel$run()
nIter <- nSteps*timeIndex
    resize(mvEWSamples, nIter)
    resize(mvWSamples, nIter)
    essVals <- 0

    for(iNode in seq_along(impSampINLAstepFnx)){
     # if(iNode == 1 | iNode == length(impSampINLAstepFnx)){
    #
    # impSampINLAstepFnx[[iNode]]$run(m, meanBeta = meanBeta, sigmaBeta = sigmaBeta, prevSamp= 0)
    # sigmaBeta <<- impSampINLAstepFnx[[t]]$updateBetaSigma()
    # meanBeta <<-   impSampINLAstepFnx[[t]]$updateBetaMean()
    # }else{
        #prevSamp = 1
      essVals <- essVals + impSampINLAstepFnx[[iNode]]$run(meanBeta = meanBeta, sigmaBeta=sigmaBeta, prevSamp = pp)
        meanBeta <<-   impSampINLAstepFnx[[iNode]]$updateBetaMean()
        sigmaBeta <<- impSampINLAstepFnx[[iNode]]$updateBetaSigma()
         #impSampINLAstepFnx[[iNode]]$returnESS()
        pp <<- 1
        #}
}
    return(essVals)
  },
  methods = list(
    getLastmeanBeta = function() {
      return(meanBeta)
      returnType(double(1))
    },
    getLastSigmaBeta = function(){
      return(sigmaBeta)
      returnType(double(2))
    }
  )
)
