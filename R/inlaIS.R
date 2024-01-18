# Adaptive multiple Importance Sampling with INLA

##  Contains code to IS
##  We have a build function (buildBootstrapFilter),
##  and step function.
importanceSamplingStepVirtual <- nimbleFunctionVirtual(
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
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
                  dfTdist,
                  adaptive
  ) {

    #Note
    # iNode = Nt
    #timeIndex = number of samples at each it = m
    m <- timeIndex

    # setting up parameters
    N <- length(fixedVals) #length of INLA parameters
    ess <- 0

    #store simulated values of beta
    betaNames <- model$expandNodeNames(nodes = beta)
    nBetaSims <- length(betaNames) #dimeension of beta
    betaVals <- rep(0, length = length(betaNames)) #vector to store beta values

    # store updated mean and sigma
    mu <- rep(0, length = nBetaSims)
    sigma <- matrix(0, nrow=nBetaSims, ncol = nBetaSims)

    # store weights and gamma
    wts <- numeric(m)
    gamma <- numeric(m)

    # Create a cummulative indexing to calculate nugget
    nugget <- seq(1, iNode, 1) #Cummulative index for N's
    nNugget <- length(nugget) #length of N's
    nugs <- numeric(nNugget)

    # Estimate increasing indices
    indInc <- (iNode - 1)*m #sum from N1 to Nt-1
    sumNt <- m + indInc

    #beta estimates at iteration t
    betaEsts <- matrix(0, nrow = m, ncol = nBetaSims)

    # Get beta estimates from steps 1 to Nt
    betaEstsUpd <- matrix(0, nrow = sumNt, ncol = (nBetaSims+1))
    rts <- matrix(0, nrow = sumNt, ncol = nBetaSims)
    rtsUpd <- matrix(0, nrow = sumNt, ncol = nBetaSims)
    muEsts <- matrix(0, nrow = sumNt, ncol = nBetaSims)

  },
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
                 #t = integer(0),
                 prevSamp= integer(0)) {
    returnType(double(0))

    k <- indInc + 1
    m <<- timeIndex


    # simulate samples
    for(i in 1:m){

      # For now, simulate beta's from proposal distribution
      if(proposal == "normal"){
        betaVals <<- rmnorm_chol(1, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE)
      }
      if(proposal == "studentT"){
        betaVals <<- rmvt_chol(1, mu = meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE)
      }
      if(proposal == "prior"){
        model$simulate(nodes = betaNames)
        betaVals <<- values(model, betaNames)
      }

        betaEsts[i, 1:nBetaSims] <<- betaVals
    }

      res <- nimbleINLA(x, y, beta= betaEsts, fixedVals,  family = fam, nCores = nCores)


        #Calulate weights
      for(i in 1:m){
        betaVals <<- betaEsts[i, 1:nBetaSims]
        nn <- 0 # start indexing to estimate gamma
      for(j in 1:iNode){
        if(proposal == "normal"){
         nn <- nn +  timeIndex * dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE)
        }
        if(proposal == "studentT"){
          nn <- nn +  timeIndex * dmvt_chol(betaVals, mu = meanBeta[j,], chol(sigmaBeta[,,j]), df= dfTdist, prec_param = FALSE, log = FALSE)
        }
        #nugs[j] <<- m * dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE)
      }
        gamma[i] <<- nn
        #gamma[i] <<- sum(nugs[1:iNode])

        # calculate numerator of wts in paper
      values(model, beta) <<- betaVals
      loglike <- res[i,1] + model$calculate(beta)


      #calculate weights
      if(iNode > 1){
      wts[i] <<- loglike - log(gamma[i]/sumNt)
      }else{
        if(proposal == "normal"){
  wts[i] <<- loglike - dmnorm_chol(betaVals, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = TRUE)
        }
        if(proposal == "studentT"){
          wts[i] <<- loglike - dmvt_chol(betaVals, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = TRUE)
        }
        if(proposal == "prior"){
          wts[i] <<- loglike - model$calculate(betaVals)
        }
  }
      # save the numerator of the weights for the updating step
      mvWSamples["logLike",i][iNode] <<- loglike
      mvEWSamples["logLike",i][iNode] <<- loglike

      # save INLA parameter samples for current time
      saveResults(fixedVals, res, ind = i)
      nimCopy(model, mvEWSamples, nodes = beta, nodesTo = beta, row=1, rowTo = k)
      nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)

      #save wts and gamma at Nt, since it won't be updates
       mvEWSamples["gamma",i][iNode] <<- gamma[i]
      mvEWSamples["wts", i][iNode] <<- wts[i]

#Saving output for updated mean and sd
betaEstsUpd[k,1:nBetaSims] <<- betaVals
betaEstsUpd[k,(nBetaSims+1)] <<- wts[i]
k <- k + 1
      }

      #Normalised weights for Effective Sample size
      maxWts <- max(wts)
      nWeights <- exp(wts - maxWts)
      nWeights <- nWeights/sum(nWeights)
      ess <<- 1/sum(nWeights^2)

#update wts and gamma for iNode > 1

if(iNode > 1 & adaptive == TRUE){
for(i in 1:m){
      for(t in 1:(iNode-1)){
        #index of weights to update
          indx <- i + (t-1)*m

          nimCopy(mvEWSamples, model, nodes = beta, nodesTo = beta, row = indx, rowTo = 1)
          betaValsNew <- values(model, beta)
          #if(prevSamp == 1){
            if(proposal == "normal"){
              priorDist <- dmnorm_chol(betaValsNew, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = FALSE)
            }
            if(proposal == "studentT"){
              priorDist <- dmvt_chol(betaValsNew, mu = meanBeta[iNode, ], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE)
            }
            if(proposal == "prior"){
              priorDist <- exp(model$calculate(beta))
              #model$simulate(nodes = betaNames)
              #betaVals <<- values(model, betaNames)
            }

          #estimate gamma update
        ret <-  mvEWSamples["gamma",i][t] #trying to trick nimble
        gammaUpd <- ret + (m * priorDist) #note that Nt = m
        mvEWSamples["gamma",i][t] <<-  gammaUpd
        mvEWSamples["wts",i][t] <<- mvEWSamples["logLike",i][t] - log(gammaUpd/sumNt)
          betaEstsUpd[indx, 1:nBetaSims] <<- betaValsNew
          betaEstsUpd[indx,(nBetaSims+1)] <<- mvEWSamples["wts",i][t]
        }
    }

}

      # calculate updated weights and gamma
  if(adaptive ==TRUE){
wtsUpd <- betaEstsUpd[1:sumNt,(nBetaSims+1)]
maxWtsUpd <- max(wts)
nWeightsUpd <- exp(wtsUpd - maxWtsUpd)
#sq <-sum(nWeightsUpd)
#print(sq)

    for(j in 1:nBetaSims){
    for(i in 1:sumNt){
      muEsts[i, j] <<- nWeightsUpd[i] * betaEstsUpd[i, j]
    }
      mu[j] <<- sum(muEsts[1:sumNt, j])/sum(nWeightsUpd[1:sumNt])
    }
    #for(i in 1:nBetaSims){
   ## mu[i] <<- #(t(nWeightsUpd[1:K]) %*% betaEstsUpd[1:K,1:nBetaSims])[1, 1:nBetaSims]
#}
#print(10000)
#for(j in 1:nBetaSims){
     for(i in 1:nBetaSims){
       indx <- i
       for(j in indx:nBetaSims){
     sigma[i,j] <<-  sum(nWeightsUpd[1:sumNt]*(betaEstsUpd[1:sumNt, i] - mu[i]) * (betaEstsUpd[1:sumNt, j] - mu[j]))/sum(nWeightsUpd[1:sumNt])
      sigma[j,i] <<-  sum(nWeightsUpd[1:sumNt]*(betaEstsUpd[1:sumNt, i] - mu[i]) * (betaEstsUpd[1:sumNt, j] - mu[j]))/sum(nWeightsUpd[1:sumNt])
       }
     }
}else{
       mu <<- meanBeta[iNode,]
       sigma <<- sigmaBeta[,,iNode]
     }
       #rr <- nWeightsUpd[i] * (betaEstsUpd[i, 1:nBetaSims] - mu[1:nBetaSims])
    #rts[i, j] <<-  nWeightsUpd[i] * (betaEstsUpd[i, j] - mu[j])
    #rtsUpd[i,j] <<- (betaEstsUpd[i, j] - mu[j])
   #  }
      # sum(weight[1:i_tot]*(eta[1:i_tot,i]-theta[[1]][i])*
      #       (eta[1:i_tot,j]-theta[[1]][j]))/(sum(weight[1:i_tot]))
#sigma <<- sigma + nWeightsUpd[i] * (t(t(betaEstsUpd[i, 1:nBetaSims] - mu[1:nBetaSims]))%*%t(betaEstsUpd[i, 1:nBetaSims] - mu[1:nBetaSims]))

    #sigma <<-  t(rts)%*%rtsUpd
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
    adaptive <- extractControlElement(control, 'adaptive',  TRUE)


    #Note
    # In the original paper:
    ## timeIndex = Nt
    ## nSteps = t

    nTarget <- length(model$expandNodeNames(nodes = target))
    if(is.null(initMean)) initMean <- rep(0, nTarget)
    if(is.null(initCov)) initCov <- diag(nTarget)
    if(is.null(nCores)) nCores <- 1
    if(is.null(dfTdist)) dfTdist <- 1

    if(!proposal %in% c("normal", "studentT", "prior")) stop("Proposal distribution must be either 'normal', 'student T' and 'prior' distributions.")


    if(proposal == "prior") adaptive <- FALSE
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

    if("gamma" %in% names) stop("change the variable name of gamma.")
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

    # create storage for mean and covariance matrices
    muStep <- matrix(0, nrow = nSteps+1, ncol = nTarget)
    sigmaStep <- array(0, dim = c(nTarget, nTarget, nSteps+1))

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
                                               dfTdist,
                                               adaptive
)
    }


    muStep[1,] <- initMean
    sigmaStep[,,1] <- initCov
    m <- timeIndex

    # prevSamp

  },
  run = function( ) {
    returnType(integer(0))
    declare(essVals, double())
    if(initModel) my_initializeModel$run()
nIter <- nSteps*timeIndex
    resize(mvEWSamples, nIter)
    resize(mvWSamples, nIter)
    essVals <- 0
    pp <- 0

    # runs for each Nt
    for(iNode in seq_along(impSampINLAstepFnx)){
      essVals <- essVals + impSampINLAstepFnx[[iNode]]$run(meanBeta = muStep, sigmaBeta=sigmaStep, prevSamp = pp)
      #if(iNode < nSteps){
        muStep[iNode+1,] <<-   impSampINLAstepFnx[[iNode]]$updateBetaMean()
        sigmaStep[ , ,iNode+1] <<- impSampINLAstepFnx[[iNode]]$updateBetaSigma()
      #}
         #impSampINLAstepFnx[[iNode]]$returnESS()
        pp <- 1
        #}
}
    return(essVals)
  },
  methods = list(
    getLastmeanBeta = function() {
      return(muStep)
      returnType(double(2))
    },
    getLastSigmaBeta = function(){
      return(sigmaStep)
      returnType(double(3))
    }
  )
)
