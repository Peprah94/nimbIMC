

# nimbleconvertToMatrix <- nimble::nimbleRcall(
#   prototype = function(
    # mvEWSamples = character(),
# target = character()
#   ) {},
#   returnType = double(2), # outcome is a vector
#   Rfun = 'convertToMatrix'
# )

# Adaptive multiple Importance Sampling with INLA

##  Contains code to IS
##  We have a build function (buildBootstrapFilter),
##  and step function.
importanceSamplingStepVirtualMultiple <- nimbleFunctionVirtual(
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
                 meanBetaVar2 = double(2),
                 sigmaBetaVar2 = double(3),
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
    },
    updateBetaMeanVar2=function(){
      returnType(double(1))
      return(muObsParams)
    },
    updateBetaSigmaVar2=function(){
      returnType(double(2))
      return(sigmaObsParams)
    }
  )
)

# Sample values and run
impSampINLAstepMultiple <- nimbleFunction(
  name = 'impSampINLAstep',
  contains = importanceSamplingStepVirtualMultiple,
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
                   isNotDiscreteTargetVar1, #indicates whether t = 0 or not
                   isNotDiscreteTargetVar2,
                   discreteTarget,
                   timeIndex,
                   vars,
                   nCores,
                   dfTdist,
                   adaptive,
                   additionalPars,
                   dataVar,
                   latentIsDependent,
                   returnLinearPred,
                   linearPred
  ) {

    #Note
    # iNode = Nt
    #timeIndex = number of samples at each it = m
    m <- timeIndex

    #check if it is the first index
    isFirstNode <- iNode == 1

    #set continuous vars to beta and discrete vars to disc
    beta <- isNotDiscreteTargetVar1
    obsParams <- isNotDiscreteTargetVar2
    disc <- discreteTarget

    # setting up parameters
    N <- length(fixedVals) #length of INLA parameters
    ess <- 0 # initialise effective sample size

    #store simulated values of beta (ecological parameters)
    betaNames <- model$expandNodeNames(nodes = beta)
    nBetaSims <- length(betaNames) #dimeension of beta
    betaVals <- rep(0, length = length(betaNames)) #vector to store beta values
    depsbetaNames <- model$getDependencies(betaNames, self = FALSE, determOnly = TRUE)

    #store simulated values of discrete (latent) target
    discTarNames <- model$expandNodeNames(nodes = discreteTarget)
    nDiscTarNames <- length(discTarNames)
    discTarVals <- rep(0, length = nDiscTarNames)




    #check if the discrete random variable is binary
    isLatentBinary <- all(model$isBinary(discTarNames) == TRUE)
    #print(isLatentBinary)
   # if(isLatentBinary == TRUE){
      #select those that are data
    binNodesToSimulate <- discTarNames[!model$isData(discTarNames)]
    nBinNodesToSimulate <- length(binNodesToSimulate)
    binNodesToFix <- discTarNames[model$isData(discTarNames)]
    nBinNodesToFix <- length(binNodesToFix)
    binNodesToSimulateVals <- rep(0, length = nBinNodesToSimulate)
    binNodesToFixVals <- rep(1, length = nBinNodesToFix)
   # } else {
    #   binNodesToSimulate <- NULL
    #   nBinNodesToSimulate <- NULL
    #   binNodesToFix <- NULL
    #   nBinNodesToFix <- NULL
    #   binNodesToSimulateVals <- NULL
    #   binNodesToFixVals <- NULL
    # }

    #get dependencies of both parameters. Will be used to calculate weights at each time
    allISparameters <- c(beta, discreteTarget)
    allISparametersDeps <-   model$getDependencies(allISparameters)

    # store updated mean and sigma
    mu <- rep(0, length = nBetaSims)
    sigma <- matrix(0, nrow=nBetaSims, ncol = nBetaSims)

    # store weights and gamma for ecological process parameters
    wts <- numeric(m)
    gamma <- numeric(m)
    betaWts <- numeric(m)


    # store weights for latent variables
    wtsLatent <- numeric(m)


    # Create a cummulative indexing to calculate nugget
    nugget <- seq(1, iNode, 1) #Cummulative index for N's
    nNugget <- length(nugget) #length of N's
    nugs <- numeric(nNugget)

    # Estimate increasing indices
    indInc <- (iNode - 1)*m #sum from N1 to Nt-1
    sumNt <- m + indInc

    #beta estimates at iteration t
    betaEsts <- matrix(0, nrow = m, ncol = nBetaSims)
    discTarEsts <- matrix(0, nrow = m, ncol = nDiscTarNames)

    # Get beta estimates from steps 1 to Nt
    betaEstsUpd <- matrix(0, nrow = sumNt, ncol = (nBetaSims+1))
    rts <- matrix(0, nrow = sumNt, ncol = nBetaSims)
    rtsUpd <- matrix(0, nrow = sumNt, ncol = nBetaSims)
    muEsts <- matrix(0, nrow = sumNt, ncol = nBetaSims)



if (!is.null(obsParams)){
  # store simulated values for observation parameters (obsPars)
  obsParamsNames <- model$expandNodeNames(nodes = obsParams)
  nObsParamsSims <- length(obsParamsNames) #dimeension of beta
  obsParamsVals <- rep(0, length = length(obsParamsNames)) #vector to store beta values
  depsObsParamsNames <- model$getDependencies(obsParamsNames, self = FALSE, determOnly = TRUE)

  # store weights and gamma for observation process parameters
  wtsObsParams <- numeric(m)
  gammaObsParams <- numeric(m)
  obsParamsWts <- numeric(m)

  # Observation Params at iteration t
  obsParamsEsts <- matrix(0, nrow = m, ncol = nObsParamsSims)

  # Get observation Params estimates from steps 1 to Nt
  obsParamsEstsUpd <- matrix(0, nrow = sumNt, ncol = (nObsParamsSims))
  rtsObsParams <- matrix(0, nrow = sumNt, ncol = nObsParamsSims)
  rtsObsParamsUpd <- matrix(0, nrow = sumNt, ncol = nObsParamsSims)
  muObsParamsEsts <- matrix(0, nrow = sumNt, ncol = nObsParamsSims)

  # store updated mean and sigma
  muObsParams <- rep(0, length = nObsParamsSims)
  sigmaObsParams <- matrix(0, nrow=nObsParamsSims, ncol = nObsParamsSims)
}

# Set the proposal distribution
if(length(proposal) == 1){
  betaProposal <- proposal
  obsParamsProposal <- proposal
} else {
  betaProposal <- proposal[1]
  obsParamsProposal <- proposal[2]
}


    #save the weights of beta and latent variables


    #check if we have additional Parameters (especially derived quantities) to save
    isNullAdditionalPars <- is.null(additionalPars)
    isNullObsParams <- is.null(obsParams)

    print(paste0("Seperate ecological and observation process:", isNullObsParams))
  fixedValsStoreMatrix <- matrix(0, nrow = timeIndex, ncol = length(fixedVals))
  #linearPredMatrix <- matrix(0, nrow = timeIndex, ncol = nDiscTarNames)
 # linearPreds <- model$getDependencies(fixedVals, determOnly = TRUE)

  lenLinPred <- length(linearPred)
  if(returnLinearPred) linearPredVals <- numeric(length(linearPred))


  },
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
                 meanBetaVar2 = double(2),
                 sigmaBetaVar2 = double(3),
                 #t = integer(0),
                 prevSamp= integer(0)) {
    returnType(double(0))

    k <- indInc + 1
    m <<- timeIndex
    linPredIndUpper <- 0
    linPredIndLower <- 0

    if(isNullObsParams){

      proposal <<- betaProposal

    if(isFirstNode & returnLinearPred) linearPredVals <<- values(model, linearPred)
#print(meanBeta[iNode,])
#print(sigmaBeta[,,iNode])
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
      # If we need linear Predictor to simulate latent nodes
    #then we fit the model z_{i-1}|beta_sim
    if(returnLinearPred)  res <- nimbleINLA(x, y, beta= betaEsts, extraVars = discTarEsts,fixedVals,  family = fam, nCores = nCores)

      #simulate discrete rvs
    for(i in 1:m){
      values(model, betaNames) <<- betaEsts[i, 1:nBetaSims]

      #for independence in latent variable and alternative MCMC sampler

      if(!latentIsDependent){
        if(iNode > 1){
        values(model, fixedVals) <<- fixedValsStoreMatrix[i, ]
       # values(model, linearPreds) <<- linearPredMatrix[i, ]
        }
      }

    if(latentIsDependent) model$calculate(depsbetaNames)
      #model$calculate()

      #if(!latentIsDependent) values(model, fixedVals) <<- fixedValsStoreMatrix[i, ]
      if(returnLinearPred){
        linPredIndLower <- (lenLinPred * (i -1)) +1
        linPredIndUpper <- (lenLinPred * i)
        linearPredVals <<- res[linPredIndLower: linPredIndUpper, N+2]
        values(model, linearPred) <<- linearPredVals
        }
     # model$calculate(depsDiscTarNames)
      model$calculate()

#print(isLatentBinary)
      if(isLatentBinary == FALSE){
      model$simulate(nodes = discTarNames)
      } else {
        model$simulate(nodes = binNodesToSimulate)
        values(model, binNodesToFix) <<- binNodesToFixVals
      }


      discTarEsts[i, 1:nDiscTarNames] <<- values(model, discTarNames)
      #print(discTarEsts[i, 1:nDiscTarNames])
      }

    if(!returnLinearPred) res <- nimbleINLA(x, y, beta= betaEsts, extraVars = discTarEsts,fixedVals,  family = fam, nCores = nCores)


    #Calulate weights
    for(i in 1:m){
      betaVals <<- betaEsts[i, 1:nBetaSims]
      values(model, beta) <<- betaVals
      values(model, discTarNames) <<- discTarEsts[i, 1:nDiscTarNames]
      #update dependencies of all parameters
      model$calculate(allISparametersDeps)

      # call the calculate of dependencies of all parameters to update logprobability of all nodes
     logProbAllPars <-  model$calculate(allISparameters)


     #####################
     # Update past importance weights
     ##########################
      nn <- 0
      for(j in 1:iNode){
        if(proposal == "normal"){
          #note that getLogProb returns the log of the log density. To get the density, we take exp of the logProb
          nn <- nn +  timeIndex * (dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE)) #+ exp(model$getLogProb(discTarNames)))
          }
        if(proposal == "studentT"){
          nn <- nn +  timeIndex * (dmvt_chol(betaVals, mu = meanBeta[j,], chol(sigmaBeta[,,j]), df= dfTdist, prec_param = FALSE, log = FALSE)) #+ exp(model$getLogProb(discTarNames)))
        }
        #nugs[j] <<- m * dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE)
      }
      gamma[i] <<- nn
      #gamma[i] <<- sum(nugs[1:iNode])

      ###############################
      # Estimating weights
      ################################

      # calculate numerator of wts in paper

      #values(model, discTarNames) <<- discTarEsts[i, 1:nDiscTarNames]

# log likelihood should include the contribution from beta and latent variables

    if(returnLinearPred){
 #N = length of fixed Vals, and the fist colum returns the marginal likelihood
      mld <- res[linPredIndLower, 1]
    } else {
        mld <- res[i, 1]
      }





      # log likelihood should include the contribution from beta and latent variables
     # ii <- (i-1)*nSites + 1

      if(latentIsDependent){
      loglike <- mld + model$calculate(beta)
      }else{
  loglike <- model$calculate(c(dataVar, beta))
      }
      #linearPredMatrix[i, ] <- res[,2]

      #calculate weights
      if(iNode > 1){
        wts[i] <<- loglike - log(gamma[i]/sumNt)
      }else{
        if(proposal == "normal"){
          wts[i] <<- loglike - dmnorm_chol(betaVals, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = TRUE) #- model$getLogProb(discTarNames)
        }
        if(proposal == "studentT"){
          wts[i] <<- loglike - dmvt_chol(betaVals, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = TRUE)# -  model$getLogProb(discTarNames)
        }
        if(proposal == "prior"){
          wts[i] <<- loglike - model$calculate(betaVals)
        }
      }

      # Calculate weights of latent
      # For latent variable, we estime weights at next level as
      #wts[i] <- p(yt +1|N)*wts[i-1]
      if(iNode > 1){
      wtsLatent[i] <<- model$calculate(discTarNames) #+ mvEWSamples["wtsLatent", i][iNode-1]
      }else{
        wtsLatent[i] <<- model$calculate(discTarNames)
      }
      # save the numerator of the weights for the updating step
      mvWSamples["logLike",i][iNode] <<- loglike
      mvEWSamples["logLike",i][iNode] <<- loglike

      # save INLA parameter samples for current time
      if(returnLinearPred){
        #N = length of fixed Vals, and the fist colum returns the marginal likelihood
        linPredIndLower <- (lenLinPred * (i -1)) +1
        saveResults(fixedVals, res, ind = linPredIndLower)
      } else {
        saveResults(fixedVals, res, ind = i)
      }

      if(!latentIsDependent){
        fixedValsStoreMatrix[i, ] <<- values(model, fixedVals)
      }

      #values(model, betaWts) <<- wts[i]
      #values(model, latentWts) <<- wtsLatent[i]

      if(isNullAdditionalPars){
        nimCopy(model, mvEWSamples, nodes = beta, nodesTo = beta, row=1, rowTo = k)
        nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)
        nimCopy(model, mvEWSamples, nodes = discTarNames, nodesTo = discTarNames, row = 1, rowTo = k)
        #nimCopy(model, mvEWSamples, nodes = latentWts, nodesTo = latentWts, row = 1, rowTo = k)
        #nimCopy(model, mvEWSamples, nodes = betaWts, nodesTo = betaWts, row = 1, rowTo = k)
      }else{
        nimCopy(model, mvEWSamples, nodes = beta, nodesTo = beta, row=1, rowTo = k)
        nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)
        nimCopy(model, mvEWSamples, nodes = discTarNames, nodesTo = discTarNames, row = 1, rowTo = k)
        nimCopy(model, mvEWSamples, nodes = additionalPars, nodesTo = additionalPars, row = 1, rowTo = k)
        #nimCopy(model, mvEWSamples, nodes = latentWts, nodesTo = latentWts, row = 1, rowTo = k)
        #nimCopy(model, mvEWSamples, nodes = betaWts, nodesTo = betaWts, row = 1, rowTo = k)
        }
      #save wts and gamma at Nt, since it won't be updates
      mvEWSamples["gamma",i][iNode] <<- gamma[i]
      mvEWSamples["wts", i][iNode] <<- wts[i]
      mvEWSamples["wtsLatent", i][iNode] <<- wtsLatent[i]

      #Saving output for updated mean and sd
      betaEstsUpd[k,1:nBetaSims] <<- betaVals
      betaEstsUpd[k,(nBetaSims+1)] <<- wts[i]
      print(k)

      #wtsLatVal <- wtsLatent[i]
     mvEWSamples["latentWeights", k][1] <<- wtsLatent[i]  #mvEWSamples["wtsLatent", i][iNode]
      mvEWSamples["betaWeights", k][1] <<- wts[i]

      k <- k + 1
    }

    #Normalised weights for Effective Sample size
    maxWts <- max(wts)
    nWeights <- exp(wts - maxWts)
    nWeights <- nWeights/sum(nWeights)
    ess <<- 1/sum(nWeights^2)

    ######################################
    #update past importance weights and gamma for iNode > 1
    ###################################

    if(iNode > 1 & adaptive == TRUE){
      for(i in 1:m){
        for(t in 1:(iNode-1)){
          #index of weights to update
          indx <- i + (t-1)*m

          nimCopy(mvEWSamples, model, nodes = allISparameters, nodesTo = allISparameters, row = indx, rowTo = 1)
          betaValsNew <- values(model, beta)

          # Update the model logProb with copied values
          model$calculate(allISparametersDeps)
          allIsparslike <- model$calculate(allISparameters)
          #print(allIsparslike)
          #if(prevSamp == 1){
          if(proposal == "normal"){
            priorDist <- dmnorm_chol(betaValsNew, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
          }
          if(proposal == "studentT"){
            priorDist <- dmvt_chol(betaValsNew, mu = meanBeta[iNode, ], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
          }
          if(proposal == "prior"){
            #priorDist <- exp(allIsparslike) #exp(model$calculate(beta))
            priorDist <- exp(model$calculate(beta))
            #model$simulate(nodes = betaNames)
            #betaVals <<- values(model, betaNames)
          }

          #estimate gamma update
          ret <-  mvEWSamples["gamma",i][t] #trying to trick nimble
          gammaUpd <- ret + (m * priorDist) #note that Nt = m
          #print(gammaUpd)
          mvEWSamples["gamma",i][t] <<-  gammaUpd
          mvEWSamples["wts",i][t] <<- mvEWSamples["logLike",i][t] - log(gammaUpd/sumNt)
          betaEstsUpd[indx, 1:nBetaSims] <<- betaValsNew
          betaEstsUpd[indx,(nBetaSims+1)] <<- mvEWSamples["wts",i][t]
          #mvEWSamples["betaWts", indx] <- mvEWSamples["wts",i][t]
        }
      }

    }

    ##########################################
    # calculate updated mean and covariance matrix for prior distribution
    #########################################

    if(adaptive ==TRUE){
      wtsUpd <- betaEstsUpd[1:sumNt,(nBetaSims+1)]
      maxWtsUpd <- max(wts)
      nWeightsUpd <- exp(wtsUpd - maxWtsUpd)
      betaWts <<- nWeightsUpd
      #print(nWeightsUpd)
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
      betaWts <<- nWeights
    }

    print(mu)
    print(sigma)
    print(muObsParams)
    print(sigmaObsParams)

    } else {

      if(isFirstNode & returnLinearPred) linearPredVals <<- values(model, linearPred)

      # simulate samples for ecological and process models
      # It's possible that the ecological and observation params are from different distributions
      # so we need to simulate each seperately

      for(i in 1:m){

        # For now, simulate beta's from proposal distribution
        if(betaProposal == "normal"){
          betaVals <<- rmnorm_chol(1, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE)
        } else if(betaProposal == "studentT"){
          betaVals <<- rmvt_chol(1, mu = meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE)
        } else if(betaProposal == "prior"){
          model$simulate(nodes = betaNames)
          betaVals <<- values(model, betaNames)
        }

        betaEsts[i, 1:nBetaSims] <<- betaVals

        # Now simulate the observation Parameters
        if(obsParamsProposal == "normal"){
          obsParamsVals <<- rmnorm_chol(1, meanBetaVar2[iNode,], chol(sigmaBetaVar2[,,iNode]), prec_param = FALSE)
        } else if(obsParamsProposal == "studentT"){
          obsParamsVals <<- rmvt_chol(1, mu = meanBetaVar2[iNode,], chol(sigmaBetaVar2[,,iNode]), df= dfTdist,prec_param = FALSE)
        } else if(obsParamsProposal == "prior"){
          model$simulate(nodes = obsParamsNames)
          obsParamsVals <<- values(model, obsParamsNames)
        }

        obsParamsEsts[i, 1: nObsParamsSims] <<- obsParamsVals

      }
      # If we need linear Predictor to simulate latent nodes
      #then we fit the model z_{i-1}|beta_sim
      if(returnLinearPred)  res <- nimbleINLA(x, y, beta= betaEsts, extraVars = discTarEsts, fixedVals,  family = fam, nCores = nCores)

      #simulate discrete rvs
      for(i in 1:m){
        values(model, betaNames) <<- betaEsts[i, 1:nBetaSims]
        values(model, obsParamsNames) <<- obsParamsEsts[i, 1: nObsParamsSims]

        #for independence in latent variable and alternative MCMC sampler

        if(!latentIsDependent){
          if(iNode > 1){
            values(model, fixedVals) <<- fixedValsStoreMatrix[i, ]
            # values(model, linearPreds) <<- linearPredMatrix[i, ]
          }
        }

        if(latentIsDependent) model$calculate(depsbetaNames)
        #model$calculate()

        #if(!latentIsDependent) values(model, fixedVals) <<- fixedValsStoreMatrix[i, ]
        if(returnLinearPred){
          linPredIndLower <- (lenLinPred * (i -1)) +1
          linPredIndUpper <- (lenLinPred * i)
          linearPredVals <<- res[linPredIndLower: linPredIndUpper, N+2]
          values(model, linearPred) <<- linearPredVals
        }
        # model$calculate(depsDiscTarNames)
        model$calculate()

        #print(isLatentBinary)
        if(isLatentBinary == FALSE){
          model$simulate(nodes = discTarNames)
        } else {
          model$simulate(nodes = binNodesToSimulate)
          values(model, binNodesToFix) <<- binNodesToFixVals
        }


        discTarEsts[i, 1:nDiscTarNames] <<- values(model, discTarNames)
        #print(discTarEsts[i, 1:nDiscTarNames])
      }

      if(!returnLinearPred) res <- nimbleINLA(x, y, beta= betaEsts, extraVars = discTarEsts,fixedVals,  family = fam, nCores = nCores)


      # Calulate weights
      for(i in 1:m){

        # retrieve valaues of ecological parameters
        betaVals <<- betaEsts[i, 1:nBetaSims]
        values(model, beta) <<- betaVals

        # retrieve values for observation process params
        obsParamsVals <<- obsParamsEsts[i, 1: nObsParamsSims]
        values(model, obsParamsNames) <<- obsParamsVals

        # retrieve values of latent state
        values(model, discTarNames) <<- discTarEsts[i, 1:nDiscTarNames]

        #update dependencies of all parameters
        model$calculate(allISparametersDeps)

        # call the calculate of dependencies of all parameters to update log probability of all nodes
        logProbAllPars <-  model$calculate(allISparameters)


        #####################
        # Update past importance weights
        ##########################
        nn <- 0
        nnObsParams <- 0

        # since each observation process parameters can have different proposals
        # we need to estimate their weights seperately
        for(j in 1:iNode){

          # estimating gamma components for ecological process
          if(betaProposal == "normal"){
            #note that getLogProb returns the log of the log density. To get the density, we take exp of the logProb
            nn <- nn +  timeIndex * (dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE)) #+ exp(model$getLogProb(discTarNames)))
          }
          if(betaProposal == "studentT"){
            nn <- nn +  timeIndex * (dmvt_chol(betaVals, mu = meanBeta[j,], chol(sigmaBeta[,,j]), df= dfTdist, prec_param = FALSE, log = FALSE)) #+ exp(model$getLogProb(discTarNames)))
          }

          # estimating gamma components for observation process
          if(obsParamsProposal == "normal"){
            #note that getLogProb returns the log of the log density. To get the density, we take exp of the logProb
            nnObsParams <- nnObsParams +  timeIndex * (dmnorm_chol(obsParamsVals, meanBetaVar2[j,], chol(sigmaBetaVar2[,,j]), prec_param = FALSE, log = FALSE)) #+ exp(model$getLogProb(discTarNames)))
          } else if(obsParamsProposal == "studentT"){
            nnObsParams <- nnObsParams +  timeIndex * (dmvt_chol(obsParamsVals, mu = meanBetaVar2[j,], chol(sigmaBetaVar2[,,j]), df= dfTdist, prec_param = FALSE, log = FALSE)) #+ exp(model$getLogProb(discTarNames)))
          }

          #nugs[j] <<- m * dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE)
        }

        gamma[i] <<- nn
        gammaObsParams[i] <<- nnObsParams
        #print(nnObsParams)
        #gamma[i] <<- sum(nugs[1:iNode])

        ###############################
        # Estimating weights
        ################################

        # calculate numerator of wts in paper

        #values(model, discTarNames) <<- discTarEsts[i, 1:nDiscTarNames]

        # log likelihood should include the contribution from beta and latent variables

        if(returnLinearPred){
          #N = length of fixed Vals, and the fist colum returns the marginal likelihood
          mld <- res[linPredIndLower, 1]
        } else {
          mld <- res[i, 1]
        }





        # log likelihood should include the contribution from beta and latent variables
        # ii <- (i-1)*nSites + 1

        if(latentIsDependent){
          loglike <- mld + model$calculate(beta)
        }else{
          loglike <- model$calculate(c(dataVar, beta))
        }
        #linearPredMatrix[i, ] <- res[,2]

        #calculate weights
        if(iNode > 1){
          wts[i] <<- loglike - log(gamma[i]/sumNt)
          wtsObsParams[i] <<- model$calculate(obsParamsNames) - log(gammaObsParams[i]/sumNt)
        }else{
          # weights of ecological parameters
          if(betaProposal == "normal"){
            wts[i] <<- loglike - dmnorm_chol(betaVals, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = TRUE) #- model$getLogProb(discTarNames)
          } else if(betaProposal == "studentT"){
            wts[i] <<- loglike - dmvt_chol(betaVals, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = TRUE)# -  model$getLogProb(discTarNames)
          } else if(betaProposal == "prior"){
            wts[i] <<- loglike - model$calculate(beta)
          }

          # weights of observation process parameters
          logLikeObsParams <-  model$calculate()
          if(obsParamsProposal == "normal"){
            wtsObsParams[i] <<- logLikeObsParams - dmnorm_chol(obsParamsVals, meanBetaVar2[iNode,], chol(sigmaBetaVar2[,,iNode]), prec_param = FALSE, log = TRUE) #- model$getLogProb(discTarNames)
          } else if(obsParamsProposal == "studentT"){
            wtsObsParams[i] <<- logLikeObsParams - dmvt_chol(obsParamsVals, meanBetaVar2[iNode,], chol(sigmaBetaVar2[,,iNode]), df= dfTdist,prec_param = FALSE, log = TRUE)# -  model$getLogProb(discTarNames)
          } else if(obsParamsProposal == "prior"){
            wtsObsParams[i] <<- logLikeObsParams - model$calculate(obsParamsNames)
          }

        }

        # Calculate weights of latent
        # For latent variable, we estime weights at next level as
        #wts[i] <- p(yt +1|N)*wts[i-1]
        if(iNode > 1){
          wtsLatent[i] <<- model$calculate(discTarNames) #+ mvEWSamples["wtsLatent", i][iNode-1]
        }else{
          wtsLatent[i] <<- model$calculate(discTarNames)
        }
        # save the numerator of the weights for the updating step
        mvWSamples["logLike",i][iNode] <<- loglike
        mvEWSamples["logLike",i][iNode] <<- loglike

        # save the numerator of the weights for the updating step
        mvWSamples["logLikeObsParams",i][iNode] <<- logLikeObsParams
        mvEWSamples["logLikeObsParams",i][iNode] <<- logLikeObsParams

        # save INLA parameter samples for current time
        if(returnLinearPred){
          #N = length of fixed Vals, and the fist colum returns the marginal likelihood
          linPredIndLower <- (lenLinPred * (i -1)) +1
          saveResults(fixedVals, res, ind = linPredIndLower)
        } else {
          saveResults(fixedVals, res, ind = i)
        }


        if(!latentIsDependent){
          fixedValsStoreMatrix[i, ] <<- values(model, fixedVals)
        }

        #values(model, betaWts) <<- wts[i]
        #values(model, latentWts) <<- wtsLatent[i]

        if(isNullAdditionalPars){
          nimCopy(model, mvEWSamples, nodes = beta, nodesTo = beta, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = obsParams, nodesTo = obsParams, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = discTarNames, nodesTo = discTarNames, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = latentWts, nodesTo = latentWts, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = betaWts, nodesTo = betaWts, row = 1, rowTo = k)
        }else{
          nimCopy(model, mvEWSamples, nodes = beta, nodesTo = beta, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = obsParams, nodesTo = obsParams, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = discTarNames, nodesTo = discTarNames, row = 1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = additionalPars, nodesTo = additionalPars, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = latentWts, nodesTo = latentWts, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = betaWts, nodesTo = betaWts, row = 1, rowTo = k)
        }

        #save wts and gamma at Nt, since it won't be updated
        mvEWSamples["gamma",i][iNode] <<- gamma[i]
        mvEWSamples["wts", i][iNode] <<- wts[i]
        mvEWSamples["gammaObsParams",i][iNode] <<- gammaObsParams[i]
        mvEWSamples["wtsObsParams", i][iNode] <<- wtsObsParams[i]
        mvEWSamples["wtsLatent", i][iNode] <<- wtsLatent[i]

        #Saving output for updated mean and sd
        betaEstsUpd[k,1:nBetaSims] <<- betaVals
        betaEstsUpd[k,(nBetaSims+1)] <<- wts[i]
        obsParamsEstsUpd[k,1: nObsParamsSims] <<- obsParamsVals
        obsParamsEstsUpd[k,1: (nObsParamsSims + 1)] <<- wtsObsParams[i]

        print(k)

        #wtsLatVal <- wtsLatent[i]
        mvEWSamples["latentWeights", k][1] <<- wtsLatent[i]  #mvEWSamples["wtsLatent", i][iNode]
        mvEWSamples["betaWeights", k][1] <<- wts[i]
        mvEWSamples["obsParamsWeights", k][1] <<- wtsObsParams[i]

        k <- k + 1
      }

      #Normalised weights for Effective Sample size
      maxWts <- max(wts)
      nWeights <- exp(wts - maxWts)
      nWeights <- nWeights/sum(nWeights)
      ess <<- 1/sum(nWeights^2)

      ######################################
      #update past importance weights and gamma for iNode > 1
      ###################################
# remember we are going to update the nodes for both the ecological and observation process
      if(iNode > 1 & adaptive == TRUE){
        for(i in 1:m){
          for(t in 1:(iNode-1)){
            #index of weights to update
            indx <- i + (t-1)*m

            nimCopy(mvEWSamples, model, nodes = allISparameters, nodesTo = allISparameters, row = indx, rowTo = 1)
            betaValsNew <- values(model, beta)
            obsParamsValsNew <- values(model, obsParams)

            # Update the model logProb with copied values
            model$calculate(allISparametersDeps)
            allIsparslike <- model$calculate(allISparameters)

            # calculate prior distribution for ecological params
            if(betaProposal == "normal"){
              priorDist <- dmnorm_chol(betaValsNew, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
            } else if(betaProposal == "studentT"){
              priorDist <- dmvt_chol(betaValsNew, mu = meanBeta[iNode, ], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
            } else if(betaProposal == "prior"){
              priorDist <- exp(model$calculate(beta))
            }

            # calculate prior distribution for ecological params
            if(obsParamsProposal == "normal"){
              priorObsParamsDist <- dmnorm_chol(obsParamsValsNew, meanBetaVar2[iNode,], chol(sigmaBetaVar2[,,iNode]), prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
            } else if(obsParamsProposal == "studentT"){
              priorObsParamsDist <- dmvt_chol(obsParamsValsNew, mu = meanBetaVar2[iNode, ], chol(sigmaBetaVar2[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
            } else if(obsParamsProposal == "prior"){
              priorObsParamsDist <- exp(model$calculate(obsParams))
            }

            #estimate gamma update for ecological process parameters
            ret <-  mvEWSamples["gamma",i][t] #trying to trick nimble
            gammaUpd <- ret + (m * priorDist) #note that Nt = m
            #print(gammaUpd)
            mvEWSamples["gamma",i][t] <<-  gammaUpd
            mvEWSamples["wts",i][t] <<- mvEWSamples["logLike",i][t] - log(gammaUpd/sumNt)
            betaEstsUpd[indx, 1:nBetaSims] <<- betaValsNew
            betaEstsUpd[indx,(nBetaSims+1)] <<- mvEWSamples["wts",i][t]

            #estimate gamma update for observation process parameters
            retObsParams <-  mvEWSamples["gammaObsParams",i][t] #trying to trick nimble
            gammaObsParamsUpd <- retObsParams + (m * priorObsParamsDist) #note that Nt = m
            mvEWSamples["gammaObsParams",i][t] <<-  gammaObsParamsUpd
            mvEWSamples["wtsObsParams",i][t] <<- mvEWSamples["logLikeObsParams",i][t] - log(gammaObsParamsUpd/sumNt)
            obsParamsEstsUpd[indx, 1:nObsParamsSims] <<-obsParamsValsNew
            obsParamsEstsUpd[indx,(nObsParamsSims+1)] <<- mvEWSamples["wtsObsParams",i][t]

          }
        }

      }

      ##########################################
      # calculate updated mean and covariance matrix for prior distribution
      #########################################

      # remember we are going to update the nodes for both the ecological and observation process
      if(adaptive == TRUE){
        # ecological process
        wtsUpd <- betaEstsUpd[1:sumNt,(nBetaSims+1)]
        maxWtsUpd <- max(wts)
        nWeightsUpd <- exp(wtsUpd - maxWtsUpd)
        betaWts <<- nWeightsUpd

        # observation process
        wtsObsParamsUpd <- obsParamsEstsUpd[1:sumNt,(nObsParamsSims+1)]  #betaEstsUpd[1:sumNt,(nBetaSims+1)]
        maxWtsObsParamsUpd <- max(wtsObsParamsUpd)
        nWeightsObsPramsUpd <- exp(wtsObsParamsUpd - maxWtsObsParamsUpd)
        obsParamsWts <<- nWeightsObsPramsUpd

        # Updating Ecological parameters mean
        for(j in 1:nBetaSims){
          for(i in 1:sumNt){
            muEsts[i, j] <<- nWeightsUpd[i] * betaEstsUpd[i, j]
          }
          mu[j] <<- sum(muEsts[1:sumNt, j])/sum(nWeightsUpd[1:sumNt])
        }

        # Updating Ecological parameters covariance matrix
        for(i in 1:nBetaSims){
          indx <- i
          for(j in indx:nBetaSims){
            sigma[i,j] <<-  sum(nWeightsUpd[1:sumNt]*(betaEstsUpd[1:sumNt, i] - mu[i]) * (betaEstsUpd[1:sumNt, j] - mu[j]))/sum(nWeightsUpd[1:sumNt])
            sigma[j,i] <<-  sum(nWeightsUpd[1:sumNt]*(betaEstsUpd[1:sumNt, i] - mu[i]) * (betaEstsUpd[1:sumNt, j] - mu[j]))/sum(nWeightsUpd[1:sumNt])
          }
        }

        # Updating observation process parameters mean
        for(j in 1:nObsParamsSims){
          for(i in 1:sumNt){
            muObsParamsEsts[i, j] <<- nWeightsObsPramsUpd[i] * obsParamsEstsUpd[i, j]
          }
          muObsParams[j] <<- sum(muObsParamsEsts[1:sumNt, j])/sum(nWeightsObsPramsUpd[1:sumNt])
        }

        # Updating observation parameters covariance matrix
        for(i in 1:nObsParamsSims){
          indx <- i
          for(j in indx:nObsParamsSims){
            sigmaObsParams[i,j] <<-  sum(nWeightsObsPramsUpd[1:sumNt]*(obsParamsEstsUpd[1:sumNt, i] - muObsParams[i]) * (obsParamsEstsUpd[1:sumNt, j] - muObsParams[j]))/sum(nWeightsObsPramsUpd[1:sumNt])
            sigmaObsParams[j,i] <<-  sum(nWeightsObsPramsUpd[1:sumNt]*(obsParamsEstsUpd[1:sumNt, i] - muObsParams[i]) * (obsParamsEstsUpd[1:sumNt, j] - muObsParams[j]))/sum(nWeightsObsPramsUpd[1:sumNt])
          }
        }
      }else{
        mu <<- meanBeta[iNode, ]
        sigma <<- sigmaBeta[ , , iNode]
        muObsParams <<- meanBetaVar2[iNode, ]
        sigmaObsParams <<- sigmaBetaVar2[ , , iNode]
        betaWts <<- nWeights
      }



      print(mu)
      print(sigma)
      print(muObsParams)
      print(sigmaObsParams)
}

    lll <- res[1,1]
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
    },
    updateBetaMeanVar2=function(){
      returnType(double(1))
      return(muObsParams)
    },
    updateBetaSigmaVar2=function(){
      returnType(double(2))
      return(sigmaObsParams)
    },
    returnBetaWts = function(){
      returnType(double(1))
      return(betaWts)
    },
    returnLatentWts = function(){
      returnType(double(1))
      return(wtsLatent)
    }
  )
)



#

inlaISmultiple <- nimbleFunction(
  name = 'inlaISmultiple',
  setup = function(model,
                   fam, x,
                   y,
                   target, #beta
                   control = list()) {

    #control list extraction
    nimbleINLA <- extractControlElement(control, 'nimbleINLA',  NULL)
    fixedVals <- extractControlElement(control, 'fixedVals',  double())
    proposal <- extractControlElement(control, 'proposal',  character())
    initMean <- extractControlElement(control, 'initMean',  NULL)
    initCov <- extractControlElement(control, 'initCov',  NULL)
    initModel <- extractControlElement(control, 'initModel',  FALSE)
    timeIndex <- extractControlElement(control, 'timeIndex',  double()) #Nt
    nSteps <- extractControlElement(control, 'nSteps',  integer()) #Number of steps at each direction
    nCores <- extractControlElement(control, 'nCores',  NULL)
    dfTdist <- extractControlElement(control, 'dfTdist',  NULL)
    adaptive <- extractControlElement(control, 'adaptive',  TRUE)
    additionalPars <- extractControlElement(control, 'additionalPars',  NULL)
    latentIsDependent <- extractControlElement(control, 'latentIsDependent',  TRUE)

    betaWts <- extractControlElement(control, 'betaWts',  NULL)
    latentWts <- extractControlElement(control, 'latentWts',  NULL)
    linearPred <- extractControlElement(control, 'linearPred',  NULL)
    returnLinearPred <- extractControlElement(control, 'returnLinearPred',  FALSE)


    #latentIsDependent must be logical
    if(!is.logical(latentIsDependent)) stop("latentIsDependent must be either TRUE or FALSE, indicating whether the latent variable is dependent on the other MCMC parameters.")
    if(!is.logical(adaptive)) stop("adaptive must be logical: either TRUE or FALSE, indicating whether the parameters of the proposal distribution are adapted or not.")

    #returnLinear must be logical
    if(!is.logical(returnLinearPred)) stop("return Linear predictor must either be TRUE or FALSE")
    # if returnLinear is TRUE, linear pred must not be false
    if(isTRUE(returnLinearPred) & is.null(linearPred)) stop("return Linear predictor cannot be TRUE and the Linear predictor variables be NULL")


     #Note
    # In the original paper:
    ## timeIndex = Nt
    ## nSteps = t
    isDic <- sapply(target, function(x){
      all(model$isDiscrete(nodes = x))==TRUE
    })
    nTarget <- length(target)

    dataVar <- y

    if(!is.character(y)) stop("'y' must be the node nae for the data variable")
    y <- model[[y]]




    if(nTarget < 2) stop("Function only works for more than two target variables")
    discreteTarget <- target[isDic]
    isNotDiscreteTarget <- model$expandNodeNames(nodes = target[!isDic], returnScalarComponents = TRUE)




   #Parameterisations for continuous vars

    if(is.null(nCores)) nCores <- 1
    if(is.null(dfTdist)) dfTdist <- 1

    if(!proposal %in% c("normal", "studentT", "prior")) stop("Proposal distribution must be either 'normal', 'student T' and 'prior' distributions.")


    if(proposal == "prior") adaptive <- FALSE
    #yExpand <- model$expandNodeNames(y, returnScalarComponents = TRUE)
    #y <- model[[y]]
    #y <- c(nimble::values(model, yExpand))
    my_initializeModel <- initializeModel(model, silent = TRUE)
    #save posterior samples
    #modelVals = modelValues(model, m = 1)
    if(is.null(additionalPars)){
    vars <- model$getVarNames(nodes = c(fixedVals, target))
    }else{
      vars <- model$getVarNames(nodes = c(fixedVals, target, additionalPars))
    }
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
    if("gamma" %in% names) stop("change the variable name of gamma.")
    # Add names and dimensions for wts and gamma
    names <- c(names,
               "wts", "wtsObsParams",
               "gamma", "gammaObsParams",
               "logLike", "logLikeObsParams",
               "wtsLatent",
               "betaWeights", "latentWeights", "obsParamsWeights")
    type <- c(type, rep("double", 10))
    size$wts <- nSteps
    size$wtsObsParams <- nSteps
    size$wtsLatent <- nSteps
    size$gamma <- nSteps
    size$gammaObsParams <- nSteps
    #size$gamma2 <- nSteps
    size$logLike <- nSteps
    size$logLikeObsParams <- nSteps
    size$betaWeights <- 1
    size$latentWeights <- 1
    size$obsParamsWeights <- 1
    #size$logLike2 <- nSteps
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
    # check if any of the continuous variables is needed to simulate the latent state
    # This is used for the cased when we want to simulate different parameters seperately
    # Here vars1 are for the parameters needed to simulate latent states
    # Vars2 are for the paramters affecting the observation process model

    contVars <- sapply(isNotDiscreteTarget, function(x){
      deps <- model$getVarNames(nodes = model$getDependencies(nodes = x, stochOnly = TRUE))
      any(deps %in% discreteTarget)
    }

    )

    if(!all(contVars == TRUE)){
      # all all the continuous variables are not used to simulate,
      # so we seperate them
      # Remember Var1 is for the ecological process parameters
      # Var2 is for the observation process

      isNotDiscreteTargetVars1 <- isNotDiscreteTarget[contVars]
      isNotDiscreteTargetVars2 <- isNotDiscreteTarget[!contVars]

      nVars1 <- length(isNotDiscreteTargetVars1)
      nVars2 <- length(isNotDiscreteTargetVars2 )

      # Set initial mean and covariance
      initMean <- rep(0, nVars1)
      initCov <- diag(nVars1)
      initMeanVar2 <- rep(0, nVars2)
      initCovVar2 <- diag(nVars2)

      #create matrix to store updated mean and covariance matrix
      muStep <- matrix(0, nrow = nSteps+1, ncol = nVars1)
      sigmaStep <- array(0, dim = c(nVars1, nVars1, nSteps+1))
      muStepVar2 <- matrix(0, nrow = nSteps+1, ncol = nVars2)
      sigmaStepVar2 <- array(0, dim = c(nVars2, nVars2, nSteps+1))

      # initialize the storage matrix
      muStep[1,] <- initMean
      sigmaStep[,,1] <- initCov
      muStepVar2[1,] <- initMeanVar2
      sigmaStepVar2[,,1] <- initCovVar2
    } else {
      isNotDiscreteTargetVars1 <- isNotDiscreteTarget[contVars]
      isNotDiscreteTargetVars2 <- NULL
      nTargetCont <- length(isNotDiscreteTargetVars1)

      if(is.null(initMean)) initMean <- rep(0, nTargetCont)
      if(is.null(initCov)) initCov <- diag(nTargetCont)
      initMeanVar2 <- rep(0, nTargetCont)
      initCovVar2 <- diag(nTargetCont)

      muStep <- matrix(0, nrow = nSteps+1, ncol = nTargetCont)
      sigmaStep <- array(0, dim = c(nTargetCont, nTargetCont, nSteps+1))
      muStepVar2 <- matrix(0, nrow = nSteps+1, ncol = nTargetCont)
      sigmaStepVar2 <- array(0, dim = c(nTargetCont, nTargetCont, nSteps+1))

      muStep[1,] <- initMean
      sigmaStep[,,1] <- initCov
      muStepVar2[1,] <- initMeanVar2
      sigmaStepVar2[,,1] <- initCovVar2
    }



    impSampINLAstepFnx <- nimbleFunctionList(importanceSamplingStepVirtualMultiple)
    #for(iNode in seq_along(nodes)){
   # beta <- target
    for(iNode in 1:nSteps){
      impSampINLAstepFnx[[iNode]] <- impSampINLAstepMultiple(model,
                                                     mvEWSamples,
                                                     mvWSamples,
                                                     fixedVals,
                                                     iNode,
                                                     x, #covariates
                                                     y, #response variable
                                                     # interInModel,
                                                     fam,
                                                     proposal, #proposal distribution
                                                     isNotDiscreteTargetVars1, # continuous variables which we specify normal or T-distribution for
                                                     isNotDiscreteTargetVars2,
                                                     discreteTarget, # Discrete random variables
                                                     timeIndex,
                                                     vars,
                                                     nCores,
                                                     dfTdist,
                                                     adaptive,
                                                     additionalPars,
                                                     dataVar,
                                                     latentIsDependent,
                                                     returnLinearPred,
                                                     linearPred
      )
    }



    m <- timeIndex


    # save the weights of beta and latent states
    wtsMatrix <- matrix(0, ncol = 2, nrow = nSteps*timeIndex)
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
    #indxs <- 1:timeIndex
    # runs for each Nt
    for(iNode in seq_along(impSampINLAstepFnx)){
      essVals <- essVals + impSampINLAstepFnx[[iNode]]$run(meanBeta = muStep, sigmaBeta = sigmaStep, meanBetaVar2 = muStepVar2, sigmaBetaVar2 = sigmaStepVar2, prevSamp = pp)
      #if(iNode < nSteps){
      muStep[iNode+1,] <<-   impSampINLAstepFnx[[iNode]]$updateBetaMean()
      sigmaStep[ , ,iNode+1] <<- impSampINLAstepFnx[[iNode]]$updateBetaSigma()

      muStepVar2[iNode+1,] <<-   impSampINLAstepFnx[[iNode]]$updateBetaMeanVar2()
      sigmaStepVar2[ , ,iNode+1] <<- impSampINLAstepFnx[[iNode]]$updateBetaSigmaVar2()
      #}
      #wtsMatrix[indxs,1] <<- impSampINLAstepFnx[[iNode]]$returnBetaWts()
     # wtsMatrix[indxs,2] <<- impSampINLAstepFnx[[iNode]]$returnLatentWts()

      pp <- 1

      #indxs <<- indxs + (iNode)*timeIndex
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
