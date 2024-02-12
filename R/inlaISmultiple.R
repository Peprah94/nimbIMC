

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
                 meanDisc = double(1),
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
    updateMuEcoPars=function(){
      returnType(double(1))
      return(muEcoPars)
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
                   obsParams,
                   ecoParams,
                   discreteTarget,
                   timeIndex,
                   vars,
                   nCores,
                   dfTdist,
                   adaptive,
                   additionalPars,
                   dataVar,
                   latentIsDependent
  ) {

    #Note
    # iNode = Nt
    #timeIndex = number of samples at each it = m
    m <- timeIndex

    #check if it is the first index
    isFirstNode <- iNode == 1


    # set names of parameters
    # ecoParams refers to ecological parameters
    # obsParams refers to observation parameters
    #get the distribution of the discrete variable

    ecoParamsProposal <- unique(model$getDistribution(ecoParams))
    obsParamsProposal <- unique(model$getDistribution(obsParams))

    # setting up parameters
    N <- length(fixedVals) #length of INLA parameters
    ess <- 0 # initialise effective sample size

    #store simulated values of ecological parameters
    ecoParams <- model$expandNodeNames(nodes = ecoParams)
    obsParams <- model$expandNodeNames(nodes = obsParams)
    nEcoParams <- length(ecoParams) #dimeension of beta
    nObsParams <- length(obsParams)
    ecoParamsVals <- rep(0, length = nEcoParams) #vector to store beta values
    depsEcoParamsNames <- model$getDependencies(ecoParams, self = FALSE, determOnly = TRUE)
    obsParamsVals <- rep(0, length = nObsParams)

    # Now we check of the ecoParams is binary
    isEcoParamsBinary <- all(model$isBinary(ecoParams) == TRUE)

    # select those that are data
    binNodesToSimulate <- ecoParams[!model$isData(ecoParams)]
    nBinNodesToSimulate <- length(binNodesToSimulate)
    binNodesToFix <- ecoParams[model$isData(ecoParams)]
    nBinNodesToFix <- length(binNodesToFix)
    binNodesToSimulateVals <- rep(0, length = nBinNodesToSimulate)
    binNodesToFixVals <- rep(1, length = nBinNodesToFix)


    #get dependencies of both parameters. Will be used to calculate weights at each time
    allISparameters <- c(ecoParams, obsParams)
    allISparametersDeps <-   model$getDependencies(allISparameters)

    # store updated mean and sigma
    mu <- rep(0, length = nObsParams)
    sigma <- matrix(0, nrow=nObsParams, ncol = nObsParams)

    # store weights and gamma for ecological process parameters
    wts <- numeric(m)
    gamma <- numeric(m)
    ecoParamsWts <- numeric(m)

    # Create a cummulative indexing to calculate nugget
    nugget <- seq(1, iNode, 1) #Cummulative index for N's
    nNugget <- length(nugget) #length of N's
    nugs <- numeric(nNugget)

    # Estimate increasing indices
    indInc <- (iNode - 1)*m #sum from N1 to Nt-1
    sumNt <- m + indInc

    # estimates at iteration t
    obsParamsEst <- matrix(0, nrow = m, ncol = nObsParams)
    ecoParamsEst <- matrix(0, nrow = m, ncol = nEcoParams)


    # Get obsParams estimates from steps 1 to Nt
    obsParamsEstUpd <- matrix(0, nrow = sumNt, ncol = (nObsParams+1))
    ecoParamsEstUpd <- matrix(0, nrow = sumNt, ncol = (nEcoParams+1))
    rts <- matrix(0, nrow = sumNt, ncol = nObsParams)
    rtsUpd <- matrix(0, nrow = sumNt, ncol = nObsParams)
    muEsts <- matrix(0, nrow = sumNt, ncol = nObsParams)

    muEcoPars <- 0.5

    #check if we have additional Parameters (especially derived quantities) to save
    isNullAdditionalPars <- is.null(additionalPars)
    isNullObsParams <- is.null(obsParams)

  fixedValsStoreMatrix <- matrix(0, nrow = timeIndex, ncol = length(fixedVals))



  },
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
                 meanDisc = double(1),
                 prevSamp= integer(0)) {
    returnType(double(0))

    k <- indInc + 1
    m <<- timeIndex
    linPredIndUpper <- 0
    linPredIndLower <- 0



     # if(isFirstNode & returnLinearPred) linearPredVals <<- values(model, linearPred)

      # simulate samples for ecological and process models
      # It's possible that the ecological and observation params are from different distributions
      # so we need to simulate each seperately

      for(i in 1:m){

        # For now, simulate beta's from proposal distribution
        if(ecoParamsProposal == "dnorm"){
          ecoParamsVals <<- rmnorm_chol(1, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE)
        } else if(ecoParamsProposal == "studentT"){
          ecoParamsVals <<- rmvt_chol(1, mu = meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE)
        } else if(ecoParamsProposal == "prior"){
          model$simulate(nodes = ecoParams)
          ecoParamsVals <<- values(model, betaNames)
        } else if(ecoParamsProposal == "dbin"){
          ranEco <- rbinom(nBinNodesToSimulate, prob = meanDisc[iNode], size = 1)
          values(model, binNodesToSimulate) <<- ranEco
          values(model, binNodesToFix) <<- binNodesToFixVals
          ecoParamsVals <<- values(model, ecoParams)
        }

        ecoParamsEst[i, 1:nEcoParams] <<- ecoParamsVals

        # Now simulate the observation Parameters
        if(obsParamsProposal == "normal" | obsParamsProposal == "dnorm"){
          obsParamsVals <<- rmnorm_chol(1, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE)
        } else if(obsParamsProposal == "studentT"){
          obsParamsVals <<- rmvt_chol(1, mu = meanBeta[iNode,], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE)
        } else if(obsParamsProposal == "prior"){
          model$simulate(nodes = obsParams)
          obsParamsVals <<- values(model, obsParams)
        }

        obsParamsEst[i, 1: nObsParams] <<- obsParamsVals

      }
      # If we need linear Predictor to simulate latent nodes
      #then we fit the model z_{i-1}|beta_sim
      res <- nimbleINLA(x, y, beta = ecoParamsEst, extraVars = ecoParamsEst, fixedVals,  family = fam, nCores = nCores)

      #simulate discrete rvs
      for(i in 1:m){
        values(model, ecoParams) <<- ecoParamsEst[i, 1:nEcoParams]
        values(model, obsParams) <<- obsParamsEst[i, 1: nObsParams]

        #for independence in latent variable and alternative MCMC sampler
        # save INLA parameter samples for current time
          saveResults(fixedVals, res, ind = i)

        model$calculate()

        # retrieve valaues of ecological parameters
        ecoParamsVals <<- ecoParamsEst[i, 1:nEcoParams]
        obsParamsVals <<- obsParamsEst[i, 1: nObsParams]


        #####################
        # Update past importance weights
        ##########################
        nn <- 0
       # nnObsParams <- 0

        # since each observation process parameters can have different proposals
        # we need to estimate their weights seperately
        for(j in 1:iNode){

          # estimating gamma components for ecological process
          if(obsParamsProposal == "normal"){
            #note that getLogProb returns the log of the log density. To get the density, we take exp of the logProb
            if(ecoParamsProposal == "dbin"){
              rt <- dbinom(ecoParamsVals, size = 1, prob = meanDisc[j], log = FALSE)
            }
            nn <- nn +  timeIndex * (dmnorm_chol(betaVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE) * rt)
          } else if(obsParamsProposal == "studentT"){
            nn <- nn +  timeIndex * (dmvt_chol(betaVals, mu = meanBeta[j,], chol(sigmaBeta[,,j]), df= dfTdist, prec_param = FALSE, log = FALSE) *rt)
          }

        }

        gamma[i] <<- nn
        #gammaObsParams[i] <<- nnObsParams
        #print(nnObsParams)
        #gamma[i] <<- sum(nugs[1:iNode])

        ###############################
        # Estimating weights
        ################################

        # calculate numerator of wts in paper

        #values(model, discTarNames) <<- discTarEsts[i, 1:nDiscTarNames]

        # log likelihood should include the contribution from beta and latent variables
       mld <- res[i, 1]






        # log likelihood should include the contribution from beta and latent variables
        # ii <- (i-1)*nSites + 1

        #if(latentIsDependent){
          loglike <- mld + model$calculate(c(dataVar, obsParams))
       # }else{
       #   loglike <- model$calculate(c(dataVar, beta))
        #}
        #linearPredMatrix[i, ] <- res[,2]

        #calculate weights
        #if(iNode > 1){
          wts[i] <<- loglike - log(gamma[i]/sumNt)

        # save the numerator of the weights for the updating step
        mvWSamples["logLike",i][iNode] <<- loglike
        mvEWSamples["logLike",i][iNode] <<- loglike

        # save the numerator of the weights for the updating step
        mvWSamples["logLikeObsParams",i][iNode] <<- loglike #logLikeObsParams
        mvEWSamples["logLikeObsParams",i][iNode] <<- loglike #logLikeObsParams



        #values(model, betaWts) <<- wts[i]
        #values(model, latentWts) <<- wtsLatent[i]

        if(isNullAdditionalPars){
          nimCopy(model, mvEWSamples, nodes = ecoParams, nodesTo = beta, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = obsParams, nodesTo = obsParams, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = discTarNames, nodesTo = discTarNames, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = latentWts, nodesTo = latentWts, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = betaWts, nodesTo = betaWts, row = 1, rowTo = k)
        }else{
          nimCopy(model, mvEWSamples, nodes = ecoParams, nodesTo = beta, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = obsParams, nodesTo = obsParams, row=1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = fixedVals, nodesTo = fixedVals, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = discTarNames, nodesTo = discTarNames, row = 1, rowTo = k)
          nimCopy(model, mvEWSamples, nodes = additionalPars, nodesTo = additionalPars, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = latentWts, nodesTo = latentWts, row = 1, rowTo = k)
          #nimCopy(model, mvEWSamples, nodes = betaWts, nodesTo = betaWts, row = 1, rowTo = k)
        }

        #save wts and gamma at Nt, since it won't be updated
        mvEWSamples["gamma",i][iNode] <<- gamma[i]
        mvEWSamples["wts", i][iNode] <<- wts[i]
        mvEWSamples["gammaObsParams",i][iNode] <<- gammaObsParams[i]
        mvEWSamples["wtsObsParams", i][iNode] <<- wts[i]
        mvEWSamples["wtsLatent", i][iNode] <<- wts[i]

        #Saving output for updated mean and sd
        ecoParamsEstUpd[k,1:nEcoParams] <<- ecoParamsVals
        ecoParamsEstUpd[k,(nEcoParams+1)] <<- wts[i]
        obsParamsEstUpd[k,1:nObsParams] <<- obsParamsVals
        obsParamsEstUpd[k, (nObsParams + 1)] <<- wts[i]

       # print(k)

        #wtsLatVal <- wtsLatent[i]
        mvEWSamples["latentWeights", k][1] <<- wts[i]  #mvEWSamples["wtsLatent", i][iNode]
        mvEWSamples["betaWeights", k][1] <<- wts[i]
        mvEWSamples["obsParamsWeights", k][1] <<- wts[i]

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
            ecoParamsValsNew <- values(model, ecoParams)
            obsParamsValsNew <- values(model, obsParams)

            # Update the model logProb with copied values
            model$calculate(allISparametersDeps)
            allIsparslike <- model$calculate(allISparameters)

            # calculate prior distribution for ecological params
            if(ecoParamsProposal == "normal"){
              priorDist <- dmnorm_chol(betaValsNew, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
            } else if(ecoParamsProposal == "studentT"){
              priorDist <- dmvt_chol(betaValsNew, mu = meanBeta[iNode, ], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
            } else if(ecoParamsProposal == "prior"){
              priorDist <- exp(model$calculate(beta))
            } else if (ecoParamsProposal == "dbin"){
              priorDist <-  dbinom(ecoParamsVals, size = 1, prob = meanDisc[iNode], log = FALSE)
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
            gammaUpd <- ret + (m * (priorDist * priorObsParamsDist))#note that Nt = m
            #print(gammaUpd)
            mvEWSamples["gamma",i][t] <<-  gammaUpd
            mvEWSamples["wts",i][t] <<- mvEWSamples["logLike",i][t] - log(gammaUpd/sumNt)
            ecoParamsEstUpd[indx, 1:nEcoParams] <<- ecoParamsVals
            ecoParamsEstUpd[indx,(nEcoParams+1)] <<- mvEWSamples["wts",i][t]

            #estimate gamma update for observation process parameters
            retObsParams <-  mvEWSamples["gammaObsParams",i][t] #trying to trick nimble
            gammaObsParamsUpd <- retObsParams + (m * priorObsParamsDist) #note that Nt = m
            mvEWSamples["gammaObsParams",i][t] <<-  gammaObsParamsUpd
            mvEWSamples["wtsObsParams",i][t] <<- mvEWSamples["logLikeObsParams",i][t] - log(gammaObsParamsUpd/sumNt)
            obsParamsEstUpd[indx, 1:nObsParamsSims] <<- obsParamsValsNew
           obsParamsEstUpd[indx,(nObsParamsSims+1)] <<- mvEWSamples["wts",i][t]

          }
        }

      }

      ##########################################
      # calculate updated mean and covariance matrix for prior distribution
      #########################################

      # remember we are going to update the nodes for both the ecological and observation process
      if(adaptive == TRUE){
        # observation process process
        wtsUpd <- obsParamsEstUpd[1:sumNt,(nObsParams+1)]
        maxWtsUpd <- max(wts)
        nWeightsUpd <- exp(wtsUpd - maxWtsUpd)
        #betaWts <<- nWeightsUpd


        # Updating observ parameters mean
        for(j in 1:nObsParams){
          for(i in 1:sumNt){
            muEsts[i, j] <<- nWeightsUpd[i] * obsParamsEstUpd[i, j]
          }
          mu[j] <<- sum(muEsts[1:sumNt, j])/sum(nWeightsUpd[1:sumNt])
        }

        # Updating observ parameters covariance matrix
        for(i in 1:nObsParams){
          indx <- i
          for(j in indx:nObsParams){
            sigma[i,j] <<-  sum(nWeightsUpd[1:sumNt]*(obsParamsEstUpd[1:sumNt, i] - mu[i]) * (obsParamsEstUpd[1:sumNt, j] - mu[j]))/sum(nWeightsUpd[1:sumNt])
            sigma[j,i] <<-  sum(nWeightsUpd[1:sumNt]*(obsParamsEstUpd[1:sumNt, i] - mu[i]) * (obsParamsEstUpd[1:sumNt, j] - mu[j]))/sum(nWeightsUpd[1:sumNt])
          }
        }

        # updating ecological parameters
        for(i in 1:sumNt){
          muEcoParsEst[i] <- sum(ecoParamsEstUpd[1:nObsParams, i] * nWeightsUpd[i])
        }
        muEcoPars <<- sum(muEcoParsEst[1:sumNt])
      }else{
        mu <<- meanBeta[iNode, ]
        sigma <<- sigmaBeta[ , , iNode]
        #betaWts <<- nWeights
        muEcoPars <<- meanDisc[iNode]
      }



      print(mu)
      print(sigma)
      print(muObsParams)
      print(sigmaObsParams)
      print(muEcoPars)


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
    updateMuEcoPars=function(){
      returnType(double(1))
      return(muEcoPars)
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
    initEcoPars <- extractControlElement(control, 'initEcoPars',  NULL)
    initModel <- extractControlElement(control, 'initModel',  FALSE)
    timeIndex <- extractControlElement(control, 'timeIndex',  double()) #Nt
    nSteps <- extractControlElement(control, 'nSteps',  integer()) #Number of steps at each direction
    nCores <- extractControlElement(control, 'nCores',  NULL)
    dfTdist <- extractControlElement(control, 'dfTdist',  NULL)
    adaptive <- extractControlElement(control, 'adaptive',  TRUE)
    additionalPars <- extractControlElement(control, 'additionalPars',  NULL)
    latentIsDependent <- extractControlElement(control, 'latentIsDependent',  FALSE)

    betaWts <- extractControlElement(control, 'betaWts',  NULL)
    latentWts <- extractControlElement(control, 'latentWts',  NULL)
    #linearPred <- extractControlElement(control, 'linearPred',  NULL)
    #returnLinearPred <- extractControlElement(control, 'returnLinearPred',  FALSE)


    #latentIsDependent must be logical
    if(!is.logical(latentIsDependent)) stop("latentIsDependent must be either TRUE or FALSE, indicating whether the latent variable is dependent on the other MCMC parameters.")
    if(!is.logical(adaptive)) stop("adaptive must be logical: either TRUE or FALSE, indicating whether the parameters of the proposal distribution are adapted or not.")

    #returnLinear must be logical
    #if(!is.logical(returnLinearPred)) stop("return Linear predictor must either be TRUE or FALSE")
    # if returnLinear is TRUE, linear pred must not be false
   # if(isTRUE(returnLinearPred) & is.null(linearPred)) stop("return Linear predictor cannot be TRUE and the Linear predictor variables be NULL")


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

    # if there are discrete variables, seperate them from the continuous ones
    discreteTarget <- target[isDic]
    contTarget <- model$expandNodeNames(nodes = target[!isDic], returnScalarComponents = TRUE)


    # I need to check the dependency of these discrete and continuous variables
    # and see which one is ecological (latent) and observation parameters
    #if(all(model$getDependencies(contTarget, stochOnly = TRUE, includeData = FALSE) %in% contTarget)){
      obsParams <- contTarget
      ecoParams <- discreteTarget
   # } else {
    #  ecoParams <- contTarget
     # obsParams <- discreteTarget
   # }


   #Parameterisations for continuous vars

    if(is.null(nCores)) nCores <- 1
    if(is.null(dfTdist)) dfTdist <- 1
      if(is.null(initEcoPars)) initEcoPars <- 0.5

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


    # Set initial mean and covariance
    nContVars <- length(contTarget)
       initMean <- rep(0, nContVars)
       initCov <- diag(nContVars)

    #   #create matrix to store updated mean and covariance matrix
      muStep <- matrix(0, nrow = nSteps+1, ncol = nContVars)
       sigmaStep <- array(0, dim = c(nContVars, nContVars, nSteps+1))
       meanEcoParsStep <- rep(0, nSteps + 1)

       #   # initialize the storage matrix
        muStep[1,] <- initMean
        sigmaStep[,,1] <- initCov
        meanEcoParsStep[1] <- initEcoPars



    # contVars <- sapply(isNotDiscreteTarget, function(x){
    #   deps <- model$getVarNames(nodes = model$getDependencies(nodes = x, stochOnly = TRUE))
    #   any(deps %in% discreteTarget)
    # }
    #
    # )
    #
    # if(!all(contVars == TRUE)){
    #   # all all the continuous variables are not used to simulate,
    #   # so we seperate them
    #   # Remember Var1 is for the ecological process parameters
    #   # Var2 is for the observation process
    #
    #   isNotDiscreteTargetVars1 <- isNotDiscreteTarget[contVars]
    #   isNotDiscreteTargetVars2 <- isNotDiscreteTarget[!contVars]
    #
    #   nVars1 <- length(isNotDiscreteTargetVars1)
    #   nVars2 <- length(isNotDiscreteTargetVars2 )
    #
    #   # Set initial mean and covariance
    #   initMean <- rep(0, nVars1)
    #   initCov <- diag(nVars1)
    #   initMeanVar2 <- rep(0, nVars2)
    #   initCovVar2 <- diag(nVars2)
    #
    #   #create matrix to store updated mean and covariance matrix
    #   muStep <- matrix(0, nrow = nSteps+1, ncol = nVars1)
    #   sigmaStep <- array(0, dim = c(nVars1, nVars1, nSteps+1))
    #   muStepVar2 <- matrix(0, nrow = nSteps+1, ncol = nVars2)
    #   sigmaStepVar2 <- array(0, dim = c(nVars2, nVars2, nSteps+1))
    #
    #   # initialize the storage matrix
    #   muStep[1,] <- initMean
    #   sigmaStep[,,1] <- initCov
    #   muStepVar2[1,] <- initMeanVar2
    #   sigmaStepVar2[,,1] <- initCovVar2
    # } else {
    #   isNotDiscreteTargetVars1 <- isNotDiscreteTarget[contVars]
    #   isNotDiscreteTargetVars2 <- NULL
    #   nTargetCont <- length(isNotDiscreteTargetVars1)
    #
    #   if(is.null(initMean)) initMean <- rep(0, nTargetCont)
    #   if(is.null(initCov)) initCov <- diag(nTargetCont)
    #   initMeanVar2 <- rep(0, nTargetCont)
    #   initCovVar2 <- diag(nTargetCont)
    #
    #   muStep <- matrix(0, nrow = nSteps+1, ncol = nTargetCont)
    #   sigmaStep <- array(0, dim = c(nTargetCont, nTargetCont, nSteps+1))
    #   muStepVar2 <- matrix(0, nrow = nSteps+1, ncol = nTargetCont)
    #   sigmaStepVar2 <- array(0, dim = c(nTargetCont, nTargetCont, nSteps+1))
    #
    #   muStep[1,] <- initMean
    #   sigmaStep[,,1] <- initCov
    #   muStepVar2[1,] <- initMeanVar2
    #   sigmaStepVar2[,,1] <- initCovVar2
    # }

    # Print all the parameters for users:
    print(paste("Multiple parameters used sampled with Importance Sampling:", is.null(discreteTarget)))
    print(paste("Ecological process parameters sampled using Importance Sampling:", discreteTarget))
    print(paste("Observation process parameters sampled using Importance Sampling:", contTarget))
    #message(paste0("Latent state sampled using Importance Sampling:", discreteTarget))
    print(paste("Other parameters sampled from INLA:", fixedVals))


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
                                                     obsParams, # continuous variables which we specify normal or T-distribution for
                                                    ecoParams,
                                                     discreteTarget, # Discrete random variables
                                                     timeIndex,
                                                     vars,
                                                     nCores,
                                                     dfTdist,
                                                     adaptive,
                                                     additionalPars,
                                                     dataVar,
                                                     latentIsDependent
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
      essVals <- essVals + impSampINLAstepFnx[[iNode]]$run(meanBeta = muStep, sigmaBeta = sigmaStep,  meanDisc = meanEcoParsStep, prevSamp = pp)
      #if(iNode < nSteps){
      muStep[iNode+1,] <<-   impSampINLAstepFnx[[iNode]]$updateBetaMean()
      sigmaStep[ , ,iNode+1] <<- impSampINLAstepFnx[[iNode]]$updateBetaSigma()
      meanEcoParsStep[iNode + 1] <<- impSampINLAstepFnx[[iNode]]$updateMuEcoPars()
      #muStepVar2[iNode+1,] <<-   impSampINLAstepFnx[[iNode]]$updateBetaMeanVar2()
      #sigmaStepVar2[ , ,iNode+1] <<- impSampINLAstepFnx[[iNode]]$updateBetaSigmaVar2()
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
    },
    getLastMeanEcoPars = function(){
      return(meanEcoParsStep)
      returnType(double(1))
    }
  )
)
