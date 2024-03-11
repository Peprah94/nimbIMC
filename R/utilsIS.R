#' @import nimble
#' @import methods

# gammaEstimation <- nimbleFunction(
#   setup = function(iNode,
#                    ecoParamsProposal,
#                    obsParamsProposal,
#                    dfTdist,
#                    timeIndex){},
#   run = function(meanBeta = double(2),
#                  sigmaBeta = double(3),
#                  meanDisc = double(1),
#                  ecoParamsVals = double(1),
#                  obsParamsVals = double(1)){
#
#     nn <- 0
#     rt <- 0
#     for(j in 1:iNode){
#       if(ecoParamsProposal == "binomial"){
#         rt <- dmybinom(ecoParamsVals, size = 1, prob = meanDisc[j], log = FALSE)
#       } else if(ecoParamsProposal == "poisson"){
#         rt <- dmypois(ecoParamsVals, lambda = meanDisc[j], log = FALSE)
#         # print(rt)
#       }
#
#       # estimating gamma components for ecological process
#       if(obsParamsProposal == "normal"){
#         #note that getLogProb returns the log of the log density. To get the density, we take exp of the logProb
#
#         nn <- nn +  timeIndex * (dmnorm_chol(obsParamsVals, meanBeta[j,], chol(sigmaBeta[,,j]), prec_param = FALSE, log = FALSE) * rt)
#       } else if(obsParamsProposal == "studentT"){
#         nn <- nn +  timeIndex * (dmvt_chol(obsParamsVals, mu = meanBeta[j,], chol(sigmaBeta[,,j]), df= dfTdist, prec_param = FALSE, log = FALSE) *rt)
#       }
#
#     }
#
#     returnType(double(0))
#     return(nn)
#   }
# )


adaptFnxVirtual <- nimbleFunctionVirtual(
  name = 'adaptFnxVirtual',
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
                 meanDisc = double(1),
                 ecoParamsEst = double(2),
                 obsParamsEst = double(2)){
    returnType(double(0))
  },
  methods = list(
    updateMeanObsParams=function(){
      returnType(double(1))
      return(mu)
    },
    updateSigmaObsParams=function(){
      returnType(double(2))
      return(sigma)
    },
    updateMeanEcoParams=function(){
      returnType(double(0))
      return(muEcoPars)
    }
  )
)

adaptingFunction <- nimbleFunction(
  contains = adaptFnxVirtual,
  setup = function(model,
                   mvEWSamples,
                   m,
                   iNode,
                   dfTdist,
                   ecoParams,
                   obsParams,
                   ecoParamsProposal,
                   obsParamsProposal,
                   adaptive){

    #store simulated values of ecological parameters
    ecoParams <- model$expandNodeNames(nodes = ecoParams)
    obsParams <- model$expandNodeNames(nodes = obsParams)
    nEcoParams <- length(ecoParams) #dimeension of beta
    nObsParams <- length(obsParams)
    ecoParamsVals <- rep(0, length = nEcoParams) #vector to store beta values
    depsEcoParamsNames <- model$getDependencies(ecoParams, self = FALSE, determOnly = TRUE)
    obsParamsVals <- rep(0, length = nObsParams)
    #get dependencies of both parameters. Will be used to calculate weights at each time
    allISparameters <- c(ecoParams, obsParams)
    allISparametersDeps <-   model$getDependencies(allISparameters)
    # Now we check of the ecoParams is binary
    isEcoParamsBinary <- all(model$isBinary(ecoParams) == TRUE)

    # select those that are data
    binNodesToSimulate <- c(ecoParams[!model$isData(ecoParams)])
    nBinNodesToSimulate <- c(length(binNodesToSimulate))
    binNodesToFix <- ecoParams[model$isData(ecoParams)]
    nBinNodesToFix <- length(binNodesToFix)
    binNodesToSimulateVals <- rep(0, length = nBinNodesToSimulate)
    binNodesToFixVals <- rep(1, length = nBinNodesToFix)

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
    muEcoParsEst <- rep(0, m)

    muEcoPars <- 0.5
    # Get obsParams estimates from steps 1 to Nt
    obsParamsEstUpd <- matrix(0, nrow = sumNt, ncol = (nObsParams+1))
    ecoParamsEstUpd <- matrix(0, nrow = sumNt, ncol = (nEcoParams+1))
    rts <- matrix(0, nrow = sumNt, ncol = nObsParams)
    rtsUpd <- matrix(0, nrow = sumNt, ncol = nObsParams)
    muEsts <- matrix(0, nrow = sumNt, ncol = nObsParams)


  },
  run = function(meanBeta = double(2),
                 sigmaBeta = double(3),
                 meanDisc = double(1),
                 ecoParamsEst = double(2),
                 obsParamsEst = double(2)){
    returnType(double(0))

    declare(priorDist, double(0))
    declare(ecoParamsValsNew, double(1))
    declare(obsParamsValsNew, double(1))
    declare(allIsparslike, double(0))
    declare(priorObsParamsDist, double(0))
    declare(gammaUpd, double(0))
    declare(gammaObsParamsUpd, double(0))
    declare(wtsUpd, double(1))
    declare(maxWtsUpd, double(0))
    priorDist <- 0
    priorObsParamsDist <- 0
    maxWtsUpd <- 0

lll <- 0

ecoParamsEstUpd <<- ecoParamsEst
obsParamsEstUpd <<- obsParamsEst

    if(iNode > 1 & adaptive == TRUE){
      for(i in 1:m){
        for(t in 1:(iNode-1)){
          #index of weights to update
          indx <- i + (t-1)*m

          nimCopy(mvEWSamples, model, nodes = allISparameters, nodesTo = allISparameters, row = indx, rowTo = 1)
          if(isEcoParamsBinary){
            ecoParamsValsNew <- values(model, binNodesToSimulate)
          } else {
            ecoParamsValsNew <- values(model, ecoParams)
          }
          obsParamsValsNew <- values(model, obsParams)

          # Update the model logProb with copied values
          model$calculate(allISparametersDeps)
          allIsparslike <- model$calculate(allISparameters)

          # calculate prior distribution for ecological params
          # if(ecoParamsProposal == "normal"){
          #   priorDist <- dmnorm_chol(betaValsNew, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
          # } else if(ecoParamsProposal == "studentT"){
          #   priorDist <- dmvt_chol(betaValsNew, mu = meanBeta[iNode, ], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
          # } else
          if(ecoParamsProposal == "prior"){
            priorDist <- exp(model$calculate(ecoParams))
          } else if (ecoParamsProposal == "binomial"){
            priorDist <-  dmybinom(ecoParamsValsNew[1:nBinNodesToSimulate], size = 1, prob = meanDisc[iNode], log = FALSE)
          } else if (ecoParamsProposal == "poisson"){
            priorDist <-  dmypois(ecoParamsValsNew, lambda = meanDisc[iNode], log = FALSE)
          }

          # calculate prior distribution for ecological params
          if(obsParamsProposal == "normal"){
            priorObsParamsDist <- dmnorm_chol(obsParamsValsNew, meanBeta[iNode,], chol(sigmaBeta[,,iNode]), prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
          } else if(obsParamsProposal == "studentT"){
            priorObsParamsDist <- dmvt_chol(obsParamsValsNew, mu = meanBeta[iNode, ], chol(sigmaBeta[,,iNode]), df= dfTdist,prec_param = FALSE, log = FALSE) #+ exp(model$getLogProb(discTarNames))
          } else if(obsParamsProposal == "prior"){
            priorObsParamsDist <- exp(model$calculate(obsParams))
          }

          #estimate gamma update for ecological process parameters
          ret <-  mvEWSamples["gamma",i][t] #trying to trick nimble
          gammaUpd <- ret + (m * (priorDist * priorObsParamsDist))#note that Nt = m
          #print(gammaUpd)
          mvEWSamples["gamma",i][t] <<-  gammaUpd
          mvEWSamples["wts",i][t] <<- mvEWSamples["logLike",i][t] - log(gammaUpd/sumNt)
          if(isEcoParamsBinary){
            ecoParamsEstUpd[indx, 1:nBinNodesToSimulate] <<- ecoParamsValsNew[1:nBinNodesToSimulate]
          } else {
            ecoParamsEstUpd[indx, 1:nEcoParams] <<- ecoParamsValsNew[1:nEcoParams]
          }
          ecoParamsEstUpd[indx,(nEcoParams+1)] <<- mvEWSamples["wts",i][t]

          #estimate gamma update for observation process parameters
          retObsParams <-  mvEWSamples["gammaObsParams",i][t] #trying to trick nimble
          gammaObsParamsUpd <- retObsParams + (m * priorObsParamsDist) #note that Nt = m
          mvEWSamples["gammaObsParams",i][t] <<-  gammaObsParamsUpd
          mvEWSamples["wtsObsParams",i][t] <<- mvEWSamples["logLikeObsParams",i][t] - log(gammaObsParamsUpd/sumNt)
          obsParamsEstUpd[indx, 1:nObsParams] <<- obsParamsValsNew
          obsParamsEstUpd[indx,(nObsParams+1)] <<- mvEWSamples["wts",i][t]

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
      maxWtsUpd <- max(wtsUpd)
      nWeightsUpd <- exp(wtsUpd - maxWtsUpd)



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

        if(isEcoParamsBinary){
          muEcoParsEst[i] <<- (sum(ecoParamsEstUpd[i, 1:nBinNodesToSimulate])/nEcoParams)* nWeightsUpd[i]
        } else {
          muEcoParsEst[i] <<- (sum(ecoParamsEstUpd[i, 1:nEcoParams])/nEcoParams)* nWeightsUpd[i]
        }
      }
      muEcoPars <<- sum(muEcoParsEst[1:sumNt])/(sum(nWeightsUpd[1:sumNt]))
    }else{
      mu <<- meanBeta[iNode, ]
      sigma <<- sigmaBeta[ , , iNode]
      #betaWts <<- nWeights
      muEcoPars <<- meanDisc[iNode]
    }


    return(lll)
  },
  methods = list(
    updateMeanObsParams=function(){
      returnType(double(1))
      return(mu)
    },
    updateSigmaObsParams=function(){
      returnType(double(2))
      return(sigma)
    },
    updateMeanEcoParams=function(){
      returnType(double(0))
      return(muEcoPars)
    }
  )
)









