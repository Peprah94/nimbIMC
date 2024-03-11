

# proposalDistributionVirtual <- nimbleFunctionVirtual(
#   setup = function(model,
#                    binNodesToSimulate,
#                    binNodesToFix,
#                    binNodesToFixVals,
#                    params){},
#   run = function(meanDisc = double(0),
#                  n = integer(0),
#                  size = integer(0)){
#     returnType(double(1))
#   }
# )

proposalDistributionVirtual <- nimbleFunctionVirtual(
  name = 'proposalDistributionVirtual',
  run = function(meanDisc = double(0),
                 n = integer(0),
                 m = integer(0),
                 mean = double(1),
                 sigma = double(2)
                 ) {
    returnType(double(2))
  }
)



# Binomial Proposal distribution
binProposal <- nimbleFunction(
  contains =  proposalDistributionVirtual,
  setup = function(model,
                   binNodesToSimulate,
                   binNodesToFix,
                   binNodesToFixVals,
                   params,
                   size){},
  run = function(meanDisc = double(0),
                 n = integer(0, default = 5),
                 m = integer(0, default = 10),
                 mean = double(1, default = c(0,0)),
                 sigma = double(2, default = diag(2))){
    paramsEst <- matrix(0, nrow = m, ncol = n)
    for(i in 1:m){
    ranEco <- rmybinom(n, prob = meanDisc, size)
    values(model, binNodesToSimulate) <<- ranEco[1:n]
    values(model, binNodesToFix) <<- binNodesToFixVals
    paramsVals <- values(model, params)
    paramsEst[i, 1:n] <- paramsVals[1:n]
    }
    returnType(double(2))
    return(paramsEst)
  }
)


# poisson Proposal distribution
poisProposal <- nimbleFunction(
  contains =  proposalDistributionVirtual,
  setup = function(model,
                   params,
                   lowerBound = double(1),
                   includeLowerBound = integer(0)){},
  run = function(meanDisc = double(0),
                 n = integer(0),
                 m = integer(0, default = 10),
                 mean = double(1, default = c(0,0)),
                 sigma = double(2, default = diag(2))){
    if(n != length(params)) stop("n should be equal to number of ecological parameters")
    paramsEst <- matrix(0, nrow = m, ncol = n)
    for(i in 1:m){
     paramsVals <- rmypois(n, lambda = meanDisc, lowerBound =lowerBound, includeLowerBound = includeLowerBound)
     paramsEst[i, 1:n] <- paramsVals[1:n]
     }
    returnType(double(2))
    return(paramsEst)
  }
)

#prior distribution
priorProposal <- nimbleFunction(
  contains =  proposalDistributionVirtual,
  setup = function(model,
                   params){},
  run = function(meanDisc = double(0),
                 n = integer(0),
                 m = integer(0, m = 10),
                 mean = double(1, default = c(0,0)),
                 sigma = double(2, default = diag(2))){
    paramsEst <- matrix(0, nrow = m, ncol = n)
    for(i in 1:m){
    model$simulate(nodes = params)
    paramsVals <- values(model, params)
    paramsEst[i, 1:n] <- paramsVals[1:n]
    }
    returnType(double(2))
    return(paramsEst)
  }
)


# Normal distribution
normalProposal <- nimbleFunction(
  contains =  proposalDistributionVirtual,
  setup = function(model,
                   params){},
  run = function(meanDisc = double(0, default = 1),
                 n = integer(0, default = 5),
                 m = integer(0, default = 10),
                 mean = double(1, default = c(0,0)),
                 sigma = double(2, default = diag(2))){
    paramsEst <- matrix(0, nrow = m, ncol = n)
    for(i in 1:m){
    paramsVals <- rmnorm_chol(1, mean, chol(sigma), prec_param = FALSE)
    paramsEst[i, 1:n] <- paramsVals[1:n]
    }
    returnType(double(2))
    return(paramsEst)
  }
)


# student T distribution
studentProposal <- nimbleFunction(
  contains =  proposalDistributionVirtual,
  setup = function(model,
                   params,
                   df){},
  run = function(meanDisc = double(0, default = 1),
                 n = integer(0, default = 5),
                 m = integer(0, default = 10),
                 mean = double(1, default = c(0,0)),
                 sigma = double(2, default = diag(2))){
  paramsEst <- matrix(0, nrow = m, ncol = n)
  for(i in 1:m){
    paramsVals <- rmvt_chol(1, mu = mean, chol(sigma), df= df,prec_param = FALSE)
    paramsEst[i, 1:n] <- paramsVals[1:n]
  }
  returnType(double(2))
  return(paramsEst)
  }
)
