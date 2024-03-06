

proposalDistributionVirtual <- nimbleFunctionVirtual(
  run = function() {
    returnType(double())
  }
)


# Binomial Proposal distribution
binProposal <- nimbleFunction(
  setup = function(model,
                   binNodesToSimulate,
                   binNodesToFix,
                   binNodesToFixVals,
                   params){},
  run = function(meanDisc = double(0),
                 n = integer(0),
                 size = integer(0)){
    ranEco <- rmybinom(n, prob = meanDisc, size)
    values(model, binNodesToSimulate) <<- ranEco[1:n]
    values(model, binNodesToFix) <<- binNodesToFixVals
    paramsVals <- values(model, params)
    returnType(double(1))
    return(paramsVals)
  }
)


# poisson Proposal distribution
poisProposal <- nimbleFunction(
  setup = function(model,
                   params){},
  run = function(meanDisc = double(0),
                 n = integer(0),
                 lowerBound = double(1),
                 includeLowerBound = integer(0)){
    if(n != length(params)) stop("n should be equal to number of ecological parameters")
    paramsVals <- rmypois(n, lambda = meanDisc, lowerBound =lowerBound, includeLowerBound = includeLowerBound)
    returnType(double(1))
    return(paramsVals)
  }
)

#prior distribution
priorProposal <- nimbleFunction(
  setup = function(model,
                   params){},
  run = function( ){
    model$simulate(nodes = params)
    paramsVals <- values(model, params)
    returnType(double(1))
    return(paramsVals)
  }
)

# Normal distribution
normalProposal <- nimbleFunction(
  setup = function(model,
                   params){},
  run = function(mean = double(1),
                 sigma = double(2)){
    paramsVals <- rmnorm_chol(1, mean, chol(sigma), prec_param = FALSE)
    returnType(double(1))
    return(paramsVals)
  }
)


# student T distribution
studentProposal <- nimbleFunction(
  setup = function(model,
                   params){},
  run = function(mean = double(1),
                 sigma = double(2),
                 df = double(0)){
    paramsVals <- rmvt_chol(1, mu = mean, chol(sigma), df= df,prec_param = FALSE)
    returnType(double(1))
    return(paramsVals)
  }
)



