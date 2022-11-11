library(nimble)
library(parallel)

cnt <- 0
listout <- list()

sampler_BASE <- nimbleFunctionVirtual(
  name = 'sampler_BASE',
  methods = list(
    reset = function() { }
  )
)

myRW_dirichlet <- nimbleFunction(
  #name = 'sampler_RW_dirichlet',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scaleOriginal       <- extractControlElement(control, 'scale',               1)
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    #calcNodes <- model$getDependencies(target,includeData = TRUE, determOnly = FALSE, omit = c("prop", "lambda", "lambda_obs", "C"), stochOnly = TRUE) #stochOnly = TRUE
    #calcNodesNoSelf <- model$getDependencies(target, includeData = TRUE, determOnly = FALSE, omit = c("prop", "lambda", "lambda_obs", "C"), stochOnly = TRUE,self = FALSE)
    calcNodes <- model$getDependencies(target,includeData = TRUE) #stochOnly = TRUE
    calcNodesNoSelf <- model$getDependencies(target, includeData = TRUE, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    targets <- model$expandNodeNames(target)
    ## numeric value generation
    #d <- length(targetAsScalar)
    d <- length(model$expandNodeNames(target))
    dr <- length(targetAsScalar)/d
    thetaVec         <- rep(0, dr)
    scaleVec         <- rep(scaleOriginal, d)
    timesRan         <- 0
    timesAcceptedVec <- rep(0, d)
    timesAdapted     <- 0
    optimalAR        <- 0.44
    gamma1           <- 0
    ## checks
    # if(length(model$expandNodeNames(target)) > 1)    stop('RW_dirichlet sampler only applies to one target node')
    #if(model$getDistribution(target) != 'ddirch')    stop('can only use RW_dirichlet sampler for dirichlet distributions')
  },
  run = function() {
     extralpD <- model$calculateDiff(calcNodesNoSelf)
     for(j in 1:dr)
    for(i in 1:d){
      #targets[i] <- targets[i]

      if(thetaVec[1] == 0)   thetaVec <<- values(model, targets[i])   ## initialization

      alphaVec <- model$getParam(targets[i], 'alpha')
      #for(i in 1:d) {
      currentValue <- thetaVec
      propLogScale <- rnorm(dr, mean = 0, sd = scaleVec[1])
      propValue <- currentValue * exp(propLogScale)
      #propValue <- exp(propLogScale)
      if(all(propValue != 0)) {
        thetaVecProp <- thetaVec
        thetaVecProp <- propValue
        values(model, targets[i]) <<- thetaVecProp / sum(thetaVecProp)
        logMHR <- sum((alphaVec-1)*propValue) - sum((alphaVec-1)*currentValue) + sum(exp(-0.5 * (log(currentValue)^2 - log(propValue)^2)/scaleVec[1])*(propValue/currentValue)) + extralpD
        jump <- decide(logMHR)
      } else jump <- FALSE
      if(adaptive & jump)   timesAcceptedVec[i] <<- timesAcceptedVec[i] + 1
      if(jump) {
        thetaVec <<- thetaVecProp
        nimCopy(from = model, to = mvSaved, row = 1, nodes = targets[i], logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
      } else {
        nimCopy(from = mvSaved, to = model, row = 1, nodes = targets[i], logProb = TRUE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
      }
      model$calculate(targets[i])                                                             ## update targets[i] logProb
      nimCopy(from = model, to = mvSaved, row = 1, nodes = targets[i], logProbOnly = TRUE)    ##
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
      timesRan         <<- 0
      timesAcceptedVec <<- numeric(d, 0)
      timesAdapted     <<- 0
      gamma1           <<- 0
    }
  )
)


assign('myRW_dirichlet', myRW_dirichlet, envir = .GlobalEnv)








dmydirch <- nimbleFunction(
  run = function(x = double(1), alpha = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    logProb <- sum(lgamma(alpha)) - lgamma(sum(alpha)) +
      sum((alpha -1) * log(x))
    if(log) return(logProb)
    else return(exp(logProb))
  })

rmydirch <- nimbleFunction(
  run = function(n = integer(0), alpha = double(1)) {
    returnType(double(1))
    if(n != 1) print("rdirch only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(p)
  })

# omega
omega_fnx <- function(p11, p22){
  omega <- matrix(NA, 2, 2)
  omega[1,1] <- p11
  omega[2,2] <- p22
  omega[1,2] <- 1- p11
  omega[2,1] <- 1- p22
  return(omega)
}

nimble_omega <- nimbleRcall(
  prototype = function(
    p11=double(0),
    p22 = double(0)#x is a matrix
    # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'omega_fnx'
)

#omega <- nimb

#my_laplace <- buildLaplace()



library(nimble)
library(dirmult)
code <-nimbleCode({


  for(i in 1:nspecies){
    omega[i, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])

 }

   #p11 ~ dunif(0.2,1)
   #p22 ~ dunif(0.2,1)

#omega[1:nspecies, 1:nspecies] <- nimble_omega(p11, p22)

  r[1:nsite,1:20] <- nimbleINLADataGenerating(omega[1:nspecies,1:nspecies]) #change 38 to a constant to be specified

  for(i in 1:nsite){
    # for(j in 1:nspecies){
    log(lambda_obs[i,1]) <- r[i,5] + r[i,6]*true_cov[i] + r[i,11]+
      r[i, 7] + r[i,8]*bias_cov[i]+ r[i,12] - log(1+exp(r[i, 7] + r[i,8]*bias_cov[i]+ r[i,12])) +
      r[i, 9] + r[i,10]*det_cov[i] - log(1+exp(r[i, 9] + r[i,10]*det_cov[i]))
    # Second species
    log(lambda_obs[i,2]) <- r[i,13] + r[i,14]*true_cov[i] + r[i,19]+
      r[i, 15] + r[i,16]*bias_cov[i]+ r[i,20] - log(1+exp(r[i, 15] + r[i,16]*bias_cov[i]+ r[i,20])) +
      r[i, 17] + r[i,18]*det_cov[i] - log(1+exp(r[i, 17] + r[i,18]*det_cov[i]))

  }

  lambda[1:nsite, 1:nspecies] <- lambda_obs[1:nsite, 1:nspecies]

  #Proportion of lambdas


  # Proportion for the multinomial distribution
  for(site.tag in 1:nsite){
    for(spe.tag in 1:nspecies){
      prop[site.tag,spe.tag] <- (lambda[site.tag, spe.tag])/sum(lambda[site.tag, 1:nspecies])
    }
  }


  # True data
  for(site.tag in 1:nsite){
    C[site.tag] ~ dcat(prop[site.tag,1:nspecies])
  }

  # Reported species
  for(site.tag in 1:nsite){
    Y[site.tag] ~ dcat(omega[C[site.tag],1:nspecies])
  }

})


## Parameterising the nimble model

#Data
inla_data <- list(Y=data_df$Y,
                  C = data_df$C,
                  true_cov = data_df$eco_cov,
                  bias_cov=data_df$samp_cov,
                  det_cov= data_df$det_cov)

#Constants
const <- list(nspecies=length(unique(data_df$C)),
              nsite = length(data_df$C),
              alpha=rep(1, length(unique(data_df$C)))
)

# Initial values
  idm_inits <- function(){list(omega = matrix(c(0.99, 0.01,
                                                0.01, 0.99),
                                              nrow=2, ncol=2, byrow = TRUE)
                               )
  }

  initsList <- idm_inits()

  #Putting all together for the creating compilation
  modelInfo <- list(
    code = code,
    constants = const,
    data = inla_data,
   inits = initsList
  )

  #Create the model in nimble
  mwtc <- nimbleModel(code,
                      data = inla_data,
                      constants = const,
                     inits = initsList)
  #)
  #library(igraph)
  #plot(mwtc$modelDef$graph)

  # Create the model in C
  Cmwtc <- compileNimble(mwtc,
                         showCompilerOutput = FALSE) #Have issues compiling


  mcmcconf <- configureMCMC(Cmwtc,
                            print=TRUE,
                            useConjugacy=FALSE,
                            monitors = c("omega"))

  #mcmcconf$removeSamplers(c("omega[1,1:2]","omega[2,1:2]"))
  #mcmcconf$addSampler(c("omega[1,1:2]"), "myRW_dirichlet")
 # mcmcconf$addSampler(c("omega[2,1:2]"), "myRW_dirichlet")
  mcmcconf$removeSamplers(c("omega"))
  mcmcconf$addSampler(c("omega"), "myRW_dirichlet")



  Rmcmc <- buildMCMC(mcmcconf)
                     #enableWAIC =FALSE)

  # Compile
  cmcmc <- compileNimble(Rmcmc,
                         project = Cmwtc,
                         resetFunctions = TRUE)

# Run the MCMC
#library(pbapply)

  mcmc.out <- runMCMC(cmcmc,
                            niter = 10,
                           #nburnin = 500,
                           inits = initsList,
                           # thin =100,
                            #setSeed = x,
                            samples=TRUE,
                            samplesAsCodaMCMC = TRUE,
                            summary = TRUE,
                            WAIC = FALSE)
save(mcmc.out, file="estimates.RData")
