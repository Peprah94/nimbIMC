
#Bayesian Lasso regression

fixedVals <- c("alpha","sigma")
fit.inla <- function(x ,
                     y ,
                     beta,
                     fixedVals,
                     #interInModel,
                     family
){
  #ii  <- get("ii",envir =  parent.frame())
  #ii <- assign("ii",ii+1,envir = parent.frame())
  #print(ii)
  data <- list(y=y, x=x)
  data$oset = data$x %*% beta
  res = INLA::inla(y ~ 1 + offset(oset), data = data, family = family,
                   control.predictor = list(compute=TRUE))
  #res = INLA::inla.rerun(res)
  fitted_values = c(res$mlik[1,1])
  #fitted_values = res$summary.fitted.values[,"mean"]
  # if(interInModel == 1){
  #   intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  # }else{
  #   intercept = 1
  # }
  intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  #ret <- cbind(fitted_values, precision)

  ret <- data.frame(mld = fitted_values, intercept, precision,
                    row.names = NULL)
  colnames(ret) <- c("mld", fixedVals)

  ret <- as.matrix(ret)
  #ret <- c(ret)
  return(ret)
}



nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("alpha","sigma")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "gaussian")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)

load("/Volumes/kwakupa/INLAWithinMCMC/Bayesian Lasso/hitters_data.RData")
df <- lassoDataDF
 ret <- nimbleINLA(df$x, c(df$y), c(0,1,1,0,0), c("alpha","sigma"),  "gaussian")

code <- nimbleCode({
  #Prior for beta1 and beta2
  alpha ~ dnorm(0,0.001)
  # for(j in 1:P){
  # beta1 ~ ddexp(location = 0,rate=est_lam)
  # beta2 ~ ddexp(location = 0,rate=est_lam)
  # beta3 ~ ddexp(location = 0,rate=est_lam)
  # beta4 ~ ddexp(location = 0,rate=est_lam)
  # beta5 ~ ddexp(location = 0,rate=est_lam)
  #
  # beta[1:5] <- nimC(beta1, beta2, beta3, beta4, beta5)
  #}

  for(j in 1:P){
    beta[j] ~ ddexp(location = 0,rate=est_lam)
  }

  #Fitting the inla with the simulated parameters
  for(i in 1:N){
  linpred[i] <- inprod(beta[1:P], x[i, 1:P]) + alpha
  }

  #linpred[1:N] <- inla.res[1:100,3]
  sigma ~ dgamma(1,0.00005)
  #linpred[1:N] <-  inla.res[1:N, 2]

  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i],tau=sigma )
  }

  # tau <- sigma
  #intercept <- inla.res[1,1]
})

## Parameterising the nimble model

#Data
df = lassoDataDF
inla_data <- list(y=as.numeric(df$y),
                  x = df$x,
                  y_obs=as.numeric(df$y))

#Constants
const <- list(N = length(df$y),
              P= ncol(df$x),
              est_lam = 1/0.73
)

# Initial values
idm_inits <- function(){list(beta=rep(1,const$P),
                             sigma = 1
)
}

initsList <- idm_inits()

data = df

stdev.samp <- .25 * solve(t(data$x)%*%data$x)

ret <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = data$x,
                               code = code,
                               family = "gaussian",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = c("inlamcmc"),
                               inlaMCsampler = "RW_INLA_block",
                               samplerControl = list(propCov = stdev.samp,
                                                     #interInModel = 0,
                                                     mu = c(rep(0,5)),
                                                     scale = 1,
                                                     adaptive = FALSE),
                               parametersToMonitor = list(inla = c("alpha","sigma"),
                                                          mcmc = c("beta")),
                               mcmcConfiguration = list(n.chains = 1,
                                                        n.iterations = 100,
                                                        n.burnin = 20,
                                                        n.thin = 1,
                                                        setSeed = TRUE,
                                                        samples=TRUE,
                                                        samplesAsCodaMCMC = TRUE,
                                                        summary = TRUE,
                                                        WAIC = FALSE)
)



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

# Create the model in C
Cmwtc <- compileNimble(mwtc,
                       showCompilerOutput = FALSE) #Have issues compiling

data = df
ret <- buildINLAmodel(mwtc,
                      fam = "gaussian",
                      x = data$x, y= c(data$y),
                      control = list(fit.inla = nimbleINLA,
                                     fixedVals = fixedVals))


cwtm <- compileNimble(ret,
                      showCompilerOutput = FALSE)

cwtm$run(beta = c(1,1,1,0,0))
cwtm$mvEWSamples

mcmcconf <- nimble::configureMCMC(Cmwtc,
                                  nodes = NULL)

#mcmcconf$removeSampler(c("beta", "a", "sigma"))


mcmcconf$addSampler(target = c("beta"),
                    type = "RW_INLA_block",
                    control = list(fit.inla = nimbleINLA,
                                   x = data$x,
                                   y = c(data$y),
                                   fixedVals = fixedVals,
                                   adaptive = TRUE,
                                   fam = "gaussian",
                                   propCov = stdev.samp,
                                   #interInModel = 0,
                                   scale = 1))

mcmcconf$printSamplers()
Rmcmc <- nimble::buildMCMC(mcmcconf)

# Compile
cmcmc <- nimble::compileNimble(Rmcmc,
                               project = Cmwtc,
                               resetFunctions = TRUE)

#cmcmc$run(1000)

cmc.out <- nimble::runMCMC(cmcmc,
                           niter = 1000,
                           nchains = 1,
                           nburnin = 100,
                           #inits = initsList,
                           thin = 1,
                           setSeed = TRUE,
                           samples = TRUE,
                           samplesAsCodaMCMC = TRUE,
                           summary = TRUE,
                           WAIC = FALSE)
cmc.out$summary


