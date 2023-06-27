# Bivariate regression


#Simulation studies

library(myphdthesis)
library(INLA)
library(inlabru)
library(nimble)

load("/Volumes/kwakupa/INLAWithinMCMC/data_for_simulations.RData")

# Nimble code
code <- nimble::nimbleCode({
  #Prior for beta1 and beta2
  #beta[1:2] ~ dmnorm(mu[1:2], cov=precision_matrix[1:2,1:2])
  for(i in 1:2){
    beta[i] ~ dnorm(0, tau = 0.001)
  }

  a ~ dnorm(0, tau = 0.001)
  #Fitting the inla with the simulated parameters
 # inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)


  for(i in 1:N){
  linpred[i] <-  a + beta[1]*x[i, 1] + beta[2]*x[i, 2]
  }

  #Bivariate linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau= sigma)
  }

  #tausig <- 1/ sigma * sigma

  sigma ~ dgamma(1,0.00005)
  #intercept <- inla.res[1,1]
})


data = df
idm_data <- list(y=data$y,
                 x = data$x,
                 y_obs=data$y,
                 inter = 1)

constants = list(N = length(data$y),
                 mu = c(0,0),
                 precision_matrix = diag(5,2),
                 fam = "gaussian")

inits <-  function(){list(beta =c(1,1),
                          a = 1,
                          sigma = 1
)
}

initsList <- inits()


mwtc <- nimble::nimbleModel(code,
                            data = idm_data,
                            constants =constants,
                            inits = initsList)


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
  data$oset = data$x %*% (beta)

  # if(interInModel == 1){
  #   formula =  1 + offset(data$oset)
  # }else{
  #   formula = - 1 + offset(data$oset)
  # }

  res = INLA::inla(y ~ 1 + offset(data$oset),
                   data = data,
                   family = family,
                   verbose=FALSE,
                   control.predictor = list(compute = TRUE))
  #fitted_values = res$summary.fitted.values[,"mean"]
  fitted_values = c(res$mlik[1,1])
  #intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
#  if(interInModel == 1){
    intercept = INLA::inla.emarginal(function(x) x,res$marginals.fixed[[1]])
 # }else{
  #  intercept = 1
 # }
  precision = INLA::inla.emarginal(function(x) x,res$marginals.hyper[[1]])
  ret <- data.frame(mld = fitted_values, intercept, precision,
                    row.names = NULL)
  colnames(ret) <- c("mld", fixedVals)

  ret <- as.matrix(ret)
  return(ret)
}

nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=double(1), #y is a vector
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("a", "sigma")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "gaussian")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'fit.inla'
)


# nimbleINLA(x = data$x,
#            y = data$y,
#            beta =  c(3,-3),
#            fixedVals = c("a", "sigma")#,
#            #interInModel = 1
#            )
# x <- data$x
# y <- data$y
# interInModel <- 1
# fixedVals = c("a", "sigma")
#
#
#
# Cmwtc <- nimble::compileNimble(mwtc,
#                                showCompilerOutput = FALSE)

ret <- INLAWiNimDataGenerating(data = c("y"),
                        covariate = data$x,
                        code = code,
                        family = "gaussian",
                        modelData = idm_data,
                        modelConstants = constants,
                        modelInits = inits,
                        nimbleINLA = nimbleINLA,
                        inlaMCMC = c("inlamcmc"),
                        inlaMCsampler = "RW_INLA_block",
                        samplerControl = list(scale = 0.75^2),
                        parametersToMonitor = list(mcmc = c("beta"),
                                                   inla = c("a", "sigma")),
                        mcmcConfiguration = list(n.chains = 1,
                                                 n.iterations = 1000,
                                                 n.burnin = 200,
                                                 n.thin = 1,
                                                 setSeed = TRUE,
                                                 samples=TRUE,
                                                 samplesAsCodaMCMC = TRUE,
                                                 summary = TRUE,
                                                 WAIC = FALSE)
                        )





# ret <- buildINLAmodel(mwtc,
#                       fam = "gaussian",
#                       x = data$x, y= "y",
#                       control = list(fit.inla = nimbleINLA,
#                                       fixedVals = fixedVals))
#
# cwtm <- compileNimble(ret,
#                       showCompilerOutput = FALSE)
# cwtm$run(beta = c(3, -3))
# cwtm$mvEWSamples[["sigma"]]
#
#
# mwtc <- nimble::nimbleModel(code,
#                             data = idm_data,
#                             constants = constants,
#                             inits = initsList)
# #cwtm$run(data$x, data$y, c(0,0))
# # Create the model in C
# Cmwtc <- nimble::compileNimble(mwtc,
#                                showCompilerOutput = FALSE) #Have issues compiling
#
# mcmcconf <- nimble::configureMCMC(Cmwtc,
#                                   nodes = NULL)
#
# #mcmcconf$removeSampler(c("beta", "a", "sigma"))
#
# mcmcconf$addSampler(target = "beta",
#                     type = "RW_INLA_block",
#                     control = list(fit.inla = nimbleINLA,
#                                    x = data$x,
#                                    y = data$y,
#                                    fixedVals = fixedVals,
#                                    fam = "gaussian",
#                                    scale = 1))
#
# mcmcconf$printSamplers()
# Rmcmc <- nimble::buildMCMC(mcmcconf)
#
# # Compile
# cmcmc <- nimble::compileNimble(Rmcmc,
#                                project = Cmwtc,
#                                resetFunctions = TRUE)
#
# #cmcmc$run(1000)
#
# cmc.out <- nimble::runMCMC(cmcmc,
#                            niter = 10,
#                            nchains = 2,
#                            nburnin = 5,
#                            #inits = initsList,
#                            thin = 1,
#                            setSeed = TRUE,
#                            samples = TRUE,
#                            samplesAsCodaMCMC = TRUE,
#                            summary = TRUE,
#                            WAIC = FALSE)
#
# cmc.out$summary$all.chains
#
# ggmcmc::ggs(cmc.out$samples)%>%
#   ggs_density()
