#multinomial data
library(nimble)
library(myphdthesis)
n.species <- 2; n.sites = 1000
lambda <-  matrix(NA, nrow = n.sites, ncol = n.species)
alpha <- c(2); beta = c(3); cov = rnorm(n.sites)
for(site.tag in 1:n.sites){
  lambda[site.tag, n.species] <- 1
for(spe.tag in 1:(n.species-1)){
lambda[site.tag, spe.tag] <- exp(alpha[spe.tag] + beta[spe.tag]*cov[site.tag])
}
}

prop <- proportions(lambda, 1)

Y <- C <-  matrix(NA, nrow = n.sites, ncol = 1)
for(site.tag in 1:n.sites){
  #for(spe.tag in 1:n.species){
    Y[site.tag,1] <- extraDistr::rcat(1, prop[site.tag, 1:n.species])
}

omega <- matrix(c(0.8, 0.2,
                  0.05, 0.95), nrow = 2, ncol = 2, byrow = TRUE)
#}
for(site.tag in 1:n.sites){
  #for(spe.tag in 1:n.species){
  C[site.tag,1] <- extraDistr::rcat(1, omega[Y[site.tag,1], 1:n.species])
}


code <- nimbleCode({
  #priors
  for(spe.tag in 1:(n.species-1)){
    alpha[spe.tag] ~ dnorm(0, 0.001)
    beta[spe.tag] ~ dnorm(0, 0.001)
  }

  for(spe.tag in 1:n.species){
    for(spe.tag1 in 1:n.species){
    alpha1[spe.tag, spe.tag1]~dexp(1)
    }
  }
  #lambda
  for(site.tag in 1:n.sites){
    lambda[site.tag, n.species] <- 1
    for(spe.tag in 1:(n.species -1)){
      lambda[site.tag, spe.tag] <- exp(alpha[spe.tag] + beta[spe.tag]*cov[site.tag])
    }
  }

  for(site.tag in 1:n.sites){
    for(spe.tag in 1:n.species){
    prop[site.tag, spe.tag] <- lambda[site.tag, spe.tag]/sum(lambda[site.tag, 1:n.species])
    }
  }

  for(site.tag in 1:n.sites){
    Y[site.tag,1] ~ dcat(prop[site.tag, 1:n.species])
  }

  for(spe.tag in 1:n.species){
    omega[spe.tag, 1:n.species] ~ ddirch(alpha1[spe.tag, 1:n.species])
  }

  for(site.tag in 1:n.sites){
    C[site.tag,1] ~ dcat(omega[Y[site.tag,1], 1:n.species])
  }

})

#Data
modelData <- list(Y = Y,
                  C = C,
                  cov = cov)

#Constants
modelConstants <- list(n.species = n.species,
              n.sites = n.sites
)

# Initial values
idm_inits <- function(){list(omega = matrix(c(0.98, 0.02, 0.01,
                                              0.99),
                                            nrow=2, ncol=2, byrow = TRUE),
                             alpha = rep(0,1),
                             beta = rep(0,1),
                             alpha1 = matrix(rep(1, length(unique(C))),
                                             nrow=2, ncol = 2)
)
}

cnt = 0
assign('myRW_dirichlet', myRW_dirichlet, envir = .GlobalEnv)

initsList <- idm_inits()


#Create the model in nimble
mwtc <- nimble::nimbleModel(code,
                            data = modelData,
                            constants = modelConstants,
                            inits = initsList)

# Create the model in C
Cmwtc <- nimble::compileNimble(mwtc,
                               showCompilerOutput = FALSE) #Have issues compiling


  mcmcconf <- nimble::configureMCMC(Cmwtc,
                                    monitors = c("alpha", "beta", "omega"))

mcmcconf$removeSampler("omega")
  mcmcconf$addSampler(target = c("omega"),
                      type = "myRW_dirichlet")


Rmcmc <- nimble::buildMCMC(mcmcconf)

# Compile
cmcmc <- nimble::compileNimble(Rmcmc,
                               project = Cmwtc,
                               resetFunctions = TRUE)

#MCMC Configurations

mcmcConfiguration =  list(n.chains = 4,
                          n.iterations = 1000,
                          n.burnin = 200,
                          n.thin = 1,
                          setSeed = TRUE,
                          samples=TRUE,
                          samplesAsCodaMCMC = TRUE,
                          summary = TRUE,
                          WAIC = FALSE)

# Run the MCMC
startTime <- Sys.time()
mcmc.out <- nimble::runMCMC(cmcmc,
                            niter = mcmcConfiguration[["n.iterations"]],
                            nchains = mcmcConfiguration[["n.chains"]],
                            nburnin = mcmcConfiguration[["n.burnin"]],
                            #inits = initsList,
                            thin = mcmcConfiguration[["n.thin"]],
                            setSeed = mcmcConfiguration[["setSeed"]],
                            samples = mcmcConfiguration[["samples"]],
                            samplesAsCodaMCMC = mcmcConfiguration[["samplesAsCodaMCMC"]],
                            summary = mcmcConfiguration[["summary"]],
                            WAIC = mcmcConfiguration[["WAIC"]])

