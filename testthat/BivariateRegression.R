library(myphdthesis)


# Bivariate regression

data("data_for_simulations")

# Nimble code
code <- nimble::nimbleCode({
  #Prior for beta1 and beta2
  #beta[1:2] ~ dmnorm(mu[1:2], cov=precision_matrix[1:2,1:2])
  for(i in 1:2){
    beta[i] ~ dnorm(0, tau = 0.01)
  }

  #Fitting the inla with the simulated parameters
  inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2],inter)

  sigma <- inla.res[1,2]
  linpred[1:N] <-  inla.res[1:N, 1]

  #Bivariate linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

  tau <- sigma
  intercept <- inla.res[1,1]
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

inits <-  function(){list(beta =c(1,1)
)
}

# Fit the INLA within Nimble model
inlaNimBivariate = INLAWiNim (data = df, code = code,
                              modelData = idm_data,
                              modelConstants = constants,
                              modelInits = inits,
                              fam = "gaussian",
                              mcmcControl = list(scale = 0.75,
                                                 adaptive = TRUE,
                                                 propCov='identity',
                                                 adaptInterval=50),
                              mcmcConfiguration =  list(n.chains = 1,
                                                        n.iterations = 1000,
                                                        n.burnin = 200,
                                                        n.thin = 1,
                                                        setSeed = TRUE,
                                                        samples=TRUE,
                                                        samplesAsCodaMCMC = TRUE,
                                                        summary = TRUE,
                                                        WAIC = FALSE))

save(inlaNimBivariate, file = "BayesianRegression.RData")

