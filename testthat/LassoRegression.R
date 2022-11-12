# Bivariate regression
library(myphdthesis)

#bayesian Lasso
data("hitters_data")

code <- nimbleCode({

  alpha ~ dnorm(0,1)


  for(j in 1:P){
    beta[j] ~ ddexp(location = 0, rate=est_lam)
  }
  #Fitting the inla with the simulated parameters
  inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:P],y_obs[1:N],beta[1:P], interInModel)

  #linpred[1:N] <- inla.res[1:100,3]
  sigma <- inla.res[1,2]
  linpred[1:N] <-  inla.res[1:N, 3] + alpha

  intercept <- inla.res[1,2]

  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

  # tau <- sigma
  #intercept <- inla.res[1,1]
})

inla_data <- list(y=as.numeric(df$y),
                  x = df$x,
                  y_obs=as.numeric(df$y),
                  interInModel = 2)

#Constants
const <- list(N = length(df$y),
              P= ncol(df$x),
              est_lam = 1/0.073
)

# Initial values
idm_inits <- function(){list(alpha = 0,
                             beta=rep(0,const$P)
)
}

inlaNimLasso = INLAWiNim (data = df, code = code,
                          modelData = inla_data,
                          modelConstants = const,
                          modelInits = idm_inits,
                          fam = "gaussian",
                          parametersToMonitor = c("beta", "alpha","sigma"),
                          mcmcSamplerChange = TRUE,
                          parametersForSamplerChange = "beta",
                          newSampler = "RW_block",
                          mcmcConfiguration =  list(n.chains = 1,
                                                    n.iterations = 1000,
                                                    n.burnin = 200,
                                                    n.thin = 1,
                                                    setSeed = TRUE,
                                                    samples=TRUE,
                                                    samplesAsCodaMCMC = TRUE,
                                                    summary = TRUE,
                                                    WAIC = FALSE))

save(inlaNimLasso, file = "BayesianLassoResults.RData")
