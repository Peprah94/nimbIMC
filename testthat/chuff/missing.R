

# Missing covariates
library(myphdthesis)
library(mice)
data(nhanes2)

# data
d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi)) # finding na's
n.mis <- length(idx.mis) # number of nans
d.mis = cbind(age = as.numeric(d.mis$age),
              bmi = d.mis$bmi,
              chl = d.mis$chl)
df = list(d.mis = d.mis, idx.mis = idx.mis)


code <- nimbleCode({

  eta[1: n.idx] ~ dmnorm(muMiss[1:n.idx], cov = covMiss[1:n.idx, 1:n.idx])

  #Fitting the inla with the simulated parameters
  inla.res[1:N, 1:5] <- nimbleINLAMissingValues(x[1:N,1:3], idxMiss[1:n.idx], eta[1: n.idx])

  sigma <- inla.res[1,2]
  linpred[1:N] <-  inla.res[1:N, 5]


  # linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

})

inla_data <- list(y = as.numeric(df$d.mis[,3]),
                  x = df$d.mis,
                  idxMiss = df$idx.mis
)

#Constants
const <- list(N = nrow(df$d.mis),
              n.idx = length(df$idx.mis),
              muMiss = rep(mean(df$d.mis[,2], na.rm = T), length(df$idx.mis)),
              covMiss = diag(mean(df$d.mis[,2], na.rm = T), length(df$idx.mis))


)

# Initial values
idm_inits <- function(){list(eta = rep(0, const$n.idx)
)
}

inlaNimMissing = INLAWiNim (data = df, code = code,
                            modelData = inla_data,
                            modelConstants = const,
                            modelInits = idm_inits,
                            fam = "gaussian",
                            parametersToMonitor = c("eta"),
                            mcmcConfiguration =  list(n.chains = 1,
                                                      n.iterations = 1000,
                                                      n.burnin = 200,
                                                      n.thin = 1,
                                                      setSeed = TRUE,
                                                      samples=TRUE,
                                                      samplesAsCodaMCMC = TRUE,
                                                      summary = TRUE,
                                                      WAIC = FALSE))

save(inlaNimMissing, file = "missingCovariates.RData")
