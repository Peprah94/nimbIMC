% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inlaWithNimbleAlternative.R
\name{INLAWiNimDataGeneratingTargetDivide}
\alias{INLAWiNimDataGeneratingTargetDivide}
\title{Fitting INLA within NIMBLE for Data Generating process}
\usage{
INLAWiNimDataGeneratingTargetDivide(
  data,
  covariate,
  code,
  family,
  modelData,
  modelConstants,
  modelInits,
  nimbleINLA,
  inlaMCMC = c("inla", "mcmc", "inlamcmc"),
  inlaMCsampler = "RW_INLA_block",
  samplerControl = list(),
  parametersToMonitor = list(mcmc = c("mcmc"), mcmc2 = c("mcmc"), inla = c("inla"),
    additionalPars = NULL),
  mcmcConfiguration = list(n.chains = 1, n.iterations = 10, n.burnin = 0, n.thin = 1,
    setSeed = TRUE, samples = TRUE, samplesAsCodaMCMC = TRUE, summary = TRUE, WAIC =
    FALSE)
)
}
\arguments{
\item{data, }{code, family,n.iterations, n.chains, n.burnin}
}
\value{
MCMC output
}
\description{
This function sets the paramters in the appropriate manner to be used by the
simulation function
}
