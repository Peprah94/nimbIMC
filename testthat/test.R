input <- misclassSetPars(nspecies = 2,
                         nsites = 1000,
                         nvisit = 1,
                         nvalidate = 333,
                         nreported = 3,
                         beta0 = c(-1,2),
                         beta1 = matrix(c(-2, 1, 2, -3), 2, 2, byrow = TRUE),
                         gamma0 = matrix(c(1, 0.5, 0,
                                           1, 1.5, 0.5), nrow = 2,
                                         ncol = 3, byrow = TRUE),
                         gamma1 = matrix(c(10, -1, -2,
                                            -4, - 8, -3), nrow = 2,
                                             ncol = 3, byrow = TRUE),
                         sigma_cov = 1,
                         sigma_class = 1,
                         sigma_corr = matrix(c(1,0.8,0.8,1), nrow = 2, byrow= T),
                         ncov = 2
                         )

data_nimble <- misclassGenData(input, "correlation", 1)

nim <- misclassFitModel(data_nimble,
                        "constant",
                        TRUE,
                        10000,
                        3,
                        5000)

data("nimble_data_ML")
nimData <- misclassFitModel(dataout,
                            "variable",
                            TRUE,
                            100,
                            2,
                            50,
                            parameters_to_monitor =c("beta0","beta1",
                                                     "gamma0","gamma1",
                                                     "accuracy", "precision",
                                                     "p", "recall"))




#####################
#Distance to nearest road

gull <- read_delim("/Volumes/kwakupa-1/misclassification/new/gull_new_dataset/occurrence.txt",
                   delim = "\t", escape_double = FALSE,
                   trim_ws = TRUE)
distGBIFdata(gull)

# Bivariate regression
code <- nimble::nimbleCode({
  #Prior for beta1 and beta2
  #beta[1:2] ~ dmnorm(mu[1:2], cov=precision_matrix[1:2,1:2])
  for(i in 1:2){
    beta[i] ~ dnorm(0, tau = 0.01)
  }

  #Fitting the inla with the simulated parameters
  inla.res[1:N, 1:3] <- nimbleINLA(x[1:N,1:2],y_obs[1:N],beta[1:2], inter)

  sigma <- inla.res[1,2]
  linpred[1:N] <-  inla.res[1:N, 1]

  #Bivariate linear model specification
  for(i in 1:N){
    y[i] ~ dnorm(linpred[i], tau=sigma)
  }

  tau <- sigma
  intercept <- inla.res[1,1]
})

data("data_for_simulations")
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

inlaNimBivariate = INLAWiNim (data = df, code = code,
                     modelData = idm_data,
                     modelConstants = constants,
                     modelInits = inits,
                     fam = "gaussian",
                     mcmcConfiguration =  list(n.chains = 1,
                                               n.iterations = 10,
                                               n.burnin = 0,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = FALSE))
#Plot results
ggSamples <- ggmcmc::ggs(inlaNimBivariate$mcmc.out$samples)
ggmcmc::ggs_pairs(ggSamples, lower = list(continuous = "density"))

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

inlanim = INLAWiNim (data = df, code = code,
                     modelData = inla_data,
                     modelConstants = const,
                     modelInits = idm_inits,
                     fam = "gaussian",
                     parametersToMonitor = c("beta", "alpha","sigma"),
                     mcmcConfiguration =  list(n.chains = 1,
                                               n.iterations = 10,
                                               n.burnin = 0,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = FALSE))
ggSamples <- ggmcmc::ggs(inlanim$mcmc.out$samples)
ggmcmc::ggs_density(ggSamples)+
  ggplot2::facet_wrap(~Parameter, ncol = 3)

# Missing covariates
data("nhanes2")

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

inlanim = INLAWiNim (data = df, code = code,
                     modelData = inla_data,
                     modelConstants = const,
                     modelInits = idm_inits,
                     fam = "gaussian",
                     parametersToMonitor = c("eta"),
                     mcmcConfiguration =  list(n.chains = 1,
                                               n.iterations = 10,
                                               n.burnin = 0,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = FALSE))
ggSamples <- ggmcmc::ggs(inlanim$mcmc.out$samples)
ggmcmc::ggs_density(ggSamples)+
  ggplot2::facet_wrap(~Parameter, ncol = 3)+
  ggplot2::theme_classic()
#ggmcmc::ggs_pairs(ggSamples)



###########################
# DataGenerating process
###########################
library(spatstat)
library(maptools)

library(sp)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(maptools)
library(spatstat)
library(gdata)
library(ggplot2)
library(gridExtra)
#library(PCDSpline)
library(foreach)
library(doParallel)
library(viridis)
library(RandomFieldsUtils)
#library(devtools)
#install_github("cran/tiff")
require(devtools)
x0 <- seq(-1, 4, length = 100)
y0 <- seq(-1,4, length = 100)
gridlocs <- expand.grid(x0,y0)

# Covariates for true ecological state
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate.im <- im(gridcov, x0, y0)

#Covariate for the sampling process
gridcov_thin <- outer(x0, y0, function(x,y) cos(2*x) - sin(2*y-4))
#gridcov_thin <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))
covariate_thin.im <- im(gridcov_thin, x0, y0)

#Covariate for the detection
gridcov_det <- outer(x0, y0, function(x,y) (x/2)^2+(y/2)^2)
covariate_detect.im <- im(gridcov_det, x0, y0)

#Input for the simulations
# input <- {list(
#   ecological = list(
#     fixed.effect=list(
#       intercept = c(0.8, 2.5, -1.5, 1.2),
#       betacov = c(1.5, -0.12, 2,-0.4)
#     ),
#     hyperparameters = list(
#       sigma2 = c(0.2, 1.2, 2, 0.1),
#       range = c(1.2, 2.5, 3.2, 0.22)
#     )
#   ),
#   sampling = list(
#     fixed.effect = list(
#       intercept = c(1.3),
#       betacov = c(-1.5)
#     ),
#     hyperparameters=list(
#       sigma2 = c(0.2),
#       range = c(2.5)
#     )
#   ),
#   detection = list(
#
#     fixed.effect = list(
#       intercept=c(2,-0.3, 5, 1.2),
#       betacov = c(-2, -0.5, -2.5, 2)
#     )
#   ),
#
#   misclassification = list(
#
#     class_prob <- matrix(c(0.9, 0.02, 0.04, 0.04,
#                            0.05, 0.89, 0.04, 0.02,
#                            0.1,0.1, 0.8, 0,
#                            0, 0.05, 0.25, 0.7),
#                          nrow=4, ncol=4, byrow = TRUE)
#
#
#   ),
#   constants = list(n.species = 4,
#                    seed= 1036610620 ),
#   plot = list(ecological=TRUE,
#               detection=FALSE,
#               sampling=FALSE,
#               all= FALSE,
#               classification=FALSE),
#   cov = list(covariate.im,
#              covariate_thin.im,
#              covariate_detect.im),
#   idxs = list(eco=c(1),
#                sampling=c(2),
#                detection=c(3))
# )}

# 2 species
nspecies <- 2#nspecies: Number of species we want to simulate
#input: List with the parameters of the model that generates CS data
input <-{list(
  ecological = list(
    fixed.effect=list(
      intercept = c(0.8, 2.5 ),
      betacov = c(1.5, -0.12)
    ),
    hyperparameters = list(
      sigma2 = c(0.2, 1.2),
      range = c(1.2, 2.5)
    )
  ),
  sampling = list(
    fixed.effect = list(
      intercept = c(1.3),
      betacov = c(-1.5)
    ),
    hyperparameters=list(
      sigma2 = c(0.2),
      range = c(2.5)
    )
  ),
  detection = list(

    fixed.effect = list(
      intercept=c(2,-0.3),
      betacov = c(-2, -0.5)
    )
  ),

  misclassification = list(

    class_prob <- matrix(c(0.7, 0.3,
                           0.35, 0.65),
                         nrow=2, ncol=2, byrow = TRUE)
  ),
  constants = list(n.species = 2,
                   seed= 1036610620 ),
  plot = list(ecological= FALSE,
              detection=FALSE,
              sampling=FALSE,
              all= FALSE,
              classification=FALSE),
  cov = list(covariate.im,
             covariate_thin.im,
             covariate_detect.im),
  idxs = list(eco=c(1),
              sampling=c(2),
              detection=c(3))
)}

## Simulating the Covariates ##
simulateddata <- generateCSData(input)

nspecies = input$constants$n.species
## Covariates need to be in SpatialPixels format ##
#Covariates for true intensity
cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate.im$v))))))
r <- raster(cov1.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
cov1.spix <- as(cov1.rast,"SpatialPixelsDataFrame")

#Covariates for first thinning
cov2.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_thin.im$v))))))
r <- raster(cov2.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov2.rast <- rasterize(cov2.sp@coords,r1,cov2.sp$cov, fun=mean,na.rm=T)
cov2.spix <- as(cov2.rast,"SpatialPixelsDataFrame")

#Covariate for second thinning
cov3.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(covariate_detect.im$v))))))
r <- raster(cov3.sp)
r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
cov3.rast <- rasterize(cov3.sp@coords,r1,cov3.sp$cov, fun=mean,na.rm=T)
cov3.spix <- as(cov3.rast,"SpatialPixelsDataFrame")

## Extra information on species detection ##

### Sampling the detections from a survey
rndpts_x0 <- runif(50, 0,3)
rndpts_y0 <- runif(50, 0,3)
rndpts <- data.frame(rndpts_x0, rndpts_y0)

det_prob <- list()
for(i in 1:nspecies){
  rndpts_lin <- rndpts %>%
    mutate(linpred = input$detection$fixed.effect$intercept[i] + input$detection$fixed.effect$betacov[i]* extract(cov3.rast,rndpts))
  det_prob[[i]] <- plogis(rndpts_lin$linpred)
}
#data_det <- rbinom(length(det_prob),1,det_prob)

detection_data <- list()
for(i in 1:nspecies){
  data_det <- vector("numeric", length(det_prob[[i]]))
  for(j in 1:length(det_prob[[i]])){
    data_det[j]<- rbinom(1,1, det_prob[[i]][j])
    detection_data[[i]] <- data_det
  }
}
#Organising as spatial dataframe
data_det_spframe <- list()
for(i in 1:nspecies){
  data_det_spframe[[i]] <- SpatialPointsDataFrame(rndpts, data = data.frame(detection_data[[i]]))
  names(data_det_spframe[[i]])<-paste0("detdata",i)
}


## Fit the model using inlabru ##

## the borders of the study region
coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
poly <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))

## the mesh
mesh <- inla.mesh.2d(loc.domain = coordsmat, offset = c(0.3, 1),
                     max.edge = c(0.1, 0.5), cutoff = 0.2)

## SPDEs definition
spdes <- list()
for(i in 1: nspecies){
  spdes[[i]] <- inla.spde2.pcmatern(mesh = mesh,
                                    # PC-prior on range: P(practic.range < 0.05) = 0.01
                                    prior.range = c(input$ecological$hyperparameters$range[i], 0.5),
                                    # PC-prior on sigma: P(sigma > 1) = 0.01
                                    prior.sigma = c(sqrt(input$ecological$hyperparameters$sigma2[i]), 0.5))
}

#SPDEs for the thinning
spde2 <- inla.spde2.pcmatern(mesh = mesh,
                             # PC-prior on range: P(practic.range < 0.05) = 0.01
                             prior.range = c(input$sampling$hyperparameters$range, 0.5),
                             # PC-prior on sigma: P(sigma > 1) = 0.01
                             prior.sigma = c(sqrt(input$sampling$hyperparameters$sigma2), 0.5))

csdata = simulateddata$thirdstage
cssampdata = simulateddata$firststage$Samp_PPFinal
detdata = data_det_spframe
covslist <- list(cov1.spix,cov2.spix,cov3.spix)
spdeslist <- list(spdes=spdes,spde2=spde2)
covs = covslist
region=poly
mesh=mesh

data_df <- data.frame(
  Y = csdata$classifications$error,
  C = csdata$classifications$true_species,
  eco_cov = extract(cov1.rast,csdata$classifications),
  samp_cov= extract(cov2.rast,csdata$classifications),
  det_cov = extract(cov3.rast,csdata$classifications))
#Testing the compiled function.
#Should give the same results as fit.inla
CnimbleINLA <- compileNimble(nimbleINLADataGenerating)
#class_prob <- matrix(c(0.9, 0.02, 0.04, 0.04,
#                       0.05, 0.89, 0.04, 0.02,
#                       0.1,0.1, 0.8, 0,
 #                      0, 0.05, 0.25, 0.7),
#                     nrow=4, ncol=4, byrow = TRUE)

class_prob <- matrix(c(0.9, 0.1,
                       0.05, 0.95),
                     nrow=2, ncol=2, byrow = TRUE)
CnimbleINLA(class_prob)


# Running INLA within NIMBLE
cnt <- 0
listout <- list()
code <-nimbleCode({


  for(i in 1:nspecies){
    omega[i, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])

  }

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

cnt = 0
assign('myRW_dirichlet', myRW_dirichlet, envir = .GlobalEnv)
inlanim = INLAWiNimDataGenerating(data = data_df,
                                  code = code,
                     modelData = inla_data,
                     modelConstants = const,
                     modelInits = idm_inits,
                     fam = "nogaussian",
                     parametersToMonitor = c("omega"),
                     mcmcConfiguration =  list(n.chains = 1,
                                               n.iterations = 1000,
                                               n.burnin = 200,
                                               n.thin = 1,
                                               setSeed = TRUE,
                                               samples=TRUE,
                                               samplesAsCodaMCMC = TRUE,
                                               summary = TRUE,
                                               WAIC = FALSE))
save(inlanim, file = "inlaDGProcess.RData")

load("estimates.RData")
ggSamples <- ggmcmc::ggs(mcmc.out$samples)
ggmcmc::ggs_traceplot(ggSamples)+
  ggplot2::facet_wrap(~Parameter, ncol = 3)+
  ggplot2::theme_classic()
