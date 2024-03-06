library(AHMbook)
RNGversion("3.5.3")
library(INLA)
library(inlabru)
library(sp)
library(sf)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(ggplot2)
library(gridExtra)
library(readr)
#library(terra)
library(tidyr)
library(stringr)
data("BerneseOberland")
library(nimble)
library(nimbleHMC)


#simulate data
dataSimulated <- simOccSpatial(nsurveys =  10,
                               mean.psi = 0.9,
                               beta = c(2 ,-2),
                               mean.p = 0.4,
                               alpha = c(-1, 1),
                               sample.size =500,
                               variance.RF = 1,
                               theta.RF = 10,
                               seeds = c(10, 100),
                               show.plots = FALSE)


N <- dataSimulated$nsurveys


dataForModel <- data.frame(Longitude = rep(dataSimulated$xcoord, N),
                           Latitude = rep(dataSimulated$ycoord, N),
                           elevationS = rep(dataSimulated$elevationS, N),
                           obsY = c(dataSimulated$y),
                           i_year = as.integer(factor(rep(1:N, each = length(dataSimulated$xcoord))) ))%>%
  as.matrix()

#x0 <- seq(-1.3, 4.3, length = 50)
#y0 <- seq(-1.3,4.3, length = 50)
#gridlocs <- expand.grid(x0,y0)

x <- cbind(as.matrix(dataForModel[, c(1,2,3,5)]),
           rep(standardize(c(dataSimulated$forest)), N),
           c(dataSimulated$wind))%>%
  as.matrix()
#y <- c(dataForModel[,4])
beta = c(0.4, -1, -1)
fixedVals <- c("intercept", "beta1", "beta2")
family = "binomial"


inlabruModelFit <- function(x, #matrix
                            y, #matrix
                            beta, #parameters fitted from elsewhere, and should be a vector
                            fixedVals,
                            family){

  # Simulate data
  coordsData <- x%>%
    as.data.frame()


  #convert y to vector
  y <- c(y)

  p <- plogis(beta[1] + beta[2]* x[,5] + beta[3]* x[,6])
  #p = c(beta)
  coordsData <- cbind(coordsData, p)
  colnames(coordsData) <- c("Longitude","Latitude","elevationS","i_year" , "forest", "wind" )

  coordsDataNew <- sf::st_as_sf(coordsData,
                                  coords = c("Longitude",
                                                "Latitude"))

  #datagbifNew <- terra::vect(dataSimulatedSP)
  #datagbifNew <- terra::project(datagbifNew,
  #                             "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

  #create mesh
  #coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
  #poly <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))

  #coordsDataNew <- as.data.frame(coordsData)
  #sp::coordinates(coordsDataNew) <- c("Longitude", "Latitude")


  max.edge = diff(range(sf::st_coordinates(coordsDataNew)[,1]))/(3*5)

  mesh1 = inla.mesh.2d(loc = st_coordinates(coordsDataNew),
                       max.edge = max.edge)

  ## the mesh
  #mesh <- inla.mesh.2d(loc = coordsData@coords, offset = c(0.3, 1),
  #                  max.edge = c(0.1, 0.5), cutoff = 0.2)

  #SPDE with exponential correlation function
  spde = inla.spde2.matern(mesh = mesh1,
                             alpha = 1.5)

  #spde <- inla.spde2.matern(
  # mesh=mesh, alpha=1.5)

  elev.spix <- SpatialPixelsDataFrame(point =  st_coordinates(coordsDataNew),
                                      data = data.frame(elev = coordsDataNew$elevationS))

  #elev.spix$elev <- AHMbook::standardize(elev.spix$elev)

  elevsq.spix <- SpatialPixelsDataFrame(point =  st_coordinates(coordsDataNew),
                                        data = data.frame(elev = (coordsDataNew$elevationS)^2))
  #elevsq.spix$elevsq <- AHMbook::standardize(elevsq.spix$elevsq)

  cmp <- as.formula("obsY~ - 1 + beta0(1) + elev(main = elev.spix, model = 'linear') + elevsq(main = elevsq.spix, model = 'linear')+ site(main= i_year, model = 'iid', n = N) + w2(main = coordinates, model = spde)")

  coordsDataNew1 <- as.data.frame(coordsData)

  coordsDataNew1 <- sp::SpatialPointsDataFrame(coords = coordsDataNew1[, c("Longitude", "Latitude")],
                                               data = coordsDataNew1)
  #sp::coordinates(coordsDataNew1) <- c("Longitude", "Latitude")
  #coordsDataNew1 <- as(coordsDataNew1, "SpatialPointsDataFrame")
  coordsDataNew1$obsY <- y

  functionConstants <- function(p, intercept, elev, elevsq, visits, w){
    linearPred <- plogis(intercept + elev + elevsq + visits + w)
    p <- abs(p-0.0001)
    detProb <- qlogis(p)
    firstTerm <- log( 1- linearPred*p)
    secondTerm <- log(1 - linearPred)
    thirdTerm <- log(1 - p)
    ret <- detProb - firstTerm + secondTerm + thirdTerm
  }

  lik1 <- inlabru::like(family,
                        formula = as.formula(paste0("obsY ~ beta0 + elev + elevsq + site +w2 + functionConstants(p, beta0, elev, elevsq, site, w2)")),
                        data = coordsDataNew1,
                        #components = cmp,
                        domain = list(coordinates = mesh1)#,
                        #samplers = poly
  )



  m_bru <- inlabru::bru(cmp, lik1,
                        #family = "binomial",
                        #data = coordsData,
                        options =
                          list(
                            bru_verbose = TRUE,
                            bru_max_iter=2,
                            control.fixed = list(expand.factor.strategy = "inla", mean = 0, prec = 1 / (100 * 100)),
                            control.family = list(link = "logit"),
                            control.inla = list(int.strategy = "eb", strategy = "gaussian"),
                            control.compute=list(return.marginals.predictor=TRUE)
                          )
  )

  #indx <- mesh1$idx$loc
 # ret <- matrix(m_bru$summary.fitted.values[indx,"mean"], nrow = length(coordsData$Longitude)/N, ncol = N, byrow = FALSE)
  fittedValues = c(m_bru$mlik[1,1])
  intercept = INLA::inla.emarginal(function(x) x,m_bru$marginals.fixed[[1]])
  elev = INLA::inla.emarginal(function(x) x,m_bru$marginals.fixed[[2]])
  elevsq = INLA::inla.emarginal(function(x) x,m_bru$marginals.fixed[[3]])
  siteSD = INLA::inla.emarginal(function(x) x,1/m_bru$marginals.hyper[[1]])
  theta1 = INLA::inla.emarginal(function(x) x,1/m_bru$marginals.hyper[[2]])
  theta2 = INLA::inla.emarginal(function(x) x,1/m_bru$marginals.hyper[[3]])

  ret1 <- cbind(fittedValues, intercept, elev, elevsq, siteSD, theta1, theta2)
  colnames(ret1) <- c("mld", fixedVals, "siteSD", "theta1", "theta2")

  #ret <- as.matrix(ret)
  return(ret1)

}

nimbleINLA <- nimble::nimbleRcall(
  prototype = function(
    x=double(2), #x is a matrix
    y=integer(2), #y is a matrix
    beta=double(1), # beta is a vector
    fixedVals = character(1, default = c("intercept", "beta1", "beta2")),
    #interInModel = double(0, default = 1),
    family = character(0, default = "binomial")
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'inlabruModelFit'
)

Cmwtc <- nimble::compileNimble(nimbleINLA)
Cmwtc(x, y = dataSimulated$y, beta = beta)


#dist(c(dataForModel[,1:2]))


code <-nimbleCode({
  # Specify priors
  #beta0 <- logit(mean.psi)
  #mean.psi <- logit(beta0)
 # alpha0 <- logit(mean.p)
  alpha0 ~ dunif(0, 1)
  beta0 ~ dnorm(0, 001)
  for(v in 1:2){
    alpha[v] ~ dnorm(0, 0.001)
    beta[v] ~ dnorm(0, 0.001)
  }

  for(site.tag in 1:nsites){
    for(visit.tag in 1:nvisits){
      logit(p[site.tag, visit.tag]) <- alpha0 + alpha[1]*forest[site.tag] + alpha[2]*wind[site.tag, visit.tag]
    }
  }


  #site effects
  #eta[1:nsites] <- inlaSimCorr(theta)

###### Occupancy model
  for(site.tag in 1:nsites){
   # for(visit.tag in 1:nvisits){
      logit(psi[site.tag]) <- beta0 + beta[1]*elev[site.tag] + beta[2]*elevsq[site.tag] #+ eta[site.tag]
   # }
  }


  # Observation model
  for(site.tag in 1:nsites){
    for(visit.tag in 1:nvisits){
      y[site.tag, visit.tag] ~ dbin(size = 1, prob = psi[site.tag]*p[site.tag, visit.tag])
    }
  }

  #


})


## Parameterising the nimble model

#Data
inla_data <- list(y = dataSimulated$y,
                  forest = standardize(dataSimulated$forest),
                  wind = dataSimulated$wind,
                  elev = dataSimulated$elevationS,
                  elevsq = (dataSimulated$elevationS)^2)

#Constants
const <- list(N = dataSimulated$nsurveys,
              nvisits = dataSimulated$nsurveys,
              nsites = length(dataSimulated$xcoord)
)
#zst <- apply(inla_data$y, 1, max)
#zst[is.na(zst)] <- 1
# Initial values
idm_inits <- function(){list(p = matrix(runif(const$N * const$nsites), nrow = const$nsites, const$N),
                             alpha = c(-1,1),
                             mean.p = 0.4,
                             beta = c(1, -1),
                             alpha0 = 0.8,
                             psi = runif(const$nsites)
)
}

initsList <- idm_inits()









ret <- INLAWiNimDataGenerating(data = c("y"),
                               covariate = x,
                               code = code,
                               family = "binomial",
                               modelData = inla_data,
                               modelConstants = const,
                               modelInits = idm_inits,
                               nimbleINLA = nimbleINLA,
                               inlaMCMC = c("inlamcmc"),
                               inlaMCsampler = "AFSS_INLA_block",
                               samplerControl = list(),
                               parametersToMonitor = list(mcmc = c("alpha0","alpha"),
                                                          inla = c("beta0", "beta")),
                               mcmcConfiguration = list(n.chains = 1,
                                                        n.iterations = 10,
                                                        n.burnin = 0,
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
# #
# # #Create the model in nimble
mwtc <- nimbleModel(code,
                    data = inla_data,
                    constants = const,
                    inits = initsList)

#
Cmwtc <- nimble::compileNimble(mwtc,
                               showCompilerOutput = FALSE) #Have issues compiling
#
mcmcconf <- nimble::configureMCMC(Cmwtc,
                                  nodes = NULL)
#
# #mcmcconf$removeSampler(c("beta", "a", "sigma"))
fixedVals = c("beta0", "beta")
obs = c("y")
mcmcconf$addSampler(target = c("alpha0","alpha"),
                    type = "RW_INLA_block",
                    control = list(fit.inla = nimbleINLA,
                                   x = x,
                                   y = obs,
                                   #obsVar = c(dataSimulated$y),
                                   fixedVals = fixedVals,
                                   fam = "binomial",
                                   scale = 1))
#
# mcmcconf$printSamplers()
Rmcmc <- nimble::buildMCMC(mcmcconf)
#
# # Compile
cmcmc <- nimble::compileNimble(Rmcmc,
                               project = Cmwtc,
                               resetFunctions = TRUE)
#
# #cmcmc$run(1000)
#
# cmc.out <- nimble::runMCMC(cmcmc,
#                            niter = 5,
#                            nchains = 2,
#                            nburnin = 1,
#                            #inits = initsList,
#                            thin = 1,
#                            setSeed = TRUE,
#                            samples = TRUE,
#                            samplesAsCodaMCMC = TRUE,
#                            summary = TRUE,
#                            WAIC = FALSE)
