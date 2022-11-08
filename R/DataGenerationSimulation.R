#' Generate Citizen Science data
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param sppixels
#'
#' @return MCMC output
#' @export
generateCSData <- function(input,
                   domain=NULL,
                   colmatrix=NULL){
source("R/metaFunctionsDataGenerating.R")
  source("R/Thinstage1.R")
  source("R/Thinstage2.R")
  source("R/Thinstage3.R")
#Retrieving necessary inputs
n.species = input$constants$n.species
seed = input$constants$seed
cov = input$cov
plot = input$plot
idxs = input$idxs

  if(class(n.species)!="numeric") stop("Number of species must be a number")

  #Setting seed for the simulation
  set.seed(seed)
  RandomFields::RFoptions(seed=seed)

  # The domain where the simulations are done
  # The default is provided here
  if(is.null(domain)){
    coordsmat <- matrix(c(0,0,3,0,3,3,0,3,0,0),ncol=2,byrow=T)
    aa <- SpatialPolygons(list(Polygons(list(Polygon(coordsmat)),ID=1)))
    win <- as.owin(aa) ##Need maptools
  }

  ## Mesh for models ##
  ## Maximum distance in the extent of study area ##
  ext <- raster::extent(aa)
  pts <- SpatialPoints(coords = matrix(ext,ncol=2),
                       proj4string = crs(aa))
  max.dist <- spDists(pts)[1,2]
  mesh <- inla.mesh.2d(loc.domain = coordsmat,
                       max.edge = c(0.02*max.dist,0.10*max.dist),
                       offset = c(0.3, 1),cutoff = 0.2*max.dist)

  ## Converting the covariates provided into a raster, image and pixels ##
  if(is.list(cov)){

    classes <- unique(sapply(cov, class))
    if(length(classes)==1){

      if(classes=="im"){
        covs.im <- cov
        covs.raster <- lapply(cov,im2rast)
        covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
      }
      else{if(classes=="RasterLayer"){
        covs.im <- lapply(cov,as.im)
        covs.raster <- cov
        covs.sppixels <- lapply(covs.raster,function(x){as(x,"SpatialPixelsDataFrame")})
      }
        else{if(classes=="SpatialPixelsDataFrame"){
          covs.im <- lapply(cov,as.im)
          covs.raster <- lapply(cov, sppixels2raster)
          covs.sppixels <- cov
        }
          else{
            stop("Covariates must be of 'im', 'RasterLayer' or 'SpatialPixelsDataFrame'")
          }
        }}

    }
    else{
      stop("All the covariate must be in the same format")
    }

    ## Indexes of the covariates##
    eco_idxs <- idxs$eco
    sampling_idxs <- idxs$sampling
    detection_idxs <- idxs$detection

    eco_covs.im <- covs.im[eco_idxs]
    eco_covs.raster <- covs.raster[eco_idxs]
    eco_covs.sppixels <- covs.sppixels[eco_idxs]

    sampling_covs.im <- covs.im[sampling_idxs]
    sampling_covs.raster <- covs.raster[sampling_idxs]
    sampling_covs.sppixels <- covs.sppixels[sampling_idxs]

    detection_covs.im <- covs.im[detection_idxs]
    detection_covs.raster <- covs.raster[detection_idxs]
    detection_covs.sppixels <- covs.sppixels[detection_idxs]

    if(!is.null(input$ecological$fixed.effect$scale)){eco_scale=input$ecological$fixed.effect$scale} ##Test
    else{eco_scale=1} ##Test
    if(!is.null(input$sampling$fixed.effect$scale)){sampling_scale=input$sampling$fixed.effect$scale}
    else{sampling_scale=1}
    if(!is.null(input$detection$fixed.effect$scale)){detection_scale=input$detection$fixed.effect$scale}
    else{detection_scale=1}
    ## Generate ecological process ##

    #Ecological Process
    Eco_PP <- list()

    ## Bringing ecological covariates ##
    p_eco <- length(eco_idxs)
    eco_form <- c()
    for(j in 1:p_eco){
      eco_form[j] <- paste0("input$ecological$fixed.effect[[",j+1,"]][i]*eco_covs.im[[",j,"]]")
    }

    eco_linpred <- paste(c("input$ecological$fixed.effect[[1]][i]",eco_form),collapse="+")

    for(i in 1:n.species){
      x0 <- eco_covs.im[[1]]$xcol
      y0 <- eco_covs.im[[1]]$yrow
      Eco_PP[[i]] <- rLGCP(model="matern",
                           mu=eval(parse(text=eco_linpred)),
                           var=input$ecological$hyperparameters[[1]][i],
                           scale=input$ecological$hyperparameters[[2]][i]/sqrt(8),
                           nu=1,
                           win = win,
                           xy=list(x=x0,y=y0))

    }

    # Storing the raster of the true intensity for each species
    # Storing the raster of the true intensity for each species
    species_rast <- list()
    Scalelambda_test_rast <- list()
    for(i in 1:n.species){
      Lam <- attr(Eco_PP[[i]], 'Lambda')
      Eco_GRF  <- log(Lam$v)
      Lambda_test <- eco_scale*Lam$v ##Test
      gridlocs <- expand.grid(x0,y0) ##  Matching resolutions between gridlocs and covariates
      df.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v))))),
                                                                          scalelambda =c(anti_t(rotate(rotate(Lambda_test))))
      )) ##Test
      r <- raster(df.sp)
      #r1<-disaggregate(r, fact=res(r)/c(0.056,0.056))
      xres <-  eco_covs.im[[1]]$xstep;yres <- eco_covs.im[[1]]$ystep## Raster resolution
      r1<-disaggregate(r, fact=res(r)/c(xres,yres))


      w1.rast <- rasterize(df.sp@coords,r1,df.sp$w, fun=mean,na.rm=T)
      w1.rastaa <- crop(w1.rast,aa)
      scalelambda.rast <- rasterize(df.sp@coords,r1,df.sp$scalelambda, fun=mean,na.rm=T) ##Test
      scalelambda.rastaa <- crop(scalelambda.rast,aa) ##Test

      species_rast[[i]] <- mask(w1.rastaa,aa)
      Scalelambda_test_rast[[i]] <- mask(scalelambda.rastaa,aa) ##Test

    }

    #environment_list <- as.list(environment())
    #### First thinning stage ##
    #print(length(environment_list))

    firststage <- firstthinning(input)

    ## Second thinning stage ##
    environment_list <- as.list(environment())
    secondstage <- secondthinning(input,environment_list)

    ## Second thinning stage ##
    environment_list <- as.list(environment())
    thirdstage <- thirdthinning(input,environment_list)


    #Add the plot later
    }
  else{
    stop("Covariates input should be a list")
  }

  ## Generating lambda_obs ##
  lambda_obs_raster <- list() ##Test
  for(i in 1:n.species){ ##Test
    lambda_obs_raster[[i]] <- Scalelambda_test_rast[[i]]*(firststage$retainprobraster)*(secondstage$detectionprobraster[[i]])
  }

  return(list(trueecological=Eco_PP,
              firststage=firststage,
              secondstage=secondstage,thirdstage=thirdstage,
              species_raster=species_rast,
              lambda_obs_raster=lambda_obs_raster))

}
