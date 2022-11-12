#' Fitting INLA within NIMBLE
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param data, code, family,n.iterations, n.chains, n.burnin
#'
#' @return MCMC output
#' @export

rotate <- function(x) (apply(t(x), 2, rev))

anti_t <- function (m){
  p <- nrow(m)
  j <- matrix(ncol = p, nrow = p, data = 0)
  for (i in 1:p) {
    j[i, p - i + 1] <- 1
  }
  return(j %*% t(m) %*% j)
}

im2rast <- function(im){
  gridlocs <- expand.grid(im$xcol,im$yrow)
  cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(im$v))))))
  r <- raster(cov1.sp)
  xstp <- im$xstep
  ystp <- im$ystep
  r1<-disaggregate(r, fact=res(r)/c(xstp,ystp))
  cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
}

#' Convert pixels to raster
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param sppixels
#'
#' @return MCMC output
#' @export
sppixels2raster <- function(sppixels){
  cov1.sp <- SpatialPointsDataFrame(coords = sppixels@coords,
                                    data = sppixels@data,
                                    proj4string = crs(sppixels))
  r <- raster(cov1.sp)
  xs <- unique(sppixels@coords[,1])
  ys <- unique(sppixels@coords[,2])
  xstp <- xs[2] - xs[1]
  ystp <- ys[2] - ys[1]
  r1<-disaggregate(r, fact=res(r)/c(xstp,ystp))
  cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
}

#' Convert image to raster
#'
#' This function sets the paramters in the appropriate manner to be used by the
#' simulation function
#'
#' @param sppixels
#'
#' @return MCMC output
#' @export
im2rast <- function(im){
  gridlocs <- expand.grid(im$xcol,im$yrow)
  cov1.sp <- SpatialPointsDataFrame(coords = gridlocs,data = data.frame(cov=c(anti_t(rotate(rotate(im$v))))))
  r <- raster(cov1.sp)
  xstp <- im$xstep
  ystp <- im$ystep
  r1<-disaggregate(r, fact=res(r)/c(xstp,ystp))
  cov1.rast <- rasterize(cov1.sp@coords,r1,cov1.sp$cov, fun=mean,na.rm=T)
}



