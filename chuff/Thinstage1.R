# First stage of thinning (sampling process)

firstthinning <- function(input, ...){

  n.species = input$constants$n.species
  seed = input$constants$seed
  cov = input$cov
  plot = input$plot
  idxs = input$idxs

environment_list <- as.list(parent.frame())
p_sampling <- length(environment_list$sampling_idxs)

sampling_form <- c()
sampling_covs.im <- environment_list$sampling_covs.im
for(j in 1:p_sampling){
  sampling_form[j] <- paste0("input$sampling$fixed.effect[[",j+1,"]]*sampling_covs.im[[",j,"]]")
}
sampling_linpred <- paste(c("input$sampling$fixed.effect[[1]]",sampling_form),collapse="+")

x0 <-  environment_list$sampling_covs.im[[1]]$xcol
y0 <-  environment_list$sampling_covs.im[[1]]$yrow
Samp_PP <- rLGCP(model="matern",mu=eval(parse(text=sampling_linpred)),
                 var=input$sampling$hyperparameters[[1]],
                 scale=input$sampling$hyperparameters[[2]]/sqrt(8),
                 nu=1,
                 win = environment_list$win,
                 xy=list(x=x0,y=y0))

Lam <- attr(Samp_PP, 'Lambda')

Samp_GRF  <- log(Lam$v)

df.sp2 <- SpatialPointsDataFrame(coords = gridlocs,
                                 data = data.frame(w=c(anti_t(rotate(rotate(log(Lam$v))))),
                                            ew=c(anti_t(rotate(rotate((Lam$v)))))))

df.sp2$retprob <- psych::logistic(df.sp2$w)
r <- raster(df.sp2)

xres <-  environment_list$sampling_covs.im[[1]]$xstep;yres <-  environment_list$sampling_covs.im[[1]]$ystep## Raster resolution
r1<-disaggregate(r, fact=res(r)/c(xres,yres))
w2.rast <- rasterize(df.sp2@coords,r1,df.sp2$w, fun=mean,na.rm=T)
#ew2.rast <- rasterize(df.sp2@coords,r1,df.sp2$ew, fun=mean,na.rm=T)
prob.rast <- rasterize(df.sp2@coords,r1,df.sp2$retprob, fun=mean,na.rm=T)

# cropping the thinning probability to the study region
w2.rastaa <- crop(w2.rast,environment_list$aa)
w2.rastaa <- mask(w2.rastaa,environment_list$aa)
prob.rastaa <- crop(prob.rast,environment_list$aa)
prob.rastaa <- mask(prob.rastaa,environment_list$aa)

## Getting probability of retaining samples ##
Eco_PPFinal <- list()

for(i in 1:n.species){
  retprobslik <- extract(prob.rast,cbind( environment_list$Eco_PP[[i]]$x, environment_list$Eco_PP[[i]]$y))
  Eco_PP1 <- SpatialPointsDataFrame(coords = cbind( environment_list$Eco_PP[[i]]$x, environment_list$Eco_PP[[i]]$y),data=data.frame(retprob=retprobslik))
  Eco_PP1$retain <- apply(Eco_PP1@data, 1,function(x){rbinom(1,1,p=x[1])})

  Eco_PPFinal[[i]] <- Eco_PP1[which(Eco_PP1$retain==1),]
}

# data from where CS sample
Samp_PPFinal <- SpatialPoints(coords=cbind(Samp_PP$x,Samp_PP$y))

return(list(retainprobraster = prob.rastaa,
            Eco_PPFinal = Eco_PPFinal,
            Samp_PPFinal = Samp_PPFinal))

}

