## Third stage of thinning (misclassification) ##

thirdthinning <- function(input, ...){

  n.species = input$constants$n.species
  seed = input$constants$seed
  cov = input$cov
  plot = input$plot
  idxs = input$idxs

  environment_list <- as.list(parent.frame())
  if(is.matrix(input$misclassification[[1]])){

    class_prob <- input$misclassification[[1]]

  if(dim(class_prob)[1]!=environment_list$n.species | dim(class_prob)[2]!=environment_list$n.species) stop(paste0("Dimensions of classification probabilities should be n.species= ", environment_list$n.species," x n.species= ",environment_list$n.species))

    if(any(class_prob<0) | any(class_prob>1)) stop("Entry of classification probability must be between 0 and 1")

  if(!all(rowSums(class_prob)==1)) stop("Sum of rows of classification probability must be one")



  for(i in 1:n.species){
    if(length(environment_list$secondstage$Eco_PPFinal_detect[[i]])>0){
    environment_list$secondstage$Eco_PPFinal_detect[[i]]$error <- apply(environment_list$secondstage$Eco_PPFinal_detect[[i]]@data, 1,function(x){which(rmultinom(1,1,p=(class_prob[i,]))==1)})
    environment_list$secondstage$Eco_PPFinal_detect[[i]]$true_species <- i}
  }

  #Putting all the data together
  cit_data <- list()
  for(i in 1:n.species){
    cit_data[[i]] <- do.call("rbind",lapply(environment_list$secondstage$Eco_PPFinal_detect, function(x) x[which(x@data$error ==i),]))
  }

  allcit_data <- do.call("rbind",cit_data)

  }

  else{
    stop("Classification probabilities must be a matrix")
  }

  return(list(classifications = allcit_data))

}
