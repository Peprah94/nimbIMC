

checkdpois <- compileNimble(dmypois)
checkdpois(x = c(3), lambda = 2, log = 1)

checkrpois <- compileNimble(rmypois)

checkrpois(n = 10, lambda = 2, includeLowerBound = 0)


saveResults = function(fixedVals = character(1),
                       res = double(2),
                       ind = integer(0)){ #which index to subset
  n <- length(fixedVals)
  vals <- numeric(n)
  #r <- character(0)
  #if(n > 1){
  for(i in seq_along(fixedVals)){
    # r <- fixedVals[i]
    vals[i] <- res[ind, i + 1]
    #model[[r]] <<- vals[i]
  }
  # }
  #for(i in seq_along(fixedVals)){
  #model[[fixedVals[i]]] <<- vals[i]
  #values(model, fixedVals[i]) <<- vals[i]
  #}
  nimble::values(model, fixedVals) <<- c(vals)
  #else{
  # vals <<- res[1, 2]
  #}
  return(vals)

}


