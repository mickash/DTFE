
divergence = function (tes1,tes2,samples=10000) {
  if (!identical(tes1$domain,tes2$domain)) stop ("Domain mismatch.")
  getProbFrom=makeProbFromTessellationFunction(tes1)
  getProbTo=makeProbFromTessellationFunction(tes2)
  sampleFromDomain=makeSampleFromDomainFunction(tes1$domain)
  sampleFromClasses=makeSampleFromClassesFunction(length(tes1$class))
  estimateKLD(getProbFrom,getProbTo,sampleFromDomain,sampleFromClasses,samples)
}

# Find KL - Divergence
estimateKLD=function(getProbFrom,getProbTo,sampleFromDomain,sampleFromClasses,samples=10000) {
  sum(sapply(1:samples,function(i) {
      X=sampleFromDomain()
      Y=sampleFromClasses()
      pFrom=getProbFrom(X,Y)
      pTo=getProbTo(X,Y)
      pFrom*log(pFrom/pTo)
    })) / samples
}

#' Make sampleFromDomain function for estimating KLD
#' 
#' Make sampleFromDomain function for estimating KLD
#' @param domain A 2 rows matrix specifying domain using mins (row 1) and maxes (row 2)
#' @return sampleFromDomain function
#' @keywords internal

makeSampleFromDomainFunction=function(domain)
  function() sapply(1:ncol(domain),function(i) runif(1,domain[1,i],domain[2,i]))
makeSampleFromClassesFunction=function(classes) function () sample(1:classes,1)
makeProbFromTessellationFunction=function(tes) {
  function (X,Y) {
    
    # 1. Work out probability of class
    # Number of class / total number of data points
    p=length(tes$classdata[[Y]]/length(tes$x))
    
    # 2. Work out probability of being in simplex
    # Volume of simpelx / total volumes of simplexes
    tri=geometry::tsearchn(tes$classdata[[Y]],tes$tes[[Y]]$tri,X)$idx
    p=p*tes$simplex_cond_probabilities[[Y]][tri]
    
    # 3. Estimate probability of value given that it is in the simplex.
    # This is the density value at point / volume of simplex
    #   - Divide by volume gives a function that integrates to 1 - a pdf
    p*((tes$simplex_models[[Y]][[tri]]%*%c(1,X))/tes$hypervolumes[[Y]][tri])
    
  }
}
