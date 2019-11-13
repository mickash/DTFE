#' Demo
#' 
#' Basic visual demonstration of density plotting, probability estimation and 
#' density estimation in 3d.
#' @export
demo=function(seed=0) {
  if (!is.na(seed)) {
    set.seed(seed)
  }
  d=getBasicData(20,20,20)
  x=d[,1:2]
  y=d[,3]
  tes=tessellations(x,y,0)
  plot(tes)
  pnts=matrix(c(0,0,20,10,10,-10),byrow=T,nrow=3)
  print("Probabilities at points: [0,0],[20,10],[10,-10]")
  print(probability(tes,pnts))
  print("Densities at points: [0,0],[20,10],[10,-10]")
  print(potential(tes,pnts))
  
}
#nn=natural_regression(x,y,order=order,verbose=3)
