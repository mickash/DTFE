#' Get a test data set.
#' 
#' Get a test data set with three classes, all drawn from normal distributions. The last 
#' can be randomly multiplied by 1 or -1 using scatterC.
#' @param An The number of class A 
#' @param Bn The number of class B 
#' @param Cn The number of class C 
#' @param Am1 The mean of X1 for class A
#' @param Am2 The mean of X2 for class A
#' @param Bm1 The mean of X1 for class B
#' @param Bm2 The mean of X2 for class B
#' @param Cm1 The mean of X1 for class C
#' @param Cm2 The mean of X2 for class C
#' @param As1 The std dev of X1 for class A
#' @param As2 The std dev of X2 for class A
#' @param Bs1 The std dev of X1 for class B
#' @param Bs2 The std dev of X2 for class B
#' @param Cs1 The std dev of X1 for class C
#' @param Cs2 The std dev of X2 for class C
#' @param scatterC Whether to 'scatter' class C - see function description.
#' @return A data matrix
#' @export
getTestData=function(
  An=40,
  Bn=20,
  Cn=80,
  Am1=10,
  Am2=0,
  Bm1=0,
  Bm2=0,
  Cm1=15,
  Cm2=15,
  As1=5,
  As2=As1,
  Bs1=5,
  Bs2=Bs1,
  Cs1=5,
  Cs2=Cs1,
  scatterC=T) {
  out=  rbind(
      cbind(rnorm(An,Am1,As1)+10,rnorm(An,Am2,As2)+10,rep(1,An)),
      cbind(rnorm(Bn,Bm1,Bs1),rnorm(Bn,Bm2,Bs2),rep(2,Bn)),
      cbind(rnorm(Cn,Cm1,Cs1)*sample(c(1,-1),Cn,replace=T),rnorm(Cn,Cm2,Cs2)*sample(c(1,-1),Cn,replace=T),rep(3,Cn))
  )
  if (scatterC && Cn>0)
      out[(An+Bn+1):(An+Bn+Cn),1:2]=out[(An+Bn+1):(An+Bn+Cn),1:2]*sample(c(1,-1),Cn*2,replace=T)
  return (out)
}
#' Get a small set of data
#' 
#' Get a small set of data. This is a wrapper for getTestData using the first three arguments.
#' @param a An
#' @param b Bn
#' @param c Cn
#' @return A data matrix
#' @export
getBasicData=function(a=4,b=6,c=8) getTestData(a,b,c)

#' Get target like data (one class inside the other)
#' 
#' Get target like data (one class inside the other)
#' @param a The number of class A
#' @param b The number of class B
#' @param dim The dimensionality of the data (2 or 3)
#' @return A data matrix
#' @export
getTargetData=function(a=40,b=40,dim=2) {
  a_=3*a
  class1=NULL
  class2=NULL
  dist1=NULL
  dist2=NULL
  if (dim==2) {
    class1=cbind(rnorm(a_,0,4),rnorm(a_,0,4),rep(1,a_))
    dist1=sqrt(class1[,1]^2 + class1[,2]^2)    
  } else if (dim==3) {
    class1=cbind(abs(rnorm(a_,0,4)),abs(rnorm(a_,0,4)),abs(rnorm(a_,0,4)),rep(1,a_))
    dist1=sqrt(class1[,1]^2 + class1[,2]^2 + class1[,3]^2)    
  } else {
    stop("Dimension of target must be 2 or 3.")
  }
  class1=class1[which(dist1<4|dist1>5),]
  class1=class1[1:min(a,nrow(class1)),]
  
  b_=5*b
  if (dim==2) {
    class2=cbind(rnorm(b_,0,4),rnorm(b_,0,4),rep(2,b_))
    dist2=sqrt(class2[,1]^2 + class2[,2]^2)    
  } else if (dim==3) {
    class2=cbind(abs(rnorm(b_,0,4)),abs(rnorm(b_,0,4)),abs(rnorm(b_,0,4)),rep(2,b_))
    dist2=sqrt(class2[,1]^2 + class2[,2]^2 + class2[,3]^2)    
  } else {
    stop("Dimension of target must be 2 or 3.")
  }
  class2=class2[which(dist2>4&dist2<5),]
  class2=class2[1:min(b,nrow(class2)),]
  
  rbind(class1,class2)
}