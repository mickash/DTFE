#' Get the hyperplane coefficients for a simplex face.
#' 
#' Get the hyperplane coefficients for a simplex face.
#' @param points The vertex points for the simplex face
#' @return A vector of coefficients for the hyperplane
getHPCoefs=function(points) {
  if (qr(points)$rank<ncol(points)) {
    # We have a single linear dependence. We must find which pair of columns it is.
    for (i in 1:(ncol(points)-1)) {
      for (j in (i+1):ncol(points)) {
        if (qr(points[,c(i,j)])$rank==1) {
          px=cbind(rep(1,nrow(points)),points[,-j])
          coef_=solve(crossprod(px))%*%crossprod(px,points[,j])          
          # j != 2 since it is bigger than i. The next line gives us alphas up to j (since beta was 1).
          coef=c(coef_[1:j],-1)    
          if (j!=ncol(points)) coef=c(alphas,coef_[(j+1):length(coef_)])
          return (coef)
        }
      }
    }
    stop ("Could not find linear dependence!")
  }
  else {
    return (c(1,solve(crossprod(points))%*%crossprod(points,rep(-1,nrow(points)))))
  }
}

makeInsideFunction=function(points,internal) {
  coefs_=getHPCoefs(points)
  alpha=coefs_[1]
  coefs=coefs_[2:length(coefs_)]
  pos=(crossprod(coefs,other)+alpha)>0
  f=function(p,closed) {
    e=crossprod(coefs,p)+alpha
    on=isTRUE(all.equal(e[1,1],0))
    inside=(e[1,1]>0)==pos
    return ( (closed && on) || inside)
  }
  return (f)  
}
