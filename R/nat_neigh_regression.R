#' Natural neighbors regression 
#' 
#' Natural neighbors regression 
#' @param x Feature matrix
#' @param y Regression values
#' @return A nn_regeression object with the following fields:
#' @export
natural_regression=function (
  x,
  y,
  order=1,
  l2=.0001,
  verbose=1
  ) {
  if (class(y)=="integer") {
    y=data.frame(as.numeric(y))
  }
  else if (class(y)=="numeric") {
    y=data.frame(y)
  }
  else if (class(y)=="matrix") {
    y=as.data.frame(y)
  }
  # if (class(x)=="matrix") {
  #   x=as.data.frame(x)
  #   cat("X matrix cast to data.frame.\n")
  # }
  # if (class(y)=="data.frame") {
  #   if (length(y)>1)
  #     cat("Only single target variables implemented. Taking first column.")
  #   y=y[,1]
  # }
  # else if (class(y)=="matrix") {
  #   if (ncol(y)>1) {
  #     cat("Only single target variables implemented. Taking first column.")
  #     y=y[,1]
  #   }
  # }    
  # if (class(y)=="data.frame" || class(y)=="matrix" ) {
  #   if (ncol(y)>1)
  #     cat("Only single target variables implemented. Taking first column.")
  #   y=y[,1]
  # }

      
  dup_forward=duplicated(x)
  if (anyDuplicated(x)) {
    warning("Duplicated feature vectors. Target variables are being aggregated (mean).")
    df=aggregate(y,x,mean)
    x=df[1:ncol(x)]
    y=df[(ncol(x)+1):ncol(df)]
  }
  
  feature_types=sapply(x,class)
  factors=which(feature_types=="factor")
  numerics=which(feature_types!="factor")
  factor_manager=get_factor_manager(x[factors],factors)
  if (length(factors)>0) {
    matches=lapply(1:factor_manager$lim,function(i){
      levels=factor_manager$get_levels(i)
      apply(x[factors],1,function(row)all(row==levels))
    })
    class_x=lapply(1:factor_manager$lim,function(i){
      as.matrix(x[matches[[i]],numerics])
    })
    class_y=lapply(1:factor_manager$lim,function(i){
      y[matches[[i]],,drop=FALSE]
    })
  }
  else {
    class_x=list(as.matrix(x))
    class_y=list(y)
  }
  if (verbose>0)
    cat("Data parsed.\n")
  
  if (length(numerics)==0)
    return (natural_regression_only_class_regressors(
      x,
      y,
      factor_manager,
      matches,
      class_x,
      class_y,
      factors,
      numerics,
      l2,
      verbose
    ))
  else if (length(numerics)==1)
    return (natural_regression_univariate_numeric_regressors(
      x,
      y,
      factor_manager,
      matches,
      class_x,
      class_y,
      factors,
      numerics,
      l2,
      verbose
    ))
  else 
    return (natural_regression_multivariate_numeric_regressors(
      x,
      y,
      factor_manager,
      matches,
      class_x,
      class_y,
      factors,
      numerics,
      l2,
      verbose
    ))
}
natural_regression_multivariate_numeric_regressors=function(
  x,
  y,
  factor_manager,
  matches,
  class_x,
  class_y,
  factors,
  numerics,
  l2,
  verbose
){
  # For each class
  # 1. Create tessellation
  # 2. Calculate convex hulls
  # 3. Calculate regression coefficients
  
  # Create tessellation for each class
  tes=lapply(1:factor_manager$lim,function(i) geometry::delaunayn(class_x[[i]],full=T))
  if (verbose>0)
    cat("Tessellations generated.\n")

  # Find convex hull points for each class
  chulls=lapply(1:factor_manager$lim,function (i) geometry::convhulln(class_x[[i]]))
  if (verbose>0)
    cat("Convex hulls calculated.\n")
  
  # Create linear model coefficients for each simplex for each class
  simplex_models=lapply(1:factor_manager$lim,function (i)
    apply(tes[[i]]$tri,1,function(pnts) {
      X=cbind(rep(1,length(pnts)),class_x[[i]][pnts,])
      Y=class_y[[i]][pnts,]
      solve(crossprod(X)+diag(l2,ncol(X)))%*%crossprod(X,Y)
    } ))
  if (verbose>0)
    cat("Density coefficients for simplexes calculated.\n")
  
  # Create return object
  out=list(
    factors=factors,
    numerics=numerics,
    tes=tes,
    class_x=class_x,
    class_y=class_y,
    dim=ncol(x),
    chulls=chulls,
    factor_manager=factor_manager,
    simplex_models=simplex_models
  )
  class(out)="natural_regression_mnr"
  return (out)
}
natural_regression_univariate_numeric_regressors=function(
  x,
  y,
  factor_manager,
  matches,
  class_x,
  class_y,
  factors,
  numerics,
  l2,
  verbose
) {
  # Create 'tessellation' for each class
  tes=lapply(1:factor_manager$lim,function(i) order(class_x[[i]][,1]))
  if (verbose>0)
    cat("Tessellations generated.\n")
  
  # Find convex hull points for each class
  chulls=lapply(1:factor_manager$lim,function (i) c(min(class_x[[i]]),max(class_x[[i]])))
  if (verbose>0)
    cat("Convex hulls calculated.\n")
  
  # Create linear model coefficients for each simplex for each class
  simplex_models=lapply(1:factor_manager$lim,function(i) {
    x_ordered=class_x[[i]][tes[[i]],,drop=FALSE]
    y_ordered=as.matrix(class_y[[i]][tes[[i]],,drop=FALSE])
    sapply(1:(nrow(x_ordered)-1),function(j){
      X=cbind(rep(1,2),x_ordered[j:(j+1),,drop=FALSE])
      Y=y_ordered[j:(j+1),,drop=FALSE]
      solve(crossprod(X)+diag(l2,ncol(X)))%*%crossprod(X,Y)
    })
  })
  if (verbose>0)
    cat("Density coefficients for simplexes calculated.\n")
  
  # Create return object
  out=list(
    factors=factors,
    numerics=numerics,
    tes=tes,
    class_x=class_x,
    class_y=class_y,
    dim=ncol(x),
    chulls=chulls,
    factor_manager=factor_manager,
    simplex_models=simplex_models
  )
  class(out)="natural_regression_unr"
  return (out)
  
}
get_factor_manager=function(factors,factor_indices) {
  levels=lapply(factors,levels)
  if (length(levels)==0){
    lim=1
    reps=NA
    get_levels=function(i) {return (NULL)}
    get_row=function(f) {return (1)}
  }
  else {
    lim=prod(sapply(levels,length))
    reps=rep(1,length(levels))
    reps[length(levels)]=1
    if (length(levels)>1) {
      for (i in length(levels):2) {
        reps[i-1]=reps[i]*length(levels[i])
      }
    }
    get_levels=function(i) {
      out=rep(NA,length(levels))
      for (j in 1:length(levels)) {
        out[j]=levels[[j]][i%/%reps[j]]
      }
      return (out)
    }     
    get_row=function(features){
      out=0
      for (i in 1:length(factor_indices)) {
        out=out+which(levels[[i]]==features[factor_indices[i]])*reps[i]
      }
      return (out)
    }
  }

  return (list(
    levels=levels,
    reps=reps,
    lim=lim,
    get_levels=get_levels,
    get_row=get_row
  ))
}