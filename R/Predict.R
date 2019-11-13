#' Predict
#' 
#' Predict regression estimate at point.
#' @param model The natural neighbors regression model
#' @param X The feature matrix or vector to regress on
#' @return Numeric giving regression estimate
#' @export
predict.natural_regression_mnr=function(model,X) {
  # inner_predict=function(feature_vector){
  #   index=model$factor_manager$get_row(x)
  #   simplex=geometry::tsearchn(model$class_x[[index]],model$tes[[index]]$tri,matrix(feature_vector[,model$numerics],nrow=1))$idx
  #   if (is.na(simplex))
  #     NA
  #   else
  #     model$simplex_models[[index]][,simplex]%*%c(1,x[model$numerics])
  # }
  
  if (class(X)=="matrix" || class(X)=="data.frame") {
    apply(X,1,predict_one.natural_regression_mnr,model)
  }
  else {
    if (model$dim==1){
      sapply(X,predict_one.natural_regression_mnr,model)
    }
    else {
      predict_one.natural_regression_mnr(x,model)
    }
  }
}
predict_one.natural_regression_mnr=function(feature_vector,model){
  index=model$factor_manager$get_row(feature_vector)
  numerics=as.numeric(feature_vector[model$numerics])
  simplex=geometry::tsearchn(model$class_x[[index]],model$tes[[index]]$tri,matrix(numerics,nrow=1))$idx
  if (is.na(simplex))
    NA
  else
    model$simplex_models[[index]][,simplex]%*%c(1,numerics)
}
#' @export
predict.natural_regression_unr=function(model,X) {
  # inner_predict=function(feature_vector){
  #   index=model$factor_manager$get_row(x)
  #   simplex=geometry::tsearchn(model$class_x[[index]],model$tes[[index]]$tri,matrix(feature_vector[,model$numerics],nrow=1))$idx
  #   if (is.na(simplex))
  #     NA
  #   else
  #     model$simplex_models[[index]][,simplex]%*%c(1,x[model$numerics])
  # }
  
  if (class(X)=="matrix" || class(X)=="data.frame") {
    apply(X,1,predict_one.natural_regression_unr,model)
  }
  else {
    sapply(X,predict_one.natural_regression_unr,model)
  }
}
predict_one.natural_regression_unr=function(feature_vector,model){
  index=model$factor_manager$get_row(feature_vector)
  numerics=as.numeric(feature_vector[model$numerics])
  if (length(numerics)!=1) {
    stop("Invalid input feature for univarite numeric regressor natural regression model.")
  }
  # Check if outside the convex hull
  if (numerics[1]<model$chulls[[index]][1] || numerics[1]>model$chulls[[index]][2]){
    NA
  }
  else {
    greatest_less_than_index=max(which(model$class_x[[index]][model$tes[[index]],1]<=numerics[1]))
    if (greatest_less_than_index==length(model$class_x[[index]][model$tes[[index]],1]))
      greatest_less_than_index=greatest_less_than_index-1 # Case where it is on the max convex hull.
    model$simplex_models[[index]][,greatest_less_than_index]%*%c(1,numerics)
  }
}

#' Predict
#' 
#' Predict (the first) most probable class at point. If probability returns NAs at a point, 
#' this returns NA at that point. If probabilities of classes are equal, this returns the
#' smaller index.
#' @param tes The tessellation object
#' @param x The point
#' @return Integer giving index of most probable class
#' @export
predict.tessellations=function(tes,x) {
  p=probability(tes,x)
  apply(p,2,function(p_){
    if (any(is.na(p_))) NA
    else which.max(p_)
  })
}
#'¨Get the density for a given class at given points
#'
#'¨Get the density for a given class at given points
#' @param tes The tessellations object
#' @param index The class index
#' @param points The points
#' @return Numeric vector giving value of density at each point.
#' @export
density=function(tes,index,points) {
  if (is.null(nrow(points))) {
    points=matrix(points,nrow=1)
    #cat("points argument cast to single row matrix.")
  }
  sims=geometry::tsearchn(tes$classdata[[index]],tes$tes[[index]]$tri,points)$idx
  sapply(1:length(sims),function(i)tes$simplex_models[[index]][,sims[i]]%*%c(1,points[i,]))
}
#'¨Get the density potential at given points
#'
#'¨Get the density potential at given points
#' @param tes The tessellations object
#' @param points The points
#' @return Numeric matrix of density values, with rows corresponding to classes, 
#' and columns to points.
#' @export
potential=function(tes,points) {
  if (is.null(nrow(points))) points=matrix(points,nrow=1)
  apply(points,1,function(point) 
    sapply(1:length(tes$tes),function(index) density(tes,index,point)))
}
#'¨Get the class probabilities at given points
#'
#'¨Get the class probabilities at given points. If no class densities are non-zero at a point
#' then the result is NAs.
#' @param tes The tessellations object
#' @param points The points
#' @return Numeric matrix of probability values, with rows corresponding to classes, 
#' and columns to points.
#' @export
probability = function (tes,points) {
  p=potential(tes,points)
  p[which(is.na(p))]=0
  apply(p,2,function(probs) {
      denom=sum(probs)
      if (denom==0) rep(NA,length(probs))
      else probs/denom
    })  
}