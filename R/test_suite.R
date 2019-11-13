#' @export
test_suite=function(order=1,n=50,verbose=1,seed=0,skip=c(),load_results=T,save_results=T,filename='test_results') {
  if (load_results){
    results=load_test_suite_results(filename)
  }
  else
    results=list()
  
  if (!is.na(seed))
    set.seed(seed)
  datasets=list(Duncan,cars,iris,quakes,trees,ToothGrowth,swiss,stackloss)
  y_indices=c(4,2,1,4,3,1,1,4)
  for (i in 1:length(y_indices)) {
    if (i %in% skip) {
      if (verbose>0)
        cat("Skipping test:",i,"\n")
    }
    else if (i<=length(results) && !is.null(results[[i]])) {
      if (verbose>0)
        cat("Test results present in loaded results. Skipping test:",i,"\n")
    }
    else {
      if (verbose>0)
        cat("Beginning test:",i,"\n")
      results[[i]]=sapply(1:n,function(j) {
        test_on_data(datasets[[i]][,-y_indices[i],drop=FALSE],datasets[[i]][,y_indices[i]],order=order,verbose=verbose)
      })
    }
  }
  if (save_results){
    save_test_suite_results(filename)
  }
  return (results)
}
load_test_suite_results=function(filename="test_results"){
  load(filename)
  return (results)
}
save_test_suite_results=function(results,filename="test_results"){
  save(results,file=filename)
}

#' @export
plot_test_results=function(res,ignore_nan){
  if (ignore_nan)
    res_=lapply(res,function(r){
      bad=which(is.nan(r[1,]) | is.nan(r[2,]))
      if (length(bad)==0)
        r[1:2,]
      else
        r[1:2,-(bad)]
    })
  else
    res_=res
  diffs=lapply(res_,function(r)r[1,]-r[2,])
  std_diffs=lapply(diffs,function(d)d/sd(d))
  boxplot(std_diffs)
  abline(h=0,col="red")
  
  print(sapply(res,function(r)c(mean(r[1,],na.rm=T),mean(r[2,],na.rm=T))))
  print(sapply(res,function(r)t.test(r[1,],r[2,])))
}
test_on_data=function(x,y,order,train_prop=.8,verbose=1) {
  indices=sample(1:nrow(x))
  train_indices=indices[1:(nrow(x)*train_prop)]
  test_indices=indices[((nrow(x)*train_prop)+1):nrow(x)]
  
  cols=ncol(x)
  ks=c(cols+1:5,cols+c(10,15))
  
  knn_start=Sys.time()
  knn=tune_knn(x[train_indices,,drop=FALSE],y[train_indices],ks,order=order,verbose=verbose-1)
  p_knn=predict(knn,x[test_indices,,drop=FALSE])
  knn_end=Sys.time()
  if (verbose>1)
    cat("KNN time:",knn_end-knn_start,"\n")

  nn_start=Sys.time()
  nn=natural_regression(x[train_indices,,drop=FALSE],y[train_indices],order=order,verbose=verbose-1)
  p_nn=predict(nn,x[test_indices,,drop=FALSE])
  nn_end=Sys.time()
  if (verbose>1)
    cat("NN time:",nn_end-nn_start,"\n")
  
  predicted=which(!is.na(p_nn))
  mse_knn=mean((p_knn[predicted]-y[test_indices][predicted])^2)
  mse_nn=mean((p_nn[predicted]-y[test_indices][predicted])^2)
  if (verbose>0)
    cat("Result of test: KNN (",knn$k,") ",mse_knn,", NN",mse_nn,"\n")
  c(mse_knn,mse_nn,knn_start-knn_end,nn_start-nn_end)
}
tune_knn=function(x,y,ks,order,train_prop=.8,verbose=1){
  indices=sample(1:nrow(x))
  train_indices=indices[1:(nrow(x)*train_prop)]
  valid_indices=indices[((nrow(x)*train_prop)+1):nrow(x)]

  m=knn_regression(x[train_indices,,drop=FALSE],y[train_indices],k=NA,order=order)
  ps=predict(m,x[valid_indices,,drop=FALSE],ks=ks)
  mses=apply(ps,1,function(p){
    mean((p-y[valid_indices])^2)
  })
  best=which.min(mses)
  k=ks[best]
  if (verbose>0)
    cat("KNN model generated: K =",k,"\n")
  knn_regression(x,y,k,order)
}
knn_regression=function(
  x,
  y,
  k,
  order,
  l2=.0001
) {
  feature_types=sapply(x,class)
  factor_indices=which(feature_types=="factor")
   numeric_indices=which(feature_types!="factor")
  out=list(
    x=x,
    x_factor=x[,factor_indices,drop=FALSE],
    x_numeric=x[,numeric_indices,drop=FALSE],
    y=y,
    k=k,
    order=order,
    factor_indices=factor_indices,
    numeric_indices=numeric_indices,
    l2=l2)
  class(out)="knn"
  return (out)
}
predict.knn=function(model,x,ks=NULL){
  sapply(1:nrow(x),function(i){
    predict_one.knn(x[i,,drop=FALSE],model,ks)
  })
  
}
edist=function(x,y){
  sqrt(sum((x-y)^2))
}
predict_one.knn=function(x,model,ks=NULL){

  # Initialize k value(s)
  if (is.null(ks))
    ks=model$k
  
  # Find data rows that match factor variables
  if (length(model$factor_indices)>0)
    subset_x=apply(model$x_factor,1,function(row)all(row==x[model$factor_indices]))
  else
    subset_x=rep(T,nrow(model$x))

  # Find distances and order these  
  dist=apply(model$x_numeric[subset_x,,drop=FALSE],1,function(row)sqrt(sum((row-x[model$numeric_indices])^2)))
  # diff=
  # dist=apply(model$x_numeric[subset_x,,drop=FALSE],1,function(tr_x){
  #   edist(tr_x,x[model$numeric_indices])
  # })
  dist_order=order(dist)
  
  # Calculate predictions for all k values
  sapply(ks,function(k){
    neighbors=dist_order[1:min(k,length(dist_order))]
    if (model$order==0)
      return (mean(model$y[neighbors]))
    else if (model$order==1) {
      #df=cbind(model$x_numeric[neighbors,,drop=FALSE],y=model$y[neighbors])
      #m=lm(y~.,df)
      X=cbind(rep(1,length(neighbors)),as.matrix(model$x_numeric[neighbors,,drop=FALSE]))
      Y=model$y[neighbors]
      B=solve(crossprod(X)+diag(model$l2,ncol(X)))%*%crossprod(X,Y)
      B[,1]%*%c(1,as.matrix(x[model$numeric_indices]))
      #df=rbind(df,c(x,NA)) # Give x the right column names etc
      #predict(m,df[nrow(df),])
      #predict(m,x)
    }
    else 
      stop("Higher order local regression not implemented.")    
  })

}