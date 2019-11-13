#' Function that checks for dominance
#' 
#' Function that checks for dominance. Passed to analyzeTessellation.
#' @param pnt The point(s) to test at
#' @param tes The tessellations object
#' @return Boolean vector specifying whether the target class's probability is the maximum at
#' each point.
#' @export
dom=function(pnts,tes) {
  d=probability(tes,pnts)
  if (is.null(ncol(d))) d=matrix(d,ncol=1)
  apply(d,2,function(d_) which.max(d_)==tes$target)
}
#' Function that checks for p-dominance
#' 
#' Function that checks for p-dominance. Passed to analyzeTessellation.
#' @param pnts The point(s) to test at
#' @param tes The tessellations object
#' @param alpha The test parameter
#' @return Boolean vector specifying whether the target class's probability is higher than 
#' alpha at each point.
#' @export
pdom=function(pnts,tes,alpha) {
  d=probability(tes,pnts)
  if (is.null(ncol(d))) d=matrix(d,ncol=1)
  apply(d,2,function(d_)d_[tes$target]>=alpha)
}
#' Function that checks for pq-dominance
#' 
#' Function that checks for pq-dominance. Passed to analyzeTessellation.
#' @param pnts The point(s) to test at
#' @param tes The tessellations object
#' @param alpha One test parameter
#' @param beta The other test parameter
#' @return Boolean specifying whether (i) the target class's probability is higher than 
#' alpha at the point; and (ii) whether the target class's raw density is higher than 
#' beta at the point.
#' @export
pqdom=function(pnts,tes,alpha,beta) {
  p=potential(tes,pnts)
  if (is.null(ncol(p))) p=matrix(p,ncol=1)
  p[which(is.na(p))]=0
  apply(p,2,function(p_) {
    denom=sum(p_)
    if (denom==0) F
    else {
      p2=p_/sum(p_)
      p_[tes$target]>=beta && p2[tes$target]>=alpha          
    }
  })
}

#' Analyze tessellation
#' 
#' Analyze the tessellation of the target class
#' @param tes The tessellation object.
#' @param func The function to use to analyze the critical points of each simplex.
#' @param ... Additional arguments required for func
#' @return Vector of boolean specifying if all critical points in the each simplex pass the 
#' criterion evaluated in func.
#' @keywords internal
#' @export
analyze=function (
  tes,
  func,
  ...
) { 
  sapply(tes$critical_points_map,function(pnts) all(func(pnts,tes,...)))
}
