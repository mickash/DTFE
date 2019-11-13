#' Create tessellations 
#' 
#' Create tessellations from data
#' @param x The data point coordinates
#' @param y The classes of the data points 
#' @param target The first to target for polytopic analysis.  Targeting
#' only functions with two classes. Default is NA, 
#' which will target the first class if there are two classes. 
#' Use 0 to target no classes. A target is required to use the analyseSimplexes function on the tessellation object.
#' @param method Not used - should be set to 1.
#' @param domain The domain of the data. If not specified, the min and max values of each column are used.
#' This can be used to calculate the size of the polytope outside the convex hull of the different tessellations
#' (see outer).
#' @param outer The size of the polytope outside the convex hull of the different tessellations. NA will 
#' calculate this from the domain.
#' @param generative If TRUE the hypervolumes of the lifted polytopes are calculated. This
#' enables the tessellations to be used as a generative model.
#' @return A tessellations object with the following fields:
#' \describe{
#'  \item{"tes"}{A list containing the class tessellations. This contains the fields:
#'    \describe{
#'      \item{"tri"}{A matrix giving the data points involved in the triangles of the tessellation. These
#'      refer to rows in the corresponding class data matrix.}
#'      \item{"neighbors"}{A list giving the neighbors for each triangle. If positive, these values refer to
#'      the triangles in the trianlge matrix. If negative, then the triangle borders the exterior polytope
#'      on a given edge.}
#'      \item{"area"}{The areas of the triangles.}
#'    }
#'  }
#'  \item{classdata}{Matrices giving the data coordinates for the different classes.}
#'  \item{x}{A matrix giving the coordinates of the complete data set.}
#'  \item{y}{A vector giving the classes of the complete data set.}
#'  \item{dim}{The number of dimensions.}
#'  \item{chulls}{A list of matrices the convex hulls of each tessellation/class data set.}
#'  \item{outer}{A vector specifying the area lying outside each tessellation.}
#'  \item{domain}{A matrix specifying the domain.}
#'  \item{classes}{A vector containing class labels (numbers if all labels are numeric).}
#'  \item{densities}{A vector specifying the point density estimates of the class data points.}
#'  \item{point_mappings}{For internal use.}
#'  \item{target}{The targeted class.}
#'  \item{critical_points_map}{For internal use.}
#'  \item{simplex_models}{A list of matrices, one for each class.  The nth column gives the 
#'  linear model associated with the nth polytope in the tessellation for that class.}
#'  \item{hypervolumes}{The hypervolumes associated with the polytopes of each tessellation.}
#'  \item{simplex_cond_probabilities}{The probability that a randomly sample from a 
#'  density comes from a given polytope. The normalized hypervolumes.}
#' }
#' @export
tessellations=function (
  x,
  y=rep(1,nrow(x)),
  target=NA,
  method=1, # Critical Points Method
  domain=NULL,
  outer=Inf, # NA calculate from domain
  generative=F,
  dummy_boundary_points=T,
  scatter=NA,
  multiply_duplicates=T,
  l2=.0001
  #  edgemap=T,
  #  enclosuremap=T
  ) {
  if (class(x)=="data.frame") {
    x=as.matrix(x)
    cat("X data frame cast to matrix.\n")
  }
  
  # Scatter if desired
  if (!is.na(scatter)) {
    if (length(scatter)==ncol(x)) {
      for (i in 1:ncol(x)) {
        x[,i]=x[,i]+rnorm(nrow(x),0,scatter[i])
      }
    }
    else if (length(scatter)==1) {
      for (i in 1:ncol(x)) {
        x[,i]=x[,i]+rnorm(nrow(x),0,scatter)
      }
    }
    else {
      stop("Invalid value: scatter argument length neither 1 nor equal to number of columns in x.")
    }
    cat("Points scattered.\n")
  }
  
  # Find domain if not specified
  if (is.null(domain)) {
    domain=rbind(apply(x,2,min),apply(x,2,max))
    cat("Domain calculated from x data.\n")
  } 
  domain_bounds=NULL
  if (dummy_boundary_points)
    domain_bounds=boundary_points(domain)
  
  # Find unique classes
  Y=unique(y)

  # Set target if NA
  if (length(Y)==2 && is.na(target)) {
    target=1
    cat("Target set to class 1.\n")
  }

  # Check targeted object has only two classes
  if (length(Y)>2 && !is.na(target) && target>0) 
    stop("Targeting can only occur with binary classes.")
  
  # Section A
  # For each class
  # 1. Create tessellation
  # 2. Calculate density
  # 3. Set up density interpolation
  
  # Separate data into classes
  #   - domain boundary points are added to each class's data (if they exist)
  classdata=list()
  multiples=list()
  x_plot=NULL
  for (class in Y) {
    cls_x=x[which(y==class),]
    agg_cls_x=as.matrix(aggregate(numdup ~., data=transform(cls_x,numdup=1), length))
    classdata[[length(classdata)+1]]=rbind(agg_cls_x[,1:ncol(x)],domain_bounds)
    multiples[[length(multiples)+1]]=agg_cls_x[,ncol(agg_cls_x)]
    x_plot=rbind(x_plot,classdata[[length(classdata)]])
  }
  cat("Data sorted by class, and duplicates recorded.\n")
  
  # Create tessellation for each class
  tes=lapply(1:length(Y),function(i) geometry::delaunayn(classdata[[i]],full=T))
  cat("Tessellations created.\n")

  # Calculate outer volumes if required
  if (is.na(outer)) {
    if (dummy_boundary_points)
      warning("Outer volumn should be set to infinity if using dummy boundary points. Continuing regardless.")
    total=prod(domain[2,]-domain[1,])
    areas=sapply(1:length(Y),function(i)sum(tes[[i]]$areas))
    outer=total-areas
  }
  else {
    outer=rep(outer,length(Y))
  }
  cat("Volume of outer polygons calculated:",outer,"\n")
  
  # Find convex hull points for each class
  chulls=lapply(1:length(Y),function (i) geometry::convhulln(classdata[[i]]))
  cat("Convex hulls calculated.\n")
  
  # Calculate point density estimations for each class
  densities=lapply(1:length(Y),function (i) 
    sapply(1:nrow(classdata[[i]]),function (j) {
      denom=sum(tes[[i]]$areas[which(tes[[i]]$tri==j,arr.ind=T)[,1]])
      if (denom==0)
        warning("Simplex areas for point is 0.")
      if (j %in% chulls[[i]]) 
        denom=denom+outer[i]
      if (multiply_duplicates)
        multiples[[i]][j]/denom # Returns 0 if denom == Inf
      else  
        1/denom # Returns 0 if denom == Inf
    }))
  cat("Point density estimates calculated.\n")
  
  # Create linear model coefficients for each simplex
  simplexmodels=lapply(1:length(Y),function (i)
    apply(tes[[i]]$tri,1,function(pnts) {
      X=cbind(rep(1,length(pnts)),classdata[[i]][pnts,])
      Y=densities[[i]][pnts]
      solve(crossprod(X)+diag(l2,ncol(X)))%*%crossprod(X,Y)
    } ))
  cat("Density coefficients for polygons calculated.\n")
  
  # Get hypervolumes if model is generative
  hypervolumes=NULL
  if (generative)
    hypervolumes=lapply(1:length(Y),function(i)
      sapply(1:nrow(tes[[i]]$tri),function(j) {
        if (all(densities[[i]][tes[[i]]$tri[j,]]==0)) 0 # If all points are on convex hull
        else {
          m=rbind(
            cbind(classdata[[i]][tes[[i]]$tri[j,],],rep(0,ncol(tes[[i]]$tri))),
            cbind(classdata[[i]][tes[[i]]$tri[j,],],densities[[i]][tes[[i]]$tri[j,]]))
          geometry::convhulln(m,"FA")$vol
        }
      }))
    
  simplex_cond_probabilities=NULL
  if (generative) simplex_cond_probabilities=lapply(1:length(Y),function(i) {
        total=sum(hypervolumes[[i]])
        sapply(1:length(hypervolumes[[i]]),function(j)hypervolumes[[i]][j]/total)
      })
    
  # For each class find simplexes in other tessellations 
  # containing each item in classdata
  # point_mappings[[i]][[j]][row] are 
  # the simplexes in tessellation j that enclose the points in 
  # class data i, ordered by row
  point_mappings=NULL
  if (!is.na(target) && target>0)
    point_mappings=lapply(1:length(Y),function(i)
      lapply(1:length(Y),function(j) {
          if (i==j) return (NA)
          geometry::tsearchn(classdata[[j]],tes[[j]]$tri,classdata[[i]])
        }))
    
  # Section B
  # For simplexes in each tesselation
  # 1. Find relevant points in other tessellations for probability estimation
  # 1.1 Find the vertices of all polytopes in other tessellations that have a non-empty
  #     intersection with this simplex.
  # 2. Calculate dominance score.
  
  # Create return object
  out=list(
    tes=tes,
    classdata=classdata,
    x=x,
    x_plot=x_plot,
    y=y,
    dim=ncol(x),
    chulls=chulls,
    outer=outer,
    domain=domain,
    classes=Y,
    densities=densities,
    point_mappings=point_mappings,  # Note that this is currently only used for critical point map (see below)
    target=ifelse(is.na(target),0,target),
    critical_points_map=NULL,
    simplex_models=simplexmodels,
    hypervolumes=hypervolumes,
    simplex_cond_probabilities=simplex_cond_probabilities,
    domain_bounds=domain_bounds,
    multiples=multiples
  )
  class(out)="tessellations"

  # If we have a target class, create critical points map.
  # [[tri]] gives matrix of critical points.
  if (!is.na(target) && target>0) 
    out$critical_points_map=createCriticalPointMap(out,method,target,ifelse(target==1,2,1))
     
  return (out)
}

boundary_points=function(domain) {
  d=ncol(domain)
  pnts=matrix(ncol=d,nrow=d^2)
  cnts=rep(1,d)
  r=1
  for (row in 1:d^2) {
    for (col in 1:d) {
      pnts[row,col]=domain[cnts[col],col]
    }
    # Update counter
    for (j in d:1) {
      cnts[j]=cnts[j]+1
      if (cnts[j]==3)
        cnts[j]=1
      else
        break
    }
  }
  return(pnts)
}