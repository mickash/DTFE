#' Plot tessellations 
#' 
#' Plot function for a tessellations object.
#' @param tes The tessellation object
#' @param dens Whether the density should be plotted. 
#' If not, only the tessellation is plotted.
#' @param single If only a single class is to be plotted, then the index of that
#' class should be specified here.
#' @param triangles The triangles to shade. Used only when plotting a single class,
#' not doing a density plot and using two dimensional data.
#' @param tricol The color to use when shading triangles.
#' @export
plot.tessellations=function(tes,dens=T,single=NULL,triangles=NULL,tricol=single) {
  if (!is.null(single)) plotsingle(tes,dens,single,triangles,tricol)
  else plotmultiple(tes,dens)
}

#' Plot tessellations
#' 
#' Plot tessellations 
#' @param tes The tessellation object
#' @param dens Whether the density should be plotted
#' @keywords internal
plotmultiple=function(tes,dens) {
  c=ncol(tes$x)
  if (c==2) 
    if (dens) densityplot(tes)
    else trimeshplot(tes)    
  else if (c==3 && !dens) tetrameshplot(tes)
  else stop("Mesh plots can only be created for two or three dimensional data, and density plots can only be created for two dimensional data.")
}

#' Plot tessellations in two dimensions.
#' 
#' Plot tessellations in two dimensions.
#' @param tes The tessellations object
#' @keywords internal
trimeshplot=function(tes){
  plot(tes$x)
  for (i in 1:length(tes$classdata)) {
    points(tes$classdata[[i]],col=i)
    geometry::trimesh(tes$tes[[i]][[1]],tes$classdata[[i]],col=i,add=T)
  }
}
#' Plot tessellations in three dimensions.
#' 
#' Plot tessellations in three dimensions.
#' @param tes The tessellations object
#' @keywords internal
tetrameshplot=function(tes){
  rgl::plot3d(tes$x_unique)
  for (i in 1:length(tes$classdata)) {
    rgl::points3d(tes$classdata[[i]],col=i)
    geometry::trimesh(tes$tes[[i]][[1]],tes$classdata[[i]],col=i,clear=F)
  }
}

#' Plot the tessellation densities.
#' 
#' Plot the tessellation densities.
#' @param tes The tessellations object
#' @param alpha The alpha parameter for the density shading.
#' @keywords internal
densityplot=function(tes,alpha=1){
  #rgl::plot3d(cbind(tes$x_unique,unlist(tes$densities)),col=1) # Get right plot size
  rgl::plot3d(cbind(tes$x_plot,unlist(tes$densities)),col=1) # Get right plot size
  #d=unlist(tes$densities)
  #rgl::plot3d(x=0,y=0,z=0,xlim=tes$domain[,1],ylim=tes$domain[,2],zlim=c(0,max(d)),col=1)
  for (i in 1:length(tes$classdata)) {
    rgl::points3d(cbind(tes$classdata[[i]],tes$densities[[i]]),col=i)
    for (j in 1:nrow(tes$tes[[i]]$tri)) {
      pnts=cbind(tes$classdata[[i]][tes$tes[[i]]$tri[j,],],tes$densities[[i]][tes$tes[[i]]$tri[j,]])
      rgl::triangles3d(pnts,col=i,alpha=alpha)
    }
  }
}

#' Plot information from a single tessellation 
#' 
#' Plot information from a single tessellation 
#' @param tes The tessellation object
#' @param dens Whether the density should be plotted
#' @param single The class to plot
#' @param triangles The triangles to shade. Used only when not doing a density plot 
#' and using two dimensional data.
#' @param tricol The color to use when shading triangles
#' @keywords internal
plotsingle=function(tes,dens,single,triangles,tricol) {
  c=ncol(tes$x)
  if (c==2)
    if (dens) singledensityplot(tes,single)
    else singletrimeshplot(tes,single,triangles,tricol)    
  else if (c==3 && !dens) singletetrameshplot(tes,single)
  else stop("Mesh plots can only be created for two or three dimensional data, and density plots can only be created for two dimensional data.")
}

#' Plot the tessellation of a single class in two dimensions.
#' 
#' Plot the tessellation of a single class. Data coordinates must be two dimensional.
#' @param tes The tessellations object
#' @param single The class to plot
#' @param triangles The triangles to shade
#' @param tricol The color to use when shading triangles
#' @keywords internal
singletrimeshplot=function(tes,single,triangles,tricol=single){
  plot(tes$classdata[[single]],col=single)
  for (tri in triangles)
    plottriangle(tes$tes[[single]],tri,tes$classdata[[single]],tricol)
  points(tes$classdata[[single]],col=single)
  geometry::trimesh(tes$tes[[single]][[1]],tes$classdata[[single]],col=single,add=T)
}
#' Plot tessellations of single class in three dimensions.
#' 
#' Plot tessellations of single class in three dimensions.
#' @param tes The tessellations object
#' @param single The class to plot
#' @keywords internal
singletetrameshplot=function(tes,single){
  rgl::plot3d(tes$classdata[[single]],col=single)
  geometry::trimesh(tes$tes[[single]][[1]],tes$classdata[[single]],col=single,clear=F)
}

#' Plot the density of a single class.
#' 
#' Plot the density of a single class. Data coordinates must be two dimensional.
#' @param tes The tessellations object
#' @param single The class to plot
#' @param alpha The alpha parameter for the density shading.
#' @keywords internal
singledensityplot=function(tes,single,alpha=1) {
  rgl::plot3d(cbind(tes$classdata[[single]],tes$densities[[single]]),col=single)
  for (j in 1:nrow(tes$tes[[single]]$tri)) {
    pnts=cbind(tes$classdata[[single]][tes$tes[[single]]$tri[j,],],tes$densities[[single]][tes$tes[[single]]$tri[j,]])
    rgl::triangles3d(pnts,col=single,alpha=alpha)
  }  
}

#' Plot a triangle
#' 
#' Plot a triangle
#' @param tes The tessellation object
#' @param tri The index of the triangle to plot 
#' @param data The data points the triangle indices refer to
#' @param tricol The color to use when plotting
#' @keywords internal
plottriangle=function(tes,tri,data,tricol) {
  pnts=data[tes$tri[tri,],]
  polygon(pnts,col=tricol)
}

