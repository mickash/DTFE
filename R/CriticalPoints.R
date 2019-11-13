# Find points within polytope and points of intersection with edges leaving these points.
# Assume binary, with mapped class ordered first.

#' Find critical points map
#' 
#' Find critical points map: The points that will be checked within a polytope.
#' Method 1 (Necessary and Sufficient) - Complexity : N * D^3
#' Method 2 (Sufficient) - Complexity : N * D
#' @param tes The tessellation object.
#' @param method The method to use.
#' @return Row matrix of critical points
#' @keywords internal
createCriticalPointMap=function(
  tes,
  method,
  tes1_index,
  tes2_index
  ) {
  lapply(1:nrow(tes$tes[[1]]$tri),findCriticalPointsForPolytope,tes1_index,tes2_index,
         tes,method)  
}

#' Find Critical Points For Polytope
#' 
#' Find Critical Points For Polytope
#' @param tri The triangle to find critical points for.
#' @param tes The tessellation object.
#' @param method The method to use.
#' @return Row matrix of critical points
#' @keywords internal
findCriticalPointsForPolytope=function(
  tri,
  tes1_index,
  tes2_index,
  tes,
  method
  ) {
  # Find enclosed vertices from other
  enclosed_indices=getEnclosedVerticesIndices(tri,tes1_index,tes2_index,tes)

  additional_points_matrix=NULL
  if (length(enclosed_indices)>0) 
    if (method==1) # Necessary and sufficient: Look at exit points from enclosed vertices.
      # We follow edges from these vertices out of polytope and find critical points
      additional_points_matrix=findExitPoints(
        enclosed_indices,
        tri,
        tes1_index,
        tes2_index,
        tes,
        F # We request row matrix for efficiency (no transpose)
      )
    else if (method==2) # Sufficient: Look at other vertices of the polytopes
      # that the enclosed vertices are part of
      additional_points_matrix=findOtherVertices(enclosed_indices,tes2_index,tes)
    else stop ("Invalid method.")
  
  # Add enclosed vertices and vertices of this 
  rbind(
    additional_points_matrix,
    tes$classdata[[tes2_index]][enclosed_indices,],
    tes$classdata[[tes1_index]][tes$tes[[tes1_index]]$tri[tri,],]
    )
}

getEnclosedVerticesIndices=function(tri,tes1_index,tes2_index,tes) 
    which(tes$point_mappings[[tes2_index]][[tes1_index]]$idx==tri)

findOtherVertices=function(
  enclosed_indices,
  tes_index,
  tes
  ) {
  # Find all vertices of all polytopes each enclosed vertex
  vertices=unlist(lapply(enclosed_indices,findOtherVerticesFromVertex,tes_index,tes))
  # Remove duplicate indices 
  vertices=unique(vertices)
  # Return matrix of points.
  tes$classdata[[tes_index]][vertices,]
}

findOtherVerticesFromVertex=function(
  index,
  tes_index,
  tes
  ) {
  # Find polytopes with the indexed vertex
  polytopes=which(tes$tes[[tes_index]]$tri==index,arr.ind=T)[,1]
  # Return indices of all vertices in polytopes (c will collapse to vector from matrix)
  c(tes$tes[[tes_index]]$tri[polytopes,])
}


#' Complexity (WC): N * D^3
#' @return Matrix of exit points
findExitPoints=function(
  enclosed_indices,
  tri,
  tes1_index,
  tes2_index,
  tes,
  col
  ) {  
  # Find exit points for each enclosed vertex 
  # We request a column matrix so we can unlist it nicely into a row matrix
  exit_points_matrix_list=lapply(
    enclosed_indices,
    findExitPointsForVertex,
    enclosed_indices,
    tri,
    tes1_index,
    tes2_index,
    tes,
    T
    )
  
  # Check that at least one of the vertices involved has exit points.
  unlisted=unlist(exit_points_matrix_list)
  if (is.null(unlisted)) return (NULL)
  
  # Turn this into a nice matrix, with rows as exit points
  out=matrix(unlisted,byrow=T,ncol=tes$dim)
  # Return column or row matrix as desired
  if (col) return (t(out))
  else return (out)
}

#' Complexity (WC): N * D^3
#' @return Matrix of exit points
findExitPointsForVertex=function(
  enclosed_index, #
  enclosed_indices, #
  tri1, #
  tes1_index, #
  tes2_index, #
  tes, #
  col
  ) {
  # Find which polytopes it belongs to
  tris=which(tes$tes[[tes2_index]]$tri==enclosed_index,arr.ind=T)[,1]
  # Find exit points for each simplex in tes2 which as enclosed index as a vertex 
  # We request a column matrix so we can unlist it nicely into a row matrix
  exit_point_matrix_list=lapply(
    tris,
    findExitPointsInPolytope,
    enclosed_index,
    enclosed_indices,
    tri1,
    tes1_index,
    tes2_index,
    tes,
    T)

  # Check that at least one of the simplexes involved has exit points.
  unlisted=unlist(exit_point_matrix_list)
  if (is.null(unlisted)) return (NULL)
  
  # Turn this into a nice matrix, with rows as exit points
  out=matrix(unlisted,byrow=T,ncol=tes$dim)
  # Return column or row matrix as desired
  if (col) return (t(out))
  else return (out)
}

#' Complexity: D^3
#' @return Matrix of exit points
findExitPointsInPolytope=function(
  tri2, #
  vertex_index, #
  enclosed_vertices, #
  tri1, #
  tes1, #
  tes2, #
  tes, #
  col
  ) {
  # Find other indices in simplex vertex is in (it is part of tri2 in tes2)
  vertices_of_tri2=tes$tes[[tes2]]$tri[tri2,]
  # Find exit points for each edge in polytope from given vertex
  exit_point_list=lapply(
    vertices_of_tri2,
    findExitPointsInPolytopeOfEdge,
    vertex_index,
    enclosed_vertices,
    tri1,
    tes1,
    tes$classdata[[tes1]],
    tes$classdata[[tes2]],
    tes
  )
  # Check that at least one of the vertices has an edge that exits.
  unlisted=unlist(exit_point_list)
  if (is.null(unlisted)) return (NULL)
  
  # Turn this into a nice matrix, with rows as exit points
  out=matrix(unlisted,byrow=T,ncol=tes$dim)
  # Return column or row matrix as desired
  if (col) return (t(out))
  else return (out)
}

#' Complexity: D^2
#' @return Point of exit or NULL if no such point
findExitPointsInPolytopeOfEdge=function(
  other_vertex,
  vertex,
  enclosed_vertices,
  tri1,
  tes1_index,
  classdata1,
  classdata2,
  tes
  ) {  
  # Check if edge leaves tri
  # Note this is true when other_vertex==vertex
  if (other_vertex%in%enclosed_vertices) return (NULL)
  # If so, find point it leaves at.
  out=findIntersectionLineConvexPolytope(
    classdata2[vertex,],
    classdata2[other_vertex,],
    classdata1[tes$tes[[tes1_index]]$tri[tri1,],]) 
  # Check success
  if (class(out)!="list") { 
    warning ("Could not find point of exit that should exist.")
    return (NULL)
  }
  # Return exit point
  return (out$exit)
}
  
