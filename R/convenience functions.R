#' @export
getTrianglesWithVertex=function(tri,vertex) which(tri==vertex,arr.ind=T)[,1]

findMappedSimplex=function(classFromDataIndices,tesFrom_index,tesTo_index,tessellations) 
  tessellations$point_mappings[[tesFrom_index]][[tesTo_index]]$idx[classFromDataIndices]

getSimplexesThatInterceptEdge=function(
  tes_from,
  tes_to,
  tri,
  tes
  ) {
  vertices=tes$dim+1
  out=c()
  for (i in 1:vertices) {
    for (j in (1:vertices)[-i]) {
      out=c(out,tes$edge_mappings[[tes_from]][[tes_to]][[tri]][[i]][[j]])
    }
  }
  unique(out)
}