findExitingFace=function(p0,p1,tri,classdata) {
  temp=findIntersectWrapper(p0,p1,tri,classdata)
  if (class(temp)=="list") return (tri[temp$outFace])
  return (FALSE) 
}
findIntersectWrapper=function(p0,p1,tri,classdata) {
  vertices=classdata[tri,]
  findIntersectionLineConvexPolytope(p0,p1,vertices)
}
#' Complexity: D^2
findIntersectionLineConvexPolytope=function(p0,p1,vertices) {
  sampleVertices=c()
  normalVectors=c()
  for (i in 1:nrow(vertices)) {
    # Each vertice has an opposing face.
    # Find a vertex on opposing face
    if (i!=1) sampleVertices=cbind(sampleVertices,vertices[1,])
    else sampleVertices=cbind(sampleVertices,vertices[2,])
    # Find normal vector on opposing face
    coefs_=getHPCoefs(vertices[-i,])
    coefs=coefs_[2:length(coefs_)]
    alpha=coefs_[1]
    # Check if it is outside or inside facing
    if ((crossprod(coefs,vertices[i,])+alpha)>0)
      coefs=coefs*-1
    # Don't worry about the intersect
    normalVectors=cbind(normalVectors,coefs)
  }
  findIntersect(p0,p1,sampleVertices,normalVectors)
}

#' Find if a line intersects a convex polytope
#' 
#' Find if a line intersects a convex polytope
#' @param p0 The start of the line
#' @param p1 The end of the line
#' @param vertices A matrix, each column giving a vertex for each face, 
#' ordered as normalvectors
#' @param normalvectors A matrix, each column giving a normal vector for each face, 
#' ordered as vertices
#' @return FALSE if no intersection found. List of entry and exit points if interset is
#' found.
#' @keywords internal
# Input: a 3D segment S from point P0 to point P1
# a 3D convex polyhedron OMEGA with n faces F0,...,Fn-1 and
# Vi = a vertex for each face Fi
# ni = an outward normal vector for each face Fi
findIntersect=function(
  p0,
  p1,
  vertices,     # One vertex for each face
  normalvectors # One for each face, in same order as above
  ) {
  # Make sure p0 != p1
  if (identical(p0,p1)) stop ("Degenerate line passed to findIntersect!")
  
#   Initialize:
#     tE = 0 for the maximum entering segment parameter;
#     tL = 1 for the minimum leaving segment parameter;
#     dS = P1 - P0 is the segment direction vector;
  tE=0  
  tL=1
  dS=p1-p0
  inFace=c()
  outFace=c()

  for (i in 1:ncol(vertices)) {
#   N = - dot product of (P0 - Vi) and ni;
#   D = dot product of dS and ni;
    N = -((p0-vertices[,i])%*%normalvectors[,i])
    D = dS%*%normalvectors[,i]
    
    if (D == 0) { 
      # Then S is parallel to the face Fi 
      if (N < 0) 
        # then P0 is outside the face Fi
        return (FALSE) # since S cannot intersect OMEGA;
      else { 
        # S cannot enter or leave OMEGA across face Fi 
        # ignore face Fi and
        # continue to process the next face;
        next # continue for R!
      }
    }
  
    t = N / D;
    if (D < 0) {
      # Then segment S is entering OMEGA across face Fi 
      if (tE==t) inFace=c(inFace,i)
      else if (tE<t) inFace=i
      tE = max(tE,t)
      if (tE > tL) 
        # Then segment S enters OMEGA after leaving
        return (FALSE) # Since S cannot intersect OMEGA
    }
    else { # D > 0 
      # Then segment S is leaving OMEGA across face Fi
      if (tL==t) outFace=c(outFace,i)
      else if (tL>t) outFace=i
      tL = min(tL,t)
      if (tL < tE) 
        # Then segment S leaves OMEGA before entering
        return (FALSE) # Since S cannot intersect OMEGA
    }
  }

  # Output: [Note: to get here, one must have tE <= tL]
  # there is a valid intersection of S with OMEGA
  # from the entering point: P(tE) = P0 + tE * dS
  # to the leaving point:    P(tL) = P0 + tL * dS
  enter=p0+tE*dS
  exit=p0+tL*dS
  return (list(enter=enter,exit=exit,inFace=inFace,outFace=outFace))
}