% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot.tessellations}
\alias{plot.tessellations}
\title{Plot tessellations}
\usage{
\method{plot}{tessellations}(tes, dens = T, single = NULL,
  triangles = NULL, tricol = single)
}
\arguments{
\item{tes}{The tessellation object}

\item{dens}{Whether the density should be plotted. 
If not, only the tessellation is plotted.}

\item{single}{If only a single class is to be plotted, then the index of that
class should be specified here.}

\item{triangles}{The triangles to shade. Used only when plotting a single class,
not doing a density plot and using two dimensional data.}

\item{tricol}{The color to use when shading triangles.}
}
\description{
Plot function for a tessellations object.
}
