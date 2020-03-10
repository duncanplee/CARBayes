\name{highlight.borders}
\alias{highlight.borders}

\title{
Creates a SpatialPoints object identifying a subset of borders between neighbouring
areas.
}


\description{
Creates a SpatialPoints object identifying a subset of borders between neighbouring
areas, which allows them to be overlayed on a map. An example is given in 
the vignette accompanying this package.
}

\usage{
highlight.borders(border.locations, spdata)
}

\arguments{
  \item{border.locations}{
A K by K matrix, where K is the number of areas, containing 3 distinct values: 
NA for non-neighbouring areas; 0 for borders between neighbouring areas to be 
highlighted on a map; and 1 for borders between neighbouring areas not to be 
highlighted on a map.
}
  \item{spdata}{
The SpatialPolygonsDataFrame object used for plotting the data and creating the 
original neighbourhood matrix W.
}
}



\value{
\item{ }{A SpatialPoints object from the sp package, which contains
the vertices of all the borders to be highlighted on the map.  The mapping can 
be done using the leaflet package, see the vignette accompanying this package
for an example.
}
}

\author{
Duncan Lee
}


\examples{
## See the vignette accompanying this package for an example of its use.
}