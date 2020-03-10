\name{combine.data.shapefile}
\alias{combine.data.shapefile}

\title{
Combines a data frame with a shapefile to create a SpatialPolygonsDataFrame object.
}

\description{
This function combines a data frame with a shapefile to create a SpatialPolygonsDataFrame
object from the `sp' package. The creation of this object allows the variables in
the data frame to be mapped using the `spplot()' function, and the neighbourhood
matrix W to be created using the `poly2nb' and `nb2mat' functions. An example is
given in the vignette accompanying this package. The mapping of the data to the 
shapefile is done by matching the rownames of the data frame to the first column
in the dbf file.
}

\usage{
combine.data.shapefile(data, shp, dbf)
}

\arguments{
\item{data}{
A data frame containing the variables relating to the K areas you wish to map or 
model. The row names of this data frame must appear in the first column of the 
dbf file.
}
  \item{shp}{
The .shp part of a shapefile containing the polygons for each of the K areas that
the data relate to.
}
  \item{dbf}{
The .dbf part of the shapefile containing a lookup table whose first column 
should inlcude the K row names of the data frame. 
}
}


\value{
\item{ }{A SpatialPolygonsDataFrame object from the sp 
package containing the conbined data and shapefile object.
}
}

\author{
Duncan Lee
}

\examples{
## See the vignette accompanying this package for an example of its use.
}
