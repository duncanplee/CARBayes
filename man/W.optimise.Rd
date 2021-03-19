\name{W.optimise}
\alias{W.optimise}


\title{Creates an optimal neighbourhood matrix  for a set of spatial data using a baseline neighbourhood matrix and a  graph based optimisation algorithm. 
}


\description{
Creates an optimal neighbourhood matrix (W.est) for a given set of spatial data (spdata) from a baseline neighbourhood matrix (W) using the graph-based optimisation algorithm proposed by Lee, Meeks and Pettersson (2021). The baseline neighbourhood matrix W should be binary and based on a simple geographical specification such as the border sharing rule. The optimal neighbourhood matrix is constructed by removing neighbour relations (i.e. setting w_ij = w_ji = 0) if they are not appropriate for the data. 
}


\usage{
W.optimise(W, spdata, add=FALSE, remove=TRUE, remove_first=FALSE)
}


\arguments{
   \item{W}{A binary K by K neighbourhood matrix, where K is the number of 
   spatial units. A typical specification would have the jkth element 
   equalling one if areas (j, k) are spatially close (e.g. share a common 
   border) and  zero otherwise. Each row must contain at least one non-zero entry.
}
  \item{spdata}{
A K by 1 vector of spatial data that you wish to optimise the neighbourhood matrix
for. The kth element of this vector must correspond to the area represented
by the kth row of the W matrix.
}

  \item{add}{
Allow the optimiser to try to add edges to the graph. Defaults to FALSE.
}

  \item{remove}{
Allow the optimiser to try to remove edges from the graph. Defaults to TRUE.
}

  \item{remove_first}{
If only one of add or remove are TRUE, then this option has no effect. If both add and remove are TRUE, then this option determines whether the optimiser tries to add or remove edges (before alternating between the two phases). Defaults to FALSE, meaning that the optimiser tries to add edges first.
} 
}


\value{
\item{W.est }{An optimised K by K neighbourhood matrix.}
}


\references{
Lee, D and Meeks, K (2020), Improved inference for areal unit count data using graph-based optimisation. arXiv:2010.10893.
}


\author{
Kitty Meeks and William Pettersson
}



\examples{
####################################################
#### Run the optmiser on simulated data on a lattice
####################################################
#### Load other libraries required
library(MASS)

#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### Split the area into two groups between which there will be a boundary.
groups <-rep(1, K) 
groups[Grid$Var1>5] <- 2

#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	

#### Generate the response data
phi <- mvrnorm(n=1, mu=groups, Sigma=0.05 * exp(-0.1 * distance))
lp <- 4 + phi
Y <- rpois(n=K, lambda =exp(lp))


# Compute the spdata object on the linear predictor scale
spdata <- log(Y) - mean(log(Y))

#### Run W matrix optmiser
W.est <- W.optimise(W, spdata)
}