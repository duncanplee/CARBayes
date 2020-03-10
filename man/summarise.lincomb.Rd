\name{summarise.lincomb}
\alias{summarise.lincomb}


\title{
Compute the posterior distribution for a linear combination of the covariates
from the linear predictor.  
}

\description{
This function takes in a CARBayes model object and computes the posterior 
distribution and posterior quantiles of a linear combination of the covariates
from the linear predictor. For example, if a quadratic effect of a covariate on 
the response was specified, then this function allows you to compute the posterior
distribution of the quadratic relationship.
}
\usage{
summarise.lincomb(model, columns=NULL, quantiles=0.5, distribution=FALSE)
}


\arguments{
  \item{model}{
A CARBayes model object from fitting one of the models in this package.
}
  \item{columns}{
A vector of column numbers stating which columns in the design matrix of covariates
the posterior distribution should be computed for.
}
 \item{quantiles}{
The vector of posterior quantiles required.
}
     \item{distribution}{
A logical value stating whether the entire posterior distribution should be 
returned or just the specified quantiles.
}
}


\value{
\item{quantiles }{A 2 dimensional array containing the requied posterior quantiles.
Each row relates to a data value, and each column to a different requested quantile. 
}
\item{posterior }{A 2 dimensional array containing the requied posterior distribution.
Each column relates to a different data value.
}
}



\author{
Duncan Lee
}


\examples{
## See the vignette accompanying this package for an example of its use.
}

