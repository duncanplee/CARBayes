\name{summarise.samples}
\alias{summarise.samples}

\title{
Summarise a matrix of Markov chain Monte Carlo samples.
}

\description{
This function takes in a matrix of Markov chain Monte Carlo (MCMC) samples from 
a CARBayes model object, such as a set of parameters or fitted values, and 
calculates posterior quantiles and exceeedence probabilities. The latter are 
probabilities of the form P(quantity > c|data), where c is a threshold chosen by
the user.
}

\usage{
summarise.samples(samples, columns=NULL, quantiles=0.5, exceedences=NULL)
}

\arguments{
  \item{samples}{
A matrix of MCMC samples obtained from a CARBayes model object.
}
  \item{columns}{
A vector of column numbers stating which columns in the matrix of MCMC samples 
summaries are required for. Defaults to all columns.
}
 \item{quantiles}{
The vector of posterior quantiles required.
}
     \item{exceedences}{
The vector of threshold levels, c, that exceedence probabilities are required for.
}
}


\value{
\item{quantiles }{A 2 dimensional array containing the requied posterior quantiles. 
Each row relates to a parameter and each column to a different requested quantile.
}
\item{exceedences }{A 2 dimensional array containing the requied exceedence 
probabilities. Each row relates to a parameter, and each column to a different 
requested exceedence probability.
}
}


\author{
Duncan Lee
}




\examples{
## See the vignette accompanying this package for an example of its use.
}

