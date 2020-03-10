\name{print.CARBayes}
\alias{print.CARBayes}


\title{
Print a summary of a fitted CARBayes model to the screen.
}
\description{
This function takes a CARBayes object and returns a summary of the fitted model.
The summary includes, for selected parameters, posterior medians and 95 percent 
credible intervals, the number of samples, the acceptance rate of the Markov chains,
the effective number of independent samples and the Geweke convergence diagnostic 
in the form of a Z-score.
}

\usage{
\method{print}{CARBayes}(x, \dots)
}



\arguments{
  \item{x}{
A CARBayes fitted model object.
}
  \item{\dots}{
Ignored.s
}
}
\author{
Duncan Lee
}
