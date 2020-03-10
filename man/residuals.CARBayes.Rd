\name{residuals.CARBayes}
\alias{residuals.CARBayes}


\title{
    Extract the residuals from a model.
}
\description{
    This function takes a CARBayes object and returns a set of residuals. The 
    allowable types of residual are "response" (raw) and "pearson" (the default).
    In each case the fitted values are based on posterior means.
}

\usage{
    \method{residuals}{CARBayes}(object, type, \dots)
}



\arguments{
    \item{object}{
        A CARBayes fitted model object.
    }
    \item{type}{
        A text string and one of "response" or "pearson". If this
        argument is omitted the default is "pearson".
    }
    \item{\dots}{
        Ignored.
    }
}
\author{
    Duncan Lee
}
