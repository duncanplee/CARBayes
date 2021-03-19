\name{CARBayes-package}
\alias{CARBayes-package}
\alias{CARBayes}
\docType{package}

\title{
Spatial Generalised Linear Mixed Models for Areal Unit Data
}


\description{
Implements a class of univariate and multivariate spatial generalised linear mixed 
models for areal unit data, with inference in a Bayesian setting using Markov chain 
Monte Carlo (MCMC) simulation. The response variable can be binomial, Gaussian, 
multinomial, Poisson or zero-inflated Poisson (ZIP), and spatial autocorrelation 
is modelled by  a set of random effects that are assigned a conditional 
autoregressive (CAR) prior distribution. A number of different models are 
available for univariate spatial data, including models with no random effects as 
well as random effects modelled by different types of CAR prior, including the 
BYM model (Besag et al. (1991) <doi:10.1007/BF00116466>), the Leroux model 
(Leroux et al. (2000) <doi:10.1007/978-1-4612-1284-3_4>) and the localised model 
(Lee et al. (2015) <doi:10.1002/env.2348>). Additionally,  a multivariate CAR 
(MCAR) model for multivariate spatial data is available, as is a two-level 
hierarchical  model for modelling data relating to individuals within areas. Full 
details are given in the vignette accompanying this package. The initial creation 
of this package was supported by the Economic and Social Research Council (ESRC) 
grant RES-000-22-4256, and on-going development has been supported by the 
Engineering and Physical Science Research Council (EPSRC) grant EP/J017442/1, ESRC 
grant ES/K006460/1, Innovate UK / Natural Environment Research Council (NERC) 
grant NE/N007352/1 and the TB Alliance.
}

\details{
\tabular{ll}{
Package: \tab CARBayes\cr
Type: \tab Package\cr
Version: \tab 5.3\cr
Date: \tab 2021-03-23\cr
License: \tab GPL (>= 2)\cr
}
}




\author{
Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}


\references{
Besag, J. and  York, J and Mollie, A (1991). Bayesian image restoration with two 
applications in spatial statistics. Annals of the Institute of Statistics and 
Mathematics 43, 1-59.

Gelfand, A and Vounatsou, P (2003). Proper multivariate conditional autoregressive
models for spatial data analysis, Biostatistics, 4, 11-25.

Kavanagh, L., D. Lee, and G. Pryce (2016). Is Poverty Decentralising? Quantifying 
Uncertainty in the Decentralisation of Urban Poverty, Annals of the American 
Association of Geographers, 106, 1286-1298.

Lee, D. and Mitchell, R (2012). Boundary detection in disease mapping studies. 
Biostatistics, 13, 415-426.

Lee, D and Sarran, C  (2015). Controlling for unmeasured confounding and spatial 
misalignment in long-term air pollution and health studies, Environmetrics, 26, 
477-487.

Leroux B, Lei X, Breslow N (2000). "Estimation of Disease Rates in SmallAreas: A 
New Mixed Model for Spatial Dependence." In M Halloran, D Berry (eds.), 
\emph{Statistical Models in Epidemiology, the Environment and Clinical Trials},
pp. 179-191. Springer-Verlag, New York.

Roberts, G and Rosenthal, J (1998). Optimal scaling of discrete approximations to
the Langevin diffusions, Journal of the Royal Statistical Society Series B 60, 
255-268.
}



\examples{
## See the examples in the function specific help files and in the vignette
## accompanying this package.
}
