\name{S.CARdissimilarity}
\alias{S.CARdissimilarity}

\title{
Fit a spatial generalised linear mixed model to data, where the random effects 
have a localised conditional autoregressive prior.
}

\description{
Fit a spatial generalised linear mixed model to areal unit data, where the response
variable can be binomial, Gaussian or Poisson. The linear predictor is modelled by
known covariates and a vector of random effects. The latter are modelled by the 
localised conditional autoregressive prior proposed by Lee and Mitchell (2012), 
and further details are given in the vignette accompanying this package. Inference
is conducted in a Bayesian setting using Markov chain Monte Carlo (MCMC) simulation.
Missing (NA) values are allowed in the response, and posterior predictive distributions
are created for the missing values using data augmentation. These are saved in the
"samples" argument in the output of the function and are denoted by "Y". For a full 
model specification see the vignette accompanying this package.
}

\usage{
S.CARdissimilarity(formula, family, data=NULL,  trials=NULL, W, 
Z, W.binary=TRUE, burnin, n.sample, thin=1, prior.mean.beta=NULL, 
prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, MALA=FALSE, verbose=TRUE)
}


\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response, offset and each covariate is a vector of length K*1. The response can 
contain missing (NA) values.
}
  \item{family}{
One of either "binomial", "gaussian" or "poisson", which respectively specify a 
binomial likelihood model with a logistic link function, a Gaussian likelihood 
model with an identity link function, or a Poisson likelihood model with a 
log link function. 
}
  \item{data}{
An optional data.frame containing the  variables in the formula.
}
  \item{trials}{
A vector the same length as the response containing the total number of trials 
for each area. Only used if family="binomial". 
}
  \item{W}{
A non-negative K by K neighbourhood matrix (where K is the number of spatial units). 
Typically a binary specification is used, where the jkth element equals one if 
areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
For this model only the matrix must be binary.
}
  \item{Z}{
A list, where each element is a K by K matrix of non-negative dissimilarity metrics.
}
  \item{W.binary}{
Logical, should the estimated neighbourhood matrix have only binary (0,1) values.
}
    \item{burnin}{
The number of MCMC samples to discard as the burn-in period.
}
  \item{n.sample}{
The number of MCMC samples to generate.
}
  \item{thin}{
The level of thinning to apply to the MCMC samples to reduce their temporal 
autocorrelation. Defaults to 1 (no thinning).
}
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are 
assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors 
are assumed). Defaults to a vector with values 100000.
}  
  \item{prior.nu2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for nu2. Defaults to c(1, 0.01) and only used if family="Gaussian".   
}
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for tau2. Defaults to c(1, 0.01).  
}
\item{MALA}{
Logical, should the function use Metropolis adjusted Langevin algorithm (MALA) 
    updates (TRUE) or simple random walk (FALSE, default) updates for the random 
    effects. Not applicable if family="gaussian".   
}
  \item{verbose}{
Logical, should the function update the user on its progress.  
}
}




\value{
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{A vector of fitted values for each area.}
    \item{residuals }{A matrix with 2 columns where each column is a type of 
        residual and each row relates to an area. The types are "response" (raw), 
        and "pearson".}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w), and
the loglikelihood.}
\item{accept }{The acceptance probabilities for the parameters.}
\item{localised.structure }{A list containing two matrices: W.posterior contains 
posterior medians for each element w_kj of the K by K neighbourhood matrix W; 
W.border.prob contains posterior probabilities that each w_kj element of the K by K
neighbourhood matrix W equals zero. This corresponds to the posterior probability
of a boundary in the random effects surface. The latter is only present if 
W.binary=TRUE, otherwise it is missing (NA). In all cases W elements that
correspond to two non-neighbouring areas have NA values.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}



\references{
Lee, D. and R. Mitchell (2012). Boundary detection in disease mapping studies. 
Biostatistics, 13, 415-426.
}

\author{
Duncan Lee
}

\examples{
#################################################
#### Run the model on simulated data on a lattice
#################################################
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
phi <- mvrnorm(n=1, mu=groups, Sigma=0.2 * exp(-0.1 * distance))
logit <- phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,K)
Y <- rbinom(n=K, size=trials, prob=prob)


#### Generate a dissimilarity metric
dissimilarity <- cbind(groups) + rnorm(K, sd=0.2)
dissimilarity.matrix <- as.matrix(dist(cbind(dissimilarity, dissimilarity), 
method="manhattan", diag=TRUE, upper=TRUE)) * W/2

Z <- list(dissimilarity.matrix=dissimilarity.matrix)

#### Run the localised smoothing model
formula <- Y ~ 1
\dontrun{model <- S.CARdissimilarity(formula=formula, family="binomial",
trials=trials, W=W, Z=Z, W.binary=TRUE, burnin=20000, n.sample=100000)}

#### Toy example for checking
model <- S.CARdissimilarity(formula=formula, family="binomial",
trials=trials, W=W, Z=Z, W.binary=TRUE, burnin=10, n.sample=50)
}