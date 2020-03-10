\name{S.CARleroux}
\alias{S.CARleroux}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
Fit a spatial generalised linear mixed model to data, where the random effects 
have a Leroux conditional autoregressive prior.
}

\description{
Fit a spatial generalised linear mixed model to areal unit data, where the response
variable can be binomial, Gaussian, Poisson or zero-inflated Poisson (ZIP). The 
linear predictor is modelled by known covariates and a vector of random effects. 
The latter are modelled by the conditional autoregressive prior proposed by 
Leroux et al. (2000), and further details are given in the vignette accompanying 
this package. Independent random effects can be obtained by setting rho=0, while 
the intrinsic CAR model can be obtained by setting rho=1. Inference is conducted 
in a Bayesian setting using Markov chain Monte Carlo (MCMC) simulation. Missing 
(NA) values are allowed in the response, and posterior predictive distributions 
are created for the missing values using data augmentation. These are saved in 
the"samples" argument in the output of the function and are denoted by "Y". For 
the ZIP model covariates can be used to estimate the probability of an observation 
being a structural zero, via a logistic regression equation. For a full model 
specification see the vignette accompanying this package.
}


\usage{
S.CARleroux(formula, formula.omega=NULL, family, data=NULL,  trials=NULL, W, burnin, 
n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,
prior.nu2=NULL, prior.tau2=NULL, prior.mean.delta=NULL, prior.var.delta=NULL, 
rho=NULL, MALA=TRUE, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response, offset and each covariate is a vector of length K*1. The response can 
contain missing (NA) values.
}
  \item{formula.omega}{
A one-sided formula object with no response variable (left side of the "~") 
needed, specifying the covariates in the logistic regression model for 
modelling the probability of an observation being a structural zero. Each 
covariate (or an offset) needs to be a vector of length K*1. Only required for 
zero-inflated Poisson models.
}
  \item{family}{
One of either "binomial", "gaussian", "poisson" or "zip", which respectively 
specify a binomial likelihood model with a logistic link function, 
a Gaussian likelihood model with an identity link function, a Poisson likelihood 
model with a log link function, or a zero-inflated Poisson model with a log link 
function. 
}
  \item{data}{
An optional data.frame containing the  variables in the formula.
}
 \item{trials}{
A vector the same length as the response containing the total number of trials 
for each area. Only used if family="binomial". 
}
   \item{W}{A non-negative K by K neighbourhood matrix (where K is the number of 
   spatial units). Typically a binary specification is used, where the jkth 
   element equals one if areas (j, k) are spatially close (e.g. share a common 
   border) and is zero otherwise. The matrix can be non-binary, but each row must 
   contain at least one non-zero entry.
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
  \item{prior.mean.delta}{
A vector of prior means for the regression parameters delta (Gaussian priors are 
assumed) for the zero probability logistic regression component of the model. 
Defaults to a vector of zeros.
}
  \item{prior.var.delta}{
A vector of prior variances for the regression parameters delta (Gaussian priors 
are assumed) for the zero probability logistic regression component of the model. 
Defaults to a vector with values 100000.
}
     \item{rho}{
The value in the interval [0, 1] that the spatial dependence parameter rho is 
fixed at if it should not be estimated. If this arugment is NULL then rho is
estimated in the model.
}
    \item{MALA}{
        Logical, should the function use Metropolis adjusted Langevin algorithm 
        (MALA) updates (TRUE, default) or simple random walk (FALSE) updates for 
        the random effects. Not applicable if family="gaussian".   
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
\item{localised.structure }{NULL, for compatability with other models.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model. If family="zip" this also includes the zero probability
logistic regression formula.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Leroux B, Lei X, Breslow N (2000). "Estimation of Disease Rates in SmallAreas: A 
New Mixed Model for Spatial Dependence." In M Halloran, D Berry (eds.), 
\emph{Statistical Models in Epidemiology, the Environment and Clinical Trials},
pp. 179-191. Springer-Verlag, New York.
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

#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <- as.matrix(dist(Grid))
W <-array(0, c(K,K))
W[distance==1] <-1 	
	
#### Generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd=0.05)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=0.4 * exp(-0.1 * distance))
logit <- x1 + x2 + theta + phi
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,K)
Y <- rbinom(n=K, size=trials, prob=prob)


#### Run the Leroux model
formula <- Y ~ x1 + x2
\dontrun{model <- S.CARleroux(formula=formula, family="binomial", 
trials=trials, W=W, burnin=20000, n.sample=100000)}

#### Toy example for checking
model <- S.CARleroux(formula=formula, family="binomial", 
trials=trials, W=W, burnin=10, n.sample=50)
}
