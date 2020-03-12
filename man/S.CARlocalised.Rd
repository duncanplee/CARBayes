\name{S.CARlocalised}
\alias{S.CARlocalised}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit a spatial generalised linear mixed model to data, where a set of spatially 
smooth random effects  are augmented with a piecewise constant intercept process.
}

\description{
Fit a spatial generalised linear mixed model to areal unit data, where the response
variable can be binomial or Poisson. Note, a Gaussian likelihood is not allowed 
because of a lack of identifiability among the parameters. The linear predictor is 
modelled by known covariates, a vector of random effects and a piecewise constant 
intercept process. The random effects are modelled by an intrinsic CAR prior, while 
the piecewise constant intercept process was proposed by Lee and Sarran (2015), and 
allow neighbouring areas to have very different values. Further details are given 
in the vignette accompanying this package. Inference is conducted in a Bayesian 
setting using Markov chain Monte Carlo (MCMC) simulation. Missing (NA) values are 
not allowed in this model. For a full model specification see the vignette 
accompanying this package.
}

\usage{
S.CARlocalised(formula, family, data=NULL, G, trials=NULL, W,
burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,
prior.tau2=NULL,prior.delta=NULL, MALA=FALSE, verbose=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response, offset and each covariate is a vector of length K*1.The  response 
cannot contain missing (NA) values.
}
  \item{family}{
One of either "binomial" or "poisson", which respectively specify a 
binomial likelihood model with a logistic link function, or a Poisson likelihood 
model with a log link function. 
}
  \item{data}{
An optional data.frame containing the  variables in the formula.
}

\item{G}{
The maximum number of distinct intercept terms (groups) to allow in the model.
}
  \item{trials}{
A vector the same length as the response containing the total number of trials 
for each area. Only used if family="binomial". 
}
   \item{W}{
A non-negative K by K neighbourhood matrix (where K is the number of spatial units). 
Typically a binary specification is used, where the jkth element equals one if areas 
(j, k) are spatially close (e.g. share a common border) and is zero otherwise. The 
matrix can be non-binary, but each row must contain at least one non-zero entry.
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
     \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
prior for tau2. Defaults to c(1, 0.01).  
}
     \item{prior.delta}{
The prior maximum for the cluster smoothing parameter delta. Defaults to 10.  
}
    \item{MALA}{
        Logical, should the function use Metropolis adjusted Langevin algorithm 
        (MALA) updates (TRUE) or simple random walk (FALSE, default) updates for 
        the regression parameters and random effects.   
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
\item{localised.structure }{A vector giving the posterior median of which intercept 
group each area is in.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Lee, D and Sarran, C  (2015). Controlling for unmeasured confounding and spatial 
misalignment in long-term air pollution and health studies, Environmetrics, 26, 
477-487.
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


#### Run the localised smoothing model
formula <- Y ~ 1
\dontrun{model <- S.CARlocalised(formula=formula, family="binomial", trials=trials,
G=2, W=W,burnin=20000, n.sample=100000)}

#### Toy example for checking
model <- S.CARlocalised(formula=formula, family="binomial", trials=trials,
G=2, W=W,burnin=10, n.sample=50)
}