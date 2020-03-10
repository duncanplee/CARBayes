\name{MVS.CARleroux}
\alias{MVS.CARleroux}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
    Fit a multivariate spatial generalised linear mixed model to data, where the 
    random effects are modelled by a multivariate conditional autoregressive model.
}

\description{
    Fit a multivariate spatial generalised linear mixed model to areal unit data, 
    where the response variable can be binomial, Gaussian, multinomial or Poisson. 
    The linear predictor is modelled by known covariates and a vector of random 
    effects. The latter account for both spatial and between variable correlation, 
    via a Kronecker product formulation. Spatial correlation is captured by the 
    conditional autoregressive (CAR) prior proposed by Leroux et al. (2000), and 
    between variable correlation is captured by a between variable covariance 
    matrix with no fixed structure. This is a type of multivariate conditional 
    autoregressive (MCAR) model. Further details are given in the vignette accompanying 
    this package. Independent (over space) random effects can be obtained by 
    setting rho=0, while the intrinsic MCAR model can be obtained by setting 
    rho=1. Inference is conducted in a Bayesian setting using Markov chain Monte 
    Carlo (MCMC) simulation. Missing (NA) values are allowed in the response, and 
    posterior predictive distributions are created for the missing values using 
    data augmentation. These are saved in the "samples" argument in the output of 
    the function and are denoted by "Y". For the multinomial model the first 
    category in the multinomial data (first column of the response matrix) is 
    taken as the baseline, and the covariates are linearly related to the log of 
    the ratio (theta_j / theta_1) for j=1,...,J, where theta_j is the probability 
    of being in category j. For a full model specification see the vignette 
    accompanying this package.
}


\usage{
MVS.CARleroux(formula, family, data=NULL,  trials=NULL, W, burnin, 
n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL,
prior.Sigma.df=NULL, prior.Sigma.scale=NULL, rho=NULL, MALA=FALSE, 
verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response and the offset (if included) should be matrices of dimension K*J, where
K is the number of spatial units and J is the number of different variables 
(categories in the multinomial model). The covariates should each be a K*1 vector,
and different regression parameters are estimated for each of the J variables. 
Missing (NA) values are allowed in the response.
}
\item{family}{
One of either "binomial", "gaussian", "multinomial", or "poisson", which 
respectively specify a binomial likelihood model with a logistic link function, 
a Gaussian likelihood model with an identity link function, a multinomial 
likelihood model with a logistic link function, or a Poisson likelihood model 
with a log link function. 
}
\item{data}{
An optional data.frame containing the  variables in the formula.
}
\item{trials}{
Only used if family="binomial" or family="multinomial". For the binomial family it
is a K*J matrix matrix the same dimension as the response. A the multinomial family
it is a vector of length K.
}
\item{W}{
    A non-negative K by K neighbourhood matrix (where K is the number of spatial 
    units). Typically a binary specification is used, where the jkth element 
    equals one if areas (j, k) are spatially close (e.g. share a common border) 
    and is zero otherwise. The matrix can be non-binary, but each row must contain 
    at least one non-zero entry.
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
\item{prior.Sigma.df}{
    The prior degrees of freedom for the Inverse-Wishart prior for Sigma. 
    Defaults to J+1.   
}
\item{prior.Sigma.scale}{
    The prior J times J scale matrix for the Inverse-Wishart prior for Sigma. 
    Defaults to the identity matrix. 
}
     \item{rho}{
The value in the interval [0, 1] that the spatial dependence parameter rho is 
fixed at if it should not be estimated. If this arugment is NULL then rho is
estimated in the model.
}
\item{MALA}{
        Logical, should the function use Metropolis adjusted Langevin algorithm 
        (MALA) updates (TRUE, default) or simple random walk (FALSE, default) updates for 
        the random effects. Not applicable if family="gaussian" or "multinomial".   
}
\item{verbose}{
Logical, should the function update the user on its progress.  
}
}




\value{
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{A matrix of fitted values for each area and response 
variable.}
\item{residuals }{A list with 2 elements, where each element is a matrix of a type  
        of residuals. Each row of a matrix relates to an area and each column to a 
        response (category). The types of residual are "response" (raw), and 
        "pearson".}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w), and 
the loglikelihood.}
\item{accept }{The acceptance probabilities for the parameters.}
\item{localised.structure }{NULL, for compatability with other models.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}

\references{
Gelfand, A and Vounatsou, P (2003). Proper multivariate conditional autoregressive
models for spatial data analysis, Biostatistics, 4, 11-25.

Kavanagh, L., D. Lee, and G. Pryce (2016). Is Poverty Decentralising? Quantifying 
Uncertainty in the Decentralisation of Urban Poverty, Annals of the American 
Association of Geographers, 106, 1286-1298.

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
K <- nrow(W)

#### Generate the correlation structures
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)

Sigma <- matrix(c(1,0.5,0, 0.5,1,0.3, 0, 0.3, 1), nrow=3)
Sigma.inv <- solve(Sigma)
J <- nrow(Sigma)
N.all <- K * J

precision.phi <- kronecker(Q.W, Sigma.inv)
var.phi <- solve(precision.phi)

#### Generate the covariate component
x1 <- rnorm(K)
x2 <- rnorm(K)
XB <- cbind(0.1 * x1 - 0.1*x2, -0.1 * x1 + 0.1*x2, 0.1 * x1 - 0.1*x2)

#### Generate the random effects
phi <- mvrnorm(n=1, mu=rep(0,N.all), Sigma=var.phi)


#### Generate the response data
lp <-as.numeric(t(XB)) + phi
prob <- exp(lp) / (1 + exp(lp))
trials.vec <- rep(100,N.all)
Y.vec <- rbinom(n=N.all, size=trials.vec, prob=prob)

#### Turn the data and trials into matrices where each row is an area.
Y <- matrix(Y.vec, nrow=K, ncol=J, byrow=TRUE)
trials <- matrix(trials.vec, nrow=K, ncol=J, byrow=TRUE)

#### Run the Leroux model
formula <- Y ~ x1 + x2
\dontrun{model <- MVS.CARleroux(formula=formula, family="binomial", 
trials=trials, W=W, burnin=20000, n.sample=100000)}

#### Toy example for checking
model <- MVS.CARleroux(formula=formula, family="binomial", 
trials=trials, W=W, burnin=10, n.sample=50)
}
