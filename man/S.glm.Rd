\name{S.glm}
\alias{S.glm}
%- Also NEED an '\alias' for EACH other topic documented here.


\title{
Fit a generalised linear model to data.
}

\description{
Fit a generalised linear model to data, where the response variable 
can be binomial, Gaussian, multinomial, Poisson or zero-inflated Poisson (ZIP). 
Inference is conducted in a Bayesian setting using Markov chain Monte Carlo 
(MCMC)  simulation. Missing (NA) values are allowed in the response, and posterior 
predictive  distributions are created  for the missing values via data 
augmentation. These are  saved in the "samples" argument in the output of the 
function and are denoted by "Y".  For the multinomial model the first category 
in the multinomial data (first column of the response matrix) is taken as the 
baseline, and the covariates are linearly related to the log of the ratio 
(theta_j / theta_1) for j=1,...,J, where theta_j is the probability of being in 
category j. For the ZIP model covariates can be used to estimate the probability
of an observation being a structural zero, via a logistic regression equation.
For a full model specification see the vignette accompanying this package.
}


\usage{
S.glm(formula, formula.omega=NULL, family, data=NULL,  trials=NULL, burnin, n.sample, 
thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, 
prior.mean.delta=NULL, prior.var.delta=NULL, verbose=TRUE)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula for the covariate part of the model using the syntax of the
lm() function. Offsets can be included here using the offset() function. The 
response, offset and each covariate are vectors of length K*1. For the multinomial 
model the response and the offset (if included) should be matrices of dimension 
K*J, where K is the number of spatial units and J is the number of different 
variables (categories in the multinomial model). The covariates should each be a 
K*1 vector, and different regression parameters are estimated for each of the 
J variables. The response can contain missing (NA) values.
}
  \item{formula.omega}{
A one-sided formula object with no response variable (left side of the "~") 
needed, specifying the covariates in the logistic regression model for 
modelling the probability of an observation being a structural zero. Each 
covariate (or an offset) needs to be a vector of length K*1. Only required for 
zero-inflated Poisson models.
}
  \item{family}{
One of either "binomial", "gaussian", "multinomial", "poisson" or "zip", which 
respectively specify a binomial likelihood model with a logistic link function, 
a Gaussian likelihood model with an identity link function, a multinomial 
likelihood model with a logistic link function, a Poisson likelihood model 
with a log link function, or a zero-inflated Poisson model with a log link 
function. 
}
  \item{data}{
An optional data.frame containing the  variables in the formula.
}
 \item{trials}{
A vector the same length as the response containing the total number of trials 
for each data point. Only used if family="binomial" or family="multinomial". 
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
    \item{verbose}{
Logical, should the function update the user on its progress.  
}
}




\value{
\item{summary.results }{A summary table of the parameters.}
\item{samples }{A list containing the MCMC samples from the model.}
\item{fitted.values }{The fitted values based on posterior medians from the model. 
For the univariate data models this is a vector, while for the multivariate data 
models this is a matrix.}
\item{residuals }{If the family is "binomial", "gaussian" or "poisson",
then this is a matrix with 2 columns, where each column is a type of residual and 
each row relates to an area. The types are "response" (raw), and "pearson". If 
family is "multinomial", then this is a list with 2 elements, where 
each element is a matrix of residuals of a different type. Each row of a matrix 
relates to an area and each column to a cateogry within the multinomial response. 
The types of residual are "response" (raw), and "pearson".}
\item{modelfit }{Model fit criteria including the Deviance Information Criterion 
(DIC) and its corresponding estimated effective number of parameters (p.d), the Log 
Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion 
(WAIC) and its corresponding estimated number of effective parameters (p.w), and
the loglikelihood.}
\item{localised.structure }{NULL, for compatability with other models.}
\item{formula }{The formula (as a text string) for the response, covariate and 
offset parts of the model. If family="zip" this also includes the zero probability
logistic regression formula.}
\item{model }{A text string describing the model fit.}
\item{X }{The design matrix of covariates.}
}



\author{
Duncan Lee
}




\examples{
#################################################
#### Run the model on simulated data on a lattice
#################################################
#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)

#### Generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
logit <- x1 + x2
prob <- exp(logit) / (1 + exp(logit))
trials <- rep(50,K)
Y <- rbinom(n=K, size=trials, prob=prob)

#### Run the model
formula <- Y ~ x1 + x2
\dontrun{model <- S.glm(formula=formula, family="binomial", trials=trials, 
burnin=20000, n.sample=100000)}

#### Toy example for checking
model <- S.glm(formula=formula, family="binomial", trials=trials, 
burnin=10, n.sample=50)
}
