gaussian.glm <- function(formula, data=NULL, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "gaussian")
K <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  
Y.DA <- Y

    
#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.nu2)  

    
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)      
    
    
    
#############################
#### Initial parameter values
#############################
#### Initial parameter values
mod.glm <- lm(Y~X.standardised-1, offset=offset)
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.unscaled)) * summary(mod.glm)$sigma
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
res.temp <- Y - X.standardised %*% beta.mean - offset
nu2 <- var(as.numeric(res.temp), na.rm=TRUE)
fitted <- as.numeric(X.standardised %*% beta) + offset

    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
nu2.posterior.shape <- prior.nu2[1] + 0.5*K

    
#### Beta update quantities
data.precision.beta <- t(X.standardised) %*% X.standardised
    if(length(prior.var.beta)==1)
    {
    prior.precision.beta <- 1 / prior.var.beta
    }else
    {
    prior.precision.beta <- solve(diag(prior.var.beta))
    }
    
    
    
###########################
#### Run the Bayesian model
###########################
#### Start timer
    if(verbose)
    {
    cat("Generating", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
    }else
    {
    percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
#### Create the MCMC samples    
    for(j in 1:n.sample)
    {
    ####################################
    ## Sample from Y - data augmentation
    ####################################
        if(n.miss>0)
        {
        Y.DA[which.miss==0] <- rnorm(n=n.miss, mean=fitted[which.miss==0], sd=sqrt(nu2))    
        }else
        {}
        
        
        
    ####################
    ## Sample from beta
    ####################
    fc.precision <- prior.precision.beta + data.precision.beta / nu2
    fc.var <- solve(fc.precision)
    beta.offset <- as.numeric(Y.DA - offset)
    beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
    fc.mean <- fc.var %*% beta.offset2
    chol.var <- t(chol(fc.var))
    beta <- fc.mean + chol.var %*% rnorm(p)        
        
        
        
    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  as.numeric(X.standardised %*% beta) + offset
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y.DA - fitted.current)^2)
    nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    

        
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + offset
    loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),K), log=TRUE)
        
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.nu2[ele, ] <- nu2
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
        }else
        {}

        
        
    ################################       
    ## print progress to the console
    ################################
        if(j %in% percentage.points & verbose)
        {
        setTxtProgressBar(progressBar, j/n.sample)
        }
    }
    
    
##### end timer
    if(verbose)
    {
    cat("\nSummarising results.")
    close(progressBar)
    }else
    {}
    

    
###################################
#### Summarise and save the results 
###################################
#### Compute the acceptance rates
accept.final <- rep(100, 2)
names(accept.final) <- c("beta", "nu2")
    
    
#### Compute the fitted deviance
mean.beta <- apply(samples.beta, 2, mean)
fitted.mean <- X.standardised %*% mean.beta + offset
nu2.mean <- mean(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),K), log = TRUE), na.rm=TRUE)

    
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)      
    

#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    
    
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
summary.hyper <- array(NA, c(1 ,7))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(samples.nu2), geweke.diag(samples.nu2)$z)
summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[nrow(summary.results)] <- c("nu2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)

    
#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(nu2.mean)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    
#### Compile and return the results
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - None\n")
    if(n.miss==0) samples.Y = NA
 
samples <- list(beta=samples.beta.orig, nu2=mcmc(samples.nu2), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
    results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)
    class(results) <- "CARBayes"
    
    #### Finish by stating the time taken    
    if(verbose)
    {
        b<-proc.time()
        cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
    }else
    {}
    return(results)
}

