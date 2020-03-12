poisson.glm <- function(formula, data=NULL,  burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, MALA=FALSE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "poisson")
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


#### Check on MALA argument
if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)

    
## Compute the blocking structure for beta     
block.temp <- common.betablock(p)
beta.beg  <- block.temp[[1]]
beta.fin <- block.temp[[2]]
n.beta.block <- block.temp[[3]]
list.block <- as.list(rep(NA, n.beta.block*2))
    for(r in 1:n.beta.block)
    {
    list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
    list.block[[r+n.beta.block]] <- length(list.block[[r]])
    }
    
    
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)      
    
    
    
#############################
#### Initial parameter values
#############################
#### Initial parameter values
mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
fitted <- exp(as.numeric(X.standardised %*% beta) + offset)     
    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples   
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept.all <- rep(0,2)
accept <- accept.all
proposal.sd.beta <- 0.01

    

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
        Y.DA[which.miss==0] <- rpois(n=n.miss, lambda=fitted[which.miss==0])    
        }else
        {}
         
        
           
    ####################
    ## Sample from beta
    ####################
    offset.temp <- offset
        if(MALA)
        {
        temp <- poissonbetaupdateMALA(X.standardised, K, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(X.standardised, K, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
        
   
         
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(X.standardised %*% beta) + offset
    fitted <- exp(lp)
    loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)

        
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <-Y.DA[which.miss==0]
        }else
        {}
        
        
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Update the proposal sds
            if(p>2)
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
            }else
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
            }
            accept.all <- accept.all + accept
            accept <- rep(0,2)
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
    
    
#### end timer
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
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.final <- c(accept.beta)
names(accept.final) <- c("beta")
    
    
#### Compute the fitted deviance
mean.beta <- apply(samples.beta, 2, mean)
fitted.mean <- exp(X.standardised %*% mean.beta  + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE), na.rm=TRUE)


#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)

    
#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    
    
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
summary.results <- summary.beta
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)

    
#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    
#### Compile and return the results
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - None\n")
    if(n.miss==0) samples.Y = NA

samples <- list(beta=samples.beta.orig, fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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
