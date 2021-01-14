gaussian.multilevelCAR <- function(formula, data=NULL, W, ind.area, ind.re=NULL, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, prior.sigma2=NULL, rho=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "gaussian")
n <- frame.results$n
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

    
#### rho
    if(is.null(rho))
    {
    rho <- runif(1)
    fix.rho <- FALSE   
    }else
    {
    fix.rho <- TRUE    
    }
    if(!is.numeric(rho) ) stop("rho is fixed but is not numeric.", call.=FALSE)  
    if(rho<0 ) stop("rho is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho is outside the range [0, 1].", call.=FALSE)  

#### CAR quantities
K <- length(unique(ind.area))
W.quants <- common.Wcheckformat(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin
K <- ncol(W)


#### Create the spatial determinant     
    if(!fix.rho)
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    }else
    {}


#### Checks and formatting for ind.area
    if(!is.vector(ind.area)) stop("ind.area is not a vector.", call.=FALSE)
    if(sum(ceiling(ind.area)==floor(ind.area))!=n) stop("ind.area does not have all integer values.", call.=FALSE)    
    if(min(ind.area)!=1) stop("the minimum value in ind.area is not 1.", call.=FALSE)    
    if(max(ind.area)!=K) stop("the maximum value in ind.area is not equal to the number of spatial areal units.", call.=FALSE)    
    if(length(table(ind.area))!=K) stop("the number of unique areas in ind.area does not equal K.", call.=FALSE)    

ind.area.list <- as.list(rep(0,K))
n.individual <- rep(0,K)
    for(r in 1:K)
    {
    ind.area.list[[r]] <- which(ind.area==r)
    n.individual[r] <- length(ind.area.list[[r]])
    }


#### Checks and formatting for ind.re
    if(!is.null(ind.re))
    {
        if(!is.factor(ind.re)) stop("ind.re is not a factor.", call.=FALSE)
    ind.re.unique <- levels(ind.re)
    q <- length(ind.re.unique)
    ind.re.list <- as.list(rep(NA,q))
    ind.re.num <- rep(NA, n)
        for(r in 1:q)
        {
        which.re.ind <- which(ind.re==ind.re.unique[r])
        ind.re.list[[r]] <- which.re.ind
        ind.re.num[which.re.ind] <- r
        }
    n.re <- as.numeric(lapply(ind.re.list,length))
    }else
    {
    cat("There are no individual level effects in this model\n")    
    }


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
    if(is.null(prior.sigma2)) prior.sigma2 <- c(1, 0.01)   
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)    
common.prior.var.check(prior.nu2)  
common.prior.var.check(prior.sigma2)  


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
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
phi.extend <- phi[ind.area]
tau2 <- var(phi) / 10
nu2 <- tau2
    if(!is.null(ind.re))
    {
    psi <- rnorm(n=q, mean=rep(0,q), sd=res.sd/5)
    psi.extend <- psi[ind.re.num]
    sigma2 <- var(psi) / 10
    }else
    {
    psi.extend <- rep(0, n)
    }
fitted <- as.numeric(X.standardised %*% beta) + phi.extend + offset + psi.extend

    
 
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
    if(!is.null(ind.re)) samples.psi <- array(NA, c(n.keep, q))
    if(!is.null(ind.re)) samples.sigma2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

    
#### Metropolis quantities
accept <- rep(0,2)
proposal.sd.rho <- 0.02
tau2.posterior.shape <- prior.tau2[1] + 0.5*K
nu2.posterior.shape <- prior.nu2[1] + 0.5*n
    if(!is.null(ind.re)) sigma2.posterior.shape <- prior.sigma2[1] + 0.5 * q 
    
 
#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
    if(rho==1) tau2.posterior.shape <- prior.tau2[1] + 0.5 * (K-n.islands)       

    
    
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
    beta.offset <- as.numeric(Y.DA - offset - phi.extend - psi.extend)
    beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
    fc.mean <- fc.var %*% beta.offset2
    chol.var <- t(chol(fc.var))
    beta <- fc.mean + chol.var %*% rnorm(p)        
        

        
    ##################
    ## Sample from nu2
    ##################
    fitted.current <-  as.numeric(X.standardised %*% beta) + phi.extend + offset + psi.extend
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y.DA - fitted.current)^2)
    nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale))    
        

            
    ####################
    ## Sample from phi
    ####################
    offset.phi <- (Y.DA - as.numeric(X.standardised %*% beta) - offset - psi.extend) / nu2
    offset.phi2 <- tapply(offset.phi, ind.area, sum, na.rm=T)
    phi <- gaussiancarmultilevelupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, n_individual=n.individual, nsites=K, phi=phi, tau2=tau2, rho=rho, nu2=nu2, offset=offset.phi2)
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    phi.extend <- phi[ind.area]
        
        
        
    #############################
    ## Sample from psi and sigma2
    #############################
        if(!is.null(ind.re))
        {
        #### Update psi 
        offset.psi <- (Y.DA - as.numeric(X.standardised %*% beta) - offset - phi.extend) / nu2
        offset.psi2 <- tapply(offset.psi, ind.re, sum, na.rm=T)
        temp1b <-  gaussiancarmultilevelupdateindiv(ind_re_list=ind.re.list, n_re=n.re, q=q, psi=psi, sigma2=sigma2, nu2=nu2, offset=offset.psi2)
        psi <- temp1b - mean(temp1b)
        psi.extend <- psi[ind.re.num]   
        
        
        #### Update sigma2 
        sigma2.posterior.scale <- 0.5 * sum(psi^2) + prior.sigma2[2] 
        sigma2 <- 1 / rgamma(1, sigma2.posterior.shape, scale=(1/sigma2.posterior.scale))
        }else
        {}
    
        
        
    ##################
    ## Sample from tau2
    ##################
    temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, rho)
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))

    
        
    ##################
    ## Sample from rho
    ##################
        if(!fix.rho)
        {
        proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)  
        temp3 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, proposal.rho)
        det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
        logprob.current <- det.Q - temp2 / tau2
        logprob.proposal <- det.Q.proposal - temp3 / tau2
        hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
        
        #### Accept or reject the proposal
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q <- det.Q.proposal
            accept[1] <- accept[1] + 1           
            }else
            {}              
        accept[2] <- accept[2] + 1           
        }else
        {}
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + phi.extend + offset + psi.extend
    loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),n), log=TRUE)

        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.nu2[ele, ] <- nu2
        samples.tau2[ele, ] <- tau2
            if(!is.null(ind.re))  samples.psi[ele, ] <- psi
            if(!is.null(ind.re)) samples.sigma2[ele, ] <- sigma2
            if(!fix.rho) samples.rho[ele, ] <- rho
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
        }else
        {}
        
        
        
    #######################################
    #### Update the acceptance rate for rho
    #######################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
        #### Update the proposal sds
            if(!fix.rho)
            {
            proposal.sd.rho <- common.accceptrates2(accept[1:2], proposal.sd.rho, 40, 50, 0.5)
            }
        accept <- c(0,0)
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
    if(!fix.rho)
    {
    accept.rho <- 100 * accept[1] / accept[2]
    }else
    {
    accept.rho <- NA    
    }
    if(!is.null(ind.re))
    {
    accept.psi <- 100 * accept[7] / accept[8]
    accept.sigma2 <- 100
    }else
    {
    accept.psi <- NA
    accept.sigma2 <- NA
    }
accept.final <- c(rep(100, 2), accept.psi, 100, 100, accept.rho, accept.sigma2)
names(accept.final) <- c("beta", "phi", "psi", "nu2", "tau2", "rho", "sigma2")
  
    
#### Compute the fitted deviance
mean.beta <- apply(samples.beta, 2, mean)
mean.phi <- apply(samples.phi, 2, mean)
mean.phi.extend <-  mean.phi[ind.area]
    if(!is.null(ind.re))
    {
    mean.psi <- apply(samples.psi, 2, mean)
    mean.psi.extend <-  mean.psi[ind.re.num]
    }else
    {
    mean.psi.extend <- rep(0,n)
    }
fitted.mean <- X.standardised %*% mean.beta + mean.phi.extend + mean.psi.extend + offset
nu2.mean <- mean(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),n), log = TRUE), na.rm=TRUE)

    
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

summary.hyper <- array(NA, c(4 ,7))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(samples.nu2), geweke.diag(samples.nu2)$z)
summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
    if(!is.null(ind.re))
    {
    summary.hyper[3, 1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(samples.sigma2), geweke.diag(samples.sigma2)$z)
    }else
    {
    summary.hyper[3, ] <- rep(NA, 7)   
    }
    
    if(!fix.rho)
    {
    summary.hyper[4, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[4, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[4, 1:3] <- c(rho, rho, rho)
    summary.hyper[4, 4:7] <- rep(NA, 4)
    }

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-3):nrow(summary.results)] <- c("nu2", "tau2", "sigma2", "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    if(is.null(ind.re)) summary.results <- summary.results[-(nrow(summary.results)-1), ]   
    

#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(nu2.mean)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    
#### Compile and return the results
    if(is.null(ind.re))
    {
    model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Multilevel Leroux CAR\n")
    }else
    {
    model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Multilevel Leroux CAR with factor random effects\n")
    }
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
    if(is.null(ind.re)) samples.sigma2 = NA    
    if(is.null(ind.re)) samples.psi = NA  

samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), zeta=mcmc(samples.psi), tau2=mcmc(samples.tau2), sigma2=mcmc(samples.sigma2), nu2=mcmc(samples.nu2), rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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

