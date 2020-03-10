poisson.lerouxCAR <- function(formula, data=NULL,  W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, rho=NULL, MALA=TRUE, verbose=TRUE)
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


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)


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
mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    
log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta.mean - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
tau2 <- var(phi) / 10
fitted <- exp(as.numeric(X.standardised %*% beta) + offset + phi)     

    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    

#### Metropolis quantities
accept.all <- rep(0,6)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
    
    

##################################
#### Set up the spatial quantities
##################################
#### CAR quantities
W.quants <- common.Wcheckformat(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin
    
    
#### Create the determinant     
    if(!fix.rho)
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q <- 0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    }else
    {}   
    
    
#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
    if(rho==1) tau2.posterior.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   

        

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
    offset.temp <- phi + offset
        if(p>2)
        {
        temp <- poissonbetaupdateMALA(X.standardised, K, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- poissonbetaupdateRW(X.standardised, K, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, proposal.sd.beta)
        }
    beta <- temp[[1]]
    regression <- X.standardised %*% beta
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
 
        
        
    ####################
    ## Sample from phi
    ####################
    beta.offset <- regression + offset
        if(MALA)
        {
        temp1 <- poissoncarupdateMALA(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y.DA, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset)
        }else
        {
        temp1 <- poissoncarupdateRW(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y.DA, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset)
        }
    phi <- temp1[[1]]
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }
        else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K    
        
        
        
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
            accept[5] <- accept[5] + 1           
            }else
            {
            }              
        accept[6] <- accept[6] + 1           
        }else
        {}
        
        
    
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(regression) + phi + offset
    fitted <- exp(lp)
    loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)


    
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.tau2[ele, ] <- tau2
            if(!fix.rho) samples.rho[ele, ] <- rho
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
        }else
        {
        }
        
        
    
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
        proposal.sd.phi <- common.accceptrates1(accept[3:4], proposal.sd.phi, 40, 50)
            if(!fix.rho) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5)
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0,0,0)
        }else
        {   
        }

        
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
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
    if(!fix.rho)
    {
    accept.rho <- 100 * accept.all[5] / accept.all[6]
    }else
    {
    accept.rho <- NA    
    }
accept.tau2 <- 100
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.tau2)
names(accept.final) <- c("beta", "phi", "rho", "tau2")
    
    
#### Compute the fitted deviance
mean.beta <- apply(samples.beta, 2, mean)
mean.phi <- apply(samples.phi, 2, mean)
fitted.mean <- exp(X.standardised %*% mean.beta + mean.phi + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE), na.rm=TRUE)


#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)



#### Transform the parameters back to the origianl covariate scale.
samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, FALSE)

        
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
summary.hyper <- array(NA, c(2 ,7))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
    if(!fix.rho)
    {
    summary.hyper[2, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[2, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[2, 1:3] <- c(rho, rho, rho)
    summary.hyper[2, 4:7] <- rep(NA, 4)
    }
    
summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-1):nrow(summary.results)] <- c("tau2", "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)
    
    
#### Compile and return the results
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Leroux CAR\n")
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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
