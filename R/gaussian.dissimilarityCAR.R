gaussian.dissimilarityCAR <- function(formula, data=NULL, W, Z, W.binary=TRUE, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, verbose=TRUE)
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


#### Dissimilarity metric matrix
    if(class(Z)!="list") stop("Z is not a list object.", call.=FALSE)
    if(sum(is.na(as.numeric(lapply(Z, sum, na.rm=FALSE))))>0) stop("Z contains missing 'NA' values.", call.=FALSE)
q <- length(Z)
	if(sum(as.numeric(lapply(Z,nrow))==K) <q) stop("Z contains matrices of the wrong size.", call.=FALSE)
	if(sum(as.numeric(lapply(Z,ncol))==K) <q) stop("Z contains matrices of the wrong size.", call.=FALSE)
	if(min(as.numeric(lapply(Z,min)))<0) stop("Z contains negative values.", call.=FALSE)

if(!is.logical(W.binary)) stop("W.binary is not TRUE or FALSE.", call.=FALSE)
if(length(W.binary)!=1) stop("W.binary has the wrong length.", call.=FALSE)

    if(W.binary)
    {
    alpha.max <- rep(NA,q)
    alpha.threshold <- rep(NA,q)
        for(k in 1:q)
        {
        Z.crit <- quantile(as.numeric(Z[[k]])[as.numeric(Z[[k]])!=0], 0.5)
        alpha.max[k] <- -log(0.5) / Z.crit
        alpha.threshold[k] <- -log(0.5) / max(Z[[k]])
        }        
    }else
    {
    alpha.max <- rep(50, q)    
    }


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)
common.prior.var.check(prior.nu2)


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)  



#############################
#### Initial parameter values
#############################
mod.glm <- lm(Y~X.standardised-1, offset=offset)
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.unscaled)) * summary(mod.glm)$sigma
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

res.temp <- Y - X.standardised %*% beta.mean - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
tau2 <- var(phi) / 10
nu2 <- tau2
alpha <- runif(n=q, min=rep(0,q), max=rep(alpha.max/(2+q))) 
fitted <- as.numeric(X.standardised %*% beta) + phi + offset



###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.alpha <- array(NA, c(n.keep, q))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

     
## Metropolis quantities
accept <- c(0,0)
proposal.sd.alpha <- 0.02 * alpha.max
tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
nu2.posterior.shape <- prior.nu2[1] + 0.5*K  



##################################
#### Set up the spatial quantities
##################################
#### CAR quantities
W.quants <- common.Wcheckformat.disimilarity(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin    
spam.W <- W.quants$spam.W   
     
     
#### Create the Z triplet form
Z.triplet <- array(NA, c(n.triplet, q))
     for(i in 1:n.triplet)
     {
     row <- W.triplet[i,1]
     col <- W.triplet[i,2]
          for(j in 1:q)
          {
          Z.triplet[i,j] <- Z[[j]][row, col]     
          }     
     }

    if(W.binary)
    {
    W.triplet[ ,3] <- as.numeric(exp(-Z.triplet %*% alpha)>=0.5)        
    }else
    {
    W.triplet[ ,3] <- as.numeric(exp(-Z.triplet %*% alpha))    
    }
W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
spam.W@entries <- W.triplet[ ,3]      
spam.Wprop <- spam.W     
W.tripletprop <- W.triplet

     
#### Create the matrix form of Q
rho <- 0.99
Q <- -rho * spam.W 
diag(Q) <- rho * rowSums(spam.W) + 1-rho
det.Q <- sum(log(diag(chol.spam(Q))))     


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
	beta.offset <- as.numeric(Y.DA - offset - phi)
	beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
	fc.mean <- fc.var %*% beta.offset2
	chol.var <- t(chol(fc.var))
	beta <- fc.mean + chol.var %*% rnorm(p)  

		
		
	##################
	## Sample from nu2
	##################
    fitted.current <-  as.numeric(X.standardised %*% beta) + phi + offset
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * sum((Y.DA - fitted.current)^2)
    nu2 <- 1 / rgamma(1, nu2.posterior.shape, scale=(1/nu2.posterior.scale)) 

    
    
	####################
	## Sample from phi
	####################
    offset.phi <- (Y.DA - as.numeric(X.standardised %*% beta) - offset) / nu2    
    phi <- gaussiancarupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=K, phi=phi, tau2=tau2, rho=rho, nu2=nu2, offset=offset.phi)
    phi <- phi - mean(phi)
		

    
	##################
	## Sample from tau2
	##################
    temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, rho)
    tau2.posterior.scale <- temp2 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))		
    
 
    
    ######################
	#### Sample from alpha
	######################
    ## Propose a value
    proposal.alpha <- alpha
        for(r in 1:q)
    	{
    	proposal.alpha[r] <- rtruncnorm(n=1, a=0, b=alpha.max[r],  mean=alpha[r], sd=proposal.sd.alpha[r])
    	}
               
    ## Create the proposal values for W and Q
        if(W.binary)
        {
        W.tripletprop[ ,3] <- as.numeric(exp(-Z.triplet %*% proposal.alpha)>=0.5)        
        }else
        {
        W.tripletprop[ ,3] <- as.numeric(exp(-Z.triplet %*% proposal.alpha))    
        }
    W.triplet.sum.prop <- tapply(W.tripletprop[ ,3], W.tripletprop[ ,1], sum)
    spam.Wprop@entries <- W.tripletprop[ ,3]     
    Qprop <- -rho * spam.Wprop 
    diag(Qprop) <- rho * rowSums(spam.Wprop) + 1-rho
    det.Qprop <- sum(log(diag(chol.spam(Qprop))))     
    temp3 <- quadform(W.tripletprop, W.triplet.sum.prop, n.triplet, K, phi, phi, rho)              
               
    #### Calculate the acceptance probability
    logprob.current <- det.Q - temp2 / tau2
    logprob.proposal <- det.Qprop - temp3 / tau2
    hastings <- sum(log(dtruncnorm(x=alpha, a=rep(0,q), b=alpha.max, mean=proposal.alpha, sd=proposal.sd.alpha)) - log(dtruncnorm(x=proposal.alpha, a=rep(0,q), b=alpha.max, mean=alpha, sd=proposal.sd.alpha))) 
    prob <- exp(logprob.proposal - logprob.current + hastings)

    #### Accept or reject the proposed value
        if(prob > runif(1))
    	{
    	alpha <- proposal.alpha
    	det.Q <- det.Qprop 
        W.triplet[ ,3] <- W.tripletprop[ ,3]
        W.triplet.sum <- W.triplet.sum.prop
        accept[1] <- accept[1] + 1
    	}else
    	{}
    accept[2] <- accept[2] + 1     
		
      
             
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(X.standardised %*% beta) + phi + offset
    loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),K), log=TRUE)


    
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
        samples.alpha[ele, ] <- alpha
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
        }else
        {}

    
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
    	{
    	#### Update the proposal sds
    	proposal.sd.alpha <- common.accceptrates2(accept[1:2], proposal.sd.alpha, 40, 50, alpha.max/4)
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
accept.alpha <- 100 * accept[1] / accept[2]
accept.final <- c(100, 100, 100, 100, accept.alpha)
names(accept.final) <- c("beta", "phi", "nu2", "tau2", "alpha")
 
        
#### Compute the fitted deviance
mean.beta <- apply(samples.beta, 2, mean)
mean.phi <- apply(samples.phi, 2, mean)
fitted.mean <- X.standardised %*% mean.beta + mean.phi + offset
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

samples.alpha <- mcmc(samples.alpha)
summary.alpha <- t(apply(samples.alpha, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.alpha <- cbind(summary.alpha, rep(n.keep, q), rep(accept.alpha,q), effectiveSize(samples.alpha), geweke.diag(samples.alpha)$z)

	if(!is.null(names(Z)))
 	{
 	rownames(summary.alpha) <- names(Z)
	}else
	{
	names.Z <- rep(NA,q)
		for(j in 1:q)
		{
		names.Z[j] <- paste("Z[[",j, "]]", sep="")
		}	
	rownames(summary.alpha) <- names.Z	
	}

summary.hyper <- array(NA, c(2 ,7))
summary.hyper[1, 1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(samples.nu2), geweke.diag(samples.nu2)$z)
summary.hyper[2, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)

summary.results <- rbind(summary.beta, summary.hyper, summary.alpha)
    if(W.binary)
    {
    alpha.min <- c(rep(NA, (p+2)), alpha.threshold)
    summary.results <- cbind(summary.results, alpha.min)    
    }else
    {}
rownames(summary.results)[(p+1):(p+2)] <- c("nu2", "tau2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    if(W.binary) summary.results[ , 8] <- round(summary.results[ , 8], 4)


#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(nu2.mean)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

 
#### Create the posterior medians for the neighbourhood matrix W
W.posterior <- array(NA, c(K,K))
    if(W.binary)
    {
    W.border.prob <- array(NA, c(K,K))    
    }else
    {
    W.border.prob <- NA
    }

    for(i in 1:K)
    {
        for(j in 1:K)
        {
            if(W[i,j]==1)
            {
            z.temp <- NA
                for(k in 1:q)
                {
                z.temp <- c(z.temp, Z[[k]][i,j])
                }	
            z.temp <- z.temp[-1]
            w.temp <- exp(-samples.alpha %*% z.temp)
                if(W.binary)
                {
                w.posterior <- as.numeric(w.temp>=0.5)
                W.posterior[i,j] <- ceiling(median(w.posterior))
                W.border.prob[i,j] <- (1 - sum(w.posterior) / length(w.posterior))    
                }else
                {
                W.posterior[i,j] <- median(w.temp)
                }
            }else
            {
            }	
        }	
    }



## Compile and return the results
## Generate the dissimilarity equation
    if(q==1)
    {
    dis.eq <- rownames(summary.results)[nrow(summary.results)]    
    }else
    {
    dis.eq <- paste(rownames(summary.alpha), "+")
    len <- length(dis.eq)
    dis.eq[len] <- substr(dis.eq[len],1,nchar(dis.eq[2])-1)    
    }

if(W.binary)
    {
    model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Binary dissimilarity CAR", "\nDissimilarity metrics - ", dis.eq, "\n")     
    }else
    {
    model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Non-binary dissimilarity CAR", "\nDissimilarity metrics - ", dis.eq, "\n")     
    }

    if(n.miss==0) samples.Y = NA
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), alpha=mcmc(samples.alpha), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=list(W.posterior=W.posterior, W.border.prob=W.border.prob),  formula=formula, model=model.string, X=X)
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
