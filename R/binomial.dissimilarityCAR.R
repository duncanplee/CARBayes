binomial.dissimilarityCAR <- function(formula, data=NULL,  trials, W, Z, W.binary=TRUE, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, MALA=FALSE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "binomial")
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


#### Check and format the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- K-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
failures.DA <- trials - Y.DA
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


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
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)


#### Compute the blocking structure for beta     
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
dat <- cbind(Y, failures)
mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta.mean - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd=res.sd)
tau2 <- var(phi) / 10
alpha <- runif(n=q, min=rep(0,q), max=rep(alpha.max/(2+q)))
lp <- as.numeric(X.standardised %*% beta) + phi + offset
prob <- exp(lp)  / (1 + exp(lp))


###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.alpha <- array(NA, c(n.keep, q))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))


#### Metropolis quantities
accept <- rep(0,6)
accept.all <- accept
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.alpha <- 0.02 * alpha.max
tau2.posterior.shape <- prior.tau2[1] + 0.5 * K



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
        Y.DA[which.miss==0] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
        failures.DA <- trials - Y.DA
        }else
        {}
        
        
        
    ####################
	## Sample from beta
	####################
    offset.temp <- phi + offset
        if(p>2)
        {
        temp <- binomialbetaupdateMALA(X.standardised, K, p, beta, offset.temp, Y.DA, failures.DA, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- binomialbetaupdateRW(X.standardised, K, p, beta, offset.temp, Y.DA, failures.DA, prior.mean.beta, prior.var.beta, proposal.sd.beta)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  

               
               
    ####################
	## Sample from phi
	####################
    beta.offset <- X.standardised %*% beta + offset
        if(MALA)
        {
        temp1 <- binomialcarupdateMALA(Wtriplet=W.triplet, Wbegfin=W.begfin, Wtripletsum=W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y.DA, failures=failures.DA, trials=trials, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset)
        }else
        {
        temp1 <- binomialcarupdateRW(Wtriplet=W.triplet, Wbegfin=W.begfin, Wtripletsum=W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y.DA, failures=failures.DA, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset)
        }
    phi <- temp1[[1]]
    phi <- phi - mean(phi)
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K

 		

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
        accept[5] <- accept[5] + 1
    	}else
    	{}
    accept[6] <- accept[6] + 1    
               
               
               
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(X.standardised %*% beta) + phi + offset
    prob <- exp(lp)  / (1 + exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)



    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
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
    	proposal.sd.alpha <- common.accceptrates2(accept[5:6], proposal.sd.alpha, 40, 50, alpha.max/4)
    	accept.all <- accept.all + accept
    	accept <- c(0,0,0,0,0,0)
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
accept.beta <- 100 * accept.all[1] / accept.all[2]
accept.phi <- 100 * accept.all[3] / accept.all[4]
accept.alpha <- 100 * accept.all[5] / accept.all[6]
accept.tau2 <- 100
accept.final <- c(accept.beta, accept.phi, accept.tau2, accept.alpha)
names(accept.final) <- c("beta", "phi", "tau2", "alpha")          

     
#### Compute the fitted deviance
mean.beta <- apply(samples.beta, 2, mean)
mean.phi <- apply(samples.phi, 2, mean)
mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.phi + offset    
mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
fitted.mean <- trials * mean.prob
deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE), na.rm=TRUE)


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

summary.hyper <- array(NA, c(1 ,7))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)

summary.results <- rbind(summary.beta, summary.hyper, summary.alpha)
    if(W.binary)
    {
    alpha.min <- c(rep(NA, (p+1)), alpha.threshold)
    summary.results <- cbind(summary.results, alpha.min)    
    }else
    {}

rownames(summary.results)[(p+1)] <- c("tau2")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    if(W.binary) summary.results[ , 8] <- round(summary.results[ , 8], 4)


## Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
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



#### Compile and return the results
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
    model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Binary dissimilarity CAR", "\nDissimilarity metrics - ", dis.eq, "\n")     
    }else
    {
    model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Non-binary dissimilarity CAR", "\nDissimilarity metrics - ", dis.eq, "\n")     
    }

    if(n.miss==0) samples.Y = NA

samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), alpha=mcmc(samples.alpha), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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
