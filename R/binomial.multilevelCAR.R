binomial.multilevelCAR <- function(formula, data=NULL, trials, W, ind.area, ind.re=NULL, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, rho=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "binomial")
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

    
#### Check and format the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- n-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
failures.DA <- trials - Y.DA
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)


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
n.individual.miss <- rep(0,K)
    for(r in 1:K)
    {
    ind.area.list[[r]] <- which(ind.area==r)
    n.individual[r] <- length(ind.area.list[[r]])
    n.individual.miss[r] <- sum(which.miss[ind.area.list[[r]]])
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
    if(is.null(prior.sigma2)) prior.sigma2 <- c(1, 0.01)   
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.var.check(prior.tau2)    
common.prior.var.check(prior.sigma2)  


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
#### Initial parameter values
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
phi.extend <- phi[ind.area]
tau2 <- var(phi) / 10
    if(!is.null(ind.re))
    {
    psi <- rnorm(n=q, mean=rep(0,q), sd=res.sd/5)
    psi.extend <- psi[ind.re.num]
    sigma2 <- var(psi) / 10
    }else
    {
    psi.extend <- rep(0, n)
    }
lp <- as.numeric(X.standardised %*% beta) + phi.extend +  psi.extend + offset
prob <- exp(lp)  / (1 + exp(lp))



###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples   
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
    if(!is.null(ind.re)) samples.psi <- array(NA, c(n.keep, q))
    if(!is.null(ind.re)) samples.sigma2 <- array(NA, c(n.keep, 1))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, n))
samples.fitted <- array(NA, c(n.keep, n))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept.all <- rep(0,8)
accept <- accept.all
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
    if(!is.null(ind.re)) proposal.sd.psi <- 0.1
proposal.sd.rho <- 0.02
tau2.posterior.shape <- prior.tau2[1] + 0.5 * K
    if(!is.null(ind.re)) sigma2.posterior.shape <- prior.sigma2[1] + 0.5 * q  
    
    
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
        Y.DA[which.miss==0] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
        failures.DA <- trials - Y.DA
        }else
        {}
        
        
        
    ####################
    ## Sample from beta
    ####################
    offset.temp <- phi.extend + offset + psi.extend
        if(p>2)
        {
        temp <- binomialbetaupdateMALA(X.standardised, n, p, beta, offset.temp, Y.DA, failures.DA, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- binomialbetaupdateRW(X.standardised, n, p, beta, offset.temp, Y.DA, failures.DA, prior.mean.beta, prior.var.beta, proposal.sd.beta)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
        
 
        
    ####################
    ## Sample from phi
    ####################
    beta.offset <- X.standardised %*% beta + offset + psi.extend
    temp1 <- binomialcarmultilevelupdate(Wtriplet=W.triplet, Wbegfin=W.begfin, Wtripletsum=W.triplet.sum, ind_area_list=ind.area.list, n_individual=n.individual, nsites=K, phi=phi, tau2=tau2, y=Y.DA, failures=failures.DA, phi_tune=proposal.sd.phi, rho=rho, offset=beta.offset)
    phi <- temp1[[1]]
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + K
    phi.extend <- phi[ind.area]
        
        
        
    #############################
    ## Sample from psi and sigma2
    #############################
    if(!is.null(ind.re))
    {
    #### Update psi 
    beta.offset <- X.standardised %*% beta + offset + phi.extend
    temp1b <-  binomialcarmultilevelupdateindiv(ind_re_list=ind.re.list, n_re=n.re, q=q, psi=psi, sigma2=sigma2, y=Y.DA, failures=failures.DA, psi_tune=proposal.sd.psi, offset=beta.offset)
    psi <- temp1b[[1]] - mean(temp1b[[1]])
    accept[7] <- accept[7] + temp1b[[2]]
    accept[8] <- accept[8] + q    
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
            accept[5] <- accept[5] + 1           
            }else
            {}              
        accept[6] <- accept[6] + 1           
        }else
        {}
        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(X.standardised %*% beta) + phi.extend +  psi.extend + offset
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
            if(!is.null(ind.re)) samples.psi[ele, ] <- psi
            if(!is.null(ind.re)) samples.sigma2[ele, ] <- sigma2
            if(!fix.rho) samples.rho[ele, ] <- rho
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
            if(!fix.rho) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5)
            if(!is.null(ind.re)) proposal.sd.psi <- common.accceptrates1(accept[7:8], proposal.sd.psi, 40, 50)
        accept.all <- accept.all + accept
        accept <- rep(0,8)
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
accept.phi <- 100 * accept.all[3] / accept.all[4]
    if(!fix.rho)
    {
    accept.rho <- 100 * accept.all[5] / accept.all[6]
    }else
    {
    accept.rho <- NA    
    }
    if(!is.null(ind.re))
    {
    accept.psi <- 100 * accept.all[7] / accept.all[8]
    accept.sigma2 <- 100
    }else
    {
    accept.psi <- NA
    accept.sigma2 <- NA
    }
accept.tau2 <- 100
accept.final <- c(accept.beta, accept.phi, accept.psi, accept.rho, accept.tau2, accept.sigma2)
names(accept.final) <- c("beta", "phi","psi", "rho", "tau2", "sigma2")


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
mean.logit <- as.numeric(X.standardised %*% mean.beta) + mean.phi.extend + mean.psi.extend + offset    
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
    
summary.hyper <- array(NA, c(3 ,7))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
    if(!is.null(ind.re))
    {
    summary.hyper[2, 1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(samples.sigma2), geweke.diag(samples.sigma2)$z)
    }else
    {
    summary.hyper[2, ] <- rep(NA, 7)   
    }
if(!fix.rho)
    {
    summary.hyper[3, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[3, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[3, 1:3] <- c(rho, rho, rho)
    summary.hyper[3, 4:7] <- rep(NA, 4)
    }

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(nrow(summary.results)-2):nrow(summary.results)] <- c("tau2", "sigma2", "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    if(is.null(ind.re)) summary.results <- summary.results[-(nrow(summary.results)-1), ]      


#### Create the Fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)


#### Compile and return the results
    if(is.null(ind.re))
    {
    model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Multilevel Leroux CAR\n")
    }else
    {
    model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Multilevel Leroux CAR with factor random effects\n")
    }
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
    if(is.null(ind.re)) samples.sigma2 = NA    
    if(is.null(ind.re)) samples.psi = NA  

samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), zeta=mcmc(samples.psi), tau2=mcmc(samples.tau2), sigma2=mcmc(samples.sigma2), rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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

