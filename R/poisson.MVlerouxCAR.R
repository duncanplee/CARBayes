poisson.MVlerouxCAR <- function(formula, data=NULL,  W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, rho=NULL, MALA=TRUE, verbose=TRUE)
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
J <- ncol(Y)
N.all <- K * J
Y.DA <- Y   

 
#### Create a missing list
    if(n.miss>0)
    {
    miss.locator <- array(NA, c(n.miss, 2))
    colnames(miss.locator) <- c("row", "column")
    locations <- which(t(which.miss)==0)
    miss.locator[ ,1] <- ceiling(locations/J)
    miss.locator[ ,2] <- locations - (miss.locator[ ,1]-1) * J
    }else
    {}


#### Check on MALA argument
if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(nrow(W)!= K) stop("The number of data points divided by the number of rows in W is not a whole number.", call.=FALSE)



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
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- J+1
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- diag(rep(1,J))
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.varmat.check(prior.Sigma.scale, J)  


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
beta <- array(NA, c(p, J))
    for(i in 1:J)
    {
    mod.glm <- glm(Y[ ,i]~X.standardised-1, offset=offset[ ,i], family="quasipoisson")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta[ ,i] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
    }

log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.vec <- rnorm(n=N.all, mean=0, sd=res.sd)
phi <- matrix(phi.vec, nrow=K, byrow=TRUE)
Sigma <- cov(phi)
Sigma.inv <- solve(Sigma)
regression <- X.standardised %*% beta
fitted <- exp(regression + phi + offset)
    

###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, J*p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.Sigma <- array(NA, c(n.keep, J, J))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

    
#### Metropolis quantities
accept.all <- rep(0,4)
accept <- accept.all
accept.all.beta <- rep(0,2*J)
accept.beta <- accept.all.beta
proposal.sd.beta <- rep(0.01, J)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
Sigma.post.df <- prior.Sigma.df + K  
    
    
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
Wstar <- diag(apply(W,1,sum)) - W
Q <- rho * Wstar + diag(rep(1-rho,K))
    
    
#### Create the determinant     
    if(!fix.rho)
    {
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q <- sum(log((rho * Wstar.val + (1-rho))))    
    }else
    {} 
    
    
#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
islands.all <- rep(islands,J)
n.islands <- max(W.islands$nc)
    if(rho==1) Sigma.post.df <- prior.Sigma.df + K - n.islands   

    
#### Specify vector variants
Y.vec <- as.numeric(t(Y))


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
        Y.DA[miss.locator] <- rpois(n=n.miss, lambda=fitted[miss.locator])    
        }else
        {}

        
        
    ###################
    ## Sample from beta
    ###################
    offset.temp <- phi + offset
        for(r in 1:J)
        {
            if(p>2)
            {
            temp <- poissonbetaupdateMALA(X.standardised, K, p, beta[ ,r], offset.temp[ ,r], Y.DA[ ,r], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block)
            }else
            {
            temp <- poissonbetaupdateRW(X.standardised, K, p, beta[ ,r], offset.temp[ ,r], Y.DA[ ,r], prior.mean.beta, prior.var.beta, proposal.sd.beta[r])
            }
        beta[ ,r] <- temp[[1]]
        accept.beta[r] <- accept.beta[r] + temp[[2]]
        accept.beta[(r+J)] <- accept.beta[(r+J)] + n.beta.block  
        }
regression <- X.standardised %*% beta          

    

    ##################
    ## Sample from phi
    ##################
    den.offset <- rho * W.triplet.sum + 1 - rho
    phi.offset <- regression + offset
        if(MALA)
        {
        temp1 <- poissonmcarupdateMALA(W.triplet, W.begfin, K, J, phi, Y.DA,  phi.offset, den.offset, Sigma.inv, rho, proposal.sd.phi)      
        }else
        {
        temp1 <- poissonmcarupdateRW(W.triplet, W.begfin, K, J, phi, Y.DA,  phi.offset, den.offset, Sigma.inv, rho, proposal.sd.phi)      
        }

    phi <- temp1[[1]]
        for(r in 1:J)
        {
        phi[ ,r] <- phi[ ,r] - mean(phi[ ,r])    
        }
    accept[1] <- accept[1] + temp1[[2]]
    accept[2] <- accept[2] + K    

        
        
    ####################
    ## Sample from Sigma
    ####################
    Sigma.post.scale <- t(phi) %*% Q %*% phi + prior.Sigma.scale
    Sigma <- riwish(Sigma.post.df, Sigma.post.scale)
    Sigma.inv <- solve(Sigma)
        
    
        
    ##################
    ## Sample from rho
    ##################
        if(!fix.rho)
        {
        ## Propose a new value
        proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
        Q.prop <- proposal.rho * Wstar + diag(rep(1-proposal.rho), K)
        det.Q.prop <-  sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))    
            
        ## Compute the acceptance rate
        logprob.current <- 0.5 * J * det.Q - 0.5 * sum(diag(t(phi) %*% Q %*% phi %*% Sigma.inv))
        logprob.proposal <- 0.5 * J * det.Q.prop - 0.5 * sum(diag(t(phi) %*% Q.prop %*% phi %*% Sigma.inv))
        hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q <- det.Q.prop
            Q <- Q.prop
            accept[3] <- accept[3] + 1           
            }else
            {}              
        accept[4] <- accept[4] + 1       
        }else
        {}
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    fitted <- exp(regression + phi + offset)
    loglike <- dpois(x=as.numeric(t(Y)), lambda=as.numeric(t(fitted)), log=TRUE)
        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- as.numeric(beta)
        samples.phi[ele, ] <- as.numeric(t(phi))
        samples.Sigma[ele, , ] <- Sigma
            if(!fix.rho) samples.rho[ele, ] <- rho
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- as.numeric(t(fitted))
            if(n.miss>0) samples.Y[ele, ] <- Y.DA[miss.locator]
        }else
        {}
        
        
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
        if(ceiling(k)==floor(k))
        {
        #### Update the proposal sds
            for(r in 1:J)
            {
                if(p>2)
                {
                proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J))], proposal.sd.beta[r], 40, 50)
                }else
                {
                proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J))], proposal.sd.beta[r], 30, 40)    
                }
            }
            
        proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            if(!fix.rho)
            {
            proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
            }
        accept.all <- accept.all + accept
        accept <- c(0,0,0,0)
        accept.all.beta <- accept.all.beta + accept.beta
        accept.beta <- rep(0,2*J)
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
accept.beta <- 100 * sum(accept.all.beta[1:J]) / sum(accept.all.beta[(J+1):(2*J)])
accept.phi <- 100 * accept.all[1] / accept.all[2]
    if(!fix.rho)
    {
    accept.rho <- 100 * accept.all[3] / accept.all[4]
    }else
    {
    accept.rho <- NA    
    }
accept.Sigma <- 100
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma)
names(accept.final) <- c("beta", "phi", "rho", "Sigma")
    

#### Compute the fitted deviance
mean.beta <- matrix(apply(samples.beta, 2, mean), nrow=p, ncol=J, byrow=F)
mean.phi <- matrix(apply(samples.phi, 2, mean), nrow=K, ncol=J, byrow=T)
fitted.mean <- exp(X.standardised %*% mean.beta + mean.phi + offset)
deviance.fitted <- -2 * sum(dpois(x=as.numeric(t(Y)), lambda=as.numeric(t(fitted.mean)), log=TRUE), na.rm=TRUE)

    
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)


#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
    for(r in 1:J)
    {
    samples.beta.orig[ ,((r-1)*p+1):(r*p)] <- common.betatransform(samples.beta[ ,((r-1)*p+1):(r*p) ], X.indicator, X.mean, X.sd, p, FALSE)
    }


#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,J*p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
col.name <- rep(NA, p*(J-1))

if(is.null(colnames(Y)))
{
    for(r in 1:J)
    {
        col.name[((r-1)*p+1):(r*p)] <- paste("Variable ", r,  " - ", colnames(X), sep="")   
    }
}else
{
    for(r in 1:J)
    {
        col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],  " - ", colnames(X), sep="")   
    }
}
rownames(summary.beta) <- col.name
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    
summary.hyper <- array(NA, c((J+1) ,7))
summary.hyper[1:J, 1] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.5)))
summary.hyper[1:J, 2] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.025)))
summary.hyper[1:J, 3] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.975)))
summary.hyper[1:J, 4] <- n.keep
summary.hyper[1:J, 5] <- accept.Sigma
summary.hyper[1:J, 6] <- diag(apply(samples.Sigma, c(2,3), effectiveSize))
    for(r in 1:J)
    {
    summary.hyper[r, 7] <- geweke.diag(samples.Sigma[ ,r,r])$z    
    }
    
    if(!fix.rho)
    {
    summary.hyper[(J+1), 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[(J+1), 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
    summary.hyper[(J+1), 4:7] <- rep(NA, 4)
    }
    
summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[((J*p)+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
#### Create the fitted values and residuals
fitted.values <- matrix(apply(samples.fitted, 2, mean), nrow=K, ncol=J, byrow=T)
response.residuals <- Y - fitted.values
pearson.residuals <- response.residuals / sqrt(fitted.values)
residuals <- list(response=response.residuals, pearson=pearson.residuals)


    
#### Compile and return the results
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects model - Leroux MCAR\n")
    if(fix.rho) samples.rho=NA
    if(n.miss==0) samples.Y = NA
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), Sigma=samples.Sigma, rho=mcmc(samples.rho), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
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
