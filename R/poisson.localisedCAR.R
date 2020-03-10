poisson.localisedCAR <- function(formula, data=NULL, G, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.delta = NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame.localised(formula, data, "poisson", trials=NA)
K <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
log.Y <- log(Y)
log.Y[Y==0] <- -0.1 
regression.vec <- frame.results$regression.vec
beta <- frame.results$beta


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Format and check the number of clusters G     
    if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
    if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
    if(G<=1) stop("G is less than 2.", call.=FALSE)    
    if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
    if(floor(G/2)==ceiling(G/2))
    {
    Gstar <- G/2    
    }else
    {
    Gstar <- (G+1)/2          
    }
  
     
#### Priors
    if(!is.null(X))
    {
        if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
        if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
        common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
    }else
    {}

    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
common.prior.var.check(prior.tau2)

    if(is.null(prior.delta)) prior.delta <- 10
    if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
    if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    


#### Compute the blocking structure for beta     
    if(!is.null(X))
    {
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
    }else
    {}


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)  



#############################
#### Initial parameter values
#############################
res.temp <- log.Y - regression.vec - offset
clust <- kmeans(res.temp,G)
lambda <- clust$centers[order(clust$centers)]
Z <- rep(1, K)
    for(j in 2:G)
    {
    Z[clust$cluster==order(clust$centers)[j]] <- j    
    }
delta <- runif(1,1, min(2, prior.delta))
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=rep(0,K), sd = res.sd)
    for(i in 1:G)
    {
    phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
    }
tau2 <- var(phi) / 10

     
     
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.phi <- array(NA, c(n.keep, K))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.Z <- array(NA, c(n.keep, K))
samples.lambda <- array(NA, c(n.keep, G))
samples.delta <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K))
samples.fitted <- array(NA, c(n.keep, K))


#### Metropolis quantities    
    if(!is.null(X))
    {
    samples.beta <- array(NA, c(n.keep, p))
    accept.all <- rep(0,8)
    proposal.sd.beta <- 0.01
    }else
    {
    accept.all <- rep(0,6)    
    }

accept <- accept.all
proposal.sd.phi <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.lambda <- 0.01
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
     ###################
     ## Sample from beta
     ###################
        if(!is.null(X))
        {
        offset.temp <- phi + offset + lambda[Z]
            if(p>2)
            {
            temp <- poissonbetaupdateMALA(X.standardised, K, p, beta, offset.temp, Y, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
            }else
            {
            temp <- poissonbetaupdateRW(X.standardised, K, p, beta, offset.temp, Y, prior.mean.beta, prior.var.beta, proposal.sd.beta)
            }
        beta <- temp[[1]]
        accept[7] <- accept[7] + temp[[2]]
        accept[8] <- accept[8] + n.beta.block  
        regression.vec <- X.standardised %*% beta              
        }else
        {}
     
                 
               
     ##################
     ## Sample from phi
     ##################
     phi.offset <- regression.vec + offset  + lambda[Z]
        if(MALA)
        {
        temp1 <- poissoncarupdateMALA(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y, phi_tune=proposal.sd.phi, rho=1, offset=phi.offset)
        }else
        {
        temp1 <- poissoncarupdateRW(Wtriplet=W.triplet, Wbegfin=W.begfin, W.triplet.sum, nsites=K, phi=phi, tau2=tau2, y=Y, phi_tune=proposal.sd.phi, rho=1, offset=phi.offset)
        }
     phi <- temp1[[1]]
          for(i in 1:G)
          {
          phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
          }
     accept[1] <- accept[1] + temp1[[2]]
     accept[2] <- accept[2] + K    
         
         
     
     ###################
     ## Sample from tau2
     ###################
     temp2 <- quadform(W.triplet, W.triplet.sum, n.triplet, K, phi, phi, 1)
     tau2.posterior.scale <- temp2 + prior.tau2[2] 
     tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))
          
         
     
    #####################
    ## Sample from lambda
    #####################
     proposal.extend <- c(-1000, lambda, 1000)
     lambda.extend <- c(-1000, lambda, 1000)
     for(i in 1:G)
     {
         proposal.extend[(i+1)] <- rtruncnorm(n=1, a=proposal.extend[i], b=proposal.extend[(i+2)], mean=lambda[i], sd=proposal.sd.lambda)    
     }
     proposal <- proposal.extend[2:(G+1)]
     lp.current <- lambda[Z] + phi + regression.vec + offset
     lp.proposal <- proposal[Z] + phi + regression.vec + offset
     prob1 <- sum((exp(lp.current) - exp(lp.proposal)))
     prob2 <- sum(Y * (lp.proposal - lp.current))
     prob <- exp(prob1 + prob2)
          if(prob > runif(1))
          {
          lambda <- proposal
          accept[3] <- accept[3] + 1  
          }else
          {
          }
     accept[4] <- accept[4] + 1       

     
    
    ################
    ## Sample from Z
    ################
    Z.offset <- phi + offset + regression.vec
    Z.proposal <- sample(1:G, size=K, replace=TRUE)
    prior <- delta * ((Z - Gstar)^2 - (Z.proposal-Gstar)^2)    
    like <- exp(Z.offset) * (exp(lambda[Z]) - exp(lambda[Z.proposal])) + Y * (lambda[Z.proposal] - lambda[Z])         
    prob <- exp(like + prior)   
    test <- prob> runif(K)         
    Z[test] <- Z.proposal[test]         

         
         
    ####################
    ## Sample from delta
    ####################
    proposal.delta <-  rtruncnorm(n=1, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)    
    prob1 <- sum((Z-Gstar)^2) * (delta - proposal.delta)        
    prob2 <- K * log(sum(exp(-delta *(1:G - Gstar)^2))) - K * log(sum(exp(-proposal.delta *(1:G - Gstar)^2)))
    hastings <- log(dtruncnorm(x=delta, a=1, b=prior.delta, mean=proposal.delta, sd=proposal.sd.delta)) - log(dtruncnorm(x=proposal.delta, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)) 
    prob <- exp(prob1 + prob2 + hastings)    
          if(prob > runif(1))
          {
          delta <- proposal.delta
          accept[5] <- accept[5] + 1  
          }else
          {
          }
     accept[6] <- accept[6] + 1       

         
         
    #########################
    ## Calculate the deviance
    #########################
    lp <- lambda[Z] + phi + regression.vec + offset
    fitted <- exp(lp)
    loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)


         
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.phi[ele, ] <- phi
        samples.lambda[ele, ] <- lambda
        samples.tau2[ele, ] <- tau2
        samples.Z[ele, ] <- Z
        samples.delta[ele, ] <- delta
        samples.loglike[ele, ] <- loglike
        samples.fitted[ele, ] <- fitted
            if(!is.null(X)) samples.beta[ele, ] <- beta    
        }else
        {
        }


    
     ########################################
     ## Self tune the acceptance probabilties
     ########################################
     k <- j/100
        if(ceiling(k)==floor(k))
        {
            if(!is.null(X))
            {
                if(p>2)
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 40, 50)
                }else
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 30, 40)    
                }
            proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates1(accept[3:4], proposal.sd.lambda, 20, 40)
            proposal.sd.delta <- common.accceptrates2(accept[5:6], proposal.sd.delta, 40, 50, prior.delta/6)
            accept.all <- accept.all + accept
            accept <- rep(0,8)
            }else
            {
            proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates1(accept[3:4], proposal.sd.lambda, 20, 40)
            proposal.sd.delta <- common.accceptrates2(accept[5:6], proposal.sd.delta, 40, 50, prior.delta/6)
            accept.all <- accept.all + accept
            accept <- rep(0,6)     
            }
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
accept.phi <- 100 * accept.all[1] / accept.all[2]
accept.lambda <- 100 * accept.all[3] / accept.all[4]
accept.delta <- 100 * accept.all[5] / accept.all[6]
accept.tau2 <- 100

    if(!is.null(X))
    {
    accept.beta <- 100 * accept.all[7] / accept.all[8]   
    accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.tau2)
    names(accept.final) <- c("beta", "lambda", "delta", "phi", "tau2")   
    }else
    {
    accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.tau2)
    names(accept.final) <- c("lambda", "delta", "phi", "tau2")   
    }

     
#### Compute the fitted deviance
mean.phi <- apply(samples.phi, 2, mean)
mean.Z <- round(apply(samples.Z,2,mean),0)
mean.lambda <- apply(samples.lambda,2,mean)
    if(!is.null(X))
    {
    mean.beta <- apply(samples.beta,2,mean)
    regression.vec <- as.numeric(X.standardised %*% mean.beta)   
    }else
    {
    regression.vec <- rep(0,K)
    }
fitted.mean <- exp(regression.vec  + mean.lambda[mean.Z] + mean.phi + offset)
deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE))

     
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)

  
#### transform the parameters back to the origianl covariate scale.
    if(!is.null(X))
    {    
    samples.beta.orig <- common.betatransform(samples.beta, X.indicator, X.mean, X.sd, p, TRUE)
    }else
    {}     

     
#### Create a summary object
summary.hyper <- array(NA, c(2 ,7))
summary.hyper[1, 1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[1, 4:7] <- c(n.keep, accept.tau2, effectiveSize(samples.tau2), geweke.diag(samples.tau2)$z)
summary.hyper[2, 1:3] <- quantile(samples.delta, c(0.5, 0.025, 0.975))
summary.hyper[2, 4:7] <- c(n.keep, accept.delta, effectiveSize(samples.delta), geweke.diag(samples.delta)$z)

samples.lambda <- mcmc(samples.lambda)
summary.lambda <- t(apply(samples.lambda, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(samples.lambda), geweke.diag(samples.lambda)$z)
Z.used <- as.numeric(names(table(samples.Z)))
summary.lambda <- summary.lambda[Z.used, ]

    if(!is.null(X))
    {
    samples.beta.orig <- mcmc(samples.beta.orig)
    summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    rownames(summary.beta) <- colnames(X)
    summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)  
    row.names(summary.results)[(p+1):nrow(summary.results)] <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
    colnames(summary.results) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }else
    {
    summary.results <- rbind(summary.lambda, summary.hyper)
    row.names(summary.results) <- c(paste("lambda", Z.used, sep=""), "tau2", "delta")
    colnames(summary.results) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }


#### Create the fitted values and residuals
fitted.values <- apply(samples.fitted, 2, mean)
response.residuals <- as.numeric(Y) - fitted.values
pearson.residuals <- response.residuals /sqrt(fitted.values)
residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)


## Compile and return the results
model.string <- c("Likelihood model - Poisson (log link function)", "\nRandom effects  model - Localised CAR model\n")
if(is.null(X)) samples.beta.orig = NA

samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi), lambda=mcmc(samples.lambda), Z=mcmc(samples.Z), tau2=mcmc(samples.tau2), delta=mcmc(samples.delta), fitted=mcmc(samples.fitted), Y=NA)          
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=mean.Z,  formula=formula, model=model.string, X=X)
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
