multinomial.MVlerouxCAR <- function(formula, data=NULL, trials, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, rho=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
    
    
#### Frame object
frame.results <- common.frame(formula, data, "multinomial")
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
N.re <- K * (J-1)


#### If only one element in Y is missing then fix it as we know the total number of trials
which.miss.row <- J-apply(which.miss,1,sum)
which.miss.1 <- which(which.miss.row==1)
    if(length(length(which.miss.1))>0)
    {
        for(r in 1:length(which.miss.1))
        {
        which.miss[which.miss.1[r], is.na(Y[which.miss.1[r], ])] <- 1
        Y[which.miss.1[r], is.na(Y[which.miss.1[r], ])] <- trials[which.miss.1[r]] - sum(Y[which.miss.1[r], ], na.rm=T)    
        }
    n.miss <- sum(is.na(Y))
    which.miss.row <- J-apply(which.miss,1,sum)
    }else
    {}
Y.DA <- Y
const.like <- lfactorial(trials[which.miss.row==0]) - apply(lfactorial(Y[which.miss.row==0, ]),1,sum)
K.present <- sum(which.miss.row==0)


#### Determine which rows have missing values
    if(n.miss>0)    which.miss.row2 <- which(which.miss.row>0)   
    

#### Check and format the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- K-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
diffs <- apply(Y, 1, sum, na.rm=T) - trials
    if(max(diffs)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)

    
#### W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    if(ceiling(N.all/K)!= floor(N.all/K)) stop("The number of data points divided by the number of rows in W is not a whole number.", call.=FALSE)


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
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- J
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- diag(rep(1,J)) / 1000
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.varmat.check(prior.Sigma.scale, J-1) 
    

#### Compute the blocking structure for beta   
block.temp <- common.betablock(p, 5)
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
beta <- array(NA, c(p, (J-1)))
    for(i in 2:J)
    {
    mod.glm <- glm(cbind(Y[ ,i], trials - Y[ ,i])~X.standardised-1, offset=offset[ ,(i-1)], family="quasibinomial")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta[ ,(i-1)] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
    }
regression <- X.standardised %*% beta

theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat[ ,-1]  / theta.hat[ ,1]) - offset - regression
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.vec <- rnorm(n=N.re, mean=0, sd=res.sd)
phi <- matrix(phi.vec, nrow=K, byrow=TRUE)
Sigma <- cov(phi)
Sigma.inv <- solve(Sigma)    


###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, (J-1)*p))
samples.phi <- array(NA, c(n.keep, N.re))
samples.Sigma <- array(NA, c(n.keep, (J-1), (J-1)))
    if(!fix.rho) samples.rho <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, K.present))
samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept.beta <- rep(0,2*(J-1))
proposal.sd.beta <- rep(0.01, (J-1))
accept <- rep(0,4)
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
            for(g in 1:length(which.miss.row2))   
            {
            ## Determine which row (area) of Y to update
            row <- which.miss.row2[g]
            
            ## Compute the vector of probabilities for that row
            lp <- c(0, regression[row, ] + phi[row, ] + offset[row, ])
            prob <- exp(lp)  / sum(exp(lp))
            
            ## Do the multinomial data augmentation
                if(which.miss.row[row]==J)
                {
                ## All the Ys are missing
                Y.DA[row, ] <- as.numeric(rmultinom(n=1, size=trials[row], prob=prob))    
                }else
                {
                ## Not all the Ys are missing
                ## Re-normalise the probabilities
                prob[!is.na(Y[row, ])] <- 0
                prob <- prob / sum(prob)
                temp <- as.numeric(rmultinom(n=1, size=trials[row]-sum(Y[row, ], na.rm=T), prob=prob))    
                Y.DA[row, which.miss[row, ]==0]  <- temp[which.miss[row, ]==0]  
                }
            }
        }else
        {}

        
        
    ###################
    ## Sample from beta
    ###################
    offset.temp <- phi + offset
        for(r in 1:(J-1))
        {
        temp <- multinomialbetaupdateRW(X.standardised, K, J, p, r, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block, rep(0, K))
        beta[ ,r] <- temp[[1]][ ,r]
        accept.beta[r] <- accept.beta[r] + temp[[2]]
        accept.beta[(r+J-1)] <- accept.beta[(r+J-1)] + n.beta.block  
        }
    regression <- X.standardised %*% beta       

    
    
    ##################    
    ## Sample from phi
    ##################
    den.offset <- rho * W.triplet.sum + 1 - rho
    phi.offset <- regression + offset
    temp1 <- multinomialmcarupdateRW(W.triplet, W.begfin, K, J, phi, Y.DA, phi.offset, den.offset, Sigma.inv, rho,  proposal.sd.phi)      
    phi <- temp1[[1]]
        for(r in 1:(J-1))
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
        logprob.current <- 0.5 * (J-1) * det.Q - 0.5 * sum(diag(t(phi) %*% Q %*% phi %*% Sigma.inv))
        logprob.proposal <- 0.5 * (J-1) * det.Q.prop - 0.5 * sum(diag(t(phi) %*% Q.prop %*% phi %*% Sigma.inv))
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
    lp <- regression + phi + offset
    lp <- cbind(rep(0,K), lp)
    prob <- exp(lp)  / apply(exp(lp),1,sum)
    fitted <- prob * trials
    loglike <-  const.like + apply(Y[which.miss.row==0, ] * log(prob[which.miss.row==0, ]),1,sum)

        
        
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
            if(n.miss>0) samples.Y[ele, ] <- t(Y.DA)[is.na(t(Y))]
        }else
        {}
        
        
        
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
        #### Update the proposal sds
            for(r in 1:(J-1))
            {
                if(p>2)
                {
                proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J-1))], proposal.sd.beta[r], 40, 50)
                }else
                {
                proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J-1))], proposal.sd.beta[r], 30, 40)    
                }
            }
        
        proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
            if(!fix.rho)
            {
            proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
            }
        accept <- c(0,0,0,0)
        accept.beta <- rep(0,2*(J-1))
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
accept.beta <- 100 * sum(accept.beta[1:(J-1)]) / sum(accept.beta[(J:(2*(J-1)))])
accept.phi <- 100 * accept[1] / accept[2]
    if(!fix.rho)
    {
    accept.rho <- 100 * accept[3] / accept[4]
    }else
    {
    accept.rho <- NA    
    }
accept.Sigma <- 100
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma)
names(accept.final) <- c("beta", "phi", "rho", "Sigma")


#### Compute the fitted deviance
mean.beta <- matrix(apply(samples.beta, 2, mean), nrow=p, ncol=(J-1), byrow=F)
mean.phi <- matrix(apply(samples.phi, 2, mean), nrow=K, ncol=(J-1), byrow=T)
mean.logit <- X.standardised %*% mean.beta + mean.phi + offset
mean.logit <- cbind(rep(0,K), mean.logit)
mean.prob <- exp(mean.logit)  / apply(exp(mean.logit),1,sum)
deviance.fitted <- -2* sum(const.like + apply(Y[which.miss.row==0, ] * log(mean.prob[which.miss.row==0, ]),1,sum))


#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)


#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
    for(r in 1:(J-1))
    {
    samples.beta.orig[ ,((r-1)*p+1):(r*p)] <- common.betatransform(samples.beta[ ,((r-1)*p+1):(r*p)], X.indicator, X.mean, X.sd, p, FALSE)
    }
    
    
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,(J-1)*p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
col.name <- rep(NA, p*(J-1))

if(is.null(colnames(Y)))
{
    for(r in 1:(J-1))
    {
        col.name[((r-1)*p+1):(r*p)] <- paste("Category ", r+1,  " - ", colnames(X), sep="")   
    }
}else
{
    for(r in 1:(J-1))
    {
        col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[(r+1)],  " - ", colnames(X), sep="")   
    }
}
rownames(summary.beta) <- col.name
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

summary.hyper <- array(NA, c(J ,7))
summary.hyper[1:(J-1), 1] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.5)))
summary.hyper[1:(J-1), 2] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.025)))
summary.hyper[1:(J-1), 3] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.975)))
summary.hyper[1:(J-1), 4] <- n.keep
summary.hyper[1:(J-1), 5] <- accept.Sigma
summary.hyper[1:(J-1), 6] <- diag(apply(samples.Sigma, c(2,3), effectiveSize))
    for(r in 1:(J-1))
    {
    summary.hyper[r, 7] <- geweke.diag(samples.Sigma[ ,r,r])$z    
    }

    if(!fix.rho)
    {
    summary.hyper[J, 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
    summary.hyper[J, 4:7] <- c(n.keep, accept.rho, effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
    }else
    {
    summary.hyper[J, 1:3] <- c(rho, rho, rho)
    summary.hyper[J, 4:7] <- rep(NA, 4)
    }

summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[(((J-1)*p)+1): nrow(summary.results)] <- c(paste(rep("Sigma",(J-1)), 1:(J-1), 1:(J-1), sep=""), "rho")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)


#### Create the fitted values and residuals
fitted.values <- matrix(apply(samples.fitted, 2, mean), nrow=K, ncol=J, byrow=T)
response.residuals <- Y - fitted.values
var.y <- fitted.values * (1-fitted.values /  trials)
## Pearson is (observed - fitted) / sd
pearson.residuals <- response.residuals / sqrt(var.y)
residuals <- list(response=response.residuals, pearson=pearson.residuals)



#### Compile and return the results
model.string <- c("Likelihood model - Multinomial (logit link function)", "\nRandom effects model - Leroux MCAR\n")
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