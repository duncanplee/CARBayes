multinomial.glm <- function(formula, data=NULL, trials, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, verbose=TRUE)
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


    
#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
common.prior.beta.check(prior.mean.beta, prior.var.beta, p)

    
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
    
    
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, (J-1)*p))
samples.loglike <- array(NA, c(n.keep, K.present))
samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
    
    
#### Metropolis quantities
accept.beta <- rep(0,2*(J-1))
proposal.sd.beta <- rep(0.01, (J-1))
    
    
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
            lp <- c(0, regression[row, ] + offset[row, ])
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
        for(r in 1:(J-1))
        {
        temp <- multinomialbetaupdateRW(X.standardised, K, J, p, r, beta, offset, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block, rep(0, K))
        beta[ ,r] <- temp[[1]][ ,r]
        accept.beta[r] <- accept.beta[r] + temp[[2]]
        accept.beta[(r+J-1)] <- accept.beta[(r+J-1)] + n.beta.block  
        }
    regression <- X.standardised %*% beta       

        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- regression + offset
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
accept.final <- accept.beta
names(accept.final) <- c("beta")
    
    
#### Compute the fitted deviance
mean.beta <- matrix(apply(samples.beta, 2, mean), nrow=p, ncol=(J-1), byrow=F)
mean.logit <- X.standardised %*% mean.beta + offset
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
summary.results <- summary.beta
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
model.string <- c("Likelihood model - Multinomial (logit link function)", "\nRandom effects model - None\n")
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