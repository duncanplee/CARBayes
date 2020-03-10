#### This file has a list of common functions in alphabetical order. These functions include:

# common.acceptrates1 - update proposal variance for a MH step based on having no max limit on the proposal var.
# common.acceptrates2 - update proposal variance for a MH step based on having  a max limit on the proposal var.
# common.betablock - Create the blocking structure for beta.
# common.betatransform - back transform the regression parameters to the original scale.
# common.burnin.nsample.thin.check - check the burnin, n.sample, thin arugments.
# common.frame - check the frame argument.
# common.frame.localised - check the frame argument for the localised model.
# common.modelfit - compute the model fit criteria.s
# common.prior.beta.check - Check the prior entered for beta.
# common.prior.var.check - check the prior entered for variance parameters.
# common.prior.varmat.check - check the prior entered for variance matrix parameters.
# common.verbose - check the verbose argument.
# common.Wcheckformat - check the W matrix.
# common.Wcheckformat.disimilarity - check the W matrix for the dissimilarity model.



#### Acceptance rates - no maximum limit on the proposal sd
common.accceptrates1 <- function(accept, sd, min, max)
{
    #### Update the proposal standard deviations
    rate <- 100 * accept[1] / accept[2]
    
    if(rate > max)
    {
        sd <- sd + 0.1 * sd
    }else if(rate < min)              
    {
        sd <- sd - 0.1 * sd
    }else
    {
    }
    
    return(sd)
}



#### Acceptance rates - maximum limit on the proposal sd
 common.accceptrates2 <- function(accept, sd, min, max, sd.max)
{
    #### Update the proposal standard deviations
    rate <- 100 * accept[1] / accept[2]
    
    if(rate > max)
    {
        sd <- sd + 0.1 * sd
        sd[which(sd>sd.max)] <- sd.max
    }else if(rate < min)              
    {
        sd <- sd - 0.1 * sd
    }else
    {
    }
    
    return(sd)
}



#### Beta blocking
common.betablock <- function(p, blocksize.beta=NULL)
{
    ## Compute the blocking structure for beta     
    if(is.null(blocksize.beta)) blocksize.beta <- 10 
       
    if(blocksize.beta >= p)
    {
        n.beta.block <- 1
        beta.beg <- 1
        beta.fin <- p
    }else
    {
        n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
        remainder <- p - n.standard * blocksize.beta
        
        if(remainder==0)
        {
            beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
            beta.fin <- seq(blocksize.beta, p, blocksize.beta)
            n.beta.block <- length(beta.beg)
        }else
        {
            beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
            beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
            n.beta.block <- length(beta.beg)
        }
    }
    
    return(list(beta.beg, beta.fin, n.beta.block))
}



#### beta back transform samples
common.betatransform <- function(samples.beta, X.indicator, X.mean, X.sd, p, localised)
{
    #### Back transform the beta values
    #### Slightly different code depending on whether the localised model is used
    samples.beta.orig <- samples.beta
    number.cts <- sum(X.indicator==1)
    
    if(localised)
    {
        #### Localised model    
        if(number.cts>0)
        {
            for(r in 1:p)
            {
                if(X.indicator[r]==1)
                {
                    samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
                }else
                {
                }
            }
        }else
        {
        }
    }else
    {
        #### Not the localised model
        if(number.cts>0)
        {
            for(r in 1:p)
            {
                if(X.indicator[r]==1)
                {
                    samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
                }else if(X.indicator[r]==2 & p>1)
                {
                    X.transformed <- which(X.indicator==1)
                    samples.temp <- as.matrix(samples.beta[ ,X.transformed])
                    for(s in 1:length(X.transformed))
                    {
                        samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
                    }
                    intercept.adjustment <- apply(samples.temp, 1,sum) 
                    samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
                }else
                {
                }
            }
        }else
        {
        }
    }
    
    #### Return the transformed samples
    return(samples.beta.orig)
}



#### Check MCMC arguments
common.burnin.nsample.thin.check <- function(burnin, n.sample, thin)
{
    #### Check for valid arguments for the burnin, n.sample and thin arguments
    if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
    if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(n.sample <= thin)  stop("thin is greater than n.sample.", call.=FALSE)
    if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
    if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
    if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
    
}


#### Read in and format the frame argument
common.frame <- function(formula, data, family)
{
    #### Overall formula object
    frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
    if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
    
    
    #### Design matrix
    ## Create the matrix
    X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    #if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
    
    n <- nrow(X)
    p <- ncol(X)
    
    ## Check for linearly related columns
    cor.X <- suppressWarnings(cor(X))
    diag(cor.X) <- 0
    if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    if(p>1)
    {
        if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
    }else
    {
    }
    
    ## Standardise the matrix
    X.standardised <- X
    X.sd <- apply(X, 2, sd)
    X.mean <- apply(X, 2, mean)
    X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
    
    for(j in 1:p)
    {
        if(length(table(X[ ,j]))>2)
        {
            X.indicator[j] <- 1
            X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
        }else if(length(table(X[ ,j]))==1)
        {
            X.indicator[j] <- 2
        }else
        {
            X.indicator[j] <- 0
        }
    }
    
    

    #### Response variable
    ## Create the response
    Y <- model.response(frame)
    J <- length(Y) / n
    which.miss <- matrix(as.numeric(!is.na(Y)), nrow=n, ncol=J)
        if(J==1) which.miss <- as.numeric(which.miss)
    n.miss <- n*J - sum(which.miss)

    
    ## Check for errors
    if(family=="binomial")
    {
        if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
        int.check <- n*J - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
        if(int.check > 0) stop("the response variable has non-integer values.", call.=FALSE)
        if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    }else if(family=="gaussian")
    {
        if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)    
    }else if(family=="poisson")
    {
            if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
            int.check <- n*J - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
            if(int.check > 0) stop("the response variable has non-integer values.", call.=FALSE)
            if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    }else if(family=="multinomial")
    {
        if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
        int.check <- n*J - n.miss - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
        if(int.check > 0) stop("the response variable has non-integer values.", call.=FALSE)
        if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    }else
    {}
    

    #### Offset variable
    offset <- try(model.offset(frame), silent=TRUE)
    #if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
        if(family=="multinomial")
        {
            if(is.null(offset))  offset <- array(0,c(n, (J-1)))
        }else
        {
            if(is.null(offset))  offset <- array(0,c(n, J))
        }
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    
    #### Return the values needed
    results <- list(n=n, p=p, X=X, X.standardised=X.standardised, X.sd=X.sd, X.mean=X.mean, X.indicator=X.indicator, 
                    offset=offset, Y=Y,  which.miss=which.miss, n.miss=n.miss)
    return(results)
}



#### Read in and format the frame argument from the localised model
common.frame.localised <- function(formula, data, family, trials)
{
    #### Overall formula object
    frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
    if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
    
    
    #### Response variable
    ## Create the response
    Y <- model.response(frame)
    n <- length(Y)
    
    ## Check for errors
    if(family=="binomial")
    {
        if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
        int.check <- n - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
        if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
        if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    }else if(family=="gaussian")
    {
        if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)    
    }else
    {
        if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
        int.check <- n - sum(ceiling(Y)==floor(Y), na.rm=TRUE)
        if(int.check > 0) stop("the response variable has non-integer values.", call.=FALSE)
        if(min(Y, na.rm=TRUE)<0) stop("the response variable has negative values.", call.=FALSE)
    }
    
    
    #### Offset variable
    offset <- try(model.offset(frame), silent=TRUE)
    #if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,n)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    
    #### Design matrix - Create and then adapt to remove the intercept term
    X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    #if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
    ptemp <- ncol(X)
    
    if(ptemp==1)
    {
        X <- NULL
        X.standardised <- NULL
        X.sd <- NULL
        X.mean <- NULL
        X.indicator <- NULL
        regression.vec <- rep(0, n)
        p <- 0
        beta <- NA
    }else
    {
        ## Check for linearly related columns
        cor.X <- suppressWarnings(cor(X))
        diag(cor.X) <- 0
        if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
        if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
        if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
        
        ## Remove the intercept term
        int.which <- which(apply(X,2,sd)==0)
        colnames.X <- colnames(X)
        X <- as.matrix(X[ ,-int.which])
        colnames(X) <- colnames.X[-int.which]
        p <- ncol(X)
        
        ## Standardise X
        X.standardised <- X
        X.sd <- apply(X, 2, sd)
        X.mean <- apply(X, 2, mean)
        X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
        
        for(j in 1:p)
        {
            if(length(table(X[ ,j]))>2)
            {
                X.indicator[j] <- 1
                X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
            }else
            {
                X.indicator[j] <- 0
            }
        }
        
        ## Compute a starting value for beta
            if(family=="binomial")
            {
            failures <- trials - Y
            mod.glm <- glm(cbind(Y, failures)~X.standardised, offset=offset, family="quasibinomial")
            beta.mean <- mod.glm$coefficients[-1]
            beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))[-1]
            beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
            regression.vec <- X.standardised %*% beta    
            }else
            {
            mod.glm <- glm(Y~X.standardised, offset=offset, family="quasipoisson")
            beta.mean <- mod.glm$coefficients[-1]
            beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))[-1]
            beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
            regression.vec <- X.standardised %*% beta    
            }
    }
    
    
    #### Return the values needed
    results <- list(n=n, p=p, X=X, X.standardised=X.standardised, X.sd=X.sd, X.mean=X.mean, X.indicator=X.indicator, 
                    offset=offset, Y=Y, regression.vec=regression.vec, beta=beta)
    return(results)
}


# Compute the DIC. WAIC,LMPL and loglikelihood
common.modelfit <- function(samples.loglike, deviance.fitted)
{
    #### WAIC
    p.w <- sum(apply(samples.loglike,2, var), na.rm=TRUE)
    mean.like <- apply(exp(samples.loglike),2,mean)
    mean.min <- min(mean.like[mean.like>0])
    mean.like[mean.like==0] <- mean.min
    lppd <- sum(log(mean.like), na.rm=TRUE)
    WAIC <- -2 * (lppd - p.w)
    
    
    #### Compute the Conditional Predictive Ordinate
    CPO <- 1/apply(exp(-samples.loglike), 2, mean)
    mean.min <- min(CPO[CPO>0])
    CPO[CPO==0] <- mean.min
    LMPL <- sum(log(CPO), na.rm=TRUE)    
    
    
    #### DIC
    mean.deviance <- -2 * sum(samples.loglike, na.rm=TRUE) /   nrow(samples.loglike)
    p.d <- mean.deviance - deviance.fitted
    DIC <- deviance.fitted + 2 * p.d
    
    
    #### loglikelihood
    loglike <- -0.5 * deviance.fitted
    
    
    #### Model fit criteria
    modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike)
    names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood")
    return(modelfit)  
}


#### Check beta prior arguments
common.prior.beta.check <- function(prior.mean.beta, prior.var.beta, p)
{
    ## Checks    
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)    
    
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
}



#### Check variance prior arguments
common.prior.var.check <- function(prior.var)
{
    ## Checks   
    if(length(prior.var)!=2) stop("the prior values for a variance parameter are the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var)) stop("the prior values for a variance parameter are not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var))!=0) stop("the prior values for a variance parameter have missing values.", call.=FALSE) 
}



#### Check variance matrix prior arguments
common.prior.varmat.check <- function(prior.varmat, J)
{
    if(nrow(prior.varmat)!=J) stop("prior.Sigma.scale is the wrong dimension.", call.=FALSE)    
    if(ncol(prior.varmat)!=J) stop("prior.Sigma.scale is the wrong dimension.", call.=FALSE)    
    if(!is.numeric(prior.varmat)) stop("prior.Sigma.scale has non-numeric values.", call.=FALSE)    
    if(sum(is.na(prior.varmat))!=0) stop("prior.Sigma.scale has missing values.", call.=FALSE)    
    if(!is.positive.definite(prior.varmat)) stop("prior.Sigma.scale is not a positive definite matrix.", call.=FALSE)
    if(!is.symmetric.matrix(prior.varmat)) stop("prior.Sigma.scale is not symmetric.", call.=FALSE)
}



#### Check the verbose option
common.verbose <- function(verbose)
{
    if(is.null(verbose)) verbose=TRUE     
    if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)
    
    if(verbose)
    {
        cat("Setting up the model.\n")
        a<-proc.time()
    }else{
        a <- 1    
    }
    return(a)
}



#### Check the W matrix
common.Wcheckformat <- function(W)
{
    #### Check W is a matrix of the correct dimension
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    n <- nrow(W)
    if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)    
    
    
    #### Check validity of inputed W matrix
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
    
    
    #### Create the triplet form
    W.triplet <- c(NA, NA, NA)
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            if(W[i,j]>0)
            {
                W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
            }else{}
        }
    }
    W.triplet <- W.triplet[-1, ]     
    n.triplet <- nrow(W.triplet) 
    W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
    n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
    
    
    #### Create the start and finish points for W updating
    W.begfin <- array(NA, c(n, 2))     
    temp <- 1
    for(i in 1:n)
    {
        W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
        temp <- temp + n.neighbours[i]
    }
    
    
    #### Return the critical quantities
    results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, n=n)
    return(results)   
}


#### Check the W matrix - Dissimilarity model
common.Wcheckformat.disimilarity <- function(W)
{
    #### Check validity of inputed W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
    n <- nrow(W)
    if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)   
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(sum(as.numeric(W)==0) + sum(as.numeric(W)==1) - n^2 !=0) stop("W has non-binary elements", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
    
    
    ## Ensure the W matrix is symmetric
    Wnew <- array(0, c(n,n))
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            if(i>j)
            {
                temp <- W[i,j]
                Wnew[i,j] <- temp
                Wnew[j,i] <- temp
            }else{}
        }
    }
    W <- Wnew  
    n.neighbours <- apply(W, 2, sum)
    spam.W <- as.spam(W)

    
    #### Create the triplet form
    W.triplet <- c(NA, NA, NA)
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            if(W[i,j]>0)
            {
                W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
            }else{}
        }
    }
    W.triplet <- W.triplet[-1, ]     
    n.triplet <- nrow(W.triplet) 
    W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
    n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
    
    
    #### Create the start and finish points for W updating
    W.begfin <- array(NA, c(n, 2))     
    temp <- 1
    for(i in 1:n)
    {
        W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
        temp <- temp + n.neighbours[i]
    }
    
    
    #### Return the critical quantities
    results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, spam.W=spam.W, n=n)
    return(results)   
}
