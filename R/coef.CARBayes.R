coef.CARBayes <- function(object,...)
{
    #### Return the estimated regression coefficient
    if(is.null(nrow(object$samples$beta)))
    {
        return(NULL)  
    }else
    {
    beta <- apply(object$samples$beta, 2, median)
    names(beta) <- colnames(object$X)    
    return(beta)
    }
}