residuals.CARBayes <- function(object, type="pearson", ...)
{
residuals <- object$residuals    
 

#### The multivariate models provides lists the univariate models provide matrices

    if(class(residuals)=="list")
    {
    #### Return one of two types of residuals
        if(type=="response")
        {
        return(residuals$response)
        }else if(type=="pearson")
        {
            return(residuals$pearson)
        }else
        {
            return("Error. That is not one of the allowable residual types.")   
        }
        
    }else
    {
    #### Return one of two types of residuals
        if(type=="response")
        {
        return(residuals$response)
        }else if(type=="pearson")
        {
        return(residuals$pearson)
        }else
        {
        return("Error. That is not one of the allowable residual types.")   
        }
    }
}