print.CARBayes <- function(x,...)
{
#### Check for missingness
    if(length(x$samples$Y)==1)
    {
    n.miss <- 0    
    }else
    {
    n.miss <- ncol(x$samples$Y)
    }
    
    if(class(x$localised.structure)=="list")
    {
        #### Print out the model fitted
        cat("\n#################\n")
        cat("#### Model fitted\n")
        cat("#################\n")
        cat(x$model)
        cat("Regression equation - ")
        print(x$formula)
        cat("Number of missing observations - ")
        cat(n.miss)
        cat("\n")
        
        #### Print out the results
        cat("\n############\n")
        cat("#### Results\n")
        cat("############\n")
        cat("Posterior quantities and DIC\n\n")
        print(x$summary.results[ ,-c(4,5)])
        cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", round(x$modelfit[5],2),"\n")

            if(length(x$localised.structure[[2]])>1)
            {
            cat("\nThe number of stepchanges identified in the random effect surface\n")
            temp <- x$localised.structure[[1]][!is.na(x$localised.structure[[1]])]
            tab <- array(NA, c(1,2))
            tab[1, ] <- c(sum(temp)/2, length(temp)/2- sum(temp)/2)
            colnames(tab) <- c("no stepchange", "stepchange")
            print(tab)
            }else
            {}
    }else if(class(x$localised.structure)=="numeric")
    {
        #### Print out the model fitted
        cat("\n#################\n")
        cat("#### Model fitted\n")
        cat("#################\n")
        cat(x$model)
        cat("Regression equation - ")
        print(x$formula)
        cat("Number of missing observations - ")
        cat(n.miss)
        cat("\n")
        
        #### Print out the results
        cat("\n############\n")
        cat("#### Results\n")
        cat("############\n")
        cat("Posterior quantities and DIC\n\n")
        print(x$summary.results[ ,-c(4,5)])
        cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", round(x$modelfit[5],2),"\n")
        cat("\nNumber of clusters with the number of data points in each one\n")
        print(table(paste("group", x$localised.structure, sep="")))
        
    }else
    {
        #### Print out the model fitted
        cat("\n#################\n")
        cat("#### Model fitted\n")
        cat("#################\n")
        cat(x$model)
            if(class(x$formula)=="formula")
            {
            cat("Regression equation - ")
            print(x$formula)    
            }else
            {
            cat("Regression equation - ")
            print(x$formula[[1]])    
            cat("Zero probability equation - ")
            print(x$formula[[2]])   
            }
        
        cat("Number of missing observations - ")
        cat(n.miss)
        cat("\n")
        
        #### Print out the results
        cat("\n############\n")
        cat("#### Results\n")
        cat("############\n")
        cat("Posterior quantities and DIC\n\n")
        print(x$summary.results[ ,-c(4,5)])
        cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "     ", "LMPL = ", round(x$modelfit[5],2),"\n")
     }
return(invisible(x))        
}



