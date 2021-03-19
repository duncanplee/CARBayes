W.optimise <- function(W, spdata, add=FALSE, remove=TRUE, remove_first=FALSE)
{
#### This function is a wrapper for the C++ function that optimises the neighbourhood matrix W to a given data set.

    
###################################################
#### Perform checks to prevent inappropriate inputs
###################################################
    if(class(W)[1]!="matrix") stop("W is not a matrix.", call.=FALSE)
    if(sum(is.na(W))>0) stop("W contains NA values.", call.=FALSE)  
    if(nrow(W)!=ncol(W)) stop("W is not a square matrix.", call.=FALSE)   
    if(length(names(table(W)))!=2) stop("W contains more than two distinct values.", call.=FALSE) 
    if(sum(as.numeric(names(table(W))) != c(0,1))>0) stop("W contains elements other than 0s and 1s.", call.=FALSE)     
    if(sum(diag(W))>0) stop("The diagonal of W contains non-zero elements.", call.=FALSE)   

    if(length(spdata)!=nrow(W)) stop("W and spdata do not have matching dimensions.", call.=FALSE)
    if(!is.numeric(spdata)) stop("spdata is not a numeric vector.", call.=FALSE)   
    if(sum(is.na(spdata))>0) stop("spdata contains NA values.", call.=FALSE)  
    if(sum(spdata==Inf)>0) stop("spdata contains Inf values.", call.=FALSE)  
    if(sum(spdata==-Inf)>0) stop("spdata contains -Inf values.", call.=FALSE)      
 
       
######################################    
#### Optimise the neighbourhood matrix
###################################### 
#### Optimise W
K <- nrow(W)
W.temp <- optimise_graph(W, spdata, add=add, remove=remove, remove_first=remove_first)


#### Create the matrix version
W.est <- array(NA, c(K,K))
    for(k in 1:K)
    {
    W.est[k, ] <- W.temp[[k]]  
    }


#### Return the result
return(W.est)
}