#include <Rcpp.h>
using namespace Rcpp;

// This file contains the following functions:

// linpredcompute - computing the linear predictor for covariates.
// quadform - computing quadratic forms phi %*% Q %*% theta.
// binomialbetaupdateMALA - update regression parameters in the binomial model using MALA
// binomialbetaupdateRW - update regression parameters in the binomial model using RW
// binomialcarupdateMALA - update random effects in the binomial model using MALA
// binomialcarupdateRW - update random effects in the binomial model using RW
// binomialindepupdateMALA - update the independent effects in the binomial model using MALA
// binomialindepupdateRW - update the independent effects in the binomial model using RW
// poissonbetaupdateMALA - update regression parameters in the poisson model using MALA
// poissonbetaupdateRW - update regression parameters in the poisson model using RW
// poissoncarupdateMALA - update random effects in the poisson model using MALA
// poissoncarupdateRW - update random effects in the poisson model using RW
// poissonindepupdateMALA - update the independent effects in the poisson model using MALA
// poissonindepupdateRW - update the independent effects in the poisson model using RW
// zipcarupdateRW - update the random effects in the zip models using RW
// zipcarupdateMALA - update the random effects in the zip models using MALA
// zipindepupdateRW - update the independent random effects in the zip models using RW
// zipindepupdateMALA - update the independent random effects in the zip models using MALA
// gaussiancarupdate - update random effects in the Gaussian model
// binomialmcarupdateMALA - update random effects in the binomial MCAR model using MALA
// binomialmcarupdateRW - update random effects in the binomial MCAR model using RW
// poissonmcarupdateMALA - update random effects in the poisson MCAR model using MALA
// poissonmcarupdateRW - update random effects in the poisson MCAR model using RW
// gaussianmcarupdateRW - update random effecs in the Gaussian MCAR model usng RW
// gaussianmcarupdateMALA - update random effects in the Gaussian MCAR model using MALA
// multinomialbetaupdateRW - update beta in the multinomial model using RW
// poissoncarmultilevelupdate - Poisson spatial random effects updates
// binomialcarmultilevelupdate - binomial spatial random effects updates
// gaussiancarmultilevelupdate - gaussian spatial random effects updates
// gaussiancarmultilevelupdateindiv - Gaussian indep random effect updates
// poissoncarmultilevelupdateindiv - Poisson indep random effect updates
// binomialcarmultilevelupdateindiv - binomial indep random effect updates


// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, 
                          NumericVector beta, NumericVector offset)
{
//Create new objects
// Compute the linear predictor
NumericVector linpred(nsites);
double temp; 


//  Compute the linear predictor via a double for loop
     for(int j = 0; j < nsites; j++)
     {
     temp = 0;
      
          for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];     
          
     linpred[j] = temp + offset[j];  
     }


// Return the result
return linpred;
}




// [[Rcpp::export]]
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                    NumericVector phi, NumericVector theta, double rho)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;
int row, col;
   
   
// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_triplet; l++)
     {
     row = Wtriplet(l,0) - 1;
     col = Wtriplet(l,1) - 1;
     tau2_quadform = tau2_quadform + phi[(Wtriplet(l,0) - 1)] * theta[(Wtriplet(l,1) - 1)] * Wtriplet(l,2); 
     }
 
 
 // Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     tau2_phisq = tau2_phisq + phi[l] * theta[l] * (rho * Wtripletsum[l] + 1 - rho);    
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}



// [[Rcpp::export]]
List binomialcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                       NumericVector Wtripletsum,const int nsites, NumericVector phi, double tau2, 
                       const NumericVector y, const NumericVector failures, NumericVector trials, const double phi_tune, 
                       double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0, rowstart=0, rowend=0;
    double acceptance, acceptance1, acceptance2,  sumphi, mala_old, mala_new;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, pold, pnew, proposal_var;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;

        for(int j = 0; j < nsites; j++)
        {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom;  
        

        // propose a value
        proposal_var = priorvar * phi_tune;
        mala_old = phinew[j] + 0.5 * proposal_var * (y[j] - (trials[j] * exp(phinew[j] + offset[j])) / (1 + exp(phinew[j] + offset[j])) - (phinew[j] - priormean) /priorvar);
        propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
            
        // Accept or reject it
        // Full conditional ratio
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        pold = exp(offset[j] + phinew[j]) / (1 + exp(offset[j] + phinew[j]));
        pnew = exp(offset[j] + propphi) / (1 + exp(offset[j] + propphi));
        oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
        newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
        acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

        // Proposal distribution ratio
        mala_new = propphi + 0.5 * proposal_var * (y[j] - (trials[j] * exp(propphi + offset[j])) / (1 + exp(propphi + offset[j])) - (propphi - priormean) /priorvar);
        acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew[j] - mala_new),2) - pow((propphi-mala_old),2)));
        acceptance = acceptance1 * acceptance2;
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                phinew[j] = propphi;
                accept = accept + 1;
            }
            else
            { 
            }    
    }

// Return the results
List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
List binomialcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                       NumericVector Wtripletsum,const int nsites, NumericVector phi, double tau2, 
                       const NumericVector y, const NumericVector failures, const double phi_tune, 
                       double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0, rowstart=0, rowend=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, pold, pnew, proposal_var;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom;  
        
        // propose a value
        proposal_var = priorvar * phi_tune;
        propphi = rnorm(1, phinew[j], sqrt(proposal_var))[0];
            
        // Accept or reject it
        // Full conditional ratio
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        pold = exp(offset[j] + phinew[j]) / (1 + exp(offset[j] + phinew[j]));
        pnew = exp(offset[j] + propphi) / (1 + exp(offset[j] + propphi));
        oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
        newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            

            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                phinew[j] = propphi;
                accept = accept + 1;
            }
            else
            { 
            }    
        }
    
    // Return the results
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List binomialbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                            NumericVector offset, NumericVector y,  NumericVector failures,
                            NumericVector trials, NumericVector prior_meanbeta, 
                            NumericVector prior_varbeta, const int nblock,
                            double beta_tune, List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = y -  trials * exp(lp_current) / (1 + exp(lp_current));
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals  
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
            p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
            oldlikebit = oldlikebit + y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
            newlikebit = newlikebit + y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
        }
        likebit = newlikebit - oldlikebit;
        
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y -  trials * exp(lp_proposal) / (1 + exp(lp_proposal));
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject hte proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}


// [[Rcpp::export]]
List binomialbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                          NumericVector offset, NumericVector y,  NumericVector failures,
                          NumericVector prior_meanbeta, NumericVector prior_varbeta, 
                          const int nblock,double beta_tune, List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }
  
  // Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    for(int g=0; g<len; g++)
    {
      beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }
    
    
    // Compute the acceptance ratio - full conditionals  
    oldlikebit = 0;
    newlikebit=0;
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
      p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
      p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
      oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
      newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }
    
    
    // Accept or reject hte proposal      
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }
  
  
  // Compute the acceptance probability and return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}








// [[Rcpp::export]]
List binomialindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
               const NumericVector failures, const NumericVector trials, const double theta_tune,  NumericVector offset)
{
// Update the independent random effects 
//Create new objects
int accept=0;
double acceptance, acceptance1, acceptance2, mala_old, mala_new;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double proptheta, pold, pnew;
NumericVector thetanew(nsites);
 
   
//  Update each random effect in turn
thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
    // propose a value
    mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - (trials[j] * exp(thetanew[j] + offset[j])) / (1 + exp(thetanew[j] + offset[j])) - thetanew[j] / sigma2);
    proptheta = rnorm(1, mala_old, theta_tune)[0];
        
    // Accept or reject it
    // Full conditional ratio
    newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
    oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
        
    pold = exp(offset[j] + thetanew[j]) / (1 + exp(offset[j] + thetanew[j]));
    pnew = exp(offset[j] + proptheta) / (1 + exp(offset[j] + proptheta));
    oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
    newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
    acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

    // Proposal distribution ratio
    mala_new = proptheta + 0.5 * pow(theta_tune, 2) * (y[j] - (trials[j] * exp(proptheta + offset[j])) / (1 + exp(proptheta + offset[j])) - proptheta / sigma2);
    acceptance2 = exp(-(0.5 / pow(theta_tune, 2)) * (pow((thetanew[j] - mala_new),2) - pow((proptheta-mala_old),2)));
    acceptance = acceptance1 * acceptance2;
        
    // Acceptace or reject the proposal
        if(runif(1)[0] <= acceptance) 
        {
        thetanew[j] = proptheta;
        accept = accept + 1;
        }
        else
        { 
        }    
    }

List out(2);
out[0] = thetanew;
out[1] = accept;
return out;
}


// [[Rcpp::export]]
List binomialindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                         const NumericVector failures, const double theta_tune,  NumericVector offset)
{
    // Update the independent random effects 
    //Create new objects
    int accept=0;
    double acceptance;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, pold, pnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
    // propose a value
    proptheta = rnorm(1, thetanew[j], theta_tune)[0];
            
    // Accept or reject it
    // Full conditional ratio
    newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
    oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            
    pold = exp(offset[j] + thetanew[j]) / (1 + exp(offset[j] + thetanew[j]));
    pnew = exp(offset[j] + proptheta) / (1 + exp(offset[j] + proptheta));
    oldlikebit = y[j] * log(pold) + failures[j] * log((1-pold));
    newlikebit = y[j] * log(pnew) + failures[j] * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            

    // Acceptace or reject the proposal
        if(runif(1)[0] <= acceptance) 
        {
        thetanew[j] = proptheta;
        accept = accept + 1;
        }
        else
        { 
        }    
    }
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List poissonindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
               const double theta_tune,  NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance, acceptance1, acceptance2, mala_old, mala_new;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double proptheta, lpold, lpnew;
NumericVector thetanew(nsites);
 
   
//  Update each random effect in turn
thetanew = theta;
     for(int j = 0; j < nsites; j++)
     {
    // propose a value
     mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - exp(thetanew[j] + offset[j]) - thetanew[j] / sigma2);
     proptheta = rnorm(1, mala_old, theta_tune)[0];
             
     // Accept or reject it
     // Full conditional ratio
     newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
     oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
     lpold = offset[j] + thetanew[j];
     lpnew = offset[j] + proptheta;
     oldlikebit = y[j] * lpold - exp(lpold);
     newlikebit = y[j] * lpnew - exp(lpnew);
     acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
             
     // Proposal distribution ratio
     mala_new = proptheta + 0.5 * pow(theta_tune, 2) * (y[j] - exp(proptheta + offset[j]) - proptheta / sigma2);
     acceptance2 = exp(-(0.5 / pow(theta_tune, 2)) * (pow((thetanew[j] - mala_new),2) - pow((proptheta-mala_old),2)));
     acceptance = acceptance1 * acceptance2;
             
     // Acceptace or reject the proposal
         if(runif(1)[0] <= acceptance) 
         {
         thetanew[j] = proptheta;
         accept = accept + 1;
         }
         else
         { 
         }    
    }

List out(2);
out[0] = thetanew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
List poissonindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                            const double theta_tune,  NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0;
    double acceptance;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, lpold, lpnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
    // propose a value
    proptheta = rnorm(1, thetanew[j], theta_tune)[0];
    
    // Accept or reject it
    // Full conditional ratio
    newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
    oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
    lpold = offset[j] + thetanew[j];
    lpnew = offset[j] + proptheta;
    oldlikebit = y[j] * lpold - exp(lpold);
    newlikebit = y[j] * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);

    // Acceptace or reject the proposal
        if(runif(1)[0] <= acceptance) 
        {
        thetanew[j] = proptheta;
        accept = accept + 1;
        }
        else
        { 
        }    
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List poissonbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                       NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                       NumericVector prior_varbeta, const int nblock,
                       double beta_tune, List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = y - exp(lp_current);
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals 
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
            newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
        }
        likebit = newlikebit - oldlikebit;
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y - exp(lp_proposal);
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject hte proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}



// [[Rcpp::export]]
List poissonbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                         NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                         NumericVector prior_varbeta, const int nblock, double beta_tune, 
                         List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }

  
// Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    for(int g=0; g<len; g++)
    {
      beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }

    
    // Compute the acceptance ratio - likelihood part
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    oldlikebit = 0;
    newlikebit=0;
    for(int j = 0; j < nsites; j++)     
    {
      oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
      newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
    }
    likebit = newlikebit - oldlikebit;
    
    
    // Compute the acceptance ratio - prior part
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }

    
    // Accept or reject the proposal    
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }


  // Return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}




// [[Rcpp::export]]
List poissoncarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                      NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                      double tau2, const NumericVector y, const double phi_tune, 
                      double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, acceptance1, acceptance2, sumphi, proposal_var, mala_old, mala_new;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value
        proposal_var = priorvar * phi_tune;
        mala_old = phinew[j] + 0.5 * proposal_var * (y[j] - exp(phinew[j] + offset[j]) - (phinew[j] - priormean) /priorvar);
        propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
            
        // Accept or reject it
        // Full conditional ratio
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        lpold = offset[j] + phinew[j];
        lpnew = offset[j] + propphi;
        oldlikebit = y[j] * lpold - exp(lpold);
        newlikebit = y[j] * lpnew - exp(lpnew);
        acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
        // Proposal distribution ratio
        mala_new = propphi + 0.5 * proposal_var * (y[j] - exp(propphi + offset[j]) - (propphi - priormean) /priorvar);
        acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew[j] - mala_new),2) - pow((propphi-mala_old),2)));
        acceptance = acceptance1 * acceptance2;
            
        // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
            phinew[j] = propphi;
            accept = accept + 1;
            }
            else
            { 
            }    
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List poissoncarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                          NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                          double tau2, const NumericVector y, const double phi_tune, 
                          double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi, proposal_var;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value
        proposal_var = priorvar * phi_tune;
        propphi = rnorm(1, phinew[j], sqrt(proposal_var))[0];
            
        // Accept or reject it
        // Full conditional ratio
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        lpold = offset[j] + phinew[j];
        lpnew = offset[j] + propphi;
        oldlikebit = y[j] * lpold - exp(lpold);
        newlikebit = y[j] * lpnew - exp(lpnew);
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
        // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
            phinew[j] = propphi;
            accept = accept + 1;
            }
            else
            { 
            }    
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List zipcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                    NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                    double tau2, const NumericVector y, const double phi_tune, 
                    double rho, NumericVector offset, NumericVector poiind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi, proposal_var;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // Different updates depending on whether the y[j] is missing or not.
        if(poiind[j]==1)
        {
            // propose a value
            proposal_var = priorvar * phi_tune;
            propphi = rnorm(1, phinew[j], sqrt(proposal_var))[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
            lpold = offset[j] + phinew[j];
            lpnew = offset[j] + propphi;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
            acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                phinew[j] = propphi;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            phinew[j] = rnorm(1, priormean, sqrt(priorvar))[0];    
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List zipcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                      NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                      double tau2, const NumericVector y, const double phi_tune, 
                      double rho, NumericVector offset, NumericVector poiind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, acceptance1, acceptance2, sumphi, proposal_var, mala_old, mala_new;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // Different updates depending on whether the y[j] is missing or not.
        if(poiind[j]==1)
        {
            // propose a value
            proposal_var = priorvar * phi_tune;
            mala_old = phinew[j] + 0.5 * proposal_var * (y[j] - exp(phinew[j] + offset[j]) - (phinew[j] - priormean) /priorvar);
            propphi = rnorm(1, mala_old, sqrt(proposal_var))[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
            oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
            lpold = offset[j] + phinew[j];
            lpnew = offset[j] + propphi;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Proposal distribution ratio
            mala_new = propphi + 0.5 * proposal_var * (y[j] - exp(propphi + offset[j]) - (propphi - priormean) /priorvar);
            acceptance2 = exp(-(0.5 / proposal_var) * (pow((phinew[j] - mala_new),2) - pow((propphi-mala_old),2)));
            acceptance = acceptance1 * acceptance2;
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                phinew[j] = propphi;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            phinew[j] = rnorm(1, priormean, sqrt(priorvar))[0];    
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List zipindepupdateRW(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                      const double theta_tune, NumericVector offset, NumericVector poiind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0;
    double acceptance;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, lpold, lpnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    
    for(int j = 0; j < nsites; j++)
    {
        // Different updates depending on whether the y[j] is missing or not.
        if(poiind[j]==1)
        {
            // propose a value
            proptheta = rnorm(1, thetanew[j], theta_tune)[0];
            
            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            lpold = offset[j] + thetanew[j];
            lpnew = offset[j] + proptheta;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
            acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                thetanew[j] = proptheta;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            thetanew[j] = rnorm(1, 0, sqrt(sigma2))[0];    
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List zipindepupdateMALA(const int nsites, NumericVector theta, double sigma2, const NumericVector y, 
                       const double theta_tune, NumericVector offset, NumericVector poiind)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0;
    double acceptance, acceptance1, acceptance2, mala_old, mala_new;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double proptheta, lpold, lpnew;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    
    for(int j = 0; j < nsites; j++)
    {
    // Different updates depending on whether the y[j] is missing or not.
        if(poiind[j]==1)
        {
            // propose a value
            mala_old = thetanew[j] + 0.5 * pow(theta_tune, 2) * (y[j] - exp(thetanew[j] + offset[j]) - thetanew[j] / sigma2);
            proptheta = rnorm(1, mala_old, theta_tune)[0];

            // Accept or reject it
            // Full conditional ratio
            newpriorbit = (0.5/sigma2) * pow(proptheta, 2); 
            oldpriorbit = (0.5/sigma2) * pow(thetanew[j], 2);
            lpold = offset[j] + thetanew[j];
            lpnew = offset[j] + proptheta;
            oldlikebit = y[j] * lpold - exp(lpold);
            newlikebit = y[j] * lpnew - exp(lpnew);
            acceptance1 = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
            
            // Proposal distribution ratio
            mala_new = proptheta + 0.5 * pow(theta_tune, 2) * (y[j] - exp(proptheta + offset[j]) - proptheta / sigma2);
            acceptance2 = exp(-(0.5 / pow(theta_tune, 2)) * (pow((thetanew[j] - mala_new),2) - pow((proptheta-mala_old),2)));
            acceptance = acceptance1 * acceptance2;
            
            // Acceptace or reject the proposal
            if(runif(1)[0] <= acceptance) 
            {
                thetanew[j] = proptheta;
                accept = accept + 1;
            }
            else
            { 
            }    
        }else
        {
            thetanew[j] = rnorm(1, 0, sqrt(sigma2))[0];    
        }
    }
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}





// [[Rcpp::export]]
NumericVector gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, double tau2, 
     double rho, double nu2, NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int rowstart=0, rowend=0;
double sumphi;
double fcprecision, fcsd, fcmean;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);


//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
      
      // propose a value  
      fcprecision = (1/nu2) + (1/priorvar);
      fcsd = pow((1/fcprecision),0.5);
      fcmean = (priormean / priorvar + offset[j]) / fcprecision;
      phinew[j] = rnorm(1, fcmean, fcsd)[0];      
      }

return phinew;
}



// [[Rcpp::export]]
List binomialmcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                        const int nsites,  const int nvar, NumericMatrix phi, 
                        NumericMatrix Y, NumericMatrix failures, NumericMatrix trials,
                        NumericMatrix phioffset, NumericVector denoffset,  
                        NumericMatrix Sigmainv, double rho, double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar), mala1(nvar), mala2(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector pold(nvar), pnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance, hastings;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {      
    // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
    fcmean = rho * sumphi / denoffset[j]; 
        
        
    // Generate the proposal distribution mean and propose a value
        for(int r=0; r<nvar; r++)
        {
        mala1[r] =  phi(j,r) + 0.5 * pow(phi_tune, 2) * (Y(j,r) - (trials(j,r) * exp(phi(j,r) + phioffset(j,r))) / (1 + exp(phi(j,r) + phioffset(j,r))) -sum(fcprec(r,_) * (phi(j,_) - fcmean)));
        propphi[r] = rnorm(1, mala1[r], phi_tune)[0];
        }

     // Generate the mala mean in reverse
        for(int r=0; r<nvar; r++)
        {
        mala2[r] =  propphi[r] + 0.5 * pow(phi_tune, 2) * (Y(j,r) - (trials(j,r) * exp(propphi[r] + phioffset(j,r))) / (1 + exp(propphi[r] + phioffset(j,r))) -sum(fcprec(r,_) * (propphi - fcmean)));
        }  
            
    // Compute the prior ratio
    diffcurrent = phi(j,_) - fcmean;
    diffprop = propphi - fcmean;
        for(int r=0; r<nvar; r++)
        {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
              
    // Likelihood ratio
    pold = exp(phioffset(j,_) + phi(j,_)) / (1 + exp(phioffset(j,_) + phi(j,_)));
    pnew = exp(phioffset(j,_) + propphi) / (1 + exp(phioffset(j,_) + propphi));
    oldlikebit = sum(Y(j,_) * log(pold) + failures(j,_) * log(1 - pold));
    newlikebit = sum(Y(j,_) * log(pnew) + failures(j,_) * log(1 - pnew));

    // Hastings ratio   
    hastings = - (sum(pow(phi(j,_) - mala2,2)) - sum(pow(propphi - mala1,2))) / (2*pow(phi_tune,2));
     
     
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit + hastings);
    if(runif(1)[0] <= acceptance) 
        {
        phi(j,_) = propphi;
        accept = accept + 1;
        }
        else
        { 
        }
    }     
            
            
    // Return the results
    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}





// [[Rcpp::export]]
List binomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                        const int nsites,  const int nvar, NumericMatrix phi, 
                        NumericMatrix Y, NumericMatrix failures,
                        NumericMatrix phioffset, NumericVector denoffset,  
                        NumericMatrix Sigmainv, double rho, double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector pold(nvar), pnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {      
        // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        
        // Generate the proposal distribution mean and propose a value
        for(int r=0; r<nvar; r++)
        {
            propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
        }
        

        // Compute the prior ratio
        diffcurrent = phi(j,_) - fcmean;
        diffprop = propphi - fcmean;
        for(int r=0; r<nvar; r++)
        {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        pold = exp(phioffset(j,_) + phi(j,_)) / (1 + exp(phioffset(j,_) + phi(j,_)));
        pnew = exp(phioffset(j,_) + propphi) / (1 + exp(phioffset(j,_) + propphi));
        oldlikebit = sum(Y(j,_) * log(pold) + failures(j,_) * log(1 - pold));
        newlikebit = sum(Y(j,_) * log(pnew) + failures(j,_) * log(1 - pnew));
        

        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phi(j,_) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    
    // Return the results
    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List poissonmcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                       const int nsites,  const int nvar, NumericMatrix phi, 
                       NumericMatrix Y, NumericMatrix phioffset, 
                       NumericVector denoffset, NumericMatrix Sigmainv, double rho, 
                       double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar), mala1(nvar), mala2(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance, hastings;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {     
    // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
    // Generate the proposal distribution mean and propose a value
        for(int r=0; r<nvar; r++)
        {
        mala1[r] =  phi(j,r) + 0.5 * pow(phi_tune, 2) * (Y(j,r) - exp(phi(j,r) + phioffset(j,r)) -sum(fcprec(r,_) * (phi(j,_) - fcmean)));
        propphi[r] = rnorm(1, mala1[r], phi_tune)[0];
        }
        
    // Generate the mala mean in reverse
        for(int r=0; r<nvar; r++)
        {
        mala2[r] =  propphi[r] + 0.5 * pow(phi_tune, 2) * (Y(j,r) - exp(propphi[r] + phioffset(j,r)) - sum(fcprec(r,_) * (propphi - fcmean)));
        }       
        
         
    // Compute the prior ratio
    diffcurrent = phi(j,_) - fcmean;
    diffprop = propphi - fcmean;
        for(int r=0; r<nvar; r++)
        {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
    // Likelihood ratio
    lpold = phioffset(j,_) + phi(j,_);
    lpnew = phioffset(j,_) + propphi;
    oldlikebit = sum(Y(j,_) * lpold - exp(lpold));
    newlikebit = sum(Y(j,_) * lpnew - exp(lpnew));
    
    // Hastings ratio   
    hastings = - (sum(pow(phi(j,_) - mala2,2)) - sum(pow(propphi - mala1,2))) / (2*pow(phi_tune,2));

    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit + hastings);
        if(runif(1)[0] <= acceptance) 
        {
        phi(j,_) = propphi;
        accept = accept + 1;
        }
        else
        { 
        }
    }     

    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List poissonmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                           const int nsites,  const int nvar, NumericMatrix phi, 
                           NumericMatrix Y, NumericMatrix phioffset, 
                           NumericVector denoffset, NumericMatrix Sigmainv, double rho, 
                           double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {     
        // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        // Generate the proposal distribution mean and propose a value
        for(int r=0; r<nvar; r++)
        {
            propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
        }
        
        
        // Compute the prior ratio
        diffcurrent = phi(j,_) - fcmean;
        diffprop = propphi - fcmean;
        for(int r=0; r<nvar; r++)
        {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        lpold = phioffset(j,_) + phi(j,_);
        lpnew = phioffset(j,_) + propphi;
        oldlikebit = sum(Y(j,_) * lpold - exp(lpold));
        newlikebit = sum(Y(j,_) * lpnew - exp(lpnew));
        

        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phi(j,_) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
List gaussianmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                         const int nsites,  const int nvar, NumericMatrix phi, 
                         NumericMatrix phioffset, NumericVector denoffset, 
                         NumericMatrix Sigmainv, double rho, NumericVector nu2,
                         double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {     
        // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        // Generate a proposal value
        for(int r=0; r<nvar; r++)
        {
            propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
        }
        
        
        // Compute the prior ratio
        diffcurrent = phi(j,_) - fcmean;
        diffprop = propphi - fcmean;
        for(int r=0; r<nvar; r++)
        {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        lpold = pow((phioffset(j,_) - phi(j,_)),2);
        lpnew = pow((phioffset(j,_) - propphi),2);
        oldlikebit = 0.5 * sum(lpold / nu2);
        newlikebit = 0.5 * sum(lpnew / nu2);                             

        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit + oldlikebit - newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phi(j,_) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List gaussianmcarupdateMALA(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                          const int nsites,  const int nvar, NumericMatrix phi, 
                          NumericMatrix phioffset, NumericVector denoffset, 
                          NumericMatrix Sigmainv, double rho, NumericVector nu2,
                          double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec(nvar, nvar);
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar), mala1(nvar), mala2(nvar);
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance, hastings;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {     
        // Calculate the prior precision and mean
        for(int r=0; r<nvar; r++)
        {
            fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = rep(0,nvar);
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        
        // Generate the proposal distribution mean and propose a value
            for(int r=0; r<nvar; r++)
            {
            mala1[r] =  phi(j,r) + 0.5 * pow(phi_tune, 2) * (phioffset(j,r) - phi(j,r) - sum(fcprec(r,_) * (phi(j,_) - fcmean)));
            propphi[r] = rnorm(1, mala1[r], phi_tune)[0];
            }
        
        // Generate the mala mean in reverse
        for(int r=0; r<nvar; r++)
        {
            mala2[r] =  propphi[r] + 0.5 * pow(phi_tune, 2) * (phioffset(j,r) - propphi[r] - sum(fcprec(r,_) * (propphi - fcmean)));
        }  
        

        // Compute the prior ratio
        diffcurrent = phi(j,_) - fcmean;
        diffprop = propphi - fcmean;
        for(int r=0; r<nvar; r++)
        {
            quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
            quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
        oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
        newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
        // Likelihood ratio
        lpold = pow((phioffset(j,_) - phi(j,_)),2);
        lpnew = pow((phioffset(j,_) - propphi),2);
        oldlikebit = 0.5 * sum(lpold / nu2);
        newlikebit = 0.5 * sum(lpnew / nu2);                             
        
        // Hastings ratio   
        hastings = - (sum(pow(phi(j,_) - mala2,2)) - sum(pow(propphi - mala1,2))) / (2*pow(phi_tune,2));
        
        // Accept or reject the value
        acceptance = exp(oldpriorbit - newpriorbit + oldlikebit - newlikebit + hastings);
        if(runif(1)[0] <= acceptance) 
        {
            phi(j,_) = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}







// [[Rcpp::export]]
List multinomialbetaupdateRW(NumericMatrix X, const int nsites, const int J, const int p, 
                             const int col, NumericMatrix beta, NumericMatrix offset, NumericMatrix y,                            
                             NumericVector prior_meanbeta, NumericVector prior_varbeta, 
                             const int nblock, double beta_tune, List block_list, NumericVector zeros)
{
    // Compute the acceptance probability for beta
    //Create new objects
    NumericMatrix lp_current(nsites, J), lp_proposal(nsites, J);
    NumericVector p_current(nsites), p_proposal(nsites);
    int accept=0;
    double oldlikebit=0, newlikebit=0, priorbit=0;
    double acceptance;
    
    
    // Create a beta old and new matrix
    NumericMatrix beta_new(p, (J-1));
    for(int j = 0; j < (J-1); j++)
    {
        beta_new(_,j) = beta(_,j);
    }
    
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        for(int g=0; g<len; g++)
        {
            beta_new(idx[g],(col-1)) = rnorm(1, beta(idx[g],(col-1)), beta_tune)[0];
        }
        
        // Compute the linear predictors
        lp_current(_, 0) = zeros;
        lp_proposal(_, 0) = zeros;   
        
        for(int g=1; g<J; g++)
        {
            lp_current(_,g) = linpredcompute(X, nsites, p, beta(_, (g-1)), offset(_, (g-1)));    
            lp_proposal(_,g) = linpredcompute(X, nsites, p, beta_new(_, (g-1)), offset(_, (g-1)));    
        }
        
        // Compute the probabilities and the likelihood component of the MH step
        oldlikebit = 0;
        newlikebit=0;
        for(int j = 0; j < nsites; j++)     
        {
            p_current = exp(lp_current(j, _)) / sum(exp(lp_current(j, _)));
            p_proposal = exp(lp_proposal(j, _)) / sum(exp(lp_proposal(j, _))); 
            oldlikebit = oldlikebit + sum(y(j, _) * log(p_current));
            newlikebit = newlikebit + sum(y(j, _) * log(p_proposal));
        }
        
        // Compute the prior component of the MH step  
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta(idx[g],(col-1))-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new(idx[g],(col-1))-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Accept or reject the proposal      
        acceptance = exp(newlikebit - oldlikebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta(idx[g], (col-1)) = beta_new(idx[g], (col-1));  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new(idx[g], (col-1)) = beta(idx[g], (col-1));  
            }   
        }
    }
    
    
    // Compute the acceptance probability and return the value
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}




// [[Rcpp::export]]
List multinomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                          const int nsites,  const int nvar, NumericMatrix phi, 
                          NumericMatrix Y, NumericMatrix phioffset, 
                          NumericVector denoffset,  NumericMatrix Sigmainv, 
                          double rho, double phi_tune)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericMatrix fcprec((nvar-1), (nvar-1));
    int rowstart=0, rowend=0, accept=0;
    NumericVector sumphi((nvar-1)), fcmean((nvar-1)), propphi((nvar-1));
    NumericVector diffcurrent(nvar), diffprop(nvar);        
    NumericVector quadcurrent(nvar), quadprop(nvar);  
    NumericVector lpold(nvar), lpnew(nvar), pold(nvar), pnew(nvar);
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
    
    //  Update each random effect in turn
    for(int j = 0; j < nsites; j++)
    {
    // Calculate the prior precision and mean
        for(int r=0; r<(nvar-1); r++)
        {
        fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
        }
        
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0,(nvar-1));
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
        fcmean = rho * sumphi / denoffset[j]; 
        
        
    // Propose a possible value
        for(int r=0; r<(nvar-1); r++)
        {
        propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
        }
        
        
    // Compute the prior ratio
    diffcurrent = phi(j,_) - fcmean;
    diffprop = propphi - fcmean;
        for(int r=0; r<(nvar-1); r++)
        {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
        }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
        
    // Likelihood ratio
    // compute the linear predictor
    lpold[0] = 0;
    lpnew[0] = 0;   
        for(int g=1; g<nvar; g++)
        {
        lpold[g] =  phi(j, (g-1))  + phioffset(j, (g-1));
        lpnew[g] =  propphi[(g-1)]  + phioffset(j, (g-1));
        }
    
    // Compute the probabilities and the likelihood component of the MH step
    pold = exp(lpold) / sum(exp(lpold));
    pnew = exp(lpnew) / sum(exp(lpnew));; 
    oldlikebit = sum(Y(j, _) * log(pold));
    newlikebit = sum(Y(j, _) * log(pnew));

    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
        phi(j,_) = propphi;
        accept = accept + 1;
        }
        else
        { 
        }
    }     
    
    
    // Return the results
    List out(2);
    out[0] = phi;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
NumericVector gaussiancarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                          NumericVector Wtripletsum, NumericVector n_individual,
                                          const int nsites, NumericVector phi, double tau2, 
                                          double rho, double nu2, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int rowstart=0, rowend=0;
    double sumphi;
    double fcprecision, fcsd, fcmean;
    double priorvardenom, priormean, priorvar;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        fcprecision = (n_individual[j]/nu2) + (1/priorvar);
        fcsd = pow((1/fcprecision),0.5);
        fcmean = (priormean / priorvar + offset[j]) / fcprecision;
        phinew[j] = rnorm(1, fcmean, fcsd)[0];      
    }
    
    
    return phinew;
}




// [[Rcpp::export]]
List binomialcarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                 NumericVector Wtripletsum, List ind_area_list, NumericVector n_individual,
                                 const int nsites, NumericVector phi, double tau2, 
                                 const NumericVector y, const NumericVector failures, const double phi_tune, 
                                 double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0, rowstart=0, rowend=0, n_current=0, datapoint=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, likebittotal;
    double priorvardenom, priormean, priorvar;
    double propphi, pold, pnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
        
        
        // Accept or reject it
        // Prior part
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        // Likelihood part
        // Determine the set of data points relating to area j
        n_current = n_individual[j];
        NumericVector individuals(n_current);
        individuals = ind_area_list[j];
        
        // Compute the data likelihood
        likebittotal = 0;
        for(int r = 0; r < n_current; r++)
        {
            datapoint = individuals[r] - 1;
            pold = exp(offset[datapoint] + phinew[j]) / (1 + exp(offset[datapoint] + phinew[j]));
            pnew = exp(offset[datapoint] + propphi) / (1 + exp(offset[datapoint] + propphi)); 
            oldlikebit = y[datapoint] * log(pold) + failures[datapoint] * log((1-pold));
            newlikebit = y[datapoint] * log(pnew) + failures[datapoint] * log((1-pnew));
            likebittotal = likebittotal + newlikebit - oldlikebit;   
        }
        
        // Compute the acceptance probability and accept or reject the proposal    
        acceptance = exp(oldpriorbit - newpriorbit + likebittotal);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    // Return the results
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}





// [[Rcpp::export]]
List poissoncarmultilevelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                NumericVector Wtripletsum, List ind_area_list, NumericVector n_individual,
                                const int nsites, NumericVector phi, 
                                double tau2, const NumericVector y, const double phi_tune, 
                                double rho, NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0, n_current=0, datapoint=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit, likebittotal;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
        
        // Accept or reject it
        // Prior part
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        // Likelihood part
        // Determine the set of data points relating to area j
        n_current = n_individual[j];
        NumericVector individuals(n_current);
        individuals = ind_area_list[j];
        
        // Compute the data likelihood
        likebittotal = 0;
        for(int r = 0; r < n_current; r++)
        {
            datapoint = individuals[r] - 1;
            lpold = offset[datapoint] + phinew[j];
            lpnew = offset[datapoint] + propphi; 
            oldlikebit = y[datapoint] * lpold - exp(lpold);
            newlikebit = y[datapoint] * lpnew - exp(lpnew);
            likebittotal = likebittotal + newlikebit - oldlikebit;   
        }
        
        // Compute the acceptance probability and accept or reject the proposal    
        acceptance = exp(oldpriorbit - newpriorbit + likebittotal);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List poissoncarmultilevelupdateindiv(List ind_re_list, NumericVector n_re,
                                     const int q, NumericVector psi, 
                                     double sigma2, const NumericVector y, const double psi_tune, 
                                     NumericVector offset)
{
    // Update the independent individual (or near individual) random effects 
    //Create new objects
    NumericVector psinew(q);
    double proppsi, lpold, lpnew, acceptance, priorbit, likebittotal, oldlikebit, newlikebit;
    int n_current, datapoint, accept=0;    
    
    
    //  Update each random effect in turn
    psinew = psi;
    
    for(int j = 0; j < q; j++)
    {
        // propose a value  
        proppsi = rnorm(1, psinew[j], sqrt(psi_tune))[0];
        
        // Prior part
        priorbit = (0.5 / sigma2) * (pow(psinew[j], 2) - pow(proppsi, 2)); 
        
        // Likelihood part
        // Determine the set of data points relating to area j
        n_current = n_re[j];
        NumericVector individuals(n_current);
        individuals = ind_re_list[j];
        
        // Compute the data likelihood
        likebittotal = 0;
        for(int r = 0; r < n_current; r++)
        {
            datapoint = individuals[r] - 1;
            lpold = offset[datapoint] + psinew[j];
            lpnew = offset[datapoint] + proppsi; 
            oldlikebit = y[datapoint] * lpold - exp(lpold);
            newlikebit = y[datapoint] * lpnew - exp(lpnew);
            likebittotal = likebittotal + newlikebit - oldlikebit;   
        }
        
        // Compute the acceptance probability and accept or reject the proposal    
        acceptance = exp(priorbit + likebittotal);
        if(runif(1)[0] <= acceptance) 
        {
            psinew[j] = proppsi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = psinew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List binomialcarmultilevelupdateindiv(List ind_re_list, NumericVector n_re,
                                      const int q, NumericVector psi, 
                                      double sigma2, const NumericVector y, const NumericVector failures,
                                      const double psi_tune, NumericVector offset)
{
    // Update the independent individual (or near individual) random effects 
    //Create new objects
    NumericVector psinew(q);
    double proppsi, pold, pnew, acceptance, priorbit, likebittotal, oldlikebit, newlikebit;
    int n_current, datapoint, accept=0;    
    
    
    //  Update each random effect in turn
    psinew = psi;
    
    for(int j = 0; j < q; j++)
    {
        // propose a value  
        proppsi = rnorm(1, psinew[j], sqrt(psi_tune))[0];
        
        // Prior part
        priorbit = (0.5 / sigma2) * (pow(psinew[j], 2) - pow(proppsi, 2)); 
        
        // Likelihood part
        // Determine the set of data points relating to area j
        n_current = n_re[j];
        NumericVector individuals(n_current);
        individuals = ind_re_list[j];
        
        // Compute the data likelihood
        likebittotal = 0;
        for(int r = 0; r < n_current; r++)
        {
            datapoint = individuals[r] - 1;
            pold = exp(offset[datapoint] + psinew[j]) / (1 + exp(offset[datapoint] + psinew[j]));
            pnew = exp(offset[datapoint] + proppsi) / (1 + exp(offset[datapoint] + proppsi)); 
            oldlikebit = y[datapoint] * log(pold) + failures[datapoint] * log((1-pold));
            newlikebit = y[datapoint] * log(pnew) + failures[datapoint] * log((1-pnew));
            likebittotal = likebittotal + newlikebit - oldlikebit;   
        }
        
        // Compute the acceptance probability and accept or reject the proposal    
        acceptance = exp(priorbit + likebittotal);
        if(runif(1)[0] <= acceptance) 
        {
            psinew[j] = proppsi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = psinew;
    out[1] = accept;
    return out;
}



// [[Rcpp::export]]
NumericVector gaussiancarmultilevelupdateindiv(List ind_re_list, NumericVector n_re,
                                               const int q, NumericVector psi, double sigma2, double nu2, 
                                               NumericVector offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    NumericVector psinew(q);
    double fcprecision, fcsd, fcmean;
    
    //  Update each random effect in turn
    psinew = psi;
    for(int j = 0; j < q; j++)
    {
        // Calculate the precision, sd and mean
        fcprecision = (n_re[j] / nu2) + (1 / sigma2);
        fcsd = pow(( 1 / fcprecision), 0.5);    
        fcmean = offset[j] / fcprecision;
        
        // propose a value
        psinew[j] = rnorm(1, fcmean, fcsd)[0]; 
    }
    
    return psinew;
}



// This file contains the following functions:
// std::vector<std::vector<int>> optimise_graph(std::vector<std::vector<int>> adj, std::vector<double> data,
//    bool add=false, bool remove=true, bool remove_first=false) {
// Note that this function appears right at the end of this file.

// Needed for set_difference
#include <algorithm>
// shared_ptr
#include <memory>

// Create anonymous namespace to hold graph classes. Anything in an anonymous
// namespace is only available in this file.
namespace {


// Iterators to nicely generate powersets with minimal data copying
template<typename T>
class PowersetIterator {
private:
  std::vector<bool> counter_;
  const std::vector<T> & elts;
  bool end_;
public:
  typedef std::shared_ptr<const std::set<T>> value_type;
  typedef void difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef std::input_iterator_tag iterator_category;
  
  explicit PowersetIterator(const std::vector<T> & orig, bool end=false) : counter_(orig.size(), false), elts(orig), end_(end) { }
  std::shared_ptr<const std::set<T>> operator*() const {
    std::set<T> set;
    for(size_t i = 0; i < elts.size(); ++i) {
      if (counter_[i] == true) {
        set.insert(elts[i]);
      }
    }
    return std::make_shared<const std::set<T>>(set);
  }
  
  bool operator==(const PowersetIterator& other) const {
    // counter_ will be all false for both begin() and end(), so we use the
    // end_ member to note that we actually are at the end
    if ((end_ == true) && (other.end_ == true)) {
      return true;
    }
    if ((end_ == true) && (other.end_ == false)) {
      return false;
    }
    if ((end_ == false) && (other.end_ == true)) {
      return false;
    }
    return counter_ == other.counter_;
  }
  
  bool operator!=(const PowersetIterator& other) const { return !(*this == other); }
  
  PowersetIterator& operator++() {
    for(size_t i = 0; i < counter_.size(); ++i) {
      counter_[i] = !counter_[i];
      if (counter_[i] == true) {
        break;
      }
      // If we've flipped all bits, and not yet flipped anything
      // false->true, we must have just been at the point where everything
      // is true aka the last subset.
      if (i == (counter_.size() - 1)) {
        end_ = true;
      }
    }
    // If we are finding the powerset of an empty set, we find one empty
    // set and then reach the end.
    if (counter_.size() == 0) {
      end_ = true;
    }
    return *this;
  }
};

template<typename T>
class Powerset {
private:
  std::vector<T> elts;
public:
  Powerset(std::set<T> elts_) : elts(elts_.begin(), elts_.end()) { }
  PowersetIterator<T> begin() { return PowersetIterator<T>(elts); }
  PowersetIterator<T> end() { return PowersetIterator<T>(elts, true); }
};



// An edge in a graph
typedef std::tuple<int,int> Edge;

// A graph, which internally stores the adjacency matrix, and the current
// neighbourhood of each vertex.
class Graph {
public:
  Graph(std::vector<std::vector<int>> adj) {
    this->adj_ = adj;
    this->size_ = adj.size();
    this->nbs_.resize(this->size_);
    int row_id = 0;
    for(const auto row: this->adj_) {
      int col_id = 0;
      for(const auto entry: row) {
        if (entry == 1) {
          this->nbs_[row_id].insert(col_id);
        }
        col_id += 1;
      }
      row_id += 1;
    }
  }
  
  Graph(int size) : adj_(size, std::vector<int>(size)) {
    this->size_ = size;
    this->nbs_.resize(this->size_);
  }
  
  std::vector<std::vector<int>> adjMatrix() const {
    return std::vector<std::vector<int>>(this->adj_);
  }
  
  int size() const {
    return this->size_;
  }
  
  bool has_edge(int u, int v) {
    return this->adj_[u][v] == 1;
  }
  
  void add_edge(const Edge& edge) {
    int v0, v1;
    std::tie(v0, v1) = edge;
    this->adj_[v0][v1] = 1;
    this->adj_[v1][v0] = 1;
    this->nbs_[v0].insert(v1);
    this->nbs_[v1].insert(v0);
  }
  
  void add_edges(const std::list<Edge> & edges) {
    for(const auto edge: edges) {
      this->add_edge(edge);
    }
  }
  
  void remove_edge(const Edge& edge) {
    int v0, v1;
    std::tie(v0, v1) = edge;
    this->adj_[v0][v1] = 0;
    this->adj_[v1][v0] = 0;
    this->nbs_[v0].erase(v1);
    this->nbs_[v1].erase(v0);
  }
  void remove_edges(std::list<Edge> to_remove) {
    for(auto edge: to_remove) {
      this->remove_edge(edge);
    }
  }
  
  int degree(int v) const {
    return this->nbs_[v].size();
  }
  
  const std::set<int> nbs(int v) const {
    return this->nbs_[v];
  }
  
  const std::list<Edge> edges() const {
    std::list<Edge> result;
    int row_id = 0;
    for(const auto row: this->adj_) {
      int col_id = 0;
      for(const auto entry: row) {
        if (entry == 1) {
          result.push_back(std::make_tuple(col_id, row_id));
        }
        // Only check and return where col_id <= row_id
        if (col_id == row_id) {
          break;
        }
        col_id += 1;
      }
      row_id += 1;
    }
    return result;
  }
  
private:
  size_t size_;
  std::vector<std::vector<int>> adj_;
  std::vector<std::set<int>> nbs_;
};

// An "EmptyGraph" has a set of potential edges. This way we can restrict
// what edges we might add.
class EmptyGraph : public Graph {
  
public:
  EmptyGraph(std::vector<std::vector<int>> adj) : Graph(adj.size()), _possible(adj) {
  }
  
  const std::set<int> possNbs(int v) const {
    return this->_possible.nbs(v);
  }
  
  const std::list<Edge> origEdges() const {
    return this->_possible.edges();
  }
  
private:
  Graph _possible;
};

// A class to do our optimising. Using a class means we can store one central
// copy of the graph and the data, and still access them quickly
class Optimiser {
public:
  Optimiser(std::vector<std::vector<int>> adj, std::vector<double> data) : _graph(adj), _data(data) { }
  
  std::vector<std::vector<int>> adjMatrix() const { return this->_graph.adjMatrix(); }
  
  void iterative_opt(bool remove, bool add, bool remove_first) {
    if ((!add) || (remove && remove_first)) {
      this->_graph.add_edges(this->_graph.origEdges());
    } else {
      // Building from a mostly empty graph
      for(int v = 0; v < this->_graph.size(); ++v) {
        if (this->_graph.degree(v) >= 1) {
          continue;
        }
        signed int bestNb = -1;
        double minDiff = 0;
        for (auto u: this->_graph.possNbs(v)) {
          double diff = fabs(this->_data[v] - this->_data[u]);
          if ((bestNb == -1) || (diff < minDiff)) {
            bestNb = u;
            minDiff = diff;
          }
        }
        //Rprintf("Setup: Adding edge (%u, %u)\n", v, bestNb);
        this->_graph.add_edge(std::make_pair(v, bestNb));
      }
    }
    bool changed = true;
    while (changed) {
      changed = false;
      if (add && (! (remove_first && remove))) {
        double avsum = this->avsum();
        double oldscore = this->score();
        auto lastAdded = this->greedy_opt_add(avsum);
        double newscore = this->score();
        //Rprintf("adding\toldscore = %f\tnewscore = %f\tchange size is %u\n", oldscore, newscore, lastAdded.size());
        if (!lastAdded.empty()) {
          if (newscore > oldscore) {
            // Good change
            changed = true;
          } else {
            // Bad change, put the edges back
            this->_graph.remove_edges(lastAdded);
          }
        }
      }
      if (remove) {
        // Flip flag off so we add in next iteration.
        remove_first = false;
        double avsum = this->avsum();
        double oldscore = this->score();
        auto lastRemoved = this->greedy_opt_remove(avsum);
        double newscore = this->score();
        //Rprintf("removing\toldscore = %f\tnewscore = %f\tchange size is %u\n", oldscore, newscore, lastRemoved.size());
        if (!lastRemoved.empty()) {
          if (newscore > oldscore) {
            // Good change
            changed = true;
          } else {
            // Bad change, put the edges back
            this->_graph.add_edges(lastRemoved);
          }
        }
      }
    }
  }
private:
  EmptyGraph _graph;
  std::vector<double> _data;
  
  double disc(int v) const {
    double sum = 0;
    for(auto nb: this->_graph.nbs(v)) {
      sum += this->_data[nb];
    }
    return this->_graph.degree(v)*std::pow(this->_data[v] - sum/this->_graph.degree(v), 2);
  }
  
  double score() const {
    double first = 0;
    double second_inner = 0;
    for(int v = 0; v < this->_graph.size(); ++v) {
      first += std::log(this->_graph.degree(v));
      second_inner += this->disc(v);
    }
    return first/2.0f - (this->_graph.size()/2.0f)*std::log(second_inner);
  }
  
  double avsum() const {
    double sum = 0;
    for(int v = 0; v < this->_graph.size(); ++v) {
      sum += this->disc(v);
    }
    return sum;
  }
  
  double vertex_val(int v, const std::set<int> & fixed_nbs, const std::set<int> & opt_nbs, double avsum) const {
    double sum = 0;
    for(auto nb: fixed_nbs) {
      sum += this->_data[nb];
    }
    for(auto nb: opt_nbs) {
      sum += this->_data[nb];
    }
    int degree = fixed_nbs.size() + opt_nbs.size();
    double vx_av = degree * std::pow(this->_data[v] - sum/degree, 2);
    if (vx_av > avsum) {
      // The discrepancy at this vertex is already worse than the total
      // discrepancy in the old graph. This means this value should be bad,
      // but also breaks the calculations so we just return something very
      // negative.
      return -1e20;
    }
    return std::log(degree)/2.0f - (this->_graph.size()/2.0f)*std::log1p(vx_av / (avsum - vx_av));
  }
  
  // BEWARE: I define a < b to be true if a[1] > b[1] to get a reversed
  // multiset.
  class CondSortedEntry : public std::tuple<std::shared_ptr<const std::set<int>>, double> {
  public:
    using BaseType = std::tuple<std::shared_ptr<const std::set<int>>, double>;
    CondSortedEntry(const BaseType& other) : BaseType(other) { }
    bool operator<(const CondSortedEntry other) const {
      return std::get<1>(*this) > std::get<1>(other);
    }
  };
  
  typedef std::multiset<CondSortedEntry> CondSortedList;
  
  CondSortedList cond_sorted_nbs_remove(int v, const std::set<int> fixed_nbs, const std::set<int> options, double avsum) const {
    CondSortedList result;
    Powerset<int> pow(options);
    for(const auto option: pow) {
      // Must have degree >= 1
      if (fixed_nbs.size() + option->size() == 0) {
        continue;
      }
      result.insert(std::make_tuple(option, this->vertex_val(v, fixed_nbs, *option, avsum)));
    }
    return result;
  }
  
  CondSortedList cond_sorted_nbs_add(int v, const std::set<int> fixed, const std::set<int> options, double avsum) const {
    CondSortedList result;
    Powerset<int> pow(options);
    for(const auto option: pow) {
      // Must have degree >= 1
      if (fixed.size() + option->size() == 0) {
        continue;
      }
      result.insert(std::make_tuple(option, this->vertex_val(v, fixed, *option, avsum)));
    }
    return result;
  }
  
  double find_best_cond(const CondSortedList & list, signed int mustHave, signed int cantHave) const {
    for(const CondSortedEntry entry: list) {
      const std::set<int> & opts = * std::get<0>(entry);
      if ((mustHave != -1) && (opts.find((int)mustHave) == opts.end())) {
        continue;
      }
      if ((cantHave != -1) && (opts.find((int)cantHave) != opts.end())) {
        continue;
      }
      return std::get<1>(entry);
    }
    // TODO Error if we reach here?
    Rprintf("ERROR ---------------- ERROR\n");
    return 1;
  }
  
  std::list<Edge> greedy_opt_add(double avsum) {
    std::list<Edge> added;
    for(int v = 0; v < this->_graph.size(); v++) {
      std::set<int> nb_options;
      std::set<int> nb_fixed;
      for(int u: this->_graph.possNbs(v)) {
        // If we already have this edge, it is fixed for now.
        if (this->_graph.has_edge(u,v)) {
          nb_fixed.insert(u);
          continue;
        }
        // Only consider edges {v,u} where v < u, so this edge has not yet
        // been added, but has been considered.
        if (u < v) {
          continue;
        }
        nb_options.insert(u);
      }
      CondSortedList v_list = this->cond_sorted_nbs_add(v, nb_fixed, nb_options, avsum);
      std::list<Edge> add_here;
      for (auto u: nb_options) {
        // -1 indicates that we don't care.
        double vx_gain = this->find_best_cond(v_list, u, -1) - this->find_best_cond(v_list, -1, u);
        std::set<int> u_options;
        std::set<int> u_fixed;
        for(int u_nb: this->_graph.possNbs(u)) {
          if (this->_graph.has_edge(u, u_nb)) {
            u_fixed.insert(u_nb);
            continue;
          }
          if (u_nb < v) {
            continue;
          }
          u_options.insert(u_nb);
        }
        CondSortedList nbr_list = this->cond_sorted_nbs_add(u, u_fixed, u_options, avsum);
        double nb_gain = this->find_best_cond(nbr_list, v, -1) - this->find_best_cond(nbr_list, -1, v);
        //Rprintf("Adding (%u, %u), vx_gain is %f and nb_gain is %f\n", v, u, vx_gain, nb_gain);
        if (vx_gain + nb_gain > 0 ) {
          add_here.push_back(std::make_pair(v, u));
        }
      }
      this->_graph.add_edges(add_here);
      for(auto edge: add_here) {
        added.push_back(edge);
      }
    }
    return added;
  }
  
  std::list<Edge> greedy_opt_remove(double avsum) {
    std::list<Edge> removed;
    for(int v = 0; v < this->_graph.size(); v++) {
      if (this->_graph.degree(v) == 1) {
        continue;
      }
      std::set<int> nb_options;
      std::set<int> nb_fixed;
      for(int u: this->_graph.nbs(v)) {
        // Only consider edges {v,u} where v < u
        if ((u < v) || (this->_graph.degree(u) == 1)) {
          nb_fixed.insert(u);
        } else {
          nb_options.insert(u);
        }
      }
      CondSortedList v_list = this->cond_sorted_nbs_remove(v, nb_fixed, nb_options, avsum);
      std::list<Edge> remove_here;
      for (auto u: nb_options) {
        // Don't remove the last edge, just break
        if ((this->_graph.degree(v) - remove_here.size()) == 1) {
          break;
        }
        // -1 indicates that we don't care.
        double vx_gain = this->find_best_cond(v_list, u, -1) - this->find_best_cond(v_list, -1, u);
        std::set<int> u_fixed;
        std::set<int> u_options;
        for(int u_nb: this->_graph.nbs(u)) {
          if ((u_nb < v) || (this->_graph.degree(u_nb) == 1)) {
            u_fixed.insert(u_nb);
          } else {
            u_options.insert(u_nb);
          }
        }
        CondSortedList nbr_list = this->cond_sorted_nbs_remove(u, u_fixed, u_options, avsum);
        double nb_gain = this->find_best_cond(nbr_list, v, -1) - this->find_best_cond(nbr_list, -1, v);
        //Rprintf("Dropping (%u, %u), vx_gain is %f and nb_gain is %f\n", v, u, vx_gain, nb_gain);
        // the gains are from the edge being _in_ the graph, so we remove
        // if the gain is negative!
        if (vx_gain + nb_gain < 0 ) {
          remove_here.push_back(std::make_pair(v, u));
        }
      }
      this->_graph.remove_edges(remove_here);
      for(auto edge: remove_here) {
        removed.push_back(edge);
      }
    }
    return removed;
  }
};
} // namespace

// [[Rcpp::export]]
std::vector<std::vector<int>> optimise_graph(const IntegerMatrix& adj, const NumericVector& data,
                                             bool add=false, bool remove=true, bool remove_first=false) {
  //std::vector<double> dd = as<std::vector<double>>(data);
  std::vector<std::vector<int>> adj_;
  for(int r = 0; r < adj.nrow(); ++r) {
    std::vector<int> row;
    for(int c = 0; c < adj.ncol(); ++c) {
      row.push_back(adj(r,c));
    }
    adj_.push_back(row);
  }
  //Optimiser opt(as<std::vector<std::vector<int>>>(adj), dd);
  //Rcout << "adj is " << adj << "\n";
  //Rcout << "data is " << data << "\n";
  Optimiser opt(adj_, as<std::vector<double>>(data));
  opt.iterative_opt(remove, add, remove_first);
  return opt.adjMatrix();
}





