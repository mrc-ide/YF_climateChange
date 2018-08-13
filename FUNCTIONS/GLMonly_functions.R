#GLM specific functions

##################################################################################
### simple proposal ###
##################################################################################
proposal = function(param, chain_cov, adapt) {
  
  no_param = length(param)
  
  if (adapt) {
    param_prop = rmvnorm(
      n = 1,
      mean = param,
      sigma = (2.38 ^ 2) * chain_cov / no_param
    ) #'optimal' scaling of chain covariance
    
  } else {
    sigma = (1e-2) ^ 2 * diag(no_param) / no_param #this is an inital proposal covariance, see [Mckinley et al 2014]
    param_prop = rmvnorm(n = 1,
                         mean = param,
                         sigma = sigma)
    
  }
  

  return(param_prop[1,])
}

##################################################################################
### GLM PRIOR ###
##################################################################################
fun_GLMprior = function(params) {
  Prior = rep(0,3)
  
  #GLM
  jj = grep("^log.adm05", names(params)) # select country parameters (parameter_type=2) the sd.prior=2 is from Kevin's original code create status
  sd.prior = 2
  
  Prior[2] =  - 0.5 * sum((params[jj] / sd.prior) ^ 2) # adjustment for reduced variation between countries?
  
  Prior[3] =  sum(dnorm(
    params[grepl("^log.adm05", names(params)) == F],
    mean = 0,
    sd = 30,
    log = TRUE
  ))
  # this term is for normally distributed non-country parameters : normal distrib with high sd (flat distrib)
  
  
  out = as.numeric( Prior )
  return( out )
}

##################################################################################
### GLM mcmc step ###
##################################################################################
GLM_MCMC_step = function(param,
                     x,
                     y,
                     chain_cov,
                     adapt, 
                     accCurrent) {
  
  ### propose new param ###
  param_prop = proposal(param, chain_cov, adapt)
 
  ### priors ###
  prior_prop = sum( fun_GLMprior(param_prop) )
  
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    like_prop = fun_glm_lnL(param_prop, x, y)
    
    ### accept/ reject ###
    accProp = like_prop + prior_prop

    p_accept= accProp - accCurrent
    if(is.na(p_accept) ){ p_accept = -Inf}
    
  } else {
    p_accept = log(0)
  }
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1)  ) { # accept:
    param = param_prop
    accCurrent = accProp
    accept = 1
  } else { # reject:
    
    accept = 0
  }
  
  return(list(param = param,  accCurrent = accCurrent, accept = accept)  )
}



