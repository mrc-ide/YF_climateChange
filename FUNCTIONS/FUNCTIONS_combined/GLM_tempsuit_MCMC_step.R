
##################################################################################
### GLM tempsuit mcmc step ###
##################################################################################
GLM_tempsuit_MCMC_step = function(param,
                                  dat_full,
                                  dat_bite,
                                  dat_mort,
                                  dat_EIP,
                                  chain_cov,
                                  adapt, 
                                  accCurrent) {
  
  ### propose new param ###
  param_prop = YFestimation::GLMproposal(param, chain_cov, adapt)
  
  ### priors ###
  prior_prop = sum( YFestimation::GLMprior(param_prop[1:19]) ) + 
                    fun_tempsuitPrior( param_prop[20:28])
  
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    
    ### TEMP SUITABILITY ###
    dat_full_temp = cbind(dat_full, 
                          temp_suitability(dat_full[,"ERAday.mean"] , 
                                           param_prop[20:28]))
    names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
    dat_full_temp$temp_suitability = dat_full_temp$temp_suitability    # this just dodges really zeros
    envdat = launch_env_dat(dat_full_temp,c34)  # see c34 is a global var!
    
    ### GET x ###
    modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability" 
    object_glm = YFestimation::fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
    x = object_glm[[2]]
    y = object_glm[[3]]
    
    ### LIKE ###
    like_prop = YFestimation::GLMlike(param_prop[1:19], x, y) + 
                fun_tempsuitLike(dat_bite, dat_mort, dat_EIP, param_prop[20:28]) 
    
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



