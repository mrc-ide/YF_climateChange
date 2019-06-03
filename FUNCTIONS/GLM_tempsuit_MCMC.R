
##################################################################################
### GLM tempsuit mcmc step ###
##################################################################################
GLM_tempsuit_MCMC_step = function(modelVec,
                                  param,
                                  dat_full,
                                  temp_type,
                                  dat_bite,
                                  dat_mort,
                                  dat_EIP,
                                  c34,
                                  chain_cov,
                                  adapt, 
                                  accCurrent) {
  
  ### propose new param ###
  param_prop = YFestimation::GLMproposal(param, chain_cov, adapt)
  
  ### priors ###
  prior_prop = sum( GLMprior_ts(param_prop[1:21]) ) + 
    fun_tempsuitPrior( param_prop[22:30])
  
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    
    ### TEMP SUITABILITY ###
    dat_full_temp = cbind(dat_full, 
                          temp_suitability(dat_full[,temp_type], 
                                           param_prop[22:30]))
    names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
    
    envdat = launch_env_dat(filepath = NA, 
                            dat_full = dat_full_temp, 
                            c34 = c34)  
    
    ### GET x ###
    
    object_glm = YFestimation::fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
    x = object_glm[[2]]
    y = object_glm[[3]]
    
    ### LIKE ###
    if(!anyNA(x) & is.finite(x)){
      like_prop = YFestimation::GLMlike(param_prop[1:21], x, y) + 
        fun_tempsuitLike(dat_bite, dat_mort, dat_EIP, param_prop[22:30]) 
    } else {
      like_prop = -Inf
    }
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




##################################################################################
### GLM tempsuit mcmc ###
##################################################################################
GLM_tempsuit_MCMC = function(Niter,
                             name_dir,
                             modelVec,
                             pars_ini,
                             dat_full,
                             temp_type,
                             dat_bite,
                             dat_mort,
                             dat_EIP,
                             c34,
                             run_id = NA,
                             plot_chain = TRUE){
  if(is.na(run_id)){
    run_id = 1
  }
  
  #set random seed
  set.seed(run_id)
  
  
  ### find posterior probability at start ###
  out = GLM_tempsuit_MCMC_step(modelVec,
                               pars_ini,
                               dat_full,
                               temp_type,
                               dat_bite,
                               dat_mort,
                               dat_EIP,
                               c34,
                               chain_cov=1,
                               adapt=0, 
                               accCurrent=-Inf)
  
  
  
  
  #mcmc setup
  chain = posteriorProb = acceptRate = NULL
  burnin = min(2*length(pars_ini), Niter)
  fileIndex = 0
  iter = 1
  
  for (iter in iter:Niter){
    #save current step
    param = out$param
    names(param) = names(pars_ini)
    accCurrent = out$accCurrent
    accept = out$accept
    
    chain = rbind(chain, param)
    posteriorProb = rbind(posteriorProb, accCurrent)
    acceptRate = rbind(acceptRate, accept)
    
    
    if(iter>100 & plot_chain ) plot(chain[,1] , type= "l") 
    
    if (iter %% 1000 == 0){
      
      colnames(acceptRate) = "acceptRate"
      colnames(posteriorProb) = "posteriorProb"
      if (iter %% 10000 == 0){
        fileIndex  = iter/10000
      }
      write.csv(cbind(chain,  posteriorProb, acceptRate)[min((fileIndex * 10000+1),iter):iter,], 
                paste0(name_dir,"/","GLM_tempsuit_rain_chain",run_id,"_output_", fileIndex,".csv") ) 
    }
    
    
    #adapt?
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }
    
    #new step
    out = GLM_tempsuit_MCMC_step(modelVec,
                                 param,
                                 dat_full,
                                 temp_type,
                                 dat_bite,
                                 dat_mort,
                                 dat_EIP,
                                 c34,
                                 chain_cov,
                                 adapt, 
                                 accCurrent)
    
  }
  
  
}


##################################################################################
### GLM prior ###
##################################################################################


GLMprior_ts = function(param) {
  
  Prior = rep(0,3)
  
  #GLM
  jj = grep("^adm05", names(param)) 
  sd.prior = 0.5 ##changed this
  
  Prior[1] =  - 0.5 * sum((param[jj] / sd.prior) ^ 2) # adjustment for reduced variation between countries?
  
  Prior[2] =  sum(dnorm(param[grepl("^adm05", names(param)) == FALSE & grepl("temp_suitability", names(param)) == FALSE],
                        mean = 0,
                        sd = 30,
                        log = TRUE))
  
  Prior[3] =  log( dtrunc(param[names(param) == "temp_suitability"], 
                          "norm", 
                          a = 0, 
                          b = Inf,  
                          mean = 0, 
                          sd = 30) ) 
  
  out = as.numeric( Prior )
  return( out )
}
