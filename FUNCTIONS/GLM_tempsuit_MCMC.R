
##################################################################################
### GLM tempsuit mcmc step ###
##################################################################################
GLM_tempsuit_MCMC_step = function(param,
                                  dat_full,
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
  prior_prop = sum( YFestimation::GLMprior(param_prop[1:20]) ) + 
                    fun_tempsuitPrior( param_prop[21:29])
  
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    
    ### TEMP SUITABILITY ###
    dat_full_temp = cbind(dat_full, 
                          temp_suitability(dat_full[,"ERAday.mean"] , 
                                           param_prop[21:29]))
    names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"

    envdat = YFestimation::launch_env_dat(filepath = NA, 
                                          dat_full = dat_full_temp, 
                                          c34 = c34)  
    
    ### GET x ###
    modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability+RFE.mean" 
    object_glm = YFestimation::fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
    x = object_glm[[2]]
    y = object_glm[[3]]
    
    ### LIKE ###
    like_prop = YFestimation::GLMlike(param_prop[1:20], x, y) + 
                fun_tempsuitLike(dat_bite, dat_mort, dat_EIP, param_prop[21:29]) 
    
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
                             pars_ini,
                             dat_full,
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
  
  #get in the right order
  pars_ini = pars_ini[c("Intercept","log.surv.qual.adm0","adm05AGO",
                        "adm05BDI","adm05ERI","adm05ETH","adm05GNB",
                        "adm05KEN","adm05MRT","adm05RWA" ,
                        "adm05SDN" ,"adm05SOM" ,"adm05SSD" ,"adm05TZA",
                        "adm05UGA","adm05ZMB" , "lon" ,"logpop","temp_suitability",
                        "RFE.mean",
                        "a_T0" ,"a_Tm" ,"a_c", "mu_T0","mu_Tm",
                        "mu_c" ,"PDR_T0","PDR_Tm","PDR_c") ]
  
  ### find posterior probability at start ###
  out = GLM_tempsuit_MCMC_step(pars_ini,
                               dat_full,
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
                paste0(name_dir,"/","GLM_tempsuit_rain_chain",fileIndex,"_output_", run_id,".csv") ) 
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
    out = GLM_tempsuit_MCMC_step(param,
                                 dat_full,
                                 dat_bite,
                                 dat_mort,
                                 dat_EIP,
                                 c34,
                                 chain_cov,
                                 adapt, 
                                 accCurrent)
    
  }
  
  
}