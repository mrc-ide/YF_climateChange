##################################################################################
### FUNCTIONS FOR TRANSDIMENSIONAL MCMC ###
##################################################################################

##################################################################################
### PROPOSAL ###
##################################################################################
fun_proposal = function(param, chain_cov, adapt, model_type) {
  #this adapts the proposal distribution covariance when adapt = 1
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
  
  names(param_prop) = names(param)
  
  if (runif(1)<0.5){
    model_type_prop = "Foi"
  } else {
    model_type_prop = "R0"
  }
  return(list(param_prop = param_prop, model_type_prop = model_type_prop) )
}


##################################################################################
### PSEUDO PRIORS ###
##################################################################################

pseudo_prior = function(param,
                        posterior_distributions) {
  
  if(length(grep("R0", names(param)))){ 
    lower_bound = 0
    param = param - 1
  } else { 
    lower_bound = 0
  }
  
  prob = rep(NA, length(param))
  for (i in 1:length(param)) {
    if (posterior_distributions[i, 3] == "beta") {
      prob[i] = dtrunc(param[i], "beta", a=lower_bound, b=Inf, shape1 = posterior_distributions[i, 1], shape2 = posterior_distributions[i, 2])
    } 
    else if (posterior_distributions[i, 3] == "gamma") {
      prob[i] = dtrunc(param[i], "gamma", a=lower_bound, b=Inf, shape = posterior_distributions[i, 1], rate = posterior_distributions[i, 2])
    }
    #else {
    #   prob[i] = dtrunc(param[i], "norm", a=lower_bound, b=Inf, mean = posterior_distributions[i, 1], sd = posterior_distributions[i, 2])
    # }
  }
  


  pp = sum(log(prob) )
  

  return(pp)
}


##################################################################################
### PRIOR ###
##################################################################################
fun_prior = function(params,  parameter_type, model_type, prob_R0, R0_posterior_distributions, FOI_posterior_distributions ) {
  #recall: vacc_eff,R0 and vc.factor are log transformed, parameter_types=1,3,4
  correction = sum(params[parameter_type == 1 |
                            parameter_type == 3 |
                            parameter_type == 4]) #correction to Jacobian
  

  #ADJUST FOR LOG TRANSFORM
  params[parameter_type == 1 |
           parameter_type == 3 |
           parameter_type == 4] = exp(params[parameter_type == 1 |
                                               parameter_type == 3 |
                                               parameter_type == 4])
  
  Prior = rep(0,8)
  
  #VACC EFF
  Prior[1]= log( dtrunc(params[ parameter_type == 1],"norm",a=0, b=1, mean = 0.975, sd = 0.05) )  #TRUNCATED NORMAL, PULLING THE EFFICACY UP TOWARDS THOSE FROM KEVINS PAPER
  
  
  #GLM
  jj = grep("^log.adm05", names(params)) # select country parameters (parameter_type=2) the sd.prior=2 is from Kevin's original code create status
  sd.prior = 0.5
  
  Prior[2] =  - 0.5 * sum((params[jj] / sd.prior) ^ 2) # adjustment for reduced variation between countries?
  
  Prior[3] =  sum(dnorm(
    params[parameter_type == 2 &
             grepl("^log.adm05", names(params)) == F],
    mean = 0,
    sd = 30,
    log = TRUE
  ))
  # this term is for normally distributed non-country parameters : normal distrib with high sd (flat distrib)
  
  # R0/FOI
  if (model_type == "Foi") {
    Prior[4] =  sum(dexp(params[grep("Foi", names(params))] , rate = 0.001, log =TRUE))   #normal prior for FOI

    Prior[5] =  pseudo_prior( params[grep("R0", names(params))], R0_posterior_distributions)
  } else if (model_type == "R0") {
    Prior[4] =  pseudo_prior( params[grep("Foi", names(params))], FOI_posterior_distributions)

    Prior[5] =  sum(dexp(params[grep("R0", names(params))] - 1, rate = 0.001, log =TRUE))
  }

  print(c(model_type,Prior[4:5]))
  #vc.factor
  Prior[6] =  dunif(params[parameter_type == 4],
                        min = 0,
                        max = 1,
                        log = TRUE) #  vc.factor
  
  #Model priors
  if(model_type == "R0"){
    Prior[7] = log(prob_R0) 
  } else if((model_type == "Foi")){
    Prior[7] =  log((1-prob_R0))
  }
  
  Prior[8] = correction
  
  out = as.numeric( Prior )
  return( out )
}

##################################################################################
### LIKELIHOOD ###
##################################################################################
fun_likelihood = function(param,
                          ii,
                          x,
                          y,
                          survey_dat,
                          no_sero_surveys,
                          foi_const_surv,
                          vc_factor,
                          age_min,
                          study_years,
                          vc_agg3d,
                          pop_agg3d,
                          pop_moments_agg,
                          t0_vac,
                          dim_year,
                          dim_age,
                          p_at_survey,
                          P_tot_survey,
                          inc_v3d_agg,
                          model_type) {
  #recall: vacc_eff,R0 and vc.factor are log transformed, parameter_types=1,3,4
  
  #sort out parameters
  
  vac_eff = exp(param[1]) # ADJUST FOR LOG TRANSFORM
  vcfac = exp(param[length(param)]) # ADJUST FOR LOG TRANSFORM
  
  ### declare lnL
  lnL = rep(NA, 41)
  
  ## GLM LIKELIHOOD
  # select GLM parameters # ii indexes the GLM parameters
  lnL[1] = fun_glm_lnL(param[ii], x, y) # lnL  GLM
  
  if (model_type == "Foi") {
    for (surveyIndex in 1:no_sero_surveys) {
      foi = exp(param[max(ii) + surveyIndex]) # ADJUST FOR LOG TRANSFORM
      
      vc_agg_ag = NULL
      ag = findInterval(0:100, age_min[[surveyIndex]])
      vc_agg_ag = aggregate(c(as.matrix(vc_agg3d[surveyIndex, paste0(study_years[surveyIndex]),]
                                        * pop_agg3d[surveyIndex, paste0(study_years[surveyIndex]),])),
                            by = list(ag = ag), sum)$x / 
                  aggregate(c(as.matrix(pop_agg3d[surveyIndex, paste0(study_years[surveyIndex]),])), 
                            by = list(ag = ag), sum)$x
      
      
      lnL[1+surveyIndex] =  fun_lnL_survey(
        foi = foi,
        obs_tot = pull(survey_dat[[surveyIndex]], samples) ,
        obs_pos = pull(survey_dat[[surveyIndex]], positives),
        pop = pop_agg3d[surveyIndex, paste0(study_years[surveyIndex]),],
        age_groups = age_min[[surveyIndex]],
        vc = vcfac * vac_eff * vc_agg_ag)
      
      
    }
  } else if (model_type == "R0") {
    # survey likelihood
    for (surveyIndex in 1:no_sero_surveys) {
      R0 = exp(param[max(ii) + no_sero_surveys + surveyIndex]) # ADJUST FOR LOG TRANSFORM
      
      lnL[1+surveyIndex] = fun_herd_lnL_survey(
        index_survey = surveyIndex,
        R0 = R0,
        foi_const = foi_const_surv[surveyIndex],
        t0_vac,
        dim_year,
        dim_age,
        study_years,
        p_at_survey,
        P_tot_survey,
        inc_v3d_agg,
        obs_tot = pull(survey_dat[[surveyIndex]], samples) ,
        obs_pos = pull(survey_dat[[surveyIndex]], positives),
        age_groups = age_min[[surveyIndex]],
        pop_moments_agg = pop_moments_agg,
        vcfac_arg = vcfac,
        vac_eff_arg = vac_eff
      )
      
    }
  }
  
  return(sum(lnL))
}


##################################################################################
### MCMC_STEP ###
##################################################################################

MCMC_step = function(param,
                     ii,
                     x,
                     y,
                     survey_dat,
                     no_sero_surveys,
                     foi_const_surv,
                     vc_factor,
                     age_min,
                     study_years,
                     vc_agg3d,
                     pop_agg3d,
                     pop_moments_agg,
                     t0_vac,
                     dim_year,
                     dim_age,
                     p_at_survey,
                     P_tot_survey,
                     inc_v3d_agg,
                     chain_cov,
                     adapt, 
                     model_type,
                     parameter_type,
                     accCurrent,
                     prob_R0, 
                     R0_posterior_distributions,
                     FOI_posterior_distributions) {
  
  ### propose new param ###
  out = fun_proposal(param, chain_cov, adapt, model_type)
  
  param_prop = out$param_prop
  names(param_prop) = names(param)
  model_type_prop = out$model_type_prop
  

  ### priors ###
  prior_prop = sum( fun_prior(param_prop,  parameter_type, model_type_prop, prob_R0, R0_posterior_distributions, FOI_posterior_distributions) )

  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    like_prop = fun_likelihood( param_prop,
                                ii,
                                x,
                                y,
                                survey_dat,
                                no_sero_surveys,
                                foi_const_surv,
                                vc_factor,
                                age_min,
                                study_years,
                                vc_agg3d,
                                pop_agg3d,
                                pop_moments_agg,
                                t0_vac,
                                dim_year,
                                dim_age,
                                p_at_survey,
                                P_tot_survey,
                                inc_v3d_agg,
                                model_type_prop)
    
    ### accept/ reject ###
    accProp = like_prop + prior_prop
    #print(c(model_type_prop, like_prop, prior_prop))
    p_accept= accProp - accCurrent
  } else {
    p_accept = log(0)
  }
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1) ) { # accept:
    param = param_prop
    accCurrent = accProp
    model_type = model_type_prop
    accept = 1
  } else { # reject:
    
    accept = 0
  }
  
  return(list(param = param, model_type = model_type, accCurrent = accCurrent, accept = accept)  )
}


##################################################################################
### Gibbs sampler ###
##################################################################################

Gibbs_MCMC_step = function(param,
                              ii,
                              x,
                              y,
                              survey_dat,
                              no_sero_surveys,
                              foi_const_surv,
                              vc_factor,
                              age_min,
                              study_years,
                              vc_agg3d,
                              pop_agg3d,
                              pop_moments_agg,
                              t0_vac,
                              dim_year,
                              dim_age,
                              p_at_survey,
                              P_tot_survey,
                              inc_v3d_agg,
                              chain_cov,
                              adapt, 
                              model_type,
                              parameter_type,
                              accCurrent,
                              prob_R0,
                              R0_posterior_distributions,
                              FOI_posterior_distributions) {
  
  ### propose new parameters ###
  out = fun_proposal(param, chain_cov, adapt, model_type)
  
  ### Accept/ reject `off' parameters`
  param_prop = param
  
  #only change "off" parameters
  if(model_type == "Foi") param_prop[grep("R0", names(param))] =  out$param_prop[grep("R0", names(param))] 
  if(model_type == "R0") param_prop[grep("Foi", names(param))] =  out$param_prop[grep("Foi", names(param))] 

  ### priors ###
  prior_prop = fun_prior(param_prop,  parameter_type, model_type, prob_R0, R0_posterior_distributions, FOI_posterior_distributions)
  prior = fun_prior(param,  parameter_type, model_type, prob_R0, R0_posterior_distributions, FOI_posterior_distributions)
  
  #accept /reject
  tmp = (runif(1))
  if(!is.finite( sum(prior_prop) ) & tmp< min( exp(sum(prior_prop[c(4,5,8)]) -sum(prior[c(4,5,8)]) ), 1) ) { # accept:
    param = param_prop
  } 


  
  ### Accept/ reject `on' parameters
  param_prop = param
  
  if(model_type == "Foi") param_prop[!grep("R0", names(param))] =  out$param_prop[!grep("R0", names(param))] 
  if(model_type == "R0") param_prop[!grep("Foi", names(param))] =  out$param_prop[!grep("Foi", names(param))] 
  
  ### priors ###
  prior_prop = sum( fun_prior(param_prop,  parameter_type, model_type, prob_R0, R0_posterior_distributions, FOI_posterior_distributions) )
  prior = sum( fun_prior(param,  parameter_type, model_type, prob_R0, R0_posterior_distributions, FOI_posterior_distributions) )

  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    like_prop = fun_likelihood( param_prop,
                                ii,
                                x,
                                y,
                                survey_dat,
                                no_sero_surveys,
                                foi_const_surv,
                                vc_factor,
                                age_min,
                                study_years,
                                vc_agg3d,
                                pop_agg3d,
                                pop_moments_agg,
                                t0_vac,
                                dim_year,
                                dim_age,
                                p_at_survey,
                                P_tot_survey,
                                inc_v3d_agg,
                                model_type)
    
    like = fun_likelihood( param,
                                ii,
                                x,
                                y,
                                survey_dat,
                                no_sero_surveys,
                                foi_const_surv,
                                vc_factor,
                                age_min,
                                study_years,
                                vc_agg3d,
                                pop_agg3d,
                                pop_moments_agg,
                                t0_vac,
                                dim_year,
                                dim_age,
                                p_at_survey,
                                P_tot_survey,
                                inc_v3d_agg,
                                model_type)
    
    ### accept/ reject ###
    accProp = like_prop + prior_prop
    accCurrent = like + prior
    p_accept= accProp - accCurrent
  } else {
    p_accept = log(0)
  }
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1) ) { # accept:
    param = param_prop
    accCurrent = accProp
  } 
  
  
  ## Propose new model ##
  model_type_prop = out$model_type_prop
  
  
  ### priors ###
  prior_prop = sum( fun_prior(param,  parameter_type, model_type_prop, prob_R0, R0_posterior_distributions, FOI_posterior_distributions) )
  
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    like_prop = fun_likelihood( param,
                                ii,
                                x,
                                y,
                                survey_dat,
                                no_sero_surveys,
                                foi_const_surv,
                                vc_factor,
                                age_min,
                                study_years,
                                vc_agg3d,
                                pop_agg3d,
                                pop_moments_agg,
                                t0_vac,
                                dim_year,
                                dim_age,
                                p_at_survey,
                                P_tot_survey,
                                inc_v3d_agg,
                                model_type_prop)
    
    ### accept/ reject ###
    accProp = like_prop + prior_prop
    #print(c(model_type, accProp))
    p_accept= accProp - accCurrent
  } else {
    p_accept = log(0)
  }
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1) ) { # accept:
    model_type = model_type_prop
    accCurrent = accProp
    accept = 1
  } else { # reject:
    
    accept = 0
  }
  
  return(list(param = param, model_type = model_type, accCurrent = accCurrent, accept = accept)  )
}


