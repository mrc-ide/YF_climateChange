
calc_Foi_tempchange = function(dat_full,
                               adjusted_params,
                               temp_param,
                               scenarios,
                               c34,
                               adm1s,
                               t0_vac_africa,
                               sero_studies,
                               dim_year,
                               dim_age,
                               p_prop_3d,
                               P_tot_2d,
                               inc_v3d,
                               pop_moments_whole,
                               varsin_nc,
                               pop1,
                               vc2d) {
  
  FOI_median_samples = NULL
  
  
  ######BASELINE#######

  ### TEMP SUITABILITY ###
  dat_full_temp = cbind(dat_full, temp_suitability(dat_full[,"ERAday.mean"] , temp_param ))
  names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
  
  envdat = launch_env_dat(dat_full_temp,c34)  
  

  
  ### GET x ###
  modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability" 
  object_glm = fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
  x = object_glm[[2]]
  

  
  ### p detect ###
  p_detect =  fun_calc_pdetect_multi_both(x=x, 
                                          adm1s= adm1s, 
                                          params=adjusted_params, 
                                          dat = envdat$dat, 
                                          t0_vac_africa, 
                                          sero_studies, 
                                          dim_year, 
                                          dim_age, 
                                          p_prop_3d,
                                          P_tot_2d, 
                                          inc_v3d, 
                                          pop_moments_whole, 
                                          varsin_nc,
                                          "Foi")
  p_detect_link = mean(p_detect)
  
  ########SCENARIOS########
  for(s in scenarios){
    
    dat_full_ts = dat_full
    dat_full_ts[,"ERAday.mean"] = dat_full_ts[,"ERAday.mean"] + s 

    ### TEMP SUITABILITY ###
    dat_full_temp_ts = cbind(dat_full_ts, temp_suitability(dat_full_ts[,"ERAday.mean"] , temp_param ) )
    names(dat_full_temp_ts)[ncol(dat_full_temp_ts)] = "temp_suitability"
    envdat_ts = launch_env_dat(dat_full_temp_ts,c34) 
    
    ## sort variance
    envdat_ts$dat$temp_suitability = envdat_ts$dat$temp_suitability * (max(envdat$dat$temp_suitability)/max(envdat_ts$dat$temp_suitability))
    
    
    ### GET x ###
    object_glm_ts = fit_glm(dat =envdat_ts$dat, depi = envdat_ts$depi, modelVec ) 
    x_ts = object_glm_ts[[2]]
    
    mypreds_ts  = fun_calcPred(coefs = as.numeric(adjusted_params)[ii], 
                               newdata=x_ts, 
                               type="link",
                               varsin=varsin_nc)
    
    
    
    Ninf_whole = exp( mypreds_ts - p_detect_link)
    
    # calculate transmission
    pop_vc_moments = create_pop30.agg_vc30.agg(pop1,vc2d)$pop.vc.moments
    
    z = -Ninf_whole
    polydeg = 5
    
    if(polydeg>0) for(i in 1:polydeg) {
      z = cbind(z,(-1)^(i+1)*pop_vc_moments[,i+1]/factorial(i-1))
    }
    
    transmission_whole = sapply(1:nrow(x_ts), function(i) polyroot(z[i,]))
    transmission_whole[abs(Arg(transmission_whole))<=1e-10] = Re(transmission_whole)[abs(Arg(transmission_whole))<=1e-10]
    transmission_whole[abs(Arg(transmission_whole))>1e-10] = NA
    
    dt = dim(transmission_whole)
    transmission_whole = as.numeric(transmission_whole)
    dim(transmission_whole) = dt
    transmission_whole = apply(transmission_whole,2,min,na.rm=T)

    
    Foi_tmp = cbind(data.frame("scenario" = s), t(transmission_whole) )
    names(Foi_tmp)[2:480] = dim_adm
    
    FOI_median_samples = rbind(FOI_median_samples, Foi_tmp)
  }
  
  
  return(FOI_median_samples)
}