###
######
#########
############
#Here begins the MCMC functions
############
#########
######
###

##################################################################################
### create_status0 ### 
##################################################################################

 # function to get randomly chosen start parameters

create_status0 = function(x, beta0, sero_studies, model_type){
  
  no_sero_surveys = length(sero_studies) 
  
  ## itialising parameters for MCMC
  parnames=  c("vac_eff", paste("log",names(beta0),sep="."), paste(model_type, sero_studies, sep="_"), paste("vc_factor_CMRs", sep="_") )  
  pars_ini = rep(NA,length(parnames))
  names(pars_ini) = parnames
  
  ## filling the pars_in vector
  # vaccine efficacy KATY IS LOG TRANSFORMING THIS
  pars_ini[1] = log( rtrunc(1,"norm", a=0, b=1, mean = 0.975, sd = 0.05) ) #randomly sample from prior for start
  
  # GLM parameters
  ii = 2:(length(beta0)+1)
  pars_ini[ii] = beta0 
  
  if(model_type == "Foi"){
    # foi for each survey KATY IS LOG TRANSFORMING THESE
    pars_ini[(max(ii)+1):(max(ii)+6)] = log(rexp(6,  rate=100)) #random sample but not from prior
    pars_ini[grep("zone", parnames)] = log(rexp(no_sero_surveys-6,  rate=100)) 
  } else if(model_type == "R0"){
    # R0 for each survey KATY IS LOG TRANSFORMING THESE
    pars_ini[(max(ii)+1):(max(ii)+6)] = log(1+rexp(6,  rate=100))
    pars_ini[grep("zone", parnames)] = log(1+rexp(no_sero_surveys-6,  rate=100))
  }
  
  #vc.factor.CMRs KATY IS LOG TRANSFORMING THESE
  pars_ini[length(pars_ini)]= log( runif(1, min=0, max=1))
  
  ## declaring status list
  status0 = list(params=pars_ini,lnL=rep(-10000, no_sero_surveys+1), prior=-10000, accept=0) 
  names(status0$lnL) = c("lnL_GLM" ,paste("lnL", sero_studies, sep="_"))
  
  # indices where vaccination affects serology 
  varsin_nc = ii[-grep("adm0",colnames(x))] - 1 
  
  ## declare vector to identify different parameter types: vacc eff=1, GLM=2, Foi=3, vc.factor.CMRs =4
  parameter_type = c(1,rep(2,length(ii)),rep(3,no_sero_surveys), 4) # THIS NOW INDEXES THE DIFFERENT PARAMETER TYPES
  
  return(list(status0=status0, pars_ini=pars_ini, parameter_type=parameter_type, varsin_nc=varsin_nc, ii = ii))
}

##################################################################################
### fun_glm_lnL ### FROM fun.glm.lnL ###
##################################################################################

# Calculating lnL of glm regression:
fun_glm_lnL = function(beta, x, y) {
  # beta are the model coefficients,
  # x the independent covariates (needs one column for the Intercept),
  # and y the binary outcome_
  # model predictions pi = p, 1-pi = q
  eta = as.numeric(x %*% beta)
  logq = -exp(eta) # -exp(X*beta) (= log(1-q) )
  logp = log(1-exp(logq)) # 
  #q = exp(logq)
  #p = 1-q
  logl = sum(logp[y==1]) + sum(logq[y==0])
  return (logl)
}


##################################################################################
### fun_proposal ### 
##################################################################################
fun_proposal = function( status_new, status, chain_cov, adapt) {
  #this adapts the proposal distribution covariance when adapt = 1
  no_param=length(status$params)
  
  if (adapt ) {

    status_new_params= rmvnorm(n=1, mean=status$params, sigma = (2.38^2)*chain_cov/no_param) #'optimal' scaling of chain covariance
    
  } else {
    
    sigma= (1e-2)^2*diag(no_param)/no_param #this is an inital proposal covariance, see [Mckinley et al 2014]
    status_new_params= rmvnorm(n=1, mean=status$params, sigma = sigma)
    
  }
  status_new$params=as.numeric(status_new_params)
  names(status_new$params)=names(status$params)  
  return(status_new)
}



##################################################################################
### fun_herd_mcmc_step ### ADAPTED FROM fun.herd.mcmc.step ###
##################################################################################
fun_herd_mcmc_step = function(ii, 
                              x, 
                              y, 
                              status, 
                              chain_cov, 
                              adapt, 
                              survey_dat, 
                              no_sero_surveys, 
                              age_min, 
                              foi_const_surv, 
                              vc_factor, 
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
                              model_type,
                              parameter_type) {

  status_new = status
  names(status_new$params)=names(status$params)
  
  ############
  # Proposals calculation
  ############
  status_new=fun_proposal(status_new, status, chain_cov, adapt) 
  
  ############
  # prior calculation
  ############
  status_new$prior = fun_prior(status_new,  parameter_type, model_type) #this includes correction term for log transformed parameters
  
  
  ############
  # Likelihood calc
  ############
  if(is.finite(status_new$prior)){ # only calculates likelihood if parameters are in the right ranges
    status_new = fun_likelihood(status_new,
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
                                inc_v3d_agg)
    
    ############
    # accept/reject
    ############
    p_accept= sum(c(status_new$lnL, status_new$prior), na.rm=F) -  sum( c(status$lnL, status$prior), na.rm=F) 
    
  } else {
    ############
    # accept/reject
    ############
    p_accept= log(0)
    
  }
  
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1) ) { # accept:
    status_new$accept=1
  } else { # reject:
    status_new = status
    status_new$accept = 0
  }
  
  return(status_new)
}



##################################################################################
### run_mcmc_model ### ADAPTED FROM run.mcmc.r0.model ###
##################################################################################
run_mcmc_model = function(ii,
                          x, 
                          y,
                          adm1s,
                          chain_length, 
                          status,
                          varsin_nc,
                          dat, 
                          t0_vac_africa, 
                          sero_studies, 
                          study_years,
                          vc_agg3d,
                          survey_dat, 
                          age_min, 
                          foi_const_surv, 
                          vc_factor,
                          dim_year, 
                          dim_age, 
                          pop_agg3d,
                          p_prop_3d,
                          P_tot_2d, 
                          inc_v3d_agg, 
                          t0_vac,
                          p_at_survey,
                          P_tot_survey,
                          pop_moments_whole,
                          pop_moments_agg,
                          instance,
                          burnin,
                          learning_phase,
                          model_type,
                          parameter_type,
                          chain_cov = 1){

  #create a directory to save the output in
  date=format(Sys.time(),"%Y%m%d")
  name_dir = paste(model_type,"_MCMC_chain", "_", date, sep='')
  if(!dir.exists(name_dir)) dir.create(name_dir,  showWarnings = TRUE)
  setwd(name_dir)
  
  no_sero_surveys=length(sero_studies)
  ## DECLARE VARIABLES
  # mcmc_params
  mcmc_params = NULL
  lnLserosurvey=NULL
  for (i in 1:no_sero_surveys) {lnLserosurvey[i] = paste("lnL_serosurvey", i, sep="_")}
  # colnames(mcmc_params) = c(names(status$params),"lnL_GLM", lnLserosurvey, "Prior")
  
  # accept
  accept = NULL
  

  t= as.numeric(Sys.time())
  seed= (t - floor(t)) * 1e8 
  set.seed(seed)
  

  for(chainIndex in 1:chain_length) {
    
    if ( chainIndex < learning_phase | runif(1)<0.1 ) { 
      adapt=0
      
    } else { #adapt after n iterations for 90% of the time
      adapt=1
      if ( chain_cov == 1  | !(chainIndex %% 100)  ){ #PUT THIS BACK IN IF STARTING FROM COMPLETELY RANDOM PARAM
      #update covariance only every 100 timesteps after 100,000 iterations to reduce computation time
      chain_cov = cov( mcmc_params[max(c(burnin,chainIndex-1e5)):(chainIndex-1), 1:length(status$params)] ) 
      #calculate chain covariance up to current iteration for inferred parameters
      }
    }
    
    status = fun_herd_mcmc_step(ii, 
                                x, 
                                y, 
                                status, 
                                chain_cov, 
                                adapt, 
                                survey_dat, 
                                no_sero_surveys, 
                                age_min, 
                                foi_const_surv, 
                                vc_factor, 
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
                                model_type,
                                parameter_type)
    
    
    accept = rbind(accept,status$accept)
    
    mcmc_params = rbind(mcmc_params, c(status$params,status$lnL,status$prior))
    
    
    #if( chainIndex>50 & instance<100) plot(rowSums(mcmc_params[50:chainIndex,(length(status$params)+1):(length(status$params)+length(status$lnL)+1)]),type="l")
    if( chainIndex>50 & instance<100) plot(mcmc_params[,1])
    
    increment = 1000
    if (!(chainIndex %% increment)) { #every n iterations, output is saved, overwrites previous 

      colnames(mcmc_params) = c(names(status$params),"lnL_GLM", lnLserosurvey, "Prior")
      
      write.csv(mcmc_params[(chainIndex-increment):chainIndex,],paste("mcmc_",instance,"_",chainIndex/increment,".csv",sep=""),row.names=FALSE)
      accept_rate=(cumsum(accept)/c(1:chainIndex) )
      write.csv(accept_rate, file=paste("accept_",instance,  ".csv", sep="" ))
      
      
    }
    
  } #end chainIndex
  
  return(mcmc_params = mcmc_params)
}


##################################################################################
### fun_prior ### 
##################################################################################
# Calculating the prior probabilities
fun_prior = function(status,  parameter_type, model_type) {
  
  params=status$params
  names(params) =names(status$params)
  #recall: vacc_eff,R0 and vc.factor are log transformed, parameter_types=1,3,4
  
  correction = sum(params[parameter_type==1 | parameter_type==3 | parameter_type==4]) #correction to Jacobian
  #ADJUST FOR LOG TRANSFORM
  params[parameter_type==1 | parameter_type==3 | parameter_type==4]=exp(params[parameter_type==1 | parameter_type==3 | parameter_type==4])
  
  Prior = rep(NA, 5) # there are 4 different parameter types, vacc_eff, GLM, R0/FOI and vc.factor.CMRs (Kevin split the GLM into two cases- consider revising)
  
  #VACC EFF
  Prior[1]= log( dtrunc(params[ parameter_type == 1],"norm",a=0, b=1, mean = 0.975, sd = 0.05) )  #TRUNCATED NORMAL, PULLING THE EFFICACY UP TOWARDS THOSE FROM KEVINS PAPER
  
  
  #GLM
  jj = grep("^log.adm05",names(params)) # select country parameters (parameter_type=2) the sd.prior=2 is from Kevin's original code create status
  sd.prior= 0.5
  
  Prior[2] = -0.5*sum( (params[jj]/sd.prior)^2 ) # adjustment for reduced variation between countries? 
  
  Prior[3] = sum( dnorm(params[parameter_type == 2 & grepl("^log.adm05",names(params))==F], mean=0, sd=30, log=TRUE))
  # this term is for normally distributed non-country parameters : normal distrib with high sd (flat distrib)
  
  # R0/FOI
  if (model_type=="Foi"){
    Prior[4] = sum( dexp( params[parameter_type == 3] , rate=0.001, log=TRUE))  
  } else if (model_type=="R0"){
    Prior[4] = sum( dexp( params[parameter_type == 3] - 1, rate=0.001, log=TRUE)) 
  }
  
  #vc.factor
  Prior[5] = dunif( params[parameter_type == 4], min = 0, max = 1, log = TRUE) #  vc.factor     

  # print(Prior)
  
  return(sum(Prior)+correction)
}


