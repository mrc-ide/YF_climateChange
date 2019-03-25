#transdimensional mcmc for estimating glm AND temp suitability

run_estimation = function(run_id){
  #########################################################################################################
  ### LIBRARIES FOR PACKAGES USED ###
  #########################################################################################################
  #########################################################################################################
  ### LIBRARIES FOR PACKAGES USED ###
  #########################################################################################################
  
  library(dplyr)
  library(readr)
  library(reshape)
  library(mvtnorm)
  library(truncdist)
  library(R.utils)
  library(tibble)
  
  library(YFestimation)
  library(KsetupR)
  
  #########################################################################################################
  ### LOADING COUNTRIES ###
  #########################################################################################################
  
  Countries = read_csv(paste0("../Data/","Countries.csv"))
  c34 = Countries$c34
  
  #########################################################################################################
  ### SOURCE FUNCTIONS ###
  #########################################################################################################
  
  sourceDirectory("FUNCTIONS", modifiedOnly = FALSE)
  
  #########################################################################################################
  ### LOAD ENVIRONMENTAL DATA ###
  #########################################################################################################
  
  Env_Table_path = (paste0("../Data/","Environment/dat_with_worldclim/dat_worldclim_baseline_2019-03-22.csv")) 
  
  dat_full = read.csv(Env_Table_path, 
                      stringsAsFactors = FALSE)
  
  dat_full = dat_full %>% add_column(worldclim_temp_mid = (dat_full$worldclim_temp_min + dat_full$worldclim_temp_max)/2)
  
  dat_full = dat_full %>% add_column(worldclim_temp_range = (dat_full$worldclim_temp_max - dat_full$worldclim_temp_min))
  
  temp_type = "worldclim_temp_mid"
  
  modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability+worldclim_temp_range+RFE.mean" 
  
  #########################################################################################################
  ### LOAD TEMPSUIT DATA ###
  #########################################################################################################
  ### a
  mordecai_biterate <- read_csv("Data/mordecai_biterate.csv")
  hamlet_biterate = read_csv("Data/hamlet_biterate.csv")
  mordecai_biterate$author = "mordecai"
  hamlet_biterate$author = "hamlet"
  names(hamlet_biterate) = names(mordecai_biterate)
  
  dat_bite = rbind(mordecai_biterate, hamlet_biterate)
  
  ### mu
  dat_mort <- read_csv("Data/Survival_mortality/SurvivalData_mordecai2018.csv")
  dat_mort = filter(dat_mort, `Time (dpi)`>0)
  
  ### PDR
  dat_EIP <- read_csv("Data/davis_EIP.csv")
  
  ########################################################################################################
  ### CREATE STATUS ###
  #########################################################################################################
  
  #first reset seed for generating the pseudo random start conditions
  index_code = as.numeric(commandArgs(TRUE))
  if(identical(index_code, numeric(0))) index_code=1
  
  t= as.numeric(Sys.time())
  seed= (t - floor(t)) * 1e8 
  set.seed(seed)
  
  
  #### INITIAL PARAM ####
  
  ### OR LOAD FROM PREVIOUS RUN ###
  file_name = get_latest_file(pattern = "GLM_tempsuit_parameter_estimates")
  prev_param = read.csv(file_name)
  
  pars_ini = c(prev_param$median)
  names(pars_ini) = prev_param$Parameter

  
  #get in the right order
  pars_ini = pars_ini[c("Intercept","log.surv.qual.adm0","adm05AGO",
                        "adm05BDI","adm05ERI","adm05ETH","adm05GNB",
                        "adm05KEN","adm05MRT","adm05RWA" ,
                        "adm05SDN" ,"adm05SOM" ,"adm05SSD" ,"adm05TZA",
                        "adm05UGA","adm05ZMB" , "lon" ,"logpop","temp_suitability",
                        "worldclim_temp_range", "RFE.mean",
                        "a_T0" ,"a_Tm" ,"a_c", "mu_T0","mu_Tm",
                        "mu_c" ,"PDR_T0","PDR_Tm","PDR_c") ]
  
  pars_ini[grep("^adm05", names(pars_ini))] = -1
  pars_ini[grep("^adm05AGO", names(pars_ini))] = 1
  #########################################################################################################
  ### MCMC ###
  #########################################################################################################
  #create a directory to save the output in
  date=format(Sys.time(),"%Y%m%d")
  name_dir = paste0("GLM_tempsuit_MCMC_chain", "_", date)
  dir.create(name_dir)
  
  Niter = 500000
  
  
  GLM_tempsuit_MCMC( Niter,
                     name_dir,
                     modelVec,
                     pars_ini,
                     dat_full,
                     temp_type,
                     dat_bite,
                     dat_mort,
                     dat_EIP,
                     c34,
                     run_id,
                     plot_chain = FALSE)
  
}
