
library(dplyr)
library(magrittr)
library(tidyr)

library(YFestimation)
library(snapalette)
library(KsetupR)

library(ggplot2)

foi <- read.csv("transmission_intensity_samples_simple.csv", stringsAsFactors = FALSE)

foi %<>% dplyr::filter(year == "now", scenario == 26)

foi %<>% group_by(adm1) %>% mutate(sample = 1:n())

foi %>%
  ggplot() +
  geom_boxplot()+
  aes(x = adm0, y = runs)

# import serology fit

filepath = "Z:/MultiModelInference/multi_model_MCMC_chain_20180622" 


mcmc_out_sero = get_chains(filepath, burnin = 1, thin = 1)

mcmc_out_f = dplyr::filter(mcmc_out_sero, model_chain == 0)

mcmc_out_f %<>% select(starts_with("Foi"))

mcmc_out_f %<>% mutate(sample = 1:nrow(mcmc_out_f))

mcmc_out_f %<>% pivot_longer(names_to = "parameter", values_to = "value", -sample)

mcmc_out_f %<>% mutate(value = exp(value))

mcmc_out_f %<>% mutate(adm0 = substr(gsub("Foi_", "", parameter), 1, 3))

mcmc_out_f %<>% filter(sample < 1001)


comp <- left_join(mcmc_out_f, foi)

comp %>% 
  filter(adm0 != "CAF") %>%
  ggplot() + 
  geom_point(alpha = 0.1, colour = "blue") + 
  aes(x = value, y = FOI) + 
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0, size = 2) +
  theme_bw()
#--------------------------------------------------------------
#
mcmc_out_f %<>% dplyr::filter(sample < 101)
mcmc_out_f %<>% mutate(study = gsub("Foi_", "", parameter))

df <- NULL
for(i in 1:seroout$no_sero_surveys){
  foi_tmp <- foi %>% dplyr::filter(adm1 %in% seroout$adm1s[[i]])
  mcmc_out_tmp <- mcmc_out_f %>% filter(study == seroout$sero_studies[i])
  
  foi_tmp %<>% group_by(sample) %>% summarise(foi = median(runs)) %>% unique()
  
  df %<>% bind_rows(data.frame(study =  seroout$sero_studies[i],
                               sero = mcmc_out_tmp$value,
                               glm = foi_tmp$foi,
                               adm0 = substr(seroout$sero_studies[i], 1, 3)))

}

df %>%
  ggplot() +
  geom_point() +
  aes(x = sero, y = glm, colour = adm0) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Serology estimates",
       y = "GLM estimates",
       colour = "Country")

















# ### script to run climate change burden for all scenarios ###
# 
# ### --------------------------------------------------------------------- ###
# library(maptools)
# library(sp) 
# library(shapefiles)
# library(Hmisc)
# library(fields)
# library(dplyr)
# library(EnvStats)
# library(readr)
# library(reshape)
# library(abind)
# library(mvtnorm)
# library(RColorBrewer)
# library(truncdist)
# library(tibble)
# 
# library(ggmcmc)
# library(mcmcplots)
# library(R.utils)
# 
# library(YFestimation)
# library(snapalette)
# library(KsetupR)
# library(YFburden)
# 
# #-----------------------------------------------------------------------------
# sourceDirectory("FUNCTIONS", modifiedOnly = FALSE)
# 
# #-----------------------------------------------------------------------------
# ### load data ###
# Env_Table_path = paste0("../Data/","Environment/dat_worldclim_all_2019-04-15.csv") 
# 
# dat_full = read.csv(Env_Table_path, stringsAsFactors=FALSE)
# 
# temp_type = "worldclim_temp_mid"
# 
# modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability+worldclim_temp_range+worldclim_rainfall" 
# 
# path = "test_GLM_tempsuit_MCMC_chain_20190603"
# 
# #read countries in
# Countries = read_csv(paste0("../Data/","Countries.csv"))
# c34 = Countries$c34
# country34 = Countries$country34
# 
# #-----------------------------------------------------------------------------
# ### get_chains
# mcmc_out = get_chains(path, burnin = 2e6, thin = 100)
# 
# plot(mcmc_out$posteriorProb, type = "l", ylab = "Posterior probability", xlab = "Iteration")
# 
# ind = which.max(mcmc_out$posteriorProb)
# 
# mcmc_out = mcmc_out[, 1:(ncol(mcmc_out)-2)]
# 
# #-----------------------------------------------------------------------------
# ### set up fit ###
# 
# # serology #
# Serology = read.csv(paste0("../Data/","Serology/serology.csv"), stringsAsFactors = FALSE)
# Serology = Serology %>% filter(country_zone != "CAF") # IGNORE CAF
# seroout = process_serology(Serology)
# 
# # population data #
# path = "../Data/"
# 
# #launch env dat just for the time being
# envdat = launch_env_dat(filepath=NA, c34, dat_full)
# 
# #function to collect , get totals for each year and proportions in each year
# all_res_pop_3d = get_pop_data_3d(path = path, c_country=c34, dat=envdat$dat)
# 
# pop1 = all_res_pop_3d$pop1                                            #population import
# pop3d = all_res_pop_3d$pop3d                                      #populations in 3d array
# P_tot_2d = all_res_pop_3d$P_tot_2d                                    #total populations for each adm and year
# p_prop_3d = all_res_pop_3d$p_prop_3d                                    #proportions of population
# 
# #get names
# dim_adm  = dimnames(pop3d)[[1]]
# dim_year = as.numeric(dimnames(pop3d)[[2]])
# dim_age  = dimnames(pop3d)[[3]]
# 
# # vaccination data #
# vaccdir = paste0("../Data/", "Vaccination/")
# # latest_vaccine_csv = "vaccination_coverage_by_adm1_year_age_base_skew0_update_2016-10-16.csv"
# latest_vaccine_csv = "Outputs/adm1_old/vaccination_coverage_by_adm1_year_age_base_skew0.csv"   
# 
# vc2d = read.csv(paste0(vaccdir,latest_vaccine_csv), 
#                 stringsAsFactors = FALSE) 
# 
# names(vc2d)[names(vc2d)=="country"]= "adm0"                          #rename countries as adm0
# names(vc2d)[names(vc2d)=="adm1"]= "adm0_adm1"                        #renames adm1 as adm0_adm1
# 
# # formally "repair_vc_data" from FOI model in Kevin's folder
# for (colIndex in 3:ncol(vc2d)){                                      
#   vc2d[,colIndex] = ifelse(is.na(vc2d[,colIndex]), vc2d[,colIndex-1], vc2d[,colIndex])
# }
# # restrict to lines in dat
# vc2d = vc2d[vc2d[,"adm0_adm1"] %in% envdat$dat[,"adm0_adm1"],]
# 
# #vc3d
# vc3d = transform_into_vc3d(vc2d %>% select(-adm0),  adm="adm1")
# 
# # t0_vac_africa #
# t0_vac_africa = calc_t0_vac_africa(vc3d)
# 
# # inc_v3d #
# inc_v3d = calc_incidence_vac_general(vc3d)
# 
# # CALCULATE population moments #
# pop_moments_whole = calc_pop_moments(p_prop_3d, 
#                                      t0_vac_africa,
#                                      dim_adm,
#                                      dim_year,
#                                      dim_age)
# 
# # aggregate #
# 
# 
# list_aggregate_pop_vc =Make_aggregate_pop_vc_3d(pop1=pop1, 
#                                                 vc2d=vc2d, 
#                                                 sero_studies=seroout$sero_studies, 
#                                                 adm1s=seroout$adm1s)
# pop_agg3d = list_aggregate_pop_vc$pop_agg3d
# vc_agg3d = list_aggregate_pop_vc$vc_agg3d
# 
# #calculate aggregated incidence (same function as before)
# inc_v3d_agg = calc_incidence_vac_general(vc_agg3d);dim(inc_v3d_agg)
# 
# #calculate aggregated moments (different fucntion before)
# pop_moments_agg = calc_pop_moments_agg(pop_agg3d,
#                                        seroout$t0_vac,
#                                        dim_year,
#                                        seroout$study_years)
# 
# 
# # R0 lookup if needed #
# load(paste0("../YellowFeverModelEstimation2017/","R0_lookup_table.Rdata") )
# 
# # pop at survey #
# foi_const_surv = rep(0, seroout$no_sero_surveys)
# 
# list_pop_at_survey = create_pop_at_survey(pop_agg3d, 
#                                           seroout$sero_studies,
#                                           dim_year)
# p_at_survey = list_pop_at_survey$p_at_survey_3d
# P_tot_survey = list_pop_at_survey$P_tot_survey_2d
# 
# 
# #-----------------------------------------------------------------------------
# # import serology fit
# 
# filepath = "Z:/MultiModelInference/multi_model_MCMC_chain_20180622" 
# 
# 
# mcmc_out_sero = get_chains(filepath, burnin = 1, thin = 1)
# 
# mcmc_out_f = filter(mcmc_out_sero, model_chain == 0)
# #-----------------------------------------------------------------------------
# 
# #####
# 
# fun_sample_transmission = function(sample_ind){
#   
#   
#   ii= 2:22
#   #### SET UP ####
#   dat_tmp = prepare_climate_dat(dat_full, "now" , 26)
#   ### TEMP SUITABILITY ###
#   dat_full_temp = cbind(dat_tmp,
#                         temp_suitability(dat_tmp[,temp_type] , 
#                                          mcmc_out[sample_ind,22:30]) )
#   names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
#   
#   envdat = YFestimation::launch_env_dat(filepath = NA, dat_full= dat_full_temp , c34 = c34)  
#   
#   ### GET x ###
#   
#   object_glm = fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
#   x = object_glm[[2]]
#   
#   
#   varsin_nc=ii[-grep("adm0",colnames(x))] - 1 
#   
#   mcmc_out_f = filter(mcmc_out_sero, model_chain == 0)
#   
#   adjusted_params = c(exp(mcmc_out_f[sample_ind,1]), 
#                       as.numeric(mcmc_out[sample_ind,1:21]), 
#                       exp(as.numeric(mcmc_out_f[sample_ind,c(2:41)])),
#                       exp(mcmc_out_f[sample_ind,ncol(mcmc_out_f)]) )
#   
#   names(adjusted_params) = c("vac_eff", names(mcmc_out)[1:21], names(mcmc_out_f)[2:41], "vc_factor_CMRs")
#   
#   
#   ### P DETECT ###
#   #get aggregated vc and pop over observation period
#   aggout=create_pop30_agg_vc30_agg(pop1, vc2d)
#   
#   #glm predictions
#   mypreds_nc  = fun_calcPred(coefs = as.numeric(adjusted_params)[ii],
#                              newdata=x,
#                              type="link",
#                              varsin=varsin_nc)
#   
#   #probability of detection
#   p_detect =  fun_calc_pdetect_multi_both(x,
#                                           ii,
#                                           seroout,
#                                           adjusted_params,
#                                           envdat$dat,
#                                           t0_vac_africa,
#                                           dim_year,
#                                           dim_age,
#                                           p_prop_3d,
#                                           P_tot_2d,
#                                           inc_v3d,
#                                           pop_moments_whole,
#                                           varsin_nc,
#                                           aggout$vc30_agg,
#                                           aggout$pop30_agg,
#                                           model_type = "Foi")
#   p_detect_link = mean(p_detect)
#   
#   polydeg = 5
#   
#   ####
#   
#   runs_clim_change = NULL
#   year = "now"
#   scenario = 26
#       
#       dat_tmp = prepare_climate_dat(dat_full, year , scenario)
#       ### TEMP SUITABILITY ###
#       dat_full_temp = cbind(dat_tmp,
#                             temp_suitability(dat_tmp[,temp_type] , 
#                                              mcmc_out[sample_ind,22:30]) )
#       names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
#       
#       envdat = YFestimation::launch_env_dat(filepath = NA, dat_full= dat_full_temp , c34 = c34)  
#       
#       ### GET x ###
#       
#       object_glm = fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
#       x = object_glm[[2]]
#       
#       
#       ii= 2:22
#       
#       varsin_nc=ii[-grep("adm0",colnames(x))] - 1 
#       
#       
#       adjusted_params = c(exp(mcmc_out_f[sample_ind,1]), 
#                           as.numeric(mcmc_out[sample_ind,1:21]), 
#                           exp(as.numeric(mcmc_out_f[sample_ind,c(2:41)])),
#                           exp(mcmc_out_f[sample_ind,ncol(mcmc_out_f)]) )
#       
#       names(adjusted_params) = c("vac_eff", names(mcmc_out)[1:21], names(mcmc_out_f)[2:41], "vc_factor_CMRs")
#       
#       
#       #glm predictions
#       mypreds_nc  = fun_calcPred(coefs = as.numeric(adjusted_params)[ii],
#                                  newdata=x,
#                                  type="link",
#                                  varsin=varsin_nc)
#       
#       #calculating number of infections over the observation period for the whole region
#       Ninf_whole = exp( mypreds_nc - p_detect_link)
#       
#       pop_vc_moments = aggout$pop_vc_moments
#       
#       if(polydeg>ncol(pop_vc_moments)) error("fun_calc_transmission_Africa: invalid value for polydeg.\n")
#       
#       z = -Ninf_whole
#       
#       if(polydeg>0) for(i in 1:polydeg) {
#         z = cbind(z,(-1)^(i+1)*pop_vc_moments[,i+1]/factorial(i-1))
#       }
#       
#       transmission_whole = sapply(1:nrow(x), function(i) polyroot(z[i,]))
#       transmission_whole[abs(Arg(transmission_whole))<=1e-10] = Re(transmission_whole)[abs(Arg(transmission_whole))<=1e-10]
#       transmission_whole[abs(Arg(transmission_whole))>1e-10] = NA
#       
#       dt = dim(transmission_whole)
#       transmission_whole = as.numeric(transmission_whole)
#       dim(transmission_whole) = dt
#       transmission_whole = apply(transmission_whole,2,min,na.rm=T)
#       # ------------------------------------------------------------------------------------#
#       
#       runs = transmission_whole
#       
#       
#       names(runs) = envdat$dat$adm0_adm1
#       
#       runs = as.data.frame(cbind(runs, adm0 = envdat$dat$adm0, adm1 = envdat$dat$adm0_adm1))
#       
#       
#       runs_clim_change = rbind(runs_clim_change, 
#                                cbind(runs, 
#                                      data.frame("year" = year, 
#                                                 "scenario" = scenario, 
#                                                 sample = sample_ind)))
# 
#   return(runs_clim_change)
# }
# 
# n_samples = 100
# 
# all_runs = lapply(base::sample(1:min(nrow(mcmc_out), nrow(mcmc_out_f)), n_samples), FUN =fun_sample_transmission)
# 
# all_runs_out = bind_rows(all_runs)
# 
# write.csv(all_runs_out, "transmission_intensity_samples_simple.csv", row.names = FALSE)
# 
# beepr::beep(3)
