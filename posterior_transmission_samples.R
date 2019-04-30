#

### script to run climate change burden for all scenarios ###

### --------------------------------------------------------------------- ###
library(maptools)
library(sp) 
library(shapefiles)
library(Hmisc)
library(fields)
library(dplyr)
library(EnvStats)
library(readr)
library(reshape)
library(abind)
library(mvtnorm)
library(RColorBrewer)
library(truncdist)
library(tibble)

library(ggmcmc)
library(mcmcplots)
library(R.utils)

library(YFestimation)
library(snapalette)
library(KsetupR)
library(YFburden)

#-----------------------------------------------------------------------------
sourceDirectory("FUNCTIONS", modifiedOnly = FALSE)

#-----------------------------------------------------------------------------
### load data ###
Env_Table_path = paste0("../Data/","Environment/dat_worldclim_all_2019-04-15.csv") 

dat_full = read.csv(Env_Table_path, stringsAsFactors=FALSE)

temp_type = "worldclim_temp_mid"

modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability+worldclim_temp_range+RFE.mean" 

path = "GLM_tempsuit_MCMC_chain_20190325"

#read countries in
Countries = read_csv(paste0("../Data/","Countries.csv"))
c34 = Countries$c34
country34 = Countries$country34

#-----------------------------------------------------------------------------
### get_chains
mcmc_out = get_chains(path, burnin = 1.5e6, thin = 100)

plot(mcmc_out$posteriorProb, type = "l", ylab = "Posterior probability", xlab = "Iteration")

ind = which.max(mcmc_out$posteriorProb)

mcmc_out = mcmc_out[, 1:(ncol(mcmc_out)-2)]

#-----------------------------------------------------------------------------
### set up fit ###

# serology #
Serology = read.csv(paste0("../Data/","Serology/serology.csv"), stringsAsFactors = FALSE)
seroout = process_serology(Serology)

# population data #
path = "../Data/"

#launch env dat just for the time being
envdat = launch_env_dat(filepath=NA, c34, dat_full)

#function to collect , get totals for each year and proportions in each year
all_res_pop_3d = get_pop_data_3d(path = path, c_country=c34, dat=envdat$dat)

pop1 = all_res_pop_3d$pop1                                            #population import
pop3d = all_res_pop_3d$pop3d                                      #populations in 3d array
P_tot_2d = all_res_pop_3d$P_tot_2d                                    #total populations for each adm and year
p_prop_3d = all_res_pop_3d$p_prop_3d                                    #proportions of population

#get names
dim_adm  = dimnames(pop3d)[[1]]
dim_year = as.numeric(dimnames(pop3d)[[2]])
dim_age  = dimnames(pop3d)[[3]]

# vaccination data #
vaccdir = paste0("../Data/", "Vaccination/")
# latest_vaccine_csv = "vaccination_coverage_by_adm1_year_age_base_skew0_update_2016-10-16.csv"
latest_vaccine_csv = "Outputs/adm1_old/vaccination_coverage_by_adm1_year_age_base_skew0.csv"   

vc2d = read.csv(paste0(vaccdir,latest_vaccine_csv), 
                stringsAsFactors = FALSE) 

names(vc2d)[names(vc2d)=="country"]= "adm0"                          #rename countries as adm0
names(vc2d)[names(vc2d)=="adm1"]= "adm0_adm1"                        #renames adm1 as adm0_adm1

# formally "repair_vc_data" from FOI model in Kevin's folder
for (colIndex in 3:ncol(vc2d)){                                      
  vc2d[,colIndex] = ifelse(is.na(vc2d[,colIndex]), vc2d[,colIndex-1], vc2d[,colIndex])
}
# restrict to lines in dat
vc2d = vc2d[vc2d[,"adm0_adm1"] %in% envdat$dat[,"adm0_adm1"],]

#vc3d
vc3d = transform_into_vc3d(vc2d,  adm="adm1")

# t0_vac_africa #
t0_vac_africa = calc_t0_vac_africa(vc3d)

# inc_v3d #
inc_v3d = calc_incidence_vac_general(vc3d)

# CALCULATE population moments #
pop_moments_whole = calc_pop_moments(p_prop_3d, 
                                     t0_vac_africa,
                                     dim_adm,
                                     dim_year,
                                     dim_age)

# aggregate #


list_aggregate_pop_vc =Make_aggregate_pop_vc_3d(pop1=pop1, 
                                                vc2d=vc2d, 
                                                sero_studies=seroout$sero_studies, 
                                                adm1s=seroout$adm1s)
pop_agg3d = list_aggregate_pop_vc$pop_agg3d
vc_agg3d = list_aggregate_pop_vc$vc_agg3d

#calculate aggregated incidence (same function as before)
inc_v3d_agg = calc_incidence_vac_general(vc_agg3d);dim(inc_v3d_agg)

#calculate aggregated moments (different fucntion before)
pop_moments_agg = calc_pop_moments_agg(pop_agg3d,
                                       seroout$t0_vac,
                                       dim_year,
                                       seroout$study_years)


# R0 lookup if needed #
load(paste0("../YellowFeverModelEstimation2017/","R0_lookup_table.Rdata") )

# pop at survey #
foi_const_surv = rep(0, seroout$no_sero_surveys)

list_pop_at_survey = create_pop_at_survey(pop_agg3d, 
                                          seroout$sero_studies,
                                          dim_year)
p_at_survey = list_pop_at_survey$p_at_survey_3d
P_tot_survey = list_pop_at_survey$P_tot_survey_2d


#-----------------------------------------------------------------------------
# import serology fit

filepath = "Z:/MultiModelInference/chains/multi_model_MCMC_chain_20180510" 


mcmc_out_sero = get_chains(filepath, burnin = 1, thin = 100)

#-----------------------------------------------------------------------------

fun_sample_transmission = function(sample_ind){
  runs_clim_change = NULL
  for(year in c(2050, 2070, "now")){
    
    if(year == "now"){scenario = "now"} else {
      
      for(scenario in c(26, 45, 60, 85)){
        
        dat_tmp = prepare_climate_dat(dat_full, year , scenario)
        ### TEMP SUITABILITY ###
        dat_full_temp = cbind(dat_tmp,
                              temp_suitability(dat_tmp[,temp_type] , 
                                               mcmc_out[sample_ind,22:30]) )
        names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
        
        envdat = YFestimation::launch_env_dat(filepath = NA, dat_full= dat_full_temp , c34 = c34)  
        
        ### GET x ###
        
        object_glm = fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
        x = object_glm[[2]]
        
        
        ii= 2:22
        
        varsin_nc=ii[-grep("adm0",colnames(x))] - 1 
        
        mcmc_out_f = filter(mcmc_out_sero, model_chain == 0)
        
        adjusted_params = c(exp(median(mcmc_out_f[,1])), 
                            as.numeric(mcmc_out[sample_ind,1:21]), #apply(mcmc_out[,1:21], 2, median, na.rm = T), # 
                            exp(apply(mcmc_out_f[,c(2:41)], 2, median, na.rm = T)),
                            exp(median(mcmc_out_f[,ncol(mcmc_out_f)])) )
        
        names(adjusted_params)[c(1,length(adjusted_params))] = c("vac_eff", "vc_factor_CMRs")
        
        
        runs = YFestimation::fun_calc_transmission_Africa(x ,
                                                          ii ,
                                                          seroout ,
                                                          params = adjusted_params ,
                                                          dat = envdat$dat ,
                                                          t0_vac_africa ,
                                                          dim_year ,
                                                          dim_age ,
                                                          p_prop_3d ,
                                                          P_tot_2d ,
                                                          inc_v3d ,
                                                          pop1,
                                                          vc2d,
                                                          varsin_nc,
                                                          polydeg = 5,
                                                          R0_lookup,
                                                          model_type = "Foi")
        names(runs) = envdat$dat$adm0_adm1
        
        runs = as.data.frame(cbind(runs, adm0 = envdat$dat$adm0))
        
        runs_adm0 = runs %>% group_by(adm0) %>% summarise(FOI = mean(as.numeric(as.character(runs))))
        
        runs_clim_change = rbind(runs_clim_change, cbind(runs_adm0, data.frame("year" = year, "scenario" = scenario, sample = sample_ind)))
      }
    }
  }
  return(runs_clim_change)
}

n_samples = 1000

all_runs = lapply(base::sample(1:nrow(mcmc_out), n_samples), FUN =fun_sample_transmission)

all_runs_out = bind_rows(all_runs)
