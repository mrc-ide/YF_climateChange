#transdimensional mcmc for estimating glm AND temp suitability

#########################################################################################################
### LIBRARIES FOR PACKAGES USED ###
#########################################################################################################

library(Hmisc)
library(fields)
library(dplyr)
library(EnvStats)
library(readr)
library(reshape)
library(abind)
library(mvtnorm)
library(truncdist)

library(YFestimation)

#########################################################################################################
### LOADING COUNTRIES ###
#########################################################################################################

Countries = read_csv(paste0("../Data/","Countries.csv"))
c34 = Countries$c34
country34 = Countries$country34

#########################################################################################################
### SOURCE FUNCTIONS ###
#########################################################################################################

R.utils::sourceDirectory("FUNCTIONS/FUNCTIONS_combined", modifiedOnly = FALSE)

#########################################################################################################
### LOAD ENVIRONMENTAL DATA ###
#########################################################################################################

Env_Table_path = (paste0("../Data/","Environment/Africa_adm1_dat_2017.csv")) 

dat_full = read.csv(Env_Table_path, stringsAsFactors=F)

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


#### INITAL PARAM ####

### OR LOAD FROM PREVIOUS RUN ###
prev_param = read.csv("GLM_tempsuit_MCMC_chain_20180823_hamlet/GLM_tempsuit_parameter_estimates.csv")

pars_ini = c(prev_param$median)
names(pars_ini) = prev_param$Parameter

#get in the right order
pars_ini = pars_ini[c("Intercept","log.surv.qual.adm0","adm05AGO","adm05BDI","adm05ERI","adm05ETH","adm05GNB",
                      "adm05KEN","adm05MRT","adm05RWA" ,"adm05SDN" ,"adm05SOM" ,"adm05SSD" ,"adm05TZA",
                      "adm05UGA","adm05ZMB" , "lon" ,"logpop","temp_suitability",
                      "a_T0" ,"a_Tm" ,"a_c", "mu_T0","mu_Tm","mu_c" ,"PDR_T0","PDR_Tm","PDR_c") ]

#########################################################################################################
### MCMC ###
#########################################################################################################
#create a directory to save the output in
date=format(Sys.time(),"%Y%m%d")
name_dir = paste0("GLM_tempsuit_MCMC_chain", "_", date, "_hamlet")

Niter = 1e6

GLM_tempsuit_MCMC(Niter, name_dir, pars_ini, dat_full, dat_bite, dat_mort, dat_EIP, plot_chain = TRUE)