#transdimensional mcmc for estimating glm AND temp suitability

#########################################################################################################
### LIBRARIES FOR PACKAGES USED ###
#########################################################################################################

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
library(truncdist)

#########################################################################################################
### SETTING THE WORKING DIRECTORY ###
#########################################################################################################

shpdir = paste0("../","shapefiles/gadm2/")


#########################################################################################################
### LOADING SHAPEFILES AND COUNTRIES ###
#########################################################################################################

#These use the shape files- which will need to be updated to gamd2.8 and a csv of the YF endemic 
#countries

#read shapefiles in
shp0 = readShapePoly(paste0(shpdir, "Africa_adm0.shp")) #if gadm2
shp1 = readShapePoly(paste0(shpdir, "Africa_adm1.shp"))

#adjust titles
shp1$adm0_adm1 = paste(shp1$ISO, shp1$ID_1, sep="_")
shp1 = shp1[order(shp1$adm0_adm1),]

#read countries in
Countries = read_csv(paste0("../Data/","Countries.csv"))
c34 = Countries$c34
country34 = Countries$country34



#########################################################################################################
### SOURCE FUNCTIONS ###
#########################################################################################################

R.utils::sourceDirectory("FUNCTIONS/FUNCTIONS_combined", modifiedOnly = FALSE)
source("FUNCTIONS/GLMonly_functions.R")

#########################################################################################################
### LOAD ENVIRONMENTAL DATA ###
#########################################################################################################

Env_Table_path = (paste0("../Data/","Environment/Africa_adm1_dat_2017.csv")) #this file is adapted by hand to have latest outbreaks

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

# ### temp suit param
# # take the median fitted parameters from individual runs
# pars_ini_ts = c( 8.311805221, 40.10016485, 0.000226577,
#                  10.20297414, 38.06140501, -0.55457294117,
#                  18.61394525, 42.19607364, 0.000151139 )
# names(pars_ini_ts) = c("a_T0", "a_Tm", "a_c", 
#                        "mu_T0", "mu_Tm", "mu_c", 
#                        "PDR_T0", "PDR_Tm", "PDR_c")
# 
# 
# ### TEMP SUITABILITY ###
# dat_full_temp = cbind(dat_full, temp_suitability(dat_full[,"ERAday.mean"] , pars_ini_ts))
# names(dat_full_temp)[ncol(dat_full_temp)] = "temp_suitability"
# envdat = launch_env_dat(dat_full_temp,c34)
# 
# ### GET x ###
# modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability" 
# object_glm = fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec ) 
# beta0 = object_glm[[1]]
# x = object_glm[[2]]
# y = object_glm[[3]]
# 
# # add GLM parameters
# pars_ini = c( beta0, pars_ini_ts)

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

### find posterior probability at start ###
out = GLM_tempsuit_MCMC_step(pars_ini,
                             dat_full,
                             dat_bite,
                             dat_mort,
                             dat_EIP,
                             chain_cov=1,
                             adapt=0, 
                             accCurrent=-Inf)

#create a directory to save the output in
date=format(Sys.time(),"%Y%m%d")
name_dir = paste0("GLM_tempsuit_MCMC_chain", "_", date, "_hamlet")
if(!dir.exists(name_dir)) dir.create(name_dir,  showWarnings = TRUE)



#plain old mcmc on fixed model prob
Niter = 1e6
chain = posteriorProb = acceptRate = NULL
burnin = 50
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
  
  
  if(iter>100 ) plot(chain[,1] , type= "l") #plot(posteriorProb, type = 'l')
  
  if (iter %% 1000 == 0){
    
    colnames(acceptRate) = "acceptRate"
    colnames(posteriorProb) = "posteriorProb"
    if (iter %% 10000 == 0){
      fileIndex  = iter/10000
    }
    write.csv(cbind(chain,  posteriorProb, acceptRate)[min((fileIndex * 10000+1),iter):iter,], 
              paste0(name_dir,"/","GLM_tempsuit_chain",fileIndex,"_output",".csv") ) 
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
                               chain_cov,
                               adapt, 
                               accCurrent)
  
}




