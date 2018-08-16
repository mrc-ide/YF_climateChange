#transdimensional mcmc for estimating glm
run_id =1
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
# shp0 = readShapePoly(paste0(shpdir, "All/","gadm28_adm0.shp")) #if gadm2.8 they do not have just Africa
# shp1 = readShapePoly(paste0(shpdir, "All/","gadm28_adm1.shp"))

#adjust titles
shp1$adm0_adm1 = paste(shp1$ISO, shp1$ID_1, sep="_")
shp1 = shp1[order(shp1$adm0_adm1),]

#read countries in
Countries = read_csv(paste0("../Data/","Countries.csv"))
c34 = Countries$c34
country34 = Countries$country34


#sort age vector
maxAge = 100
ageVec = c(0:maxAge)


#########################################################################################################
### SOURCE FUNCTIONS ###
#########################################################################################################

R.utils::sourceDirectory("FUNCTIONS")

#########################################################################################################
### LOAD ENVIRONMENTAL DATA ###
#########################################################################################################

Env_Table_path = (paste0("../Data/","Environment/Africa_adm1_dat_2017.csv")) #this file is adapted by hand to have latest outbreaks

dat_full = read.csv(Env_Table_path, stringsAsFactors=F)

### TEMP SUITABILITY ###
dat_full = cbind(dat_full, temp_suitability_mordecai(dat_full[,"ERAday.mean"]))
names(dat_full)[ncol(dat_full)] = "temp_suitability"

envdat = launch_env_dat(dat_full,c34)
dat = envdat$dat

#########################################################################################################
### FIT GLM ###
#########################################################################################################

#read in models
#modelVec=  "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+ERAday.mean "                                            
modelVec = "cas.or.out~log.surv.qual.adm0+adm05+lon+logpop+temp_suitability " 

object_glm = fit_glm(dat =envdat$dat, depi = envdat$depi, modelVec )   #fit_glm from Kevin's code, adapted to take models as input
names(object_glm)

beta0 = object_glm[[1]]
x = object_glm[[2]]
y = object_glm[[3]]



########################################################################################################
### CREATE STATUS ###
#########################################################################################################

#first reset seed for generating the pseudo random start conditions
index_code = as.numeric(commandArgs(TRUE))
if(identical(index_code, numeric(0))) index_code=1

t= as.numeric(Sys.time())
seed= (t - floor(t)) * 1e8 
set.seed(seed)

StartParamFoi=read.csv(paste0("../YellowFeverModelEstimation2017/","StartParam_","Foi",".csv"),header = T)
StartParamR0=read.csv(paste0("../YellowFeverModelEstimation2017/","StartParam_","R0",".csv"),header = T)

## itialising parameters for MCMC
parnames =  paste("log", names(beta0), sep = ".")

pars_ini = rep(NA,length(parnames))
names(pars_ini) = parnames


# GLM parameters
ii = 2:(length(beta0)+1)
pars_ini = beta0


#########################################################################################################
### MCMC ###
#########################################################################################################



### find posterior probability at start ###
out = GLM_MCMC_step(pars_ini,
                    x,
                    y,
                    chain_cov=1,
                    adapt=0, 
                    accCurrent=-Inf)






#create a directory to save the output in
date=format(Sys.time(),"%Y%m%d")
name_dir = paste("GLM_MCMC_chain", "_", date, sep='')
if(!dir.exists(name_dir)) dir.create(name_dir,  showWarnings = TRUE)



#plain old mcmc on fixed model prob
Niter = 1e5
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
  
  
  if(iter>100 & run_id<100) plot(chain[,1] , type= "l") #plot(posteriorProb, type = 'l')
  
  if (iter %% 1000 == 0){
    
    colnames(acceptRate) = "acceptRate"
    colnames(posteriorProb) = "posteriorProb"
    if (iter %% 10000 == 0){
      fileIndex  = iter/10000
    }
    write.csv(cbind(chain,  posteriorProb, acceptRate)[min((fileIndex * 10000+1),iter):iter,], 
              paste0(name_dir,"/","GLM_chain",fileIndex,"_output",run_id,".csv") ) 
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
  out = GLM_MCMC_step(param,
                      x,
                      y,
                      chain_cov,
                      adapt, 
                      accCurrent)
  
  
}




