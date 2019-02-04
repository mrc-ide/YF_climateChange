
run_id = 1
### ESTIMATE MORTALITY ###

library(ggplot2)
library(dplyr)
library(readr)
library(mvtnorm)

library(YFestimation)

### DATA ###

dat <- read_csv("Z:/YF_climateChange/Data/Survival_mortality/SurvivalData_mordecai2018.csv")

dat = filter(dat, `Time (dpi)`>0)
################################################################
### FUNCTIONS ###
################################################################
source('Z:/YF_climateChange/FUNCTIONS/temp_suitability.R')

likelihood = function(param, dat){
  
  T0 = as.numeric( param[1] )
  Tm = as.numeric( param[2])
  c = as.numeric( param[3])
  
  Temp = unique(dat$Temp)
  
  lf = quad(Temp, T0, Tm , c)
  lf[lf<0] = 0
  
  mu = ifelse( lf<=0, 1, 1/ lf ) +1e-8
  
  LL = NULL
  for (i in 1:length(Temp)){
    
    dat_sub  = filter(dat, Temp == Temp[i])
    
    LL = rbind( LL, sum(  dbinom(dat_sub$Dead, 
                                 size = dat_sub$Dead+dat_sub$Alive, 
                                 prob = mu[i] , 
                                 log = TRUE) , 
                          na.rm = TRUE) )
  }
  
  return( sum(LL, na.rm = TRUE) )
}
################################################################
prior = function(param){
  
  T0 = as.numeric(param[1] )
  Tm = as.numeric(param[2])
  c = as.numeric(param[3])
  
  prior_prob = log( dnorm(T0, mean = 11, sd = 1) ) + 
    log( dnorm(Tm, mean = 37, sd = 0.8) ) +
    log( dnorm(c, mean = -3e-1, sd = 2e-2))
  
  return(prior_prob)
}
################################################################
MCMC_step = function(param, dat, chain_cov, adapt, accCurrent){
  
  #new param
  param_prop = YFestimation::GLMproposal(as.numeric(param), chain_cov, adapt)
  names(param_prop) = names(param)
  
  #prior_prop
  prior_prop = prior(  param_prop )
  
  #if fintie calc likelihood
  if(is.finite(prior_prop)){
    like_prop = likelihood(  param_prop , dat)
    
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


#########################################################################################################
### Setup ###
#########################################################################################################
pars_ini =  data.frame("T0" = 11,
                       "Tm" = 37,
                       "c" = -3e-1 )


#########################################################################################################
### MCMC ###
#########################################################################################################

### find posterior probability at start ###
out = MCMC_step(pars_ini, 
                dat,
                chain_cov=1,
                adapt=0, 
                accCurrent=-Inf)

#create a directory to save the output in
date=format(Sys.time(),"%Y%m%d")
name_dir = paste0("mortality_MCMC_chain", "_", date)
if(!dir.exists(name_dir)) dir.create(name_dir,  showWarnings = TRUE)


#plain old mcmc on fixed model prob
Niter = 5e4
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
              paste0(name_dir,"/","chain",fileIndex,"_output",".csv") ,row.names = FALSE) 
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
  out = MCMC_step(param,
                  dat,
                  chain_cov,
                  adapt, 
                  accCurrent)
  
  
}



####################################################################################
### OUTPUT ###
####################################################################################

library(ggmcmc)
library(mcmcplots)

#remove burnin
chain_sub1 = chain[10000:nrow(chain),]

chain_sub2 = convert.mcmc.list(chain_sub1)

chain_sub = ggs(chain_sub2 )

ggs_traceplot(chain_sub)

ggs_histogram(chain_sub)


hpdout = HPDinterval(chain_sub2)

param_est = data.frame("lower" = hpdout[[1]][,1] ,
                       "median" =  apply(chain_sub1, 2, median ),
                       "upper" = hpdout[[1]][,2] )
write.csv(param_est, paste0(name_dir, "/", "mortality_estimates.csv") )

#### PLOT FIT

Temp = unique(dat$Temp)

lf = quad(Temp, param_est[1,2],param_est[2,2], param_est[3,2])
lf[lf<0] = 0

mu = ifelse( lf<=0, 1, 1/ lf ) +1e-8

plot(dat$Temp, dat$Dead/(dat$Dead+dat$Alive), pch = 20, col="red", xlab = "Temperature", ylab = "Mortality rate")
lines(Temp, mu)
dev.copy(png,paste0(name_dir, "/",'fit.png') )
dev.off()
