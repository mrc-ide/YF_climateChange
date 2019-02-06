
run_id = 1
### ESTIMATE BITE RATE ###

library(ggplot2)
library(dplyr)

library(YFestimation)

### DATA ###

davis_EIP <- read_csv("Z:/YF_climateChange/Data/davis_EIP.csv")

################################################################
### FUNCTIONS ###
################################################################
source('Z:/YF_climateChange/FUNCTIONS/temp_suitability.R')

likelihood = function(param, dat){
  
  T0 = as.numeric( param[1] )
  Tm = as.numeric( param[2])
  c = as.numeric( param[3])
  
  Temp = dat$T
  
  PDR = briere(Temp, T0, Tm , c)
  PDR[PDR<=0] = 1e-4
  

    LL = sum( log( dnorm(dat$PDR, mean = PDR, sd = 0.08)) ) # sd set from fit to data

  
  return(LL)
}
################################################################
prior = function(param){
  
  T0 = as.numeric(param[1] )
  Tm = as.numeric(param[2])
  c = as.numeric(param[3])
  
  prior_prob = log( dnorm(T0, mean = 18.3, sd = 3) ) +  # these are set from Tesla for zikv
               log( dnorm(Tm, mean = 42.3, sd =1) ) +
               log( dnorm(c, mean = 1.74e-4, sd =5e-5))
  
  return(prior_prob)
}
################################################################
MCMC_step = function(param, dat, chain_cov, adapt, accCurrent){
  
  #new param
  param_prop = YFestimation::GLMproposal(as.numeric(param), 
                                         chain_cov, 
                                         adapt)
  names(param_prop) = names(param)
  
  #prior_prop
  prior_prop = prior( exp( param_prop) )
  
  #if fintie calc likelihood
  if(is.finite(prior_prop)){
    like_prop = likelihood( exp( param_prop ), dat)
    
    accProp = like_prop + prior_prop + sum(param_prop) #correction for log scale
    
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
pars_ini = log( data.frame("T0" = 18.3,
                           "Tm" = 42.3,
                           "c" = 1.74e-4) )


#########################################################################################################
### MCMC ###
#########################################################################################################

### find posterior probability at start ###
out = MCMC_step(pars_ini,
                davis_EIP,
                chain_cov=1,
                adapt=0, 
                accCurrent=-Inf)

#create a directory to save the output in
date=format(Sys.time(),"%Y%m%d")
name_dir = paste0("PDR_MCMC_chain", "_", date)
if(!dir.exists(name_dir)) dir.create(name_dir,  showWarnings = TRUE)


#plain old mcmc on fixed model prob
Niter = 1e4
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
                  davis_EIP,
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
chain_sub1 = chain[1000:nrow(chain),]

chain_sub2 = convert.mcmc.list(chain_sub1)

# chain_sub = ggs(chain_sub2 )
# 
# ggs_traceplot(chain_sub)
# 
# ggs_histogram(chain_sub)


hpdout = HPDinterval(chain_sub2)

param_est = data.frame("lower" = exp(hpdout[[1]][,1]) ,
                       "median" = exp( apply(chain_sub1, 2, median) ),
                       "upper" = exp(hpdout[[1]][,2]) )
write.csv(param_est, paste0(name_dir, "/", "PDR_estimates.csv") )

# plot fit

plot(davis_EIP$T, briere(davis_EIP$T, T0 = exp( median(chain_sub1[,1]) ),
                   Tm = exp( median(chain_sub1[,2]) ),
                   c = exp( median(chain_sub1[,3]) )  ), type = "l", ylim = c(0,0.5), lwd =1, 
     xlab = "Temperature", ylab = "PDR")

for (i in 1:1000){
  lines(davis_EIP$T, briere(davis_EIP$T, T0 = exp( chain_sub1[sample(nrow(chain_sub1), 1),1] ),
                      Tm = exp( chain_sub1[sample(nrow(chain_sub1), 1),2]) ,
                      c = exp( chain_sub1[sample(nrow(chain_sub1), 1),3])   ), col = rgb(0,0,0,alpha=0.01)  )
}
points(davis_EIP$T, davis_EIP$PDR, col = "red", pch =20)

dev.copy(png,paste0(name_dir, "/",'PDR_fit.png') )
dev.off()

