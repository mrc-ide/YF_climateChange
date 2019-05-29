###########################################################
### FUNCTIONs for temp suitability ###
###########################################################
briere = function(Temp, parm){
  out = parm[3]*Temp*(Temp-parm[1])*(parm[2]-Temp)^0.5
  return(out)
}

quad = function(Temp, parm){
  out = -parm[3]*(Temp-parm[1])*(parm[2] -Temp)
  return(out)
}
###########################################################

###########################################################
### FUNCTION to calculate temperature suitability index ###
###########################################################

temp_suitability = function(Temp, param){
  
  #bite rate
  a = briere(Temp, 
             as.numeric(param[ grep("^a_", names(param))] )) # new estimate medians
  
  # mu = m mortality
  lf = quad(Temp, 
            as.numeric(param[ grep("^mu_", names(param))]))
  mu = 1/ ifelse( lf<=0, 1, lf ) #guard against negatives and zeros
  
  
  # PDR = parasite development rate = 1/EIP  
  PDR = briere(Temp, 
               as.numeric(param[ grep("^PDR_", names(param))])) #
  PDR[PDR<0]=0
  
  a[is.na(a)] = PDR[is.na(PDR)] = 0
  
  # suitability
  Z = (a^2 * exp(-mu / PDR) ) / mu 
  
  return(Z)
}



##################################################################################
### temp suit PRIOR ###
##################################################################################
fun_tempsuitPrior = function( param ){
  
  #prior distributions for the three parameters defining bite rate, mortality rate
  #and parasite development rate
  prior_prob_a =  log( dtrunc(param[names(param) == "a_T0"] , 
                              "norm",
                              a=0,
                              b = Inf,
                              mean = 13.35, 
                              sd = 2*2)) +   #set from Mordecai
    dnorm(param[names(param) == "a_Tm"], 
          mean = 40.08, 
          sd = 0.05*2, 
          log = TRUE)  +
    log( dtrunc(param[names(param) == "a_c"], 
                "norm", 
                a = 0, 
                b = Inf,  
                mean = 2.02e-4, 
                sd = 2e-5*2) ) #keep it positive
  
  prior_prob_mu = log( dtrunc(param[names(param) == "mu_T0"], 
                              "norm",
                              a = 0,
                              b = Inf,
                              mean = 11, 
                              sd = 1*2, 
                              log = TRUE))  +    # set from Tesla
    dnorm(param[names(param) == "mu_Tm"], 
          mean = 37, 
          sd = 0.8*2, 
          log = TRUE)  +
    log( dtrunc(param[names(param) == "mu_c"], 
                "norm", 
                a=-Inf, 
                b = 0, 
                mean = -3e-1, 
                sd = 2e-2*2) ) #keep it negative
  
  prior_prob_PDR = log( dtrunc(param[names(param) == "PDR_T0"], 
                               "norm",
                               a = 0,
                               b = Inf,
                               mean = 18.3, 
                               sd = 3*2, 
                               log = TRUE)) +  # these are set from Tesla for zikv
    dnorm(param[names(param) == "PDR_Tm"], 
          mean = 42.3, 
          sd =1*2, 
          log = TRUE) +
    log( dtrunc(param[names(param) == "PDR_c"], 
                "norm", 
                a = 0, 
                b = Inf, 
                mean = 1.74e-4, 
                sd =5e-5*2) )
  
  
  return(as.numeric( prior_prob_a + prior_prob_mu + prior_prob_PDR ))
}

##################################################################################
### temp suit LIKE ###
##################################################################################
fun_tempsuitLike = function(dat_bite, dat_mort, dat_EIP, param){
  
  #bite rate
  Temp = dat_bite$T
  a = briere(Temp, 
             param[ grep("^a_", names(param))])
  a[a<=0] = 1e-4
  LL_a = sum(  dexp(dat_bite$bite_rate, rate = a, log = TRUE) )
  
  #mortality
  Temp = unique(dat_mort$Temp)
  lf = quad(Temp, 
            param[ grep("^mu_", names(param))])
  lf[lf<0] = 0
  mu = ifelse( lf<=0, 1, 1/ lf )+ 1e-8
  LL_mu = NULL
  for (i in 1:length(Temp)){
    
    dat_sub  = filter(dat_mort, Temp == Temp[i])
    
    LL_mu = rbind( LL_mu, sum(  dbinom(dat_sub$Dead, 
                                       size = dat_sub$Dead+dat_sub$Alive, 
                                       prob = mu[i] , 
                                       log = TRUE) , na.rm = TRUE) )
  }
  
  #PDR
  Temp = dat_EIP$T
  PDR = briere(Temp, 
               param[ grep("^PDR_", names(param))])
  PDR[PDR<=0] = 1e-4
  LL_PDR = sum(  dnorm(dat_EIP$PDR, 
                       mean = PDR, 
                       sd = 0.08, 
                       log = TRUE) )
  
  # put together
  LL = LL_a + sum(LL_mu, na.rm = TRUE) + LL_PDR
  
  return(LL)
}