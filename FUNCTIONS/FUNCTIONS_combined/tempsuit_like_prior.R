###########################################################
### FUNCTIONs for temp suitability estimation ###
###########################################################

##################################################################################
### temp suit PRIOR ###
##################################################################################
fun_tempsuitPrior = function( param ){
  
  a_T0 = param[1]; a_Tm = param[2]; a_c = param[3]
  mu_T0 = param[4]; mu_Tm = param[5]; mu_c = param[6]
  PDR_T0 = param[7]; PDR_Tm = param[8]; PDR_c = param[9]
  
  prior_prob_a =  dnorm(a_T0, mean = 13.35, sd = 2, log = TRUE) +   #set from Mordecai
                  dnorm(a_Tm, mean = 40.08, sd = 0.05, log = TRUE)  +
                  dnorm(a_c, mean = 2.02e-4, sd = 2e-5, log = TRUE)
  
  prior_prob_mu = dnorm(mu_T0, mean = 11, sd = 1, log = TRUE)  +    # set from Tesla
                  log( dnorm(mu_Tm, mean = 37, sd = 0.8, log = TRUE) ) +
                  log( dnorm(mu_c, mean = -3e-1, sd = 2e-2, log = TRUE))
  
  prior_prob_PDR = dnorm(PDR_T0, mean = 18.3, sd = 3, log = TRUE) +  # these are set from Tesla for zikv
                   dnorm(PDR_Tm, mean = 42.3, sd =1, log = TRUE) +
                   dnorm(PDR_c, mean = 1.74e-4, sd =5e-5, log = TRUE)
  
  prior_prob = prior_prob_a + prior_prob_mu + prior_prob_PDR
  
  return(prior_prob)
}

##################################################################################
### temp suit LIKE ###
##################################################################################
fun_tempsuitLike = function(dat_bite, dat_mort, dat_EIP, param){
  
  a_T0 = param[1]; a_Tm = param[2]; a_c = param[3]
  mu_T0 = param[4]; mu_Tm = param[5]; mu_c = param[6]
  PDR_T0 = param[7]; PDR_Tm = param[8]; PDR_c = param[9]
  
  #bite rate
  Temp = dat$T
  a = briere(Temp, a_T0, a_Tm , a_c)
  a[a<=0] = 1e-4
  LL_a = sum(  dexp(dat_bite$bite_rate, rate = a, log = TRUE) )
  
  #mortality
  Temp = unique(dat_mort$Temp)
  lf = quad(Temp, mu_T0, mu_Tm , mu_c)
  lf[lf<0] = 1e-5
  mu = ifelse( lf<=0, 1, 1/ lf )
  LL_mu = NULL
  for (i in 1:length(Temp)){
    
    dat_sub  = filter(dat_mort, Temp == Temp[i])
    
    LL_mu = rbind( LL_mu, sum(  dbinom(dat_sub$Alive, size = dat_sub$Dead+dat_sub$Alive, prob = lf[i] , log = TRUE) , na.rm = TRUE) )
  }
  
  #PDR
  Temp = dat_EIP$T
  PDR = briere(Temp, T0, Tm , c)
  PDR[PDR<=0] = 1e-4
  LL_PDR = sum(  dnorm(dat_EIP$PDR, mean = PDR, sd = 0.08, log = TRUE) )
  
  # put together
  LL = LL_a + sum(LL_mu, na.rm = TRUE) + LL_PDR
  
  return(LL)
}