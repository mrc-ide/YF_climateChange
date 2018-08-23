###########################################################
### FUNCTIONs for temp suitability ###
###########################################################
briere = function(Temp, T0, Tm, c){
  out = c*Temp*(Temp-T0)*(Tm-Temp)^0.5
  return(out)
}

quad = function(Temp, T0, Tm, c){
  out = -c*(Temp-T0)*(Tm-Temp)
  return(out)
}
###########################################################

###########################################################
### FUNCTION to calculate temperature suitability index ###
###########################################################

temp_suitability = function(Temp, param){
  
  a_T0 = param[1]; a_Tm = param[2]; a_c = param[3]
  mu_T0 = param[4]; mu_Tm = param[5]; mu_c = param[6]
  PDR_T0 = param[7]; PDR_Tm = param[8]; PDR_c = param[9]
  
  #bite rate
  a = briere(Temp, T0=a_T0, Tm=a_Tm, c=a_c) # new estimate medians
  
  # mu = m mortality
  lf = quad(Temp, T0=mu_T0, Tm=mu_Tm, c=mu_c)
  mu = 1/ ifelse( lf<=0, 1, lf ) #guard against negatives and zeros
  
  
  # PDR = parasite development rate = 1/EIP  
  PDR = briere(Temp, T0=PDR_T0, Tm=PDR_Tm, c=PDR_c) #
  PDR[PDR<0]=0
  
  a[is.na(a)] = PDR[is.na(PDR)] = 0
  
  # suitability
  Z = (a^2 * exp(-mu / PDR) ) / mu 
  
  return(Z)
}