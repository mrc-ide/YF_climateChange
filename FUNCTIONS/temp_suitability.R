###########################################################
### FUNCTIONs to for relationshipos in mordecai ###
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

temp_suitability_mordecai = function(Temp){
  
  #using mean values of mordecai
  
  #bite rate
  a = briere(Temp, T0=13.35, Tm=40.08, c=2.02e-4) 
  
  # no eggs produced in a day
  EFD = briere(Temp, T0=14.58, Tm=34.61, c=8.56e-3) 
  
  # p_EA = prob survive from egg to adult
  p_EA = quad(Temp, T0=13.56, Tm=38.29, c=-5.99e-3)
  
  # MDR = mosquito development rate
  MDR = briere(Temp, T0=11.36, Tm=39.17, c=7.86e-5)
  
  # mu = m mortality
  mu = 1/ ifelse( quad(Temp, T0=9.16, Tm=37.73, c=-1.48e-1)<=0, 1, quad(Temp, T0=9.16, Tm=37.73, c=-1.48e-1) ) #guard against negatives and zeros
  
  # b = propn infectious bites m->h
  b = briere(Temp, T0=17.05, Tm=35.83, c=8.49e-4)
  
  # c = propn infectious bites h->m
  c = briere(Temp, T0=12.22, Tm=37.46, c=4.91e-4)
  
  # PDR = parasite development rate = 1/EIP  
  PDR = briere(Temp, T0=10.68, Tm=45.90, c=6.65e-5)
  
  a[is.na(a)] = EFD[is.na(EFD)] = p_EA[is.na(p_EA)] = MDR[is.na(MDR)] = b[is.na(b)] = c[is.na(c)] = PDR[is.na(PDR)] = 0
  
  # suitability
  Z = (a^2 * b * c * exp(-mu / PDR) * EFD * p_EA * MDR) / mu^3 
  
  return(Z)
}