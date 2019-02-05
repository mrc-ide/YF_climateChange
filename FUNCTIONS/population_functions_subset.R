#########################################################################################################
### LAUNCH ENVIRONMENTAL DATA ### STRAIGHT FROM ORIGINAL CODE ###
#########################################################################################################

launch_env_dat = function(dat_full, c34, delete_surv_AGO = TRUE) {
  
  depvar = "cas.or.out" # which presence/abence outcome to consider. Alternatives : "cases" # "outbreaks" # 
  depi = match(depvar,names(dat_full)) # column number for the chosen outcome
  
  # excluding LC1,3,5,15 as they only ever cover <5% - doesn't make sense to use these...
  ex_i = match(c("LC1","LC3","LC5","LC15"),names(dat_full))
  dat_full = dat_full[,-ex_i]
  
  # adding a categorical "dominant land cover" variable:
  LC_i = grep("LC",names(dat_full))
  LC_dom = names(dat_full)[LC_i][apply(dat_full[,LC_i],1,which.max)]
  dat_full = cbind(dat_full, LC_dom) # 37 potential covariates/levels
  
  
  ###############
  dat_full$surv.qual.adm0[is.na(dat_full$surv.qual.adm0)] = 0
  dat_full$surv.qual.adm1[is.na(dat_full$surv.qual.adm1)] = 0
  
  if(delete_surv_AGO==T){# setting the surv.qual.adm0 for AGO to 0 (very few cases reported): 
    dat_full$surv.qual.adm0[dat_full$adm0=="AGO"] = 0 
  }
  
  #adding log(surv_qual_adm0)
  dat_full = cbind(dat_full, log.surv.qual.adm0 = log(dat_full$surv.qual.adm0)) # log_surv better fits data than surv
  dat_full$log.surv.qual.adm0[is.infinite(dat_full$log.surv.qual.adm0 )] = 0  # A way not to apply beta(log_surveillance) to coutries uotside of YFSD
  
  adm05 = dat_full$adm0
  adm05[dat_full$surv.qual.adm0>0] = "AFR" # a same categorical variable for all countries within the YFSD
  dat_full = cbind(dat_full, adm05=adm05) # for countries outside YFSD, 1 categorical variable per country
  
  dat = dat_full[dat_full$adm0 %in% c34,]
  dat$adm05 = as.factor(as.character(dat$adm05))
  
  
  v1 = apply(dat,2,var)
  for(i in 8:(ncol(dat))) {
    if(!is.factor(dat[,i]) & !is.character(dat[,i])) {
      dat[,i] = dat[,i]/sqrt(v1[i])
      dat_full[,i] = dat_full[,i]/sqrt(v1[i])  
    } 
    # Explanation:
    # If we fit the model on dat and if we want to project estimates on dat_full, 
    # variabble from dat_full need to be expressed on the same scale than those from dat , thus we normalize dat_full relatively to dat
    
  }
  
  dat = dat[,names(dat)!="surv_qual_adm0"]
  depi = match(depvar,names(dat))
  
  
  # I do that because pop is ordered that way and we need to match both pop and dat
  dat = dat[ order(dat$adm0_adm1), ]
  
  return(list(depi = depi, dat = dat))
} 

