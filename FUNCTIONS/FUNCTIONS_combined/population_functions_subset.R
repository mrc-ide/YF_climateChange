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

#########################################################################################################
### FIT GLM ### EDITED TO TAKE MODEL CSV ###
#########################################################################################################
fit_glm = function(dat, depi, models) {
  
  fm_list = models
  
  model=1
  model = as.numeric(model)
  fm_best = fm_list[model]
  
  bm = glm(as.formula(fm_best), data=dat, family=binomial(link="cloglog")) #bm = best model
  
  # setting up the evaluation of the likelihood:
  beta = coefficients(bm)
  vl=NULL
  for(i in 1:ncol(dat)) if(length(grep(names(dat)[i], names(beta)))>0) vl = c(vl,i) # select the variables used in the GLM
  x = cbind(Intercept=1,dat[,vl])
  j_expand = !sapply(1:ncol(x), function(i) is.numeric(x[,i]) & !is.factor(x[,i]))
  x_num = NULL
  for(j in 1:ncol(x)) {   # create indicative variables for adm05
    if(j_expand[j]==T) {
      tab = table(1:nrow(x),x[,j])
      colnames(tab) = paste(colnames(x)[j],colnames(tab),sep="")
      x_num = cbind(x_num,tab[,-1])
    } else {
      x_num = cbind(x_num,x[,j])
      colnames(x_num)[ncol(x_num)] = colnames(x)[j]
    }
  }
  x = x_num # covariate matrix
  rm(x_num)
  y = dat[,depi]  # dependant variable vector
  
  
  beta0 = coefficients(bm)
  names(beta0)[names(beta0)=="(Intercept)"] = "Intercept"
  
  mm = match(names(beta0),colnames(x))
  x = x[,mm]
  
  beta0 = beta0[match(colnames(x),names(beta0))]
  nn = names(beta0)
  
  return(list(beta0=beta0, x=x, y=y))
  
}
##################################################################################
### fun_calcPred ### TAKEN FROM KEVIN'S CODE fun.calcPred ###
##################################################################################
fun_calcPred = function(coefs,
                        newdata,
                        type = "response",
                        varsin = NA) {
  
  if(!(type %in% c("link","response"))) stop("fun_calcPred: invalid type of predictions specified_\n")
  if(is.na(varsin[1])) varsin = 1:length(coefs)
  eta = newdata[,varsin] %*% coefs[varsin]
  if(type=="link") {
    preds = eta #X*beta
  } else if(type=="response") {
    preds = 1-exp(-exp(eta)) # q = 1 - exp (- exp(X*beta))
  }
  return(preds)
}