# function file to work with model_reestimation

#########################################################################################################
### LAUNCH ENVIRONMENTAL DATA ### STRAIGHT FROM ORIGINAL CODE ###
#########################################################################################################

launch_env_dat = function(env_table, c34, delete_surv_AGO = T) {
  
  
  dat_full = read.csv(env_table, stringsAsFactors=F)
  
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
  #table(dat$adm05)
  
  #summary(dat)
  #str(dat)
  v1 = apply(dat,2,var)
  for(i in 8:(ncol(dat))) {
    if(!is.factor(dat[,i]) & !is.character(dat[,i])) {
      dat[,i] = dat[,i]/sqrt(v1[i])
      dat_full[,i] = dat_full[,i]/sqrt(v1[i])  
    } 
    # Explanation:
    # If we fit the model on dat and if we want to project estimates on dat_full, variabble from dat_full need to be expressed on the same scale than those from dat , thus we normalize dat_full relatively to dat
    
  }
  
  # dat = dat[,names(dat)!="adm0"] ligne enlevee au 04/09/15, a remettre si pb
  dat = dat[,names(dat)!="surv_qual_adm0"]
  depi = match(depvar,names(dat))
  
  
  dim(dat)
  class(dat$adm0_adm1)
  
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

#########################################################################################################
### IMPORT POPULATION DATA ### EDITED TO ONLY TAKE LS2014  ###
#########################################################################################################
import_pop_data_LS2014 = function(path, c_country, adm) {
  
  pop_year_age = NULL  
  
  for(adm0 in c_country){
    tmp = read.csv(paste0(path,"Population", "/from_LandScan2014/", adm, "/pop_year_age_", adm, "_", adm0, "_LandScan2014.csv"),h=T)
    pop_year_age = rbind(pop_year_age, tmp)
  }
  
  pop_year_age[,adm] = paste(pop_year_age[,"adm0"], pop_year_age[,adm], sep="_")
  colnames(pop_year_age)[colnames(pop_year_age)==adm] = paste("adm0", adm, sep = "_")
  return(pop_year_age)
}

#########################################################################################################
### IMPORT POPULATION DATA ### EDITED TO ONLY TAKE LS2017 at adm1 ###
#########################################################################################################
import_pop_data_LS2017 = function(path, c_country, adm) {
  
  pop_year_age = read.csv(paste0(path,"Population", "/from_LandScan2014/", "population_by_adm1_year_age_1950-2100_LandScan2014_gadm2_countries.csv"),h=T)
  #pop_year_age = read.csv(paste0(path,"Population", "/from_LandScan2017/", "population_by_adm1_year_age_1950-2100_LandScan2017_gadm2_8_countries.csv"),h=T)
  pop_year_age = pop_year_age[pop_year_age$adm0 %in% c_country,]
  
  pop_year_age[,adm] = paste(pop_year_age[,"adm0"], pop_year_age[,adm], sep="_")
  colnames(pop_year_age)[colnames(pop_year_age)==adm] = paste("adm0", adm, sep = "_")
  return(pop_year_age)
}
#########################################################################################################
### add_1940_1950 ### TAKEN FROM KEVIN'S CODE ###
#########################################################################################################
add_1940_1950 = function(pop2d){
  if(min(pop2d$year) != 1950)  stop("pop2d should start in 1950")
  
  for(y in 1949:1940) {
    pop1_early = pop2d[pop2d$year==y+1,]
    pop1_early$year = y
    pop2d = rbind(pop1_early,pop2d)
  }
  
  return(pop2d)
}

#########################################################################################################
### Make_pop3d_Ptot2d_Pprop3d ### TAKEN FROM KEVIN'S CODE transform_into_pop3d, get_P_tot_2d###
#########################################################################################################
Make_pop3d_Ptot2d_Pprop3d = function(pop2d, adm){
  
  adm0_adm = paste0("adm0_", adm)
  
  #pop3d =3d array- dim [adm, year, age]
  pop3d=rep(NA, nrow(pop2d)*(ncol(pop2d)-3))
  dim(pop3d)=c(length(table(pop2d[,adm0_adm])), length(table(pop2d$year)), ncol(pop2d)-3)
  
  dim_adm  = names(table(pop2d[,adm0_adm]))
  dim_year = as.numeric(names(table(pop2d$year))) 
  dim_age  = names(pop2d)[4:length(names(pop2d))] 
  
  length_adm=length(dim_adm);   length_year=length(dim_year);   length_age=length(dim_age)   #numbers of elements in each dimension
  
  for(ageIndex in 1:length_age) { 
    for(yearIndex in min(dim_year):max(dim_year)) { 
      mm = match(pop2d[,adm0_adm][pop2d$year==yearIndex],dim_adm) 
      pop3d[mm,yearIndex-min(dim_year)+1,ageIndex]=pop2d[pop2d$year==yearIndex,ageIndex+3]
    }
  }
  
  
  #P_tot2d = sum(pop3d(adm,year))
  P_tot_2d = rep(NA, length_adm*length_year);   dim(P_tot_2d)=c(length_adm,length_year )
  
  #p_prop_3d = pop3d/P_tot_2d
  p_prop_3d=rep(NA, length_adm*length_year*length_age);   dim(p_prop_3d)=dim(pop3d)
  
  for(admIndex in 1:length_adm){
    for (yearIndex in 1:length_year){
      
      P_tot_2d[admIndex,yearIndex] = sum(pop3d[admIndex,yearIndex,], na.rm=T)
      
      for (ageIndex in 1:length_age){
        p_prop_3d[admIndex,yearIndex,ageIndex] = pop3d[admIndex, yearIndex, ageIndex]/P_tot_2d[admIndex, yearIndex]
      }
    }
  }
  
  #names
  dimnames(p_prop_3d)[[1]] = dimnames(pop3d)[[1]] = rownames(P_tot_2d) = dim_adm
  dimnames(p_prop_3d)[[2]] = dimnames(pop3d)[[2]] = colnames(P_tot_2d) = dim_year
  dimnames(p_prop_3d)[[3]] = dimnames(pop3d)[[3]] = dim_age
  
  return(list(pop3d=pop3d, P_tot_2d=P_tot_2d,p_prop_3d=p_prop_3d))
  
}


#########################################################################################################
### GET POPULATION 3D ### EDITED TO ONLY TAKE LS2017 and ADM1 ###
#########################################################################################################
get_pop_data_3d = function(path, c_country, dat = dat) {
  
  adm="adm1"
  adm0_adm = paste0("adm0_", adm)
  
  pop1 = import_pop_data_LS2017(path=path, c_country=c_country, adm=adm)
  
  
  #restrict to lines in dat THIS IS CHANGED FOR 2017!
  pop1 = pop1[pop1[,adm0_adm] %in% dat[,adm0_adm],]
  
  #make sure data starts from 1940
  pop1 = add_1940_1950(pop1)
  
  #pop3d= 3d array- dim [adm, year, age]
  #P_tot_2d = sum(pop3d(adm,year))
  #p_prop_3d = pop3d/P_tot_2d
  out = Make_pop3d_Ptot2d_Pprop3d(pop1, adm)
  
  return(list(pop1= pop1, pop3d= out$pop3d, P_tot_2d = out$P_tot_2d, p_prop_3d = out$p_prop_3d))
}

#########################################################################################################
### transform_into_vc3d ### TAKEN FROM KEVIN'S CODE and simplified ###
#########################################################################################################
transform_into_vc3d = function(vc2d, adm){
  
  adm0_adm = paste0("adm0_", adm)
  
  vc3d=rep(NA, nrow(vc2d)*(ncol(vc2d)-3))
  dim(vc3d)=c(length(table(vc2d[,adm0_adm])), length(table(vc2d$year)), ncol(vc2d)-3)
  
  dim_adm = names(table(vc2d[,adm0_adm]))
  dim_year = as.numeric(names(table(vc2d$year))) 
  dim_age = names(vc2d)[4:length(names(vc2d))] 
  
  for(ageIndex in 1:length(dim_age)) { 
    for(yearIndex in min(dim_year):max(dim_year)) { 
      mm = match(vc2d[,adm0_adm][vc2d$year==yearIndex],dim_adm) 
      vc3d[mm,yearIndex-min(dim_year)+1,ageIndex] = vc2d[vc2d$year==yearIndex,ageIndex+3]
    }
  }
  
  dimnames(vc3d)[[1]] = dim_adm
  dimnames(vc3d)[[2]] = dim_year
  dimnames(vc3d)[[3]] = dim_age
  
  return(vc3d)
  
}

#########################################################################################################
### calc_incidence_vac ### TAKEN FROM KEVIN'S CODE and checked to generalise to aggregated ###
#########################################################################################################
calc_incidence_vac_general = function(vc3d){
  inc_v3d=rep(NA, dim(vc3d)[1]*dim(vc3d)[2]*dim(vc3d)[3] )
  dim(inc_v3d)= dim(vc3d)
  dimnames(inc_v3d)=dimnames(vc3d)
  
  for (admIndex in 1:dim(vc3d)[1] ) { #
    inc_v3d[admIndex,1,]= 0 # no vaccination in 1940
    
    for(year in 2:(dim(vc3d)[2]-1)){
      inc_v3d[admIndex,year,1] = vc3d[admIndex,year+1,1] # for a0, incidence = coverage in the a1 age group the next year
      
      for(age in 2:(dim(vc3d)[3]-1)) { 
        if(!is.na(vc3d[admIndex,year,age]) & vc3d[admIndex,year,age]==1){ #if vc is already =1, incidence =0
          inc_v3d[admIndex,year,age]=0
        } else {                             
          inc_v3d[admIndex,year,age]=  ( vc3d[admIndex,year+1,age+1] - vc3d[admIndex,year,age] ) / (1 - vc3d[admIndex,year,age] ) # among those suscpetible at y-1
          inc_v3d[admIndex,year,age]=ifelse(inc_v3d[admIndex,year,age]<10e-15,0,inc_v3d[admIndex,year,age])
        }# with rounding, some values become negative
      }
      inc_v3d[admIndex,year,dim(vc3d)[3]]=0 # incidence vaccination = 0 for age=100
    }
    
    inc_v3d[admIndex,dim(vc3d)[2],]=0# incidence vaccination = 0 for year=2050
    inc_v3d[admIndex,dim(vc3d)[2],1] = vc3d[admIndex,dim(vc3d)[2],1] # except among new borns
  }
  
  return(inc_v3d)
}

#########################################################################################################
### calc_t0_vac_Africa ### TAKEN FROM KEVIN'S CODE and simplified ###
#########################################################################################################
calc_t0_vac_africa = function(vc3d){
  t0_vac_africa = rep(NA,dim(vc3d)[1])
  for (i in 1:dim(vc3d)[1]){
    year_i = 1940
    sum_vac=0
    while(sum_vac == 0 & year_i<=2050) {
      sum_vac = sum( vc3d[i, year_i-1940+1,], na.rm=T)
      year_i = year_i + 1
    }
    t0_vac_africa[i]=year_i-2
  }
  names(t0_vac_africa)=dimnames(vc3d)[[1]]
  return(t0_vac_africa)
}

#########################################################################################################
### calc_pop_moments ### TAKEN FROM KEVIN'S CODE and generalized ###
#########################################################################################################

calc_pop_moments = function(p_prop,
                            t0_vac_,
                            dim_adm,
                            dim_year,
                            dim_age) {
  
  ageVec=c(0:100)# a matrix of age with the same dimension that pop_agg for 1 fixed year
  
  minYear=1983 #min year of data
  maxYear=2017 #last year of data
  maxDuration=maxYear-minYear + 1 #duration of data for second for loop
  
  length_adm=length(dim_adm)
  
  pop_moments_whole = rep(NA, length_adm*6)
  dim(pop_moments_whole)=c(length_adm,6)
  
  for (admIndex in 1:length_adm){
    if (t0_vac_[admIndex]<maxYear){
      
      
      vec = which(dim_year==1940):which(dim_year==t0_vac_[admIndex])
      
      pop_moments_whole_tmp = rep(NA, length(vec)*6)
      dim(pop_moments_whole_tmp) = c(length(vec), 6)
      
      for(yearIndex in vec){
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 1] = sum(p_prop[admIndex, yearIndex,], na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 2] = sum(p_prop[admIndex, yearIndex,]*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 3] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 4] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 5] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 6] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec*ageVec*ageVec*ageVec, na.rm=T)
      }
      
    } else if (t0_vac_[admIndex]>=maxYear){
      
      pop_moments_whole_tmp = rep(NA, maxDuration*6)
      dim(pop_moments_whole_tmp) = c(maxDuration, 6)
      vec = which(dim_year==minYear):which(dim_year==maxYear)
      
      for(yearIndex in vec){
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 1] = sum(p_prop[admIndex, yearIndex,], na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 2] = sum(p_prop[admIndex, yearIndex,]*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 3] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 4] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 5] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_whole_tmp[yearIndex-min(vec)+1, 6] = sum(p_prop[admIndex, yearIndex,]*ageVec*ageVec*ageVec*ageVec*ageVec, na.rm=T)
      }
    }
    pop_moments_whole[admIndex,] = apply(pop_moments_whole_tmp,2,mean)
  }
  return(pop_moments_whole)
}
#########################################################################################################
### calc_pop_moments_agg ### TAKEN FROM KEVIN'S CODE and simplified ###
#########################################################################################################

calc_pop_moments_agg = function(pop_agg3d, t0_vac_, dim_year, study_years) {
  
  ageVec=c(0:100)# a matrix of age with the same dimension that pop_agg for 1 fixed year
  
  minYear=1983 #min year of data
  maxYear=2017 #last year of data
  maxDuration=maxYear-minYear + 1 #duration of data for second for loop
  n_serosurveys=length(study_years)
  
  pop_moments_agg = rep(NA, n_serosurveys*6)
  dim(pop_moments_agg)=c(n_serosurveys,6)
  
  for (index_survey in 1:n_serosurveys){
    
    if (t0_vac_[index_survey]<study_years[[index_survey]]){
      
      vec = which(dim_year==1940):which(dim_year==t0_vac_[index_survey])
      pop_moments_agg_tmp = rep(NA, length(vec)*6)
      dim(pop_moments_agg_tmp) = c(length(vec), 6)
      
      for(year in vec){
        p_prop_agg = pop_agg3d[index_survey, year,]/(sum(pop_agg3d[index_survey, year,], na.rm=T))
        pop_moments_agg_tmp[year, 1] = sum( p_prop_agg, na.rm=T)
        pop_moments_agg_tmp[year, 2] = sum(p_prop_agg*ageVec, na.rm=T)
        pop_moments_agg_tmp[year, 3] = sum(p_prop_agg*ageVec*ageVec, na.rm=T)
        pop_moments_agg_tmp[year, 4] = sum(p_prop_agg*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_agg_tmp[year, 5] = sum(p_prop_agg*ageVec*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_agg_tmp[year, 6] = sum(p_prop_agg*ageVec*ageVec*ageVec*ageVec*ageVec, na.rm=T)
      }
      
    } else if (t0_vac_[index_survey]>=maxYear){
      
      vec = which(dim_year==minYear):which(dim_year==maxYear)
      pop_moments_agg_tmp = rep(NA, length(vec)*6)
      dim(pop_moments_agg_tmp) = c(length(vec), 6)
      
      for(year in vec){
        p_prop_agg = pop_agg3d[index_survey, year,]/(sum(pop_agg3d[index_survey, year,], na.rm=T))
        pop_moments_agg_tmp[year-min(vec)+1, 1] = sum( p_prop_agg, na.rm=T)
        pop_moments_agg_tmp[year-min(vec)+1, 2] = sum(p_prop_agg*ageVec, na.rm=T)
        pop_moments_agg_tmp[year-min(vec)+1, 3] = sum(p_prop_agg*ageVec*ageVec, na.rm=T)
        pop_moments_agg_tmp[year-min(vec)+1, 4] = sum(p_prop_agg*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_agg_tmp[year-min(vec)+1, 5] = sum(p_prop_agg*ageVec*ageVec*ageVec*ageVec, na.rm=T)
        pop_moments_agg_tmp[year-min(vec)+1, 6] = sum(p_prop_agg*ageVec*ageVec*ageVec*ageVec*ageVec, na.rm=T)
      }
    }
    pop_moments_agg[index_survey,] = apply(pop_moments_agg_tmp,2,mean)
  }  
  return(pop_moments_agg)
}

#########################################################################################################
### Make_aggregate_pop_vc ### TAKEN FROM KEVIN'S CODE create_pop.agg_vc.agg ###
#########################################################################################################
## average population structure and vaccination coverage among the same serosurvey, 
# 1 averaged value per year until the year of the survey
Make_aggregate_pop_vc = function(pop1, vc2d, sero_studies,adm1s) {
  
  no_sero_surveys=length(sero_studies)
  
  vc_tmp = pop_tmp = pop_agg = vc_agg = NULL
  
  for(yearIndex in 1940:2100) { # if needed replace 2050 by max(unlist(study_years)) 
    for(surveyIndex in 1:no_sero_surveys) {
      
      pop_tmp = pop1[pop1$year == yearIndex & pop1$adm0_adm1 %in% adm1s[[surveyIndex]],]
      vc_tmp = vc2d[vc2d$year == yearIndex & vc2d$adm0_adm1 %in% adm1s[[surveyIndex]],]
      if (nrow(pop_tmp)==1) {
        pop_tmp$adm0 = sero_studies[surveyIndex]
        vc_tmp$adm0 = sero_studies[surveyIndex]
      } else
        if(nrow(pop_tmp)>1) { # need to do some summing/averaging...
          pop_tmp[is.na(vc_tmp)] = NA
          vc_tmp = colSums(pop_tmp[,-(1:3)]*vc_tmp[,-(1:3)],na.rm=T)/colSums(pop_tmp[,-(1:3)],na.rm=T) 
          dim(vc_tmp) = c(1,length(vc_tmp))
          colnames(vc_tmp) = paste0("a",0:100)
          vc_tmp = data.frame(adm0 = sero_studies[surveyIndex], adm0_adm1 = NA, year = yearIndex, vc_tmp)
          
          pop_tmp = colSums(pop_tmp[,-(1:3)])
          dim(pop_tmp) = c(1,length(pop_tmp))
          colnames(pop_tmp) = paste0("a",0:100)
          pop_tmp = data.frame(adm0 = sero_studies[surveyIndex], adm0_adm1 = NA, year = yearIndex,pop_tmp)
        }
      pop_agg = rbind(pop_agg,pop_tmp)
      vc_agg = rbind(vc_agg,vc_tmp)
      
    }
  }
  
  pop_agg[,-(1:3)][is.na(pop_agg[,-(1:3)])] = 0
  vc_agg[,-(1:3)][is.na(vc_agg[,-(1:3)])] = 0
  
  return(list(pop_agg=pop_agg,vc_agg=vc_agg))
}

#########################################################################################################
### Make_aggregate_pop_vc_3d ### TAKEN FROM KEVIN'S CODE create_pop.agg_vc.agg ###
#########################################################################################################

Make_aggregate_pop_vc_3d = function(pop1, vc2d, sero_studies,adm1s){
  
  agg=Make_aggregate_pop_vc(pop1, vc2d, sero_studies,adm1s)
  
  pop_agg=agg$pop_agg
  vc_agg=agg$vc_agg
  
  ## pass in 3d
  pop_agg3d=rep(NA, nrow(pop_agg)*(ncol(pop_agg)-3))
  dim(pop_agg3d)=c(length(sero_studies), length(table(pop_agg$year)), ncol(pop_agg)-3)
  vc_agg3d=rep(NA, nrow(vc_agg)*(ncol(vc_agg)-3))
  dim(vc_agg3d)=c(length(table(vc_agg$adm0)), length(table(vc_agg$year)), ncol(vc_agg)-3)
  
  dim_survey = sero_studies # no. surveys
  dim_year = as.numeric(names(table(pop_agg$year))) 
  dim_age = names(pop1)[4:length(names(pop_agg))] 
  
  for(ageIndex in 1:length(dim_age)) { # dim_age 
    for(yearIndex in min(dim_year):max(dim_year)) { 
      mm = match(pop_agg$adm0[pop_agg$year==yearIndex],dim_survey) 
      pop_agg3d[mm,yearIndex-min(dim_year)+1,ageIndex]=pop_agg[pop_agg$year==yearIndex,ageIndex+3]
      vc_agg3d[mm,yearIndex-min(dim_year)+1,ageIndex]=vc_agg[vc_agg$year==yearIndex,ageIndex+3]
    }
  }
  dim(pop_agg3d)
  dimnames(pop_agg3d)[[1]] = dimnames(vc_agg3d)[[1]] = dim_survey
  dimnames(pop_agg3d)[[2]] = dimnames(vc_agg3d)[[2]] = dim_year
  dimnames(pop_agg3d)[[3]] = dimnames(vc_agg3d)[[3]] = dim_age
  
  
  return(list(pop_agg3d = pop_agg3d, vc_agg3d=vc_agg3d))
  
}

###
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