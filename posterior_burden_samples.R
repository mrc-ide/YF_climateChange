### --------------------------------------------------------------------- ###
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
library(RColorBrewer)
library(truncdist)
library(tibble)
library(montagu)
library(ggmcmc)
library(mcmcplots)
library(R.utils)
library(parallelsugar)

library(YFestimation)
library(snapalette)
library(KsetupR)
library(YFburden)

#-----------------------------------------------------------------------------
sourceDirectory("FUNCTIONS", modifiedOnly = FALSE)

#-----------------------------------------------------------------------------
transmission_proj = read.csv( "transmission_intensity_samples.csv", stringsAsFactors = FALSE)

transmission_proj = transmission_proj %>% mutate(year = as.character(year), scenario = as.character(scenario))

transmission_proj = transmission_proj[,-1]

param_samples = spread(transmission_proj, year, FOI)
param_samples[is.na(param_samples)] = 0

### interpolate for each year ###
fun_interp = function(i){
  df = param_samples[i,]
  new_df = data.frame(adm0 = df$adm0,
                      scenario = df$scenario,
                      sample = df$sample,
                      year = 2018:2070,
                      FOI = NA)
  new_df$FOI = approx(x = c(2018, 2050, 2070),
                      y = c(df$now, df$`2050`, df$`2070`),
                      xout = c(2018:2070))$y
  out = spread(new_df, year, FOI)
  return(out)
}

param_samples_interp = lapply(1:nrow(param_samples), fun_interp)
param_samples_interp = bind_rows(param_samples_interp)

#-----------------------------------------------------------------------------

montagu::montagu_server_global_default_set(
  montagu::montagu_server("production", "montagu.vaccineimpact.org"))

touchstone_version = "201710gavi-5"

### expectation id ###
exp_id = montagu_expectations("IC-Garske", touchstone_version)$id

### population ###
pop_all = montagu_demographic_data(type_code = "int_pop",
                                   touchstone_id = touchstone_version)

pop_all = pop_all %>% filter( year>=1939 & country_code %in% transmission_proj$adm0)


### fixed values ###
P_severe = 0.12
P_severeDeath = 0.47

d_acute = 17.8 / 365 #report
dw_acute = 0.172
d_conv = 28 / 365
dw_conv = 0.024

vac_eff = 0.9447782 # from median of serology estiamtion 0622

### sort out coverage ###
GAVI_vac = montagu_coverage_data("IC-Garske", touchstone_version,  "yf-preventive-gavi")

GAVI_vac_fix = filter(GAVI_vac, year< 2018)

fix_tmp = NULL
for(i in 2018:2070){
  tmp = filter(GAVI_vac, year == 2018)
  tmp$year = i
  fix_tmp = bind_rows(fix_tmp, tmp)
}

GAVI_vac_fix = bind_rows(GAVI_vac_fix, fix_tmp)


#-----------------------------------------------------------------------------

years = 2018:2070

# infections_out = NULL
# for(i in 1:nrow(param_samples_interp)){
  
fun_calc_burden = function(i){
  df = param_samples_interp[i, ]
  
  #get coverage for that country
  coverage_country_prev = filter(GAVI_vac_fix, country_code == df$adm0)
  
  #population for that country
  pop_country = filter(pop_all, country_code == df$adm0)
  
  ### want to reshape into previous format ###
  pop_new = NULL
  for (y in max(pop_country$year) : min(pop_country$year)){
    tmp = c(y, filter(pop_country, year == y)$value)
    
    if(length(tmp)< 102) {tmp = c(tmp, rep(NA, 102-length(tmp)))}
    pop_new = rbind(pop_new, tmp)
  }
  colnames(pop_new) = c("year", 0:100)
  pop_new = pop_new[order(pop_new[,1]),]
  
  
  tmp_burden = NULL
  for(y in 1:length(years)){
    
    if(y == 1){
      immunity_start = fun_immunityStart(model_type = "Foi",
                                         transmission_param = as.numeric(df[3+y]),
                                         age_max = 100,
                                         pop = pop_new,
                                         old_coverage = coverage_country_prev)
    }
    
    #calculate burden in year of interest
    
    out_prev = run_infections_unit(model_type = "Foi",
                                   transmission_param = as.numeric(df[3+y]),
                                   vac_eff,
                                   years = years[y],
                                   age_max = 100,
                                   pop = pop_new,
                                   coverage = coverage_country_prev,      #this is a subset of pop_moments_whole for adm1s of interest
                                   immunity_start)
    
    
    immunity_start = out_prev$immunity
    
    tmp_burden = bind_rows(tmp_burden, data.frame(adm0 = df$adm0,
                                                  scenario = df$scenario,
                                                  sample = df$sample,
                                                  year = years[y],
                                                  infections = sum(out_prev$new_infections)))
    
  }
  
  tmp_burden_out = spread(tmp_burden, year, infections)
  #infections_out = bind_rows(infections_out, spread(tmp_burden, year, infections))
  return(tmp_burden_out)
}

infections_out_l = mclapply(1:nrow(param_samples_interp), fun_calc_burden, nc.cores = 6)
infections_out = bind_rows(infections_out_l)
