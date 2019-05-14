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
transmission_proj = read.csv( "transmission_intensity_med.csv", stringsAsFactors = FALSE)

transmission_proj = transmission_proj %>% mutate(adm0 = substr(adm0_adm1, 1,3))

transmission_proj1 = transmission_proj %>% group_by(adm0, year, scenario) %>% summarise(median_FOI = mean(runs))

transmission_proj1 = transmission_proj1 %>% ungroup() %>% mutate(year = as.character(year), scenario = as.character(scenario))

param_samples = filter(transmission_proj1, !(scenario %in% c(45, 60, 85) & year == "now"))
param_samples$scenario[param_samples$year == "now"] = "now"

param_samples = spread(param_samples, year, median_FOI) # spread and split
param_samples_now = filter(param_samples, scenario == "now")
param_samples = filter(param_samples, scenario != "now")

param_samples$now = filter(transmission_proj1, year == "now")$median_FOI

param_samples_now$`2050` = param_samples_now$`2070` = param_samples_now$now

param_samples = bind_rows(param_samples, param_samples_now)


### interpolate for each year ###
fun_interp = function(i){
  df = param_samples[i,]
  new_df = data.frame(adm0 = df$adm0,
                      scenario = df$scenario,
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

write.csv(param_samples_interp, "transmission_intensity_median_interp.csv", row.names = FALSE)


#param_samples_interp = read.csv( "transmission_intensity_median_interp.csv", stringsAsFactors  = FALSE)

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

write.csv(pop_all, "population.csv", row.names = FALSE)



vac_eff = 0.9447782 # from median of serology estiamtion 0622

### sort out coverage ###
GAVI_vac = montagu_coverage_data("IC-Garske", touchstone_version,  "yf-preventive-gavi")

GAVI_vac_fix = filter(GAVI_vac, year< 2019)

fix_tmp = NULL
for(i in 2019:2070){
  tmp = filter(GAVI_vac, year == 2018 & activity_type == "routine")
  tmp$year = i
  fix_tmp = bind_rows(fix_tmp, tmp)
}

GAVI_vac_fix = bind_rows(GAVI_vac_fix, fix_tmp)


#-----------------------------------------------------------------------------

years = 2018:2070

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
  pop_out = NULL
  for(y in 1:length(years)){
    
    if(y == 1){
      immunity_start = fun_immunityStart(model_type = "Foi",
                                         transmission_param = as.numeric(df[2+y]),
                                         age_max = 100,
                                         pop = pop_new,
                                         old_coverage = coverage_country_prev)
      
      out_prev = run_infections_unit(model_type = "Foi",
                                     transmission_param = as.numeric(df[2+y]),
                                     vac_eff,
                                     years_in = years[y],
                                     age_max = 100,
                                     pop = pop_new,
                                     coverage = coverage_country_prev,      #this is a subset of pop_moments_whole for adm1s of interest
                                     immunityStart = immunity_start)
    } else{
      
      #calculate burden in year of interest
      
      out_prev = run_infections_unit_changing_FOI(model_type = "Foi",
                                                  transmission_param = as.numeric(df[2+y]),
                                                  vac_eff,
                                                  years_in = years[y],
                                                  age_max = 100,
                                                  pop = pop_new,
                                                  coverage = coverage_country_prev,      #this is a subset of pop_moments_whole for adm1s of interest
                                                  immunityStart = immunity_start)
    }
    
    immunity_start = out_prev$immunity
    
    tmp_burden = bind_rows(tmp_burden, data.frame(adm0 = df$adm0,
                                                  scenario = df$scenario,
                                                  year = years[y],
                                                  infections = sum(out_prev$new_infections) ))
    
  }
  
  tmp_burden_out = spread(tmp_burden, year, infections)
  return(tmp_burden_out)
}


ind = 1:nrow(param_samples_interp)
infections_out_l = mclapply(ind, fun_calc_burden, mc.cores = 2)
infections_out = bind_rows(infections_out_l)


write.csv(infections_out, paste0("infections_per_scenario_year_country_median.csv"), row.names = FALSE)
