
prepare_climate_dat = function(dat, year, scenario){
  
  dat_full = dat
  
  if(year == "now"){
  #baseline
  dat_full = dat_full %>% mutate(worldclim_temp_mid = (dat_full$X1970_1990tnnow + dat_full$X1970_1990txnow)/2)
  
  dat_full = dat_full %>% mutate(worldclim_temp_range = (dat_full$X1970_1990txnow - dat_full$X1970_1990tnnow))
  
  dat_full = dat_full %>% mutate(worldclim_rainfall = dat_full$X1970_1990prnow)
  
  } else {
    
    dat_full = dat_full %>% mutate(worldclim_temp_mid = (dat_full[, paste0("X", year, "tn", scenario)] + 
                                                                        dat_full[, paste0("X", year, "tx", scenario)])/2)
    
    dat_full = dat_full %>% mutate(worldclim_temp_range = (dat_full[, paste0("X", year, "tx", scenario)] - 
                                                                 dat_full[, paste0("X", year, "tn", scenario)]))
    
    dat_full = dat_full %>% mutate(worldclim_rainfall = dat_full[,paste0("X", year, "pr", scenario)] )
    
  }
  
  return(dat_full)
}