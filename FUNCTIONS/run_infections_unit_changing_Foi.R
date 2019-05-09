

run_infections_unit_changing_FOI = function(model_type = "Foi",
                                           transmission_param,
                                           vac_eff = 1,
                                           years_in,
                                           age_max,
                                           pop,
                                           coverage,
                                           immunityStart) {
  
  
  years = years_in
  ages = c(0:age_max)
  n_years = length(years)
  
  ### get rid of NA in pop ###
  pop[is.na(pop)] = 0
  ####
  
  
  ### declare set up
  immunity = new_infections = matrix(NA, nrow = n_years, ncol = age_max + 1)
  
  
  immunity = immunityStart  # immunity profile taking into account previous vaccination and infection
  
  
  # iterate over the years
  immunity = update_immunity(immunity)
  
  
  ## generate new infections
  
  tmp = switch(model_type,
               "Foi" = generate_infections_static(transmission_param,
                                                  pop[which( pop[,1] %in% years), 2:102],
                                                  immunity),
               "R0" = generate_infections_R0(transmission_param,
                                             pop[which( pop[,1] %in% years), 2:102],
                                             immunity))
  
  new_infections = tmp$new_infections
  immunity = tmp$immunity
  
  
  
  ## add vaccination at the end of the year: coverage_df: year, age, coverage.
  if (nrow(coverage) > 0)
    for (y in 1:nrow(coverage)) {
      if (coverage$year[y] == years) {
        
        immunity = add_vaccination(coverage$coverage[y],
                                                vac_eff = vac_eff,
                                                age_first = coverage$age_first[y],
                                                age_last = coverage$age_last[y],
                                                immunity,
                                                skew = 0)#coverage$skew[y])
      }
    }
  
  
  
  ### only want immunity for years_in
  immunity = immunity
  new_infections = new_infections
  
  return(list(immunity = immunity, new_infections = new_infections))
}
