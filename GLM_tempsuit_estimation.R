


source("run_estimation.R")

if(!parall){
  run_estimation(run_id = 1)
}

#########################################################################################################
### parallel ###
#########################################################################################################
if(parall){
  library(parallelsugar)
  
  mclapply(X = c(1:6),
           FUN = function(run_id){run_estimation(run_id)})
}
