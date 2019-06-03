#### run on cluster ####

# drat:::addRepo("mrc-ide")
# install.packages("didehpc")
library(didehpc)

didehpc_config()

context::context_log_start()

root = "contexts"

packages <- list( attached = c("dplyr", "readr", "reshape", "mvtnorm", "truncdist","R.utils", "tibble", "YFestimation", "KsetupR"))

sources <- c("FUNCTIONS/GLM_tempsuit_MCMC.R", 
             "FUNCTIONS/temp_suit.R", 
             "FUNCTIONS/calc_FOI_tempchange.R", 
             "FUNCTIONS/calc_FOI_tempprecipchange.R",
             "run_estimation.R")

ctx <- context::context_save(root, packages = packages, sources = sources,
                             package_sources = provisionr::package_sources(local = c("Z:/KsetupR", "Z:/YFestimation", "Z:/Burden_package/YFburden") ))


obj <- didehpc::queue_didehpc(ctx)


# run for a number

run_id = 100:120
grp <- obj$lapply(run_id, run_estimation)


grp$status()
grp$results()
