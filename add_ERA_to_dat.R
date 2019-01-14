

###script to add the mean temp per year to dat###

dat_old = read.csv(paste0("../Data/","Environment/Africa_adm1_dat_2018_gadm36.csv"), stringsAsFactors = FALSE) 

ERA_means = read.csv(paste0("../Data/","Environment/ERA/mean_annual_admin.csv"), stringsAsFactors = FALSE)

#get ERA means in the same format/ order

library(dplyr)

ERA_means$ID_1 = gsub("_1", "", ERA_means$ID_1)

ERA_means$ID_1 = gsub("\\.", "_", ERA_means$ID_1)

ERA_means = arrange(ERA_means, ID_1)

dat_old = arrange(dat_old, adm0_adm1)

#only want the era means that are in dat_old

ERA_means = ERA_means %>% filter(ID_1 %in% dat_old$adm0_adm1)

#combine

if(all.equal(dat_old$adm0_adm1, ERA_means$ID_1)){
  dat_new = bind_cols(dat_old, ERA_means[, c(4:ncol(ERA_means))])
}

#give the year columns sensible names
names(dat_new)[grep("^X", names(dat_new))] = gsub("X", "ERA_mean_", names(dat_new)[grep("^X", names(dat_new))])


### test and compare ###
plot(dat_new[, "ERAday.mean"], type = "l", col = "red")
for (i in 1:18){
  lines(dat_new[, 50+i], col = "green")
}
lines(dat_new[, "ERAnight.mean"],  col = "blue")


write.csv(dat_new, paste0("../Data/","Environment/Africa_adm1_dat_2018_gadm36_ERAannual.csv"), row.names = FALSE)