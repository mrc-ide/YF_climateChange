library(dplyr)
library(magrittr)
library(ggplot2)

df <- read.csv("Z:/YF_climateChange/transmission_intensity_med_rainfall_temp.csv", stringsAsFactors = FALSE)

df %<>% mutate(adm0 = substr(adm0_adm1, 1, 3))




df %>%
  filter(year != "now") %>%
  group_by(year, rainfall_temp, adm0) %>%
  summarise(mean_runs = mean(runs)) %>%
  ggplot() +
  geom_col(position = "dodge") +
  aes(x = year, y = mean_runs, fill = rainfall_temp) +
  facet_wrap(adm0 ~., scales = "free_y") +
  theme_bw() +
  scale_fill_viridis_d(option = "magma", end = 0.9)+
  labs(x = "Year", y= "Mean force of infection", fill = "Rainfall or temperature varying")
