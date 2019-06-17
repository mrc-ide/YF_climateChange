
library(dplyr)
library(tidyr)
library(ggplot2)
library(gganimate)
library(ggridges)
library(snapalette)
library(magrittr)

snapal = "Stavanger"
#------------------------------#
# transmission samples --------#
transmission_proj = read.csv( "transmission_intensity_samples.csv", stringsAsFactors = FALSE)

transmission_proj %<>% mutate(year = as.character(year), scenario = as.character(scenario))  %>% 
  mutate(year = factor(year, levels = c("now", 2050, 2070)))

tp = transmission_proj %>% filter( !(scenario %in% c(45, 60, 85) & year == "now"))
tp$scenario[tp$year == "now"] = "now"


#------------------------------#
#get files in infections directory
fil = list.files("infections_stop")

inf_df = NULL
for(i in 1:length(fil)){
  inf_df[[i]] = mutate_all( read.csv(paste0("infections_stop/", fil[i]), stringsAsFactors = FALSE), as.character)
}

inf_df = bind_rows(inf_df)

inf_df = gather(inf_df,   "Year", "Infections", -c(adm0, scenario, sample))

inf_df$Year = gsub("X", "", inf_df$Year)

inf_df %<>% mutate(Infections = as.numeric(Infections)) %>% filter(Year %in% c(2018, 2050, 2070))

#-------------#

n_samples = length(unique(inf_df$sample))

P_severe_runs = rbeta(n_samples, 6.367309,44.60736)
P_death_runs = rbeta(n_samples, 16.43466, 18.49048)

tmp_df = data.frame(sample = unique(inf_df$sample),
                    P_severe = P_severe_runs,
                    P_death = P_death_runs)

inf_df %<>% left_join( tmp_df, by = "sample")

inf_df= unique(inf_df)

inf_df %<>% mutate(Cases = inf_df$Infections * inf_df$P_severe)
inf_df %<>% mutate(Deaths = inf_df$Cases * inf_df$P_death)

rm(P_severe_runs, P_death_runs, n_samples, tmp_df, transmission_proj)

#----------#


pop_all = read.csv("population.csv", stringsAsFactors = FALSE)

pop_all %<>% filter(year %in% 2018:2070) %>% 
  group_by(country_code, year) %>% 
  summarise(total_pop = sum(value))

pop_all %<>% dplyr::rename(adm0 = country_code, Year = year) %>%
  mutate(Year = as.character(Year))

inf_df %<>% left_join( pop_all, by = c("adm0", "Year"))

rm(pop_all)

#-------------------#
# p<- ggplot(inf_df, 
#            aes(x = scenario, y=Deaths/total_pop, fill = scenario, colour = scenario)) +
#   geom_violin(show.legend = FALSE, alpha = 0.7) +
#   scale_fill_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
#   scale_colour_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
#   ylab("Deaths per capita")+
#   scale_y_log10() +
#   theme_bw()#+
#   #facet_wrap(adm0~.)
# 
# p + transition_time(as.integer(Year) ) +
#   labs(title = "Year: {frame_time}")
# 
# 
# 
# p<- ggplot(inf_df, 
#            aes(x = scenario, y=Deaths/total_pop, fill = scenario, colour = scenario)) +
#   geom_violin(show.legend = FALSE, alpha = 0.7) +
#   scale_fill_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
#   scale_colour_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
#   ylab("Deaths per capita")+
#   #xlab("Cases per capita")+
#   #scale_y_log10() +
#   #scale_x_log10()+
#   facet_wrap()
#   theme_bw()
# 
# p + transition_time(as.integer(Year) ) +
#   labs(title = "Year: {frame_time}")

#----------------#

tmp = inf_df %>% filter( scenario == "now")
tmp %<>% dplyr::rename(Deaths_now = Deaths) %>% select("adm0", "Year", "Deaths_now", "sample")

inf_df %<>% left_join( tmp, by = c("adm0", "Year", "sample"))

rm(tmp)

inf_df %<>% mutate(relative_deaths = 100* (Deaths - Deaths_now)/ Deaths_now  ) 

inf_df %<>% unique()

inf_df %<>% mutate(death_diff = Deaths - Deaths_now ) 

inf_df %<>% unique()

#--------------------#

west = c("BEN", "BFA", "GMB", "GHA", "GIN", "GNB", "CIV", "LBR", "MLI", "MRT", "NER", "NGA", "SEN", "SLE", "TGO")

central = c("TCD", "CAF", "CMR", "GAB", "COD", "COG", "AGO", "SSD", "SDN")

inf_df %<>% mutate(WE = ifelse(adm0 %in% west, "West", 
                               ifelse(adm0 %in% central, "Central and Sudan",
                                      "East")))


# p<- ggplot(inf_df, 
#            aes(x = scenario, y=Deaths/total_pop, fill = scenario, colour = scenario)) +
#   geom_violin(show.legend = FALSE, alpha = 0.7) +
#   scale_fill_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
#   scale_colour_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
#   ylab("Deaths per capita")+
#   #xlab("Cases per capita")+
#   scale_y_log10() +
#   #scale_x_log10()+
#   facet_wrap(WE~.) +
#   theme_bw()
# 
# p + transition_time(as.integer(Year) ) +
#   labs(title = "Year: {frame_time}")




p<- ggplot(filter(inf_df, scenario != "now"), 
           aes(x = scenario, y=relative_deaths, fill = scenario, colour = scenario)) +
  geom_jitter(show.legend = FALSE, alpha = 0.7, size = 2) +
  scale_fill_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
  scale_colour_manual(values = c(snapalette(snapal)[c(1:4)], "black"))+
  ylab("Percentage change in deaths")+
  facet_wrap(WE~.) +
  theme_bw() +
  theme(text = element_text(size = 20))

ap = p + transition_time(as.integer(Year) ) +
  labs(title = "Year: {frame_time}") + 
  shadow_trail()
  
animate(ap, height = 800, width = 1000)

anim_save("Percentage_change_in_deaths.gif")
