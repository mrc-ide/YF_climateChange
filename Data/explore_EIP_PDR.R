
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

#data
davis_EIP <- read_csv("Z:/YF_climateChange/Data/davis_EIP.csv")


histogram(davis_EIP$PDR)

#mordecai fit
davis_EIP$fit= briere(davis_EIP$T, T0=18.3, Tm=42.3, c=0.000174) # this is for zika

ggplot(davis_EIP) + geom_point( aes(x = T, y = PDR)) + geom_line( aes(x = T, y = fit))

##################################################################################################################
#normailty
ggdensity(davis_EIP$PDR, main = "Density plot of bite rate")

ggqqplot(davis_EIP$PDR)

shapiro.test(davis_EIP$PDR) #if p value >0.05 then the data is normally distributed

