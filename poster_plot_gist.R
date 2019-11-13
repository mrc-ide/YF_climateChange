
library(ggplot2)
library(magrittr)

df = data.frame(Year = c(rep(2050, 4), rep(2070, 4)),
                Scenario = rep(c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5"), 2),
                med = c(10.8, 16.7, 15.5, 24.9, 10.0, 19.3, 21.3, 39.6),
                low = c(-2.4, -2.4, -2.8, -2.2, -0.7, -2.8, -4.6, -2.9),
                high = c(37.9, 57.4, 51.8, 88.3, 34.1, 71.1, 77.7, 178.6))

df %<>% mutate(Year = as.character(Year))


dodge <- position_dodge(width=0.9)

ggplot(df)+
  geom_pointrange( aes(x = Scenario, ymin = low, ymax = high, y = med, colour = Year),
                   position = dodge,
                   size = 2) +
  theme_bw()+
  ylab("Percentage change in deaths")+
  scale_color_manual(values = viridis::magma(5)[c(2,3)]) +
  theme(text = element_text(size = 20))
