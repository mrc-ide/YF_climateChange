
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

#data
mordecai_biterate <- read_csv("Z:/YF_climateChange/Data/mordecai_biterate.csv")
hamlet_biterate = read_csv("Z:/YF_climateChange/Data/hamlet_biterate.csv")

mordecai_biterate$author = "mordecai"
hamlet_biterate$author = "hamlet"

names(hamlet_biterate) = names(mordecai_biterate)

biterate = rbind(mordecai_biterate, hamlet_biterate)


histogram(biterate$bite_rate)

#mordecai fit
biterate$fit = briere(biterate$T, T0=13.35, Tm=40.08, c=2.02e-4)

ggplot(biterate) + geom_point( aes(x = T, y = bite_rate, colour = author)) + geom_line( aes(x = T, y = fit))

##################################################################################################################
#normailty
ggdensity(biterate$bite_rate, main = "Density plot of bite rate")

ggqqplot(biterate$bite_rate)

shapiro.test(biterate$bite_rate) #if p value >0.05 then the data is normally distributed

## BOTH DATA TOGETHER ARE not NORMAL ##

shapiro.test(mordecai_biterate$bite_rate)

shapiro.test(hamlet_biterate$bite_rate)

## SEPARATELY THEY ARE ##

ggqqplot(biterate, x = "bite_rate", color = "author", palette = c("#00AFBB", "#E7B800"))
##################################################################################################################
# exponential

fit1 = MASS::fitdistr(biterate$bite_rate, "exponential")
ks.test(biterate$bite_rate, "pexp", fit1$estimate)

# p value> 0.05 so can accept exponentially distributed

qqplot(qexp(ppoints(length(biterate$bite_rate))), biterate$bite_rate)
abline(a = 0, b = 1)
