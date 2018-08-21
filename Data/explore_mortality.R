
## explore lifespan ##
library(dplyr)
library(readr)
library(ggplot2)

SurvivalData <- read_csv("Z:/YF_climateChange/Data/SurvivalData_mordecai2018.csv")

# we keep ZIKv infected because this is a factor we may see varying mosquito lifespan in the wild

ggplot(SurvivalData) + geom_point( aes(x = `Time (dpi)`, y = Alive/(Dead + Censored + Alive), colour = Temp))


SurvivalData$alive_propn = SurvivalData$Alive/(SurvivalData$Dead +  SurvivalData$Alive)


# plot(density(filter(SurvivalData, Temp == 16 & `Time (dpi)` == 1)$alive_propn), xlim = c(0.95,1))
# lines(density(filter(SurvivalData, Temp == 20 & `Time (dpi)` == 1)$alive_propn))
# lines(density(filter(SurvivalData, Temp == 24 & `Time (dpi)` == 1)$alive_propn))
# lines(density(filter(SurvivalData, Temp == 28 & `Time (dpi)` == 1)$alive_propn))
# lines(density(filter(SurvivalData, Temp == 32 & `Time (dpi)` == 1)$alive_propn))
# lines(density(filter(SurvivalData, Temp == 34 & `Time (dpi)` == 1)$alive_propn))
# lines(density(filter(SurvivalData, Temp == 36 & `Time (dpi)` == 1)$alive_propn))


### dead propn at each time step and temp
SurvivalData$dead_propn = SurvivalData$Dead/(SurvivalData$Dead +  SurvivalData$Alive)

hist(SurvivalData$dead_propn)

ggplot(SurvivalData) + geom_point( aes(x = Temp, y = dead_propn, colour = `Time (dpi)`) )
