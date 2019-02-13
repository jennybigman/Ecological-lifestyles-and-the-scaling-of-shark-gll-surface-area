# Code for analyses in Bigman et al. (2018) Journal of Morphology

# Please cite the following if code or data is used:
# Bigman JS, Pardo SA, Prinzing TS, Dando M, Wegner NC, Dulvy NK. 
# Ecological lifestyles and the scaling of shark gill surface area. 
# Journal of Morphology. 2018;1–9. https://doi.org/10.1002/jmor.20879

# packages
# devtools::install_github("dgrtwo/broom") 
library(tidyverse)
library(tidyr)
library(cowplot)
library(ggplot2)
library(MASS) 
library(broom)
library(lme4)
library(nlme)
library(forcats)

#### Data wrangling ####

## Mustelus gill surface area data
MusCaliGSA <- read.csv("MustelusCalifornicusGSA_Data.csv", header = TRUE, 
                       strip.white = TRUE) %>% tbl_df() # smoothhound raw data
## Gill surface area data for other 11 sharks
AllSharkData <- read.csv("Bigman_raw_data_combined.csv", header = TRUE, 
                         strip.white = TRUE) %>% tbl_df

## Mustelus data transformations

#Gill surface area
MusCaliGSA$MassG <- (MusCaliGSA$MassKG*1000)
MusCaliGSA$Log10MassG <- log10(MusCaliGSA$MassG)
MusCaliGSA$Log10GSAcm2 <- log10(MusCaliGSA$GSAcm2)

# Total filament length
MusCaliGSA$Log10FilamentLength <- log10(MusCaliGSA$TotalFilamentLength2Sides)

# Average lamellar frequency
MusCaliGSA$Log10LamFreq <- log10(MusCaliGSA$LamellarFrequency)

# Mean bilateral lamellar surface area
MusCaliGSA$Log10LamSA <- log10(MusCaliGSA$BilateralLamellarSA)


## Remaining shark data transformations
AllSharkData$Log10MassG <- log10(AllSharkData$MassG)
AllSharkData$Log10GSAcm2 <- log10(AllSharkData$GSAcm2)

AllSharkData$Log10CenterMassG <- scale(AllSharkData$Log10MassG, 
                                       center = log10(5000), scale = FALSE)


#### Analyses ####

## Mustelus gill surface area 

# Total gill surface area regression and prediction
GSA.mass.log <- lm(Log10GSAcm2 ~ Log10MassG, data = MusCaliGSA)
summary(GSA.mass.log)

predict <- predict(GSA.mass.log, interval="prediction")
MusCaliGSA <- cbind(MusCaliGSA, predict)
MusCaliGSA$Fit.NoLog <- 10^(MusCaliGSA$fit) 
MusCaliGSA$Upper.NoLog <- 10^(MusCaliGSA$upr)
MusCaliGSA$Lower.NoLog <- 10^(MusCaliGSA$lwr)


# Total filament length regression and prediction
FilamentLength.mass <- lm(Log10FilamentLength ~ Log10MassG, data = MusCaliGSA)
summary(FilamentLength.mass)

predict.FilamentLengthLog10 <- predict(FilamentLength.mass, interval="prediction")
MusCaliGSA.FilamentLength <- data.frame(MusCaliGSA$MassG, MusCaliGSA$TotalFilamentLength2Sides,
                                        MusCaliGSA$Log10FilamentLength, MusCaliGSA$Log10MassG,
                                        predict.FilamentLengthLog10)
MusCaliGSA.FilamentLength$FitFilLength.NoLog <- 10^(MusCaliGSA.FilamentLength$fit) 
MusCaliGSA.FilamentLength$UpperFilLength.NoLog <- 10^(MusCaliGSA.FilamentLength$upr)
MusCaliGSA.FilamentLength$LowerFilLength.NoLog <- 10^(MusCaliGSA.FilamentLength$lwr)


# Average lamellar frequency regression and prediction
LamFreq.mass <- lm(Log10LamFreq ~ Log10MassG, data = MusCaliGSA)
summary(LamFreq.mass)

predict.Log10LamFreq <- predict(LamFreq.mass, interval="prediction")
MusCaliGSA.LamFreq <-data.frame(MusCaliGSA$MassG, MusCaliGSA$LamellarFrequency, 
                                MusCaliGSA$Log10MassG, MusCaliGSA$Log10LamFreq,
                                predict.Log10LamFreq)
MusCaliGSA.LamFreq$FitFilLength.NoLog <- 10^(MusCaliGSA.LamFreq$fit) 
MusCaliGSA.LamFreq$UpperFilLength.NoLog <- 10^(MusCaliGSA.LamFreq$upr)
MusCaliGSA.LamFreq$LowerFilLength.NoLog <- 10^(MusCaliGSA.LamFreq$lwr)


# Mean bilateral lamellar surface area regression and prediction
LamSA.mass <- lm(Log10LamSA ~ Log10MassG, data = MusCaliGSA)
summary(LamSA.mass)

predict.Log10LamSA <- predict(LamSA.mass, interval="prediction")
MusCaliGSA.LamSA <- data.frame(MusCaliGSA$MassG, MusCaliGSA$BilateralLamellarSA,
                               MusCaliGSA$Log10LamSA, MusCaliGSA$Log10MassG,
                               predict.Log10LamSA)
MusCaliGSA.LamSA$FitFilLength.NoLog <- 10^(MusCaliGSA.LamSA$fit) 
MusCaliGSA.LamSA$UpperFilLength.NoLog <- 10^(MusCaliGSA.LamSA$upr)
MusCaliGSA.LamSA$LowerFilLength.NoLog <- 10^(MusCaliGSA.LamSA$lwr)


#### Scaling of gill surface area for all 12 sharks ####

## Regressions to estimate coefficients
fitted_models.Log10CenterMass <- AllSharkData %>% 
  group_by(Species) %>% 
  do(fits.Log10CenterMass = lm(Log10GSAcm2 ~ Log10CenterMassG, data = .))

model.output.Log10CenterMass <- fitted_models.Log10CenterMass %>%
  tidy(fits.Log10CenterMass)
fits.Log10CenterMass <- fitted_models.Log10CenterMass %>% 
  augment(fits.Log10CenterMass)

## Do intercepts and slopes differ across species? 
mod1 <- lm(Log10GSAcm2 ~ Log10CenterMassG * Species - 1, data = AllSharkData) # 
coef(mod1)
summary(mod1)

## Do ecological lifestyle traits differ across species? 

# 1. Caudal fin aspect ratio
ActivityLevelDifferences  <- lme(Log10GSAcm2 ~ Log10CenterMassG * AspectRatio, 
                                 random = ~ Log10CenterMassG | Species,
                                 data = AllSharkData, control = lmeControl(opt = "optim"))
coef(ActivityLevelDifferences)
summary(ActivityLevelDifferences)

# 2. Habitat type 
HabitatDifferences  <- lme(Log10GSAcm2 ~ Log10CenterMassG * HabitatType, 
                           random = ~ Log10CenterMassG | Species,
                           data = AllSharkData, control = lmeControl(opt = "optim"))
coef(HabitatDifferences)
summary(HabitatDifferences)

# 3. Maximum size
MaxSizeDifferences  <- lme(Log10GSAcm2 ~ Log10CenterMassG * LogMaxMassG, 
                           random = ~ Log10CenterMassG | Species,
                           data = AllSharkData, control = lmeControl(opt = "optim"))
coef(MaxSizeDifferences)
summary(MaxSizeDifferences)








