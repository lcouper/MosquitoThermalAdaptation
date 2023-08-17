##### Life history trait experiment downstream analysis #####

# Q1: Investigating variation in thermal performance of fitness
# Q2: What climatic variables associated with that variation (if it exists)?
# Q3: What life history traits have greatest impact on temperature-dependence of fitness? 

#### Q1 Variation in thermal performance of fitness #####
# Investigating extent of variation in different parameters
# And evidence for specialist-generalist vs hotter-is-better patterns

library(tidyr)
library(ggplot2)

# Set working directory
setwd("~/Documents/Current Projects/LifeHistoryTraitExp/Analysis_TraitFits/")
load("fit_meansd_inf.Rsave")
load("fit_meansd_uni.Rsave")

# Range of CTmax, CTmin, Topt, Tbreath, and Pmax (across all populations)
Mins = as.data.frame(aggregate(mean ~ param, FUN = min, data = fitdf_inf))
Maxs = as.data.frame(aggregate(mean ~ param, FUN = max, data = fitdf_inf))
Ranges = cbind(Mins, Maxs[,2])
colnames(Ranges) = c("param", "min", "max")
Ranges$Amount = Ranges$max - Ranges$min
# see minimal variation in Topt and CTmax. Greater variation in CTmin. 
# but no clear between-population variation in either of these 3 parameters

# Repeat with uniform priors
Mins = as.data.frame(aggregate(mean ~ param, FUN = min, data = fitdf_uni))
Maxs = as.data.frame(aggregate(mean ~ param, FUN = max, data = fitdf_uni))
Ranges = cbind(Mins, Maxs[,2])
colnames(Ranges) = c("param", "min", "max")
Ranges$Amount = Ranges$max - Ranges$min

# Generalist-specialist trade-off investigation #
# would appear as negative correlation between thermal breadth and max trait performance 
library(boot)
gs1 = fitdf_inf[fitdf_inf$param == "Tbreadth" | fitdf_inf$param == "Pmax", c(1,4, 5)]
gs = gs1 %>% spread(param, mean)
summary(lm(Pmax ~ Tbreadth, data = gs))
corr(as.matrix(cbind(gs$Pmax, gs$Tbreadth)))

# not finding any real relationship
plot(Pmax ~ Tbreadth, data = gs, pch = 16, cex = 1.4,
     xlab = "Thermal breadth (\u00B0C)", ylab = "Max performance",
     cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
abline(lm(Pmax ~ Tbreadth, data = gs), col = "red", cex = 2)
text(13.75, 9.8, expression(R^2 == 0.54), cex = 1.5)
text(13.75, 8.8, substitute(paste(italic('p'), " = 0.015")), cex = 1.5)

# repeat with uniform priors 
load("fit_meansd_uni.Rsave")
gsu1 = fitdf_uni[fitdf_uni$param == "Tbreadth" | fitdf_uni$param == "Pmax", c(1,4, 5)]
gsu = gsu1 %>% spread(param, mean)
summary(lm(Pmax ~ Tbreadth, data = gsu)) # Relationship does not hold
corr(as.matrix(cbind(gsu$Pmax, gsu$Tbreadth)))
# again not finding any real relationship

# Hotter is better investigation #
# would appear as positive correlation between Topt and max trait performance
hisb1 = fitdf_inf[fitdf_inf$param == "Topt" | fitdf_inf$param == "Pmax", c(1,4, 5)]
hisb = hisb1 %>% spread(param, mean)
summary(lm(Pmax ~ Topt, data = hisb)) # Not significant
corr(as.matrix(cbind(hisb$Pmax, hisb$Topt)))

# repeat with uniform priors to check if relationship holds
hisbu1 = fitdf_uni[fitdf_uni$param == "Topt" | fitdf_uni$param == "Pmax", c(1,4, 5)]
hisbu = hisbu1 %>% spread(param, mean)
summary(lm(Pmax ~ Topt, data = hisbu)) # now a significant negative correlation...
corr(as.matrix(cbind(hisbu$Pmax, hisbu$Topt)))


plot(Pmax ~ Topt, data = hisb, pch = 16, cex = 1.4,
     xlab = "Thermal optima (\u00B0C)", ylab = "Max performance",
     cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
abline(lm(Pmax ~ Topt, data = hisb), col = "red", cex = 2)
text(22.45, 9.8, expression(R^2 == 0.19), cex = 1.5)
text(22.45, 8.8, substitute(paste(italic('p'), " = 0.204")), cex = 1.5)


#### Q2. Life history traits with biggest impact on fitness ######

#### Q2 part 0: exploratory analysis of variation in fitness ####
##### Q2 part 1: variation between life history traits ######
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

load("pLS_meansd_inf.Rsave") # dataframe = pLSdf_inf
load("LDR_meansd_inf.Rsave") # dataframe = LDRdf_inf
load("pPS_meansd_inf.Rsave") # dataframe = pPSdf_inf
load("PDR_meansd_inf.Rsave") # dataframe = PDRdf_inf
load("AL_meansd_inf.Rsave")  # dataframe = ALdf_inf
load("fec_meansd_inf.Rsave")  # dataframe = fecdf_inf
load("fit_meansd_inf.Rsave") # dataframe = fitdf_inf

pLSdf_inf$trait = "Larval.Survival"
LDRdf_inf$trait = "LarvalDevRate"
pPSdf_inf$trait = "Pupal.Survival"
PDRdf_inf$trait = "PupalDevRate"
ALdf_inf$trait = "AdultLifespan"
fecdf_inf$trait = "Fecundity"
fitdf_inf$trait = "Fitness"

load("pLS_meansd_uni.Rsave") # dataframe = pLSdf_uni
load("LDR_meansd_uni.Rsave") # dataframe = LDRdf_uni
load("pPS_meansd_uni.Rsave") # dataframe = pPSdf_uni
load("PDR_meansd_uni.Rsave") # dataframe = PDRdf_uni
load("AL_meansd_uni.Rsave")  # dataframe = ALdf_uni
load("fec_meansd_uni.Rsave")  # dataframe = fecdf_uni
load("fit_meansd_uni.Rsave") # dataframe = fitdf_uni

pLSdf_uni$trait = "Larval.Survival"
LDRdf_uni$trait = "LarvalDevRate"
pPSdf_uni$trait = "Pupal.Survival"
PDRdf_uni$trait = "PupalDevRate"
ALdf_uni$trait = "AdultLifespan"
fecdf_uni$trait = "Fecundity"
fitdf_uni$trait = "Fitness"

AllTraitsDf = rbind(pLSdf_inf, LDRdf_inf, 
                     pPSdf_inf, PDRdf_inf, ALdf_inf, fecdf_inf, fitdf_inf)

AllTraitsDf2 = rbind(pLSdf_uni, LDRdf_uni, 
                    pPSdf_uni, PDRdf_uni, ALdf_uni, fecdf_uni, fitdf_uni)

# Calculate means for each trait and by parameter
# These values become Supplemental Table 4.
MeansDf = AllTraitsDf %>% group_by(param, trait) %>% 
  dplyr::summarise(mean_value = mean(mean),
                   min_value = min(mean), max_value = max(mean)) %>% as.data.frame()
# These values become Supplemental Table 5.
MeansDf2 = AllTraitsDf2 %>% group_by(param, trait) %>% 
  dplyr::summarise(mean_value = mean(mean),
                   min_value = min(mean), max_value = max(mean)) %>% as.data.frame()

# Calculate ranges across all traits for each parameter, excluding fitness
# These values referenced in results
AllTraitsDf = AllTraitsDf[AllTraitsDf$trait != "Fitness",]
RangesDf = AllTraitsDf %>% group_by(param, trait) %>% 
  dplyr::summarise(range_value = max(mean) - min(mean)) %>% as.data.frame()





##### Q2 part 2: impact of life history traits on fitness (old) #####
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Fitness_Breakdown")
df = read.csv("Fitness_SensitivityAnalysis.csv")[,-1]
library(lme4)

# scale function: converts values into z-score (i.e.: (x - xmean) / x_stdev)
# allowing us to compare coefficient estimates between traits directly
# includes a fixed effect for temperature treatment and a random intercept for population
# each life history trait then included as random effect

# Note: previously I tried re-scaling all variables between 0-1 before applying
# the z-score transformation, to ensure probability and continuous traits were
# comparable. I,e: 
#dfscale = df
#dfscale[,c(3:9)] = lapply(dfscale[,c(3:9)], function(x) rescale(x, to = c(0,1)))
# However, the coefficient estimates were identical, 
# so I proceeded to use the original data for the z-score transformation

# Outcome in first model is the average fitness per individual 
mod1 = lmer(scale(PerIndFit) ~ scale(Larval.Survival) + scale(LarvalDevRate) + 
              scale(Pupal.Survival) + scale(PupalDevRate) + 
              scale(AdultLifespan) + scale(EggCount) +
              (1 |Population) + factor(Temp.Treatment), data = df)
summary(mod1) # This output used in Table 1 of manuscript

# Outcome in second model is the percent individuals with positive fitness 
mod2 = lmer(scale(PerPosFit) ~ scale(Larval.Survival) + scale(LarvalDevRate) + scale(Pupal.Survival) + 
              scale(PupalDevRate) + scale(AdultLifespan) + scale(EggCount) +
              (1|Population) + factor(Temp.Treatment) - 1, data = dfscale)
summary(mod2) # This output used in Table 1 of manuscript

# are egg counts and ldr strongly correlated?
tester = read.csv("Fitness_SensitivityAnalysis_IndividualLevel.csv", header = T)[,-1]
t2 = tester %>% group_by(Temp.Treatment) %>% 
  summarise(corr = cor(as.matrix(LarvalDevRate, EggCount))) %>% as.data.frame()

#### Q3. Temperature drivers of variation in fitness ######

setwd("~/Documents/Current Projects/LifeHistoryTraitExp")
library(dplyr)

# Load TPC parameters (CTmin, CTmax, Topt, Pmax, Tbreadth) for fitness
load("Analysis_TraitFits/fit_meansd_inf.Rsave")

# Pull in calculated temperature variables (from PRISM data)
TempVars = read.csv("TemperatureData/PRISM/TempVarsAtTreeholes.csv", header = T)[,-1]

# Part 1: Examine correlations between temperature variables
library(corrplot)

# re-organized columns and convert to matrix
TempMat = as.matrix(TempVars[,c(2,6,3,4,5)])
colnames(TempMat) = c("Annual mean", "Jan-March mean", "Seasonal variation",
                      "Spring/Summer daily max", "Days > 35 (\u00B0C)")
TempCorr = cor(TempMat)
corrplot(TempCorr,  addCoef.col = 'black', type = "lower",
         tl.col = "black", tl.srt = 45)
# For significance-testing below: calculate the effective number of tests
library(poolr)
meff(TempCorr, method = "nyholt") 

# Part 2: Investigate correlations between TPC parameters and temperature variables
library(plyr)
df = join(fitdf_inf, TempVars, by = "Population", type = "left")

# Examine relationship between Pmax of fitness and temperature variables
dfPmax = df[df$param == "Pmax",]

cor.test(dfPmax$mean, dfPmax$avgAnnualMean) # annual mean: 0.662, sig
cor.test(dfPmax$mean, dfPmax$avgSeasonalVar) # seasonal var: -0.002
cor.test(dfPmax$mean, dfPmax$avgExtreme) # extreme: 0.449
cor.test(dfPmax$mean, dfPmax$avgNumDays) # >3C: 0.326
cor.test(dfPmax$mean, dfPmax$avgWetTemp) # mean wet temp: 0.581

# PvalAdj = p.adjust(PvalList,  method = "fdr", n = 3)

# Only other TPC parameters for individual life history traits with sig variation:
# CTmax LDR, CTmax PDR, Topt PDR
# Investigate relationships between these params and temperature variables:

# Load TPC parameters for LDR and PDR
load("Analysis_TraitFits/LDR_meansd_inf.Rsave")
load("Analysis_TraitFits/PDR_meansd_inf.Rsave")

# Combine with temperature data
df2 = join(LDRdf_inf, TempVars, by = "Population", type = "left")
df3 = join(PDRdf_inf, TempVars, by = "Population", type = "left")

# Examine relationship between CTmax of LDR and temperature variables
df2Tmax = df2[df2$param == "Tmax",]

cor.test(df2Tmax$mean, df2Tmax$avgAnnualMean) # annual mean: 0.140, non-sig
cor.test(df2Tmax$mean, df2Tmax$avgSeasonalVar) # seasonal var: 0.629, non-sig
cor.test(df2Tmax$mean, df2Tmax$avgExtreme) # extreme: 0.229, non-sig
cor.test(df2Tmax$mean, df2Tmax$avgNumDays) # >35C: 0.489, non-sig
cor.test(df2Tmax$mean, df2Tmax$avgWetTemp) # mean wet temp: -0.169, non-sig

# Examine relationship between CTmax of PDR and temperature variables
df3Tmax = df3[df3$param == "Tmax",]

cor.test(df3Tmax$mean, df3Tmax$avgAnnualMean) # annual mean: 0.637, sig
cor.test(df3Tmax$mean, df3Tmax$avgSeasonalVar) # seasonal var: 0.307, non-sig
cor.test(df3Tmax$mean, df3Tmax$avgExtreme) # extreme: 0.678, sig
cor.test(df3Tmax$mean, df3Tmax$avgNumDays) # >35 C: 0.695, sig
cor.test(df3Tmax$mean, df3Tmax$avgWetTemp) # mean wet temp: 0.386, nonsig

# see if significant after removing POW
df3TmaxNoPOW = df3Tmax[-10,]
cor.test(df3TmaxNoPOW$mean, df3TmaxNoPOW$avgAnnualMean) # now not significant
cor.test(df3TmaxNoPOW$mean, df3TmaxNoPOW$avgExtreme) # now barely not significant
cor.test(df3TmaxNoPOW$mean, df3TmaxNoPOW$days35) # still significant

# adjust p-values?
# PvalAdj = p.adjust(PvalList,  method = p.adjust.methods, n = length(PvalList))

# Examine relationship between Topt of PDR and temperature variables
df3Topt = df3[df3$param == "Topt",]

cor.test(df3Topt$mean, df3Topt$avgAnnualMean) # annual mean: 0.669
cor.test(df3Topt$mean, df3Topt$avgSeasonalVar)# seasonal var: 0.286
cor.test(df3Topt$mean, df3Topt$avgExtreme)# extreme: 0.708
cor.test(df3Topt$mean, df3Topt$avgNumDays)# >35 C: 0.703
cor.test(df3Topt$mean, df3Topt$avgWetTemp) # mean wet temp: 0.422

# see if significant after removing POW
df3ToptNoPOW = df3Topt[-10,]
cor.test(df3ToptNoPOW$mean, df3ToptNoPOW$avgAnnualMean) # now not significant
cor.test(df3ToptNoPOW$mean, df3ToptNoPOW$avgExtreme) # still significant
cor.test(df3ToptNoPOW$mean, df3ToptNoPOW$days35) # still significant


# adjust p-values?
# PvalAdj = p.adjust(PvalList,  method = p.adjust.methods, n = length(PvalList))

#### Q4: How often thermal limits exceed at present ####

setwd("~/Documents/Current Projects/LifeHistoryTraitExp/TemperatureData/PRISM/")

# pull in climate data files
over35 = read.csv("Treehole_Exceed35_2000_end2021.csv", header = T, row.names = 1)
over35$mean = rowMeans(over35)
over35s = read.csv("Treehole_Exceed35_SUBSET_2000_end2021.csv", header = T, row.names = 1)
over35s$mean = rowMeans(over35s)
over31.6 = read.csv("Treehole_DaysExceed31.6_2000_end2021.csv", header = T, row.names = 1)
over31.6$mean = rowMeans(over31.6)
over31.6s = read.csv("Treehole_DaysExceed31.6_SUBSET_2000_end2021.csv", header = T, row.names = 1)
over31.6s$mean = rowMeans(over31.6s)


##### Temperature Data Comparisons: iButtons #####

setwd("~/Documents/Current Projects/LifeHistoryTraitExp/TemperatureData/Comparisons")

prism = read.csv("PRISMdata_POW02.csv", header = T)
ibut = read.csv("iButtondata_POW02.csv", header = T)[-1]
ibut$Date.Time = as.POSIXct(ibut$Date.Time, format = "%m/%d/%y %H:%M")

# calculate average daily temperature from ibutton data (measured every 4 hours)
ibut$Day = strftime(ibut$Date.Time, format="%m/%d/%y")
ibutMean = ibut %>% dplyr::group_by(Day) %>%
  dplyr::summarize(mean_X1 = mean(Temp.C)) %>% as.data.frame()

# remove first two rows of ibutton data to match prism data
ibutMean = ibutMean[-c(1:2),]
# remove last row of prims data to match ibutton data
prism = prism[-nrow(prism),]

# Combine
TempData = cbind.data.frame(ibutMean,prism[,2])
colnames(TempData) = c("Date", "iButton", "PRISM")
TempData$Date = as.POSIXct(TempData$Date, format = "%m/%d/%y")

# Plot
plot(TempData$iButton ~ TempData$Date, type = "l", ylim = c(5,32),
     ylab = "Temperature (C)", xlab = "Day", main = "Mean Daily Temperature",
     cex.axis = 1.5, lwd = 1.5)
points(TempData$PRISM ~ TempData$Date, type = "l", col = "red", lwd = 1.5)
legend("topleft", legend = c("iButton", "PRISM"), col = c("black", "red"), lty = 1, lwd = 4, cex = 1.5)
mean(TempData$PRISM - TempData$iButton) # PRISM on average 3C higher than ibutton

cor.test(TempData$PRISM, TempData$iButton)

# Repeat for daily max temperatures
ibutMax = ibut %>% dplyr::group_by(Day) %>%
  dplyr::summarize(max_temp = max(Temp.C)) %>% as.data.frame()
prism2= read.csv("PRISMdata_POW02_Max.csv", header = T)
ibutMax = ibutMax[-c(1:2),]
prism2 = prism2[-nrow(prism2),]
TempData2 = cbind.data.frame(ibutMax,prism2[,2])
TempData2$Day = as.POSIXct(TempData2$Day, format = "%m/%d/%y")
colnames(TempData2) = c("Date", "iButton", "PRISM")

plot(TempData2$iButton ~ TempData2$Date, type = "l", ylim = c(5,38),
     ylab = "Temperature (C)", xlab = "Day", main = "Max Daily Temperature", 
     cex.axis = 1.5, lwd = 1.5)
points(TempData2$PRISM ~ TempData2$Date, type = "l", col = "red", lwd = 1.5)
legend("topleft", legend = c("iButton", "PRISM"), col = c("black", "red"), lty = 1, lwd = 4, cex = 1.5)
mean(TempData2$PRISM - TempData2$iButton) # PRISM on average 5C higher than ibutton

cor.test(TempData2$PRISM, TempData2$iButton)


# Calculate mean daily max in spring and summer (3/19 - 9/23) 
# to compare with PRISM data source
mean(ibutMax$max_temp[52:240]) # 23.39947 (same variable from PRISM was 27.77)

# Calculate number of days where temperatures > 35
which(ibut$Temp.C > 35) # never exceeded
# maximum temperature recorded here was 34
which(ibutMax$max_temp > 31.6) # 8 days where temps > 31.6


#### Repeat for SB04 ####
prism = read.csv("PRISMdata_SB04.csv", header = T)
ibut = read.csv("iButtondata_SB04.csv", header = T)
ibut$Date.Time = as.POSIXct(ibut$Date.Time, format = "%m/%d/%y %H:%M")
ibut = ibut[-c(1:5),] # Delete first iButton data rows to match prism dataset
# remove last row of prims data to match ibutton data
prism = prism[-nrow(prism),]

# calculate average daily temperature from ibutton data (measured every 4 hours)
ibut$Day = strftime(ibut$Date.Time, format="%m/%d/%y")
ibutMean = ibut %>% dplyr::group_by(Day) %>%
  dplyr::summarize(mean_X1 = mean(Temp.C)) %>% as.data.frame()

# Combine
TempData = cbind.data.frame(ibutMean,prism[,2])
colnames(TempData) = c("Date", "iButton", "PRISM")
TempData$Date = as.POSIXct(TempData$Date, format = "%m/%d/%y")

# Plot
plot(TempData$iButton ~ TempData$Date, type = "l", ylim = c(5,36),
     ylab = "Temperature (C)", xlab = "Day", main = "Mean Daily Temperature",
     cex.axis = 1.5, lwd = 1.5)
points(TempData$PRISM ~ TempData$Date, type = "l", col = "red", lwd = 1.5)
legend("topleft", legend = c("iButton", "PRISM"), col = c("black", "red"), lty = 1, lwd = 4, cex = 1.5)
mean(TempData$PRISM - TempData$iButton) # PRISM on average 0.7C cooler than ibutton

cor.test(TempData$PRISM, TempData$iButton)


### Repeat SB04 with daily maxima ###

# calculate average daily max temperature from ibutton data (measured every 4 hours)
ibutMax = ibut %>% dplyr::group_by(Day) %>%
  dplyr::summarize(max_temp = max(Temp.C)) %>% as.data.frame()

prism2= read.csv("PRISMdata_SB04_Max.csv", header = T)
# remove last row of prims data to match ibutton data
prism2 = prism2[-nrow(prism2),]

# Combine
TempData2 = cbind.data.frame(ibutMax,prism2[,2])
colnames(TempData2) = c("Date", "iButton", "PRISM")
TempData2$Date = as.POSIXct(TempData2$Date, format = "%m/%d/%y")

# Plot
plot(TempData2$iButton ~ TempData2$Date, type = "l", ylim = c(5,50),
     ylab = "Temperature (C)", xlab = "Day", main = "Max Daily Temperature",
     cex.axis = 1.5, lwd = 1.5)
points(TempData2$PRISM ~ TempData2$Date, type = "l", col = "red", lwd = 1.5)
legend("topleft", legend = c("iButton", "PRISM"), col = c("black", "red"), lty = 1, lwd = 4, cex = 1.5)
mean(TempData2$PRISM - TempData2$iButton) # PRISM on average 5C higher than ibutton

cor.test(TempData2$PRISM, TempData2$iButton)



# Calculate mean daily max in spring and summer (3/19 - 9/23) 
# to compare with PRISM data source
mean(ibutMax$max_temp[30:218]) # 29.52 (same variable from PRISM was 27.86)

# Calculate number of days where temperatures > 35
which(ibutMax$max_temp > 35) # 45 days where temps > 35
which(ibutMax$max_temp > 31.6) # 87 days where temps > 31.6







###### Correlation Plots. Manuscript figure 6 #####
par(mar = c(4.5,4.5,2,1))
LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors = c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", 
                  "#abd9e9", "#74add1", "#4575b4", "#313695")

# Pmax and annual mean
dfPmax$Population = factor(dfPmax$Population, levels = LatOrder)
dfPmax = dfPmax[order(dfPmax$Population),]
plot(dfPmax$mean ~ dfPmax$avgAnnualMean, 
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2,
     main = "Maximum fitness", cex.main = 2,
     col = rev(LatColors),
     xlab = "Annual mean temperature (\u00B0C)", ylab = "offspring per individual")
abline(lm(dfPmax$mean ~ dfPmax$avgAnnualMean), col = "black", lwd= 2)
text("r = 0.66", x = 12.2, y = 9.7, cex = 2)

# PDR tmax and annual mean
df3Tmax$Population = factor(df3Tmax$Population, levels = LatOrder)
df3Tmax = df3Tmax[order(df3Tmax$Population),]
plot(df3Tmax$mean ~ df3Tmax$avgAnnualMean, 
     col = rev(LatColors),
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2, 
     main = "Thermal maximum for pupal dev. rate", cex.main = 2,
     xlab = "Annual mean temperature (\u00B0C)", ylab = "")
abline(lm(df3Tmax$mean ~ df3Tmax$avgAnnualMean), col = "black", lwd= 2)
text("r = 0.64", x = 12.2, y = 33.6, cex = 2)

# PDR tmax and avg # days >35 C
df3Tmax$Population = factor(df3Tmax$Population, levels = LatOrder)
df3Tmax = df3Tmax[order(df3Tmax$Population),]
plot(df3Tmax$mean ~ df3Tmax$avgNumDays, 
     col = rev(LatColors),
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2, 
     main = "Thermal maximum for pupal dev. rate", cex.main = 2,
     xlab = "# days >35 \u00B0C", ylab = "")
abline(lm(df3Tmax$mean ~ df3Tmax$avgNumDays), col = "black", lwd= 2)
text("r = 0.70", x = 8, y = 33.6, cex = 2)

# PDR tmax and avg # days >31.6C
df3Tmax$Population = factor(df3Tmax$Population, levels = LatOrder)
df3Tmax = df3Tmax[order(df3Tmax$Population),]
plot(df3Tmax$mean ~ df3Tmax$days31.6, 
     col = rev(LatColors),
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2, 
     main = "Thermal maximum for pupal dev. rate", cex.main = 2,
     xlab = "# days >31.6 \u00B0C", ylab = "")
abline(lm(df3Tmax$mean ~ df3Tmax$days31.6), col = "black", lwd= 2)
text("r = 0.47", x = 26, y = 33.6, cex = 2)




# PDR tmax and warm season maximum
df3Tmax$Population = factor(df3Tmax$Population, levels = LatOrder)
df3Tmax = df3Tmax[order(df3Tmax$Population),]
plot(df3Tmax$mean ~ df3Tmax$avgExtreme,
     col = rev(LatColors),
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2, 
     main = "Thermal maximum for pupal dev. rate", cex.main = 2,
     xlab = "Warm season max (\u00B0C)", ylab = "")
abline(lm(df3Tmax$mean ~ df3Tmax$avgExtreme), col = "black", lwd= 2)
text("r = 0.68", x = 23.3, y = 33.6, cex = 2)

# PDR topt and avg daily maximum in spring/summer
df3Topt$Population = factor(df3Topt$Population, levels = LatOrder)
df3Topt = df3Topt[order(df3Topt$Population),]
plot(df3Topt$mean ~ df3Tmax$avgExtreme, 
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2, 
     col = rev(LatColors), 
     main = "Thermal optimum for pupal dev. rate", cex.main = 2,
     xlab = "Daily max in Spring/Summer (\u00B0C)", ylab = "")
abline(lm(df3Topt$mean ~ df3Tmax$avgExtreme), col = "black", lwd= 2)
text("r = 0.71", x = 23.3, y = 27.6, cex = 2)


# PDR topt and days >35C
df3Topt$Population = factor(df3Topt$Population, levels = LatOrder)
df3Topt = df3Topt[order(df3Topt$Population),]
plot(df3Topt$mean ~ df3Tmax$days35, 
     pch =16, cex = 2.5, cex.lab = 2, cex.axis = 2, 
     col = rev(LatColors), 
     main = "Thermal optimum for pupal dev. rate", cex.main = 2,
     xlab = "# days >35 \u00B0C", ylab = "")
abline(lm(df3Topt$mean ~ df3Tmax$days35), col = "black", lwd= 2)
text("r = 0.70", x = 8, y = 27.6, cex = 2)





##### 1-D plot of days > 35c #####

LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColorsP = c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")

TempVars$Population = factor(TempVars$Population, levels = LatOrder)
TempVars = TempVars[order(TempVars$Population),]
x = TempVars$avgNumDays

ggplot(data.frame(x), aes(x=x, y=0)) +
  geom_point(size = 10)  +
  annotate("segment",x=1,xend=50, y=0, yend=0, size=2) +
  annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=2) +
  annotate("segment",x=50,xend=50, y=-0.1,yend=0.1, size=2) +
  scale_x_continuous(limits = c(1,50)) +
  scale_y_continuous(limits = c(-0.1,0.1)) +
  theme(panel.background = element_blank())

x = TempVars$days31.6
ggplot(data.frame(x), aes(x=x, y=0)) +
  geom_point(size = 10)  +
  annotate("segment",x=1,xend=100, y=0, yend=0, size=2) +
  annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=2) +
  annotate("segment",x=100,xend=100, y=-0.1,yend=0.1, size=2) +
  scale_x_continuous(limits = c(1,100)) +
  scale_y_continuous(limits = c(-0.1,0.1)) +
  theme(panel.background = element_blank())

##### Additional analyses ####
##### Is there a climate/geographical gradient in body size? ######
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")

data = read.csv("LifeHistoryTraitExp_Data.csv", header = T)[,-1]
data$Population = dplyr::recode(data$Population, "MARIN29" = 'MAR2', "MARIN35" = 'MAR1')
data$Temp.Treatment = dplyr::recode(data$Temp.Treatment, '32.18' = 32, '27.5' = 28, '23.75' = 24, '17.12' = 17, '13.35' = 13, '5.49' = 5) 

sizes = data[,c(1:3,13,15)]
# keep only with left wing length data
sizes = sizes[!is.na(sizes$Left.Wing.Length),]
sizes$Population = as.factor(sizes$Population)
sizes$Temp.Treatment = as.factor(sizes$Temp.Treatment)
sizes$Sex = as.factor(sizes$Sex)
sizes$Counter = 1

# separate males and females
sizesM = sizes[sizes$Sex == "M",]
sizesF = sizes[sizes$Sex == "F",]

LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
sizesM$Population = factor(sizesM$Population, levels = LatOrder)
sizesF$Population = factor(sizesF$Population, levels = LatOrder)

LatColors2 = rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"))

sizesM %>% ggplot(aes(x = Temp.Treatment, y = Left.Wing.Length, fill = Population)) +
  geom_boxplot() + 
  labs(x = "Temperature treatment (\u00B0C)", y = "Wing length (mm)") + 
  scale_fill_manual(values = LatColors2) + 
#  ggtitle("Adult males") + 
  theme_minimal() +
  theme(axis.text=element_text(size=15), 
      axis.title = element_text(size = 16),
      legend.text=element_text(size=15),
      legend.title = element_text(size=16),
      plot.title = element_text(hjust = 0.5, size = 24),
      legend.position= "none",
      panel.border = element_rect(colour = "black", fill = NA)) 

sizesF %>% ggplot(aes(x = Temp.Treatment, y = Left.Wing.Length, fill = Population)) +
  geom_boxplot() + 
  labs(x = "Temperature treatment (\u00B0C)", y = "Wing length (mm)") + 
  scale_fill_manual(values = LatColors2) + 
  theme_minimal() + 
#  ggtitle("Adult females") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 

##### Different metrics of fitness (prop non-0, avg per ind): Supp Figure 7 Heatmap #####
# Supplemental Table 7
library(ggplot2)

lfdf = read.csv("~/Documents/Current_Projects/LifeHistoryTraitExp/Fitness_Breakdown/Fitness_SensitivityAnalysis.csv", header= T)[-1]
lfdf$Temp.Treatment = dplyr::recode(lfdf$Temp.Treatment, '32.18' = 32, '27.5' = 28, '23.75' = 24, '17.12' = 17, '13.35' = 13, '5.49' = 5) %>% 
as.factor()
lfdf$Population = dplyr::recode(lfdf$Population, "MARIN35" = 'MAR1', "MARIN29" = 'MAR2')
lfdf$FitPerSurv = round(lfdf$FitPerSurv, 2)
lfdf$PerIndFit = round(lfdf$PerIndFit, 2)
lfdf$PerPosFit = round(lfdf$PerPosFit, 2)
lfdf$EggCount[is.na(lfdf$EggCount)] <- 0
lfdf$EggCount = round(lfdf$EggCount, 1)
lfdf$FitPerSurv[is.na(lfdf$FitPerSurv)] <- 0

LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
lfdf$Population = factor(lfdf$Population, levels = LatOrder)

ggplot(lfdf, aes(x = Temp.Treatment, y = Population, fill = PerIndFit)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "#DB4325") +
  coord_fixed() + 
  geom_text(aes(label=PerIndFit), size = 3) + 
  xlab("Temperature (\u00B0C)") + theme_minimal()

ggplot(lfdf, aes(x = Temp.Treatment, y = Population, fill = PerPosFit)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "#DB4325") +
  coord_fixed() + 
  geom_text(aes(label=PerPosFit), size = 3) + 
  xlab("Temperature (\u00B0C)") + theme_minimal()

ggplot(lfdf, aes(x = Temp.Treatment, y = Population, fill = FitPerSurv)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "#DB4325") +
  coord_fixed() + 
  geom_text(aes(label=FitPerSurv), size = 3) + 
  xlab("Temperature (\u00B0C)") + theme_minimal()
