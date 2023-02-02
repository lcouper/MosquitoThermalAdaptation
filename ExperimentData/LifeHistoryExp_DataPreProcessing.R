#### Life History Experiment & Reproduction Experiment: Data Pre-Processing ####

# Lisa Couper, Stanford University
# Code overview: 
# Pre-processing steps for main experiment:
# 1. take-in raw experimental data
# Remove any rows with issues
# Calculate traits (e.g., length larval development, lifespan, EFD)


#### Pre-processing for main experiment ####

# Bring in data
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")
data = read.csv("LifeHistoryData_PreProcessing.csv", header = T, na.strings = c(" ", "NA"))

# Remove unneeded columns & rows with any issues 
# (i.e., labeled "Remove" in Usage column, n = 11)

data = data[data$Usage != "Remove",]
data = data[data$Usage != "DontSequence",]
data = data[data$Usage != "SomeMissing",]
data = data[,c(1:14, 23,24)]
data = data[!apply(data == "", 1, all),] # remove any empty rows
# Should be 1751 rows 

# Convert date to date class 
data$Date.Started = as.Date(data$Date.Started, format = "%m/%d/%y")
data$Date.Larval.Death = as.Date(data$Date.Larval.Death, format = "%m/%d/%y")
data$Date.Pupation = as.Date(data$Date.Pupation, format = "%m/%d/%y")
data$Date.Adult.Eclosion = as.Date(data$Date.Adult.Eclosion, format = "%m/%d/%y")
data$Date.Adult.Death = as.Date(data$Date.Adult.Death, format = "%m/%d/%y")

# Convert temp treatment assignment to numeric temperature #
# Below is the avg temp recorded by ibutton during experiment
# F = 32.18, E = 27.5, D = 23.75, C = 17.12, B = 13.35, A = 5.49
data$Temp.Treatment = dplyr::recode(data$Temp.Treatment, F = 32.18, E = 27.5, D = 23.75, C = 17.12, B = 13.35, A = 5.49)

# Calculate traits 
# time = days. rates = 1/days

# Larval to pupal development time
data$Length.Larval.Dev = as.numeric(data$Date.Pupation - (data$Date.Started) + 1) 
# Subtract 1 since date started refers to day individuals were put into treatment (when they were 1 day old)

# Pupal to adult development time
data$Length.Pupal.Dev= as.numeric(data$Date.Adult.Eclosion - data$Date.Pupation)

# Larval to adult development time (aka juvenile development)
data$Length.Juvenile.Dev= as.numeric(data$Date.Adult.Eclosion - (data$Date.Started) + 1) 

# Convert times to rates
data$LarvalDevRate = 1/data$Length.Larval.Dev
data$PupalDevRate = 1/data$Length.Pupal.Dev
data$JuvenileDevRate = 1/data$Length.Juvenile.Dev

# larval survival (remove data already in column)
data$Larval.Survival = NA
data$Larval.Survival[!is.na(data$Date.Pupation)] <- "1"
data$Larval.Survival[!is.na(data$Date.Larval.Death)] <- "0"
# Individuals at 5.49 treatment with an 'NA' for larval survival = 
# survived until end of experiment (82-84 days). Denote these as a "1"
data$Larval.Survival[data$Temp.Treatment == "5.49" & is.na(data$Larval.Survival)] <- "1"

# pupal survival (remove data already in column)
data$Pupal.Survival = NA
data$Pupal.Survival[!is.na(data$Date.Adult.Eclosion)] <- "1"
data$Date.Pupal.Death[data$Date.Pupal.Death == ""] <- NA
data$Pupal.Survival[!is.na(data$Date.Pupal.Death)] <- "0"

# juvenile survival (create new column)
data$Juvenile.Survival = 0
data$Juvenile.Survival[!is.na(data$Date.Adult.Eclosion)] <- "1"

# adult lifespan
data$AdultLifespan = as.numeric(data$Date.Adult.Death - data$Date.Adult.Eclosion)

# diapause
data$Larval.Diapause[data$Larval.Diapause == "" ] = "n"

# estimated eggs laid in first clutch
# Using Washburn et al. 1989 proxy
data$Left.Wing.Length = as.numeric(data$Left.Wing.Length)
data$EggCount = (data$Left.Wing.Length * 51.33) - 87.96

# Calculate fitness 

# this function calculates fitness as follows: 
# survival to egg-laying time x eggs laid in first clutch
# survival to egg-laying time determined from reproduction validation experiment as:
# 4 (min # of days to sexual maturity from Yee & Anderson) + 
# 6 days for egg-laying for 28 and 24 C treatments, 7 days for 17 C, 13 days at 13 C (i.e., 10, 11, 17)
# survival is then multiplied by the estimated # eggs laid to calculate fitness

fitness = function(x)
{Survival = rep(NA, nrow(x)) 
#Fitness = rep(NA, nrow(x)) 
EggCount = x$EggCount
 x$AdultLifespan[is.na(x$AdultLifespan)] <- 0
 x$EggCount[is.na(x$EggCount)] <- 0
      for (i in 1:nrow(x))
        if (x$Temp.Treatment[i] == 32.18 | x$Temp.Treatment[i] == 5.49) {
                Survival[i] = 0 # adults either did not emerge or did not survive >1 day at these treatments
        } else if (x$Temp.Treatment[i] == 27.50 | x$Temp.Treatment[i] == 23.75) {
                Survival[i] = as.numeric(x$AdultLifespan[i] >= 10)
        } else if (x$Temp.Treatment[i] == 17.12) {
                 Survival[i] = as.numeric(x$AdultLifespan[i] >= 11)
        } else if (x$Temp.Treatment[i] == 13.35) {
                 Survival[i] = as.numeric(x$AdultLifespan[i] >= 17)
        } else {
                Survival[i] = 0
        }
   Fitness = Survival*(EggCount/2) # divided by 2 to account for both males and females being included
  print(Fitness)
   }

data$fitness = fitness(data)
# remove any rows where wing was broken and could not be measured 
# as these would falsely appear as '0' fitness
# these rows coded as 'yes' in 'Remove.From.Fitness.Count' column
data = data[data$Remove.From.Fitness.Count != "yes",]
data$fitness[is.na(data$fitness)] <- 0

## Calculate a sex-ratio weighted fitness 
# (e.g., weighted by sex ratio of each population/temp.treatment)

SexWeightedFitness = function(x)
{Survival = rep(NA, nrow(x)) 
EggCount = x$EggCount
x$AdultLifespan[is.na(x$AdultLifespan)] <- 0
x$EggCount[is.na(x$EggCount)] <- 0
for (i in 1:nrow(x))
  if (x$Temp.Treatment[i] == 32.18 | x$Temp.Treatment[i] == 5.49) {
    Survival[i] = 0 # adults either did not emerge or did not survive >1 day at these treatments
  } else if (x$Temp.Treatment[i] == 27.50 | x$Temp.Treatment[i] == 23.75) {
    Survival[i] = as.numeric(x$AdultLifespan[i] >= 10)
  } else if (x$Temp.Treatment[i] == 17.12) {
    Survival[i] = as.numeric(x$AdultLifespan[i] >= 11)
  } else if (x$Temp.Treatment[i] == 13.35) {
    Survival[i] = as.numeric(x$AdultLifespan[i] >= 17)
  } else {
    Survival[i] = 0
  }
x$Survival = Survival
# convert sex to 0 (Male) and 1 (Female)
x$SexNum = as.numeric(x$Sex == "F")
# calculate the sex ratio
z = aggregate(x$SexNum ~ x$Population + x$Temp.Treatment, FUN = sum, na.rm = T)
z2 = aggregate(x$SexNum ~ x$Population + x$Temp.Treatment, FUN = length)
z3 = z[,3]/z2[,3]
SexRatios = cbind(z[,1:2],z3)
colnames(SexRatios) = c("Population", "Temp.Treatment", "SexRatio")
df = join(x, SexRatios, by = c("Population", "Temp.Treatment"))
df$WeightedEggCount = df$EggCount * df$SexRatio
df$Fitness = df$Survival*df$WeightedEggCount 
print(df$Fitness)
}

data$SexWeightedFitness = SexWeightedFitness(data)

# again remove any rows where wing was broken and could not be measured 
# as these would falsely appear as '0' fitness
# these rows coded as 'yes' in 'Remove.From.Fitness.Count' column
data = data[data$Remove.From.Fitness.Count != "yes",]
data$SexWeightedFitness[is.na(data$SexWeightedFitness)] <- 0

#### Output cleaned data for main experiment #####
write.csv(data, "LifeHistoryTraitExp_Data.csv")
dataAdults = data[,c(1,2,3,13, 24,15)]
dataAdults = dataAdults[!is.na(dataAdults$AdultLifespan),]
#### Life history trait sensitivity analysis #####

library(dplyr)
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")
lfdata = read.csv("LifeHistoryTraitExp_Data.csv", header = T)[,-1]

# Want to calculate the following variables as population averages for each temp treatment:
# prop. larval survival, larval dev rate, prop. pupal survival, pupal dev rate, adult life span, fecundity
# as well as outcome measures: avg fitness per individual, prop. positive fitness 

# First, add column to indicate whether individuals had positive fitness
lfdata$PosFitness = as.numeric(lfdata$fitness > 0)

lfdf = lfdata %>% group_by(Population, Temp.Treatment) %>%
  dplyr::summarise(across(c(Larval.Survival, LarvalDevRate, Pupal.Survival, PupalDevRate, 
                     AdultLifespan, EggCount, SexWeightedFitness, PosFitness), mean, na.rm = T)) %>%
  as.data.frame()

# Convert NaNs to NAs (these are cases where no individuals survived or developed for that temp/pop)
lfdf$Larval.Survival[is.nan(lfdf$Larval.Survival)] <- NA
lfdf$LarvalDevRate[is.nan(lfdf$LarvalDevRate)] <- NA
lfdf$Pupal.Survival[is.nan(lfdf$Pupal.Survival)] <- NA
lfdf$PupalDevRate[is.nan(lfdf$PupalDevRate)] <- NA
lfdf$AdultLifespan[is.nan(lfdf$AdultLifespan)] <- NA
lfdf$EggCount[is.nan(lfdf$EggCount)] <- NA

colnames(lfdf)[c(9,10)] = c("PerIndFit", "PerPosFit")

#### Output data for sensitivity analysis ####
write.csv(lfdf, "~/Documents/Current_Projects/LifeHistoryTraitExp/Fitness_Breakdown/Fitness_SensitivityAnalysis.csv")

#### Provide individual-level trait data for Desire ####
lfdataind = lfdata[,c(1,2,3,5,9,13,20,21,24,25,27,28)]
write.csv(lfdataind, "~/Documents/Current_Projects/LifeHistoryTraitExp/Fitness_Breakdown/Fitness_SensitivityAnalysis_IndividualLevel.csv")


#### Pre-Processing for reproduction experiment ####

# Bring in data 
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/")
data2 = read.csv("Analysis_TraitFits/ReproExp_Data_PreProcessing.csv", header = T, na.strings = c(" ", "NA"))

# Convert date to date class
data2$Bloodfed.Date = as.Date(data2$Bloodfed.Date, format = "%m/%d/%y")
data2$Egg.Date = as.Date(data2$Egg.Date, format = "%m/%d/%y")
data2$Death.Date = as.Date(data2$Death.Date, format = "%m/%d/%y")

# Remove any white space in population name 
# this was causing issues later
library(stringr)
data2$Population = str_trim(data2$Population, "right")  

# Correct the temperature treatment values 
# based on temperatures recorded by iButton during experiment
# F = 32.18, E = 27.5, D = 23.75, C = 17.12, B = 13.35, A = 5.49
data2$Temperature = dplyr::recode(data2$Temperature, '28' = 27.5, '24' = 23.75, '18' = 17.12, '14' = 13.35, '4' = 5.49)

# Calculate length gonotrophic cycle and EFD
# note gonotrophic cycle = 1/biting rate
data2$CycleLength = as.numeric(data2$Egg.Date - data2$Bloodfed.Date) # length gonotrophic cycle
data2$EFD = data2$Eggs / data2$CycleLength

#### Output cleaned data for repro experiment ####
#write.csv(data2, "Analysis_TraitFits/ReproExp_Data.csv")

#### Repro Preliminary analysis : Estimate wing length - fecundity relationship #####

# from Washburn et al. 1989: eggs in first clutch = 51.33 x female wing length (mm) - 87.96)
data2wl = data2[!is.na(data2$Winglength),] # remove rows with no wing length data
data2wl = data2wl[data2wl$Eggs != 0,] # or with 0 eggs

# which temperature treatment to use? all?
temp14 = data2wl[data2wl$Temperature == "13.35",]
summary(lm(temp14$Eggs ~ temp14$Winglength))

temp17 = data2wl[data2wl$Temperature == "17.12",]
summary(lm(temp17$Eggs ~ temp17$Winglength))

temp24 = data2wl[data2wl$Temperature == "23.75",]
summary(lm(temp24$Eggs ~ temp24$Winglength))

temp28 = data2wl[data2wl$Temperature == "27.5",]
summary(lm(temp28$Eggs ~ temp28$Winglength))

# Relationship only looks good for 24 and 14 C treatments
par(mar = c(4.4,4.7,1,1))
plot(temp24$Eggs ~ temp24$Winglength, ylim = c(0,250), xlim = c(1.75, 3.75),
     pch = 16, cex = 2.5, xlab = "Female Winglength (mm)", 
     ylab = "No. of Eggs in First Clutch", axes = FALSE,
     cex.lab = 1.8)
axis(1, at=seq(1.75, 3.75, 0.5), cex.axis = 1.4)
axis(2, at=seq(0, 250, 50), cex.axis = 1.4)
abline(lm(temp24$Eggs ~ temp24$Winglength), lty = 2, lwd = 2)
abline(a = -87.96, b = 51.33, col = "red", lty = 1, lwd = 2)
box(col = "black")

plot(temp14$Eggs ~ temp14$Winglength, ylim = c(0,250), xlim = c(1.75, 3.75),
     pch = 16, cex = 2.5, xlab = "Female Winglength (mm)", 
     ylab = "No. of Eggs in First Clutch", axes = FALSE,
     cex.lab = 1.8)
axis(1, at=seq(1.75, 3.75, 0.5), cex.axis = 1.4)
axis(2, at=seq(0, 250, 50), cex.axis = 1.4)
abline(lm(temp14$Eggs ~ temp14$Winglength), lty = 2, lwd = 2)
abline(a = -87.96, b = 51.33, col = "red", lty = 1, lwd = 2)
box(col = "black")


#### Climate data Pre-processing ####

setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/TemperatureData_PRISM/")

# Pull in climate data (calculated in Google Earth Engine) 

T1 = read.csv("Treehole_AnnualAvgTemp_2000_end2021.csv", header = T)
# remove 2021 since some treeholes weren't yet sampled
T1 = subset(T1, select = -Year2021)
# calculate average across years per treehole
T1$avgAnnualMean = rowMeans(T1[,-1])

T2a = read.csv("Treehole_AvgWarmestMonth_2000_end2021.csv", header = T)
T2b = read.csv("Treehole_AvgColdestMonth_2000_end2021.csv", header = T)
# calculate difference, then average
# remove 2021 since some treeholes weren't yet sampled
T2a = subset(T2a, select = -Year2021)
T2b = subset(T2b, select = -Year2021)

T2 = T2a[,c(2:22)] - T2b[,c(2:22)]
T2$avgSeasonalVar = rowMeans(T2)

T3 = read.csv("Treehole_Extreme_2000_end2021.csv", header = T)
# remove 2021 since some treeholes not yet sampled
T3 = subset(T3, select = -Year2021)
# calculate average across years per treehole
T3$avgExtreme = rowMeans(T3[,-1])

T4 = read.csv("Treehole_Exceed35_2000_end2021.csv", header = T)
# remove 2021 since some treeholes not yet sampled
T4 = subset(T4, select = -Year2021)
# calculate average across years per treehole
T4$avgNumDays = rowMeans(T4[,-1])

T5 = read.csv("Treehole_WetMean_2000_end2021.csv", header = T)
# remove 2021 since some treeholes not yet sampled
T5 = subset(T5, select = -Year2021)
# calculate average across years per treehole
T5$avgWetTemp = rowMeans(T5[,-1])

T6 = read.csv("Treehole_DaysExceed31.6_2000_end2021.csv", header = T)
# remove 2021 since some treeholes not yet sampled
T6 = subset(T6, select = -Year2021)
# calculate average across years per treehole
T6$Days31.6 = rowMeans(T6[,-1])

# Create df
TempVars = cbind.data.frame(T1[,c(1,23)], T2$avgSeasonalVar, 
                            T3$avgExtreme, T4$avgNumDays, T5$avgWetTemp,
                            T6$Days31.6)
colnames(TempVars)[c(3,4,5,6,7)] = c("avgSeasonalVar", "avgExtreme", 
                                   "avgNumDays", "avgWetTemp", "days31.6")

#### Output cleaned climate data ####
#write.csv(TempVars, "TempVarsAtTreeholes.csv")
