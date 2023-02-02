###### Life History Trait Experiment Analysis PreProcessing & Prelim Analysis ######

# Code overview:
# 1. Compile data frames of TPC parameters for each trait 
# 2. Compile bioclim data for each population (Treehole)

##### 1. Load workspaces and packages #####

setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")
load("pLS_meansd.Rsave")
load("pPS_meansd.Rsave")
load("LDR_meansd.Rsave")
load("PDR_meansd.Rsave")
load("AL_meansd.Rsave")
#load("EFD_meansd.Rsave")
#load("BR_meansd.Rsave")
#load("pEA_meansd.Rsave")
#load("M_meansd.Rsave")

pLSdf$Trait = "Larval.Survival"
pPSdf$Trait = "Pupal.Survival"
LDRdf$Trait = "LarvalDevRate"
PDRdf$Trait = "PupalDevRate"
ALdf$Trait = "AdultLifespan"
#EFDdf$Trait = "EFD"
#BRdf$Trait = "BitingRate"
#pEAdf$Trait = "pEA"
#Mdf$Trait = "MosquitoDensity"

df_list <- list(pLSdf, pPSdf, LDRdf, PDRdf, ALdf) #  can later add:  EFDdf, BRdf, pEAdf, Mdf)
df = Reduce(function(x, y) merge(x, y, all=TRUE), df_list)  

# order dataframe based on trait
# order from coldest to hottest temp in wettest quarter
TraitOrder = c("Larval.Survival", "LarvalDevRate", "Pupal.Survival", "PupalDevRate",
               "AdultLifespan")  # can later add:  "EFD", "BitingRate", "pEA", "MosquitoDensity"
df$Trait = factor(df$Trait, levels = TraitOrder)
df = df[order(df$Trait),]

# concatenate parameter and trait columns into for ease of analysis
df$TraitParam = paste(df$Trait, df$Param)

# write.csv(df, "TPCparameters_AllTraits.csv")

##### 2. Merge with bioclim data #####

bioclim = read.csv("~/Documents/Current_Projects/LifeHistoryTraitExp/BioclimVars_ForTreeholes.csv", header = T)

# remove unneeded columns from df
df = df[,-c(2,4,5)]
# merge (by population) with TPC parameter data 
data = merge(df, bioclim, by = "Population", all.x = TRUE)

data = data[order(data$TraitParam),]

# output
write.csv(data, "TPCparameters_AllTraits_BioClimData.csv")

##### 3. Analyze trait-parameter & bioclim variable correlations ######

# Split data into individual trait-parameter combinations
dflist <- split(data , f = data$TraitParam)

# Function to calculate correlation between trait/parameter and bioclim variables
biocorr = function(x)
{x = x[,-c(1,3)]
datacor = cor(x[ , colnames(x) != "mean"],  x$mean)
return(datacor)}

traitparams = unique(data$TraitParam)

CorrsDf = as.data.frame(matrix(NA, nrow = 19, ncol = 25))
colnames(CorrsDf) = traitparams
rownames(CorrsDf) = colnames(data)[4:22]

for (i in 1:length(traitparams))
{CorrsDf[,i] = as.numeric(biocorr(dflist[[i]]))}

# sum (abs value) each row to see which bioclim variable has strongest correlations
head(CorrsDf)
CorrsDf2 = CorrsDf
CorrsDf2$Sums = rowSums(abs(CorrsDf2))

##### 4. Correlation plot ######
library(corrplot)
corrplot(as.matrix(CorrsDf), method = "number", 
         tl.col = "black", tl.srt = 45)


