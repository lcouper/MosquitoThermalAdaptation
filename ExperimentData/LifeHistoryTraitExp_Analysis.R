#### Life History Trait Experiment Analysis #####

# create df 

df1 = c(fit.HOP.out.inf$BUGSoutput$sims.matrix[,2], fit.MAR1.out.inf$BUGSoutput$sims.matrix[,2],
fit.MAR2.out.inf$BUGSoutput$sims.matrix[,2], fit.WAW.out.inf$BUGSoutput$sims.matrix[,2],
fit.EUG.out.inf$BUGSoutput$sims.matrix[,2], fit.PLA.out.inf$BUGSoutput$sims.matrix[,2],
fit.SB.out.inf$BUGSoutput$sims.matrix[,2], fit.JRA.out.inf$BUGSoutput$sims.matrix[,2],
fit.PAR.out.inf$BUGSoutput$sims.matrix[,2], fit.POW.out.inf$BUGSoutput$sims.matrix[,2])

df2 = rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 7500)

df = cbind.data.frame(df1, df2)
colnames(df) = c("Tm", "Population")
fit = aov(Tm ~ Population, data = df)





#### 0. Bring in data ####

setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")
data = read.csv("LifeHistoryTraitExp_IndividualLevelTraits.csv", header = T)
data2 = read.csv("LifeHistoryTraitExp_PopulationLevelTraits.csv", header = T)
dataAll = read.csv("TPCparameters_AllTraits_BioClimData.csv", header = T)

#### 1. Between-population variation in performance: Individual level traits ####

# Larval dev rate
traitA = data[data$Trait.Name == "LarvalDevRate",]

fitA1 = lm(value ~ Temperature, data = traitA)
fitA2 = lm(value ~ Temperature + Population, data = traitA)
fitA3 = lm(value ~ Temperature + Population + Temperature*Population, data = traitA)
anova(fit1, fit2, fit3)

# Pupal dev rate 
traitB = data[data$Trait.Name == "PupalDevRate",]

fitB1 = lm(value ~ Temperature, data = traitB)
fitB2 = lm(value ~ Temperature + Population, data = traitB)
fitB3 = lm(value ~ Temperature + Population + Temperature*Population, data = traitB)
anova(fitB1, fitB2, fitB3)

# Adult life span
traitC = data[data$Trait.Name == "AdultLifespan",]

fitC1 = lm(value ~ Temperature, data = traitC)
fitC2 = lm(value ~ Temperature + Population, data = traitC)
fitC3 = lm(value ~ Temperature + Population + Temperature*Population, data = traitC)
anova(fitC1, fitC2, fitC3)
# no significant temp*population interaction here

#### 2. Between-population variation in performance: Population level traits ####

# but issue = there's now only one data point per population & treatment...

fitD1 = lm(prob.Larval.Survival ~ Temperature, data = data2)
fitD2 = lm(prob.Larval.Survival ~ Temperature + Population, data = data2)
fitD3 = lm(prob.Larval.Survival ~ Temperature + Population + Temperature*Population, data = data2)
anova(fitD1, fitD2, fitD3)


fitE1 = lm(prob.Pupal.Survival ~ Temperature, data = data2)
fitE2 = lm(prob.Pupal.Survival ~ Temperature + Population, data = data2)
fitE3 = lm(prob.Pupal.Survival ~ Temperature + Population + Temperature*Population, data = data2)
anova(fitE1, fitE2, fitE3)



#### 3. Test if parameters are significantly associated with bioclim variables #####

# Split data into individual trait-parameter combinations
dflist <- split(dataAll , f = dataAll$TraitParam)

# Function to calculate p-value for each trait/parameter & bioclim variable combo
biopvals = function(x)
{x = x[,-c(1,3)]
fitBio1 = lm(mean ~ Bio1, data = x)
fitBio2 = lm(mean ~ Bio2, data = x)
fitBio3 = lm(mean ~ Bio3, data = x)
fitBio4 = lm(mean ~ Bio4, data = x)
fitBio5 = lm(mean ~ Bio5, data = x)
fitBio6 = lm(mean ~ Bio6, data = x)
fitBio7 = lm(mean ~ Bio7, data = x)
fitBio8 = lm(mean ~ Bio8, data = x)
fitBio9 = lm(mean ~ Bio9, data = x)
fitBio10 = lm(mean ~ Bio10, data = x)
fitBio11 = lm(mean ~ Bio11, data = x)
fitBio12 = lm(mean ~ Bio12, data = x)
fitBio13 = lm(mean ~ Bio13, data = x)
fitBio14 = lm(mean ~ Bio14, data = x)
fitBio15 = lm(mean ~ Bio15, data = x)
fitBio16 = lm(mean ~ Bio16, data = x)
fitBio17  = lm(mean ~ Bio17, data = x)
fitBio18 = lm(mean ~ Bio18, data = x)
fitBio19  = lm(mean ~ Bio19, data = x)

pvalBio1 = anova(fitBio1)$'Pr(>F)'[1]
pvalBio2 = anova(fitBio2)$'Pr(>F)'[1]
pvalBio3 = anova(fitBio3)$'Pr(>F)'[1]
pvalBio4 = anova(fitBio4)$'Pr(>F)'[1]
pvalBio5 = anova(fitBio5)$'Pr(>F)'[1]
pvalBio6 = anova(fitBio6)$'Pr(>F)'[1]
pvalBio7 = anova(fitBio7)$'Pr(>F)'[1]
pvalBio8 = anova(fitBio8)$'Pr(>F)'[1]
pvalBio9 = anova(fitBio9)$'Pr(>F)'[1]
pvalBio10 = anova(fitBio10)$'Pr(>F)'[1]
pvalBio11 = anova(fitBio11)$'Pr(>F)'[1]
pvalBio12 = anova(fitBio12)$'Pr(>F)'[1]
pvalBio13 = anova(fitBio13)$'Pr(>F)'[1]
pvalBio14 = anova(fitBio14)$'Pr(>F)'[1]
pvalBio15 = anova(fitBio15)$'Pr(>F)'[1]
pvalBio16 = anova(fitBio16)$'Pr(>F)'[1]
pvalBio17 = anova(fitBio17)$'Pr(>F)'[1]
pvalBio18 = anova(fitBio18)$'Pr(>F)'[1]
pvalBio19 = anova(fitBio19)$'Pr(>F)'[1]

pvals = c(pvalBio1, pvalBio2, pvalBio3, pvalBio4, pvalBio5, pvalBio6, pvalBio7, pvalBio8, pvalBio9,
          pvalBio10, pvalBio11, pvalBio12, pvalBio13, pvalBio14, pvalBio15, pvalBio16, pvalBio17, 
          pvalBio18, pvalBio19)
print(pvals)
}


traitparams = unique(dataAll$TraitParam)

PvalsDf = as.data.frame(matrix(NA, nrow = 19, ncol = 35))
colnames(PvalsDf) = traitparams
rownames(PvalsDf) = colnames(dataAll)[4:22]

for (i in 1:length(traitparams))
{PvalsDf[,i] = as.numeric(biopvals(dflist[[i]]))}

# Create plot and highlight significant p values
library(corrplot)

# Remove EFD and biting rate for now
PvalsDf2 = PvalsDf[,-c(6:15)]

corrplot(as.matrix(PvalsDf2), method = "color", cl.lim=c(0,1),
         col = colorRampPalette(c("brown", "darkred", "red","pink", "white"))(100),
         tl.col = "black", tl.srt = 45)

# overlay corrplot highlighting significant cells 
ids = which(PvalsDf2 < 0.05)
Dup <- as.matrix(PvalsDf2)
Dup[-ids] <- NA

corrplot(as.matrix(Dup),  method = "number", col = "black",
          na.label= " ", addgrid.col = "black", add = TRUE, cl.pos = "n",
         number.cex = 0.7, bg = "transparent", tl.col = "transparent")
         





#### 4. Plot for NSF ORCC: Kelsey Knockdown data #####

kddata = read.csv("~/Downloads/Kelsey_knockdown_data_subset.csv", header = T)[1:63,]
popll = read.csv("~/Downloads/PopLatLongs.csv", header = T)[1:394,]
popll$SampleID = toupper(popll$SampleID)
kddata$SampleID = toupper(kddata$SampleID)

kdlats = merge(kddata, popll, by = "SampleID", all.x = T)
kdlats$Lat = round(as.numeric(kdlats$Lat), digits = 1)
kdlats$SampleID = as.factor(kdlats$SampleID)
PopOrder = c("EUG.03", "WAW.12", "PLA.01", "HOP.37", "JR.A",  
             "MARIN29", "MARIN35", "POW.02", "PAR.07", "SB.04")
kdlats$SampleID = factor(kdlats$SampleID, levels = PopOrder)
kdlats = kdlats[!is.na(kdlats$SampleID),]

colsNew = c("#882255", "#05263B", 
            "#43758A", "#88CCEE", "#999933", "#F1B416", "#E56A19", "#900000")

ggplot(kdlats, aes(Lat, Decimal)) + 
  geom_point(aes(colour='black', 
                 fill = factor(SampleID)), shape=21, size = 5.5) + 
  scale_fill_manual(values=colsNew) +
  scale_color_manual(values = 'black') + theme_minimal() +
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 



