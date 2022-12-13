# Trait Fits for Pupal Survival ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for Pupal survival pPS 
# 10 Aedes sierrensis populations with uniform priors 

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit pupal survival thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to pupal survival prior thermal responses
#  5) Fit oupal survival thermal responses with data-informed priors
#  6) Plot parameters

###### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.pPS = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]

# remove all with missing data (i.e., individuals that died as larvae or in diapause)
data.pPS = data.pPS[data.pPS$Temp.Treatment != "5.49",]
data.pPS = data.pPS[!is.na(data.pPS$Pupal.Survival),] 

# Plot
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), data = data.pPS, ylab = "prob. pupal survival", xlab = "Temperature")

# only HOP, WAW, SB, POW had larvae survive the 32 C treatment 
# for all other populations, add a '0' for pupal survival at 32 
# (i.e., pupae couldn't survive if no larvae survived)

ManualMar1 = c(NA,"MARIN35", 32.18, rep(NA,5),0, rep(NA,17))
ManualMar2 = c(NA,"MARIN29", 32.18, rep(NA,5),0, rep(NA,17))
ManualEug = c(NA,"EUG", 32.18, rep(NA,5),0, rep(NA,17))
ManualPLA = c(NA,"PLA", 32.18, rep(NA,5),0, rep(NA,17))
ManualJRA = c(NA,"JRA", 32.18, rep(NA,5),0, rep(NA,17))
ManualPAR = c(NA,"PAR", 32.18, rep(NA,5),0, rep(NA,17))

data.pPS = rbind(data.pPS, ManualMar1, ManualMar2, ManualEug, 
                 ManualPLA, ManualJRA, ManualPAR)

data.pPS$Temp.Treatment = as.numeric(data.pPS$Temp.Treatment)
data.pPS$Pupal.Survival = as.integer(data.pPS$Pupal.Survival)


# Subset data by population
data.pPS.HOP <- subset(data.pPS, Population == "HOP")
data.pPS.MAR1 <- subset(data.pPS, Population == "MARIN35")
data.pPS.MAR2 <- subset(data.pPS, Population == "MARIN29")
data.pPS.WAW <- subset(data.pPS, Population == "WAW")
data.pPS.EUG <- subset(data.pPS, Population == "EUG")
data.pPS.PLA <- subset(data.pPS, Population == "PLA")
data.pPS.SB <- subset(data.pPS, Population == "SB")
data.pPS.JRA <- subset(data.pPS, Population == "JRA")
data.pPS.PAR <- subset(data.pPS, Population == "PAR")
data.pPS.POW <- subset(data.pPS, Population == "POW")


##### 2a. JAGS Models ####

# Quadratic Model with uniform priors and binomial likelihood model :

sink("quad.binom.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 42)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
      for(i in 1:N.obs)
      {trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
      trait[i] ~ dbin(ifelse(trait.mu[i] > 0.99, 0.99, trait.mu[i]), 1)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- ifelse(-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) > 1, 0.99, 
    -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]))
    }
    
    } # close model 
    ",fill=TRUE)
sink()


##### 2b. Shared settings for all models ####

# inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

# Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

# MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

#  Temp sequence for derived quantity calculations
# For actual fits
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(5, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##### 3. Fit using uniform priors #####
##### 3a. Hopland #####

# Set data
data <- data.pPS.HOP

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom1.txt
pPS.HOP.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.HOP.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.HOP, ylab = "pPS for Hopland", xlab = "Temperature")
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##### 3b. Marin35 (MAR1) ######

# Set data
data <- data.pPS.MAR1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom1.txt
pPS.MAR1.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.MAR1.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.MAR1, ylab = "pPS for MAR1land", xlab = "Temperature")
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3c. Marin29 (MAR2) #####

# Set data
data <- data.pPS.MAR2

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.MAR2.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.MAR2.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.MAR2, ylab = "pPS for MAR2land", xlab = "Temperature")
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3d. Wawona (WAW) #####

# Set data
data <- data.pPS.WAW

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.WAW.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.WAW.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.WAW, ylab = "pPS for WAW", xlab = "Temperature")
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3e. Eugene (EUG) #####

# Set data
data <- data.pPS.EUG

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.EUG.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.EUG.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.EUG, ylab = "pPS for EUG", xlab = "Temperature")
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3f. Placer (PLA) ######
# Set data
data <- data.pPS.PLA

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.PLA.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.PLA.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.PLA, ylab = "pPS for PLA", xlab = "Temperature")
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3g. Santa Barbara (SB) #####
# Set data
data <- data.pPS.SB

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.SB.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.SB.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.SB, ylab = "pPS for SB", xlab = "Temperature")
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3h. Jasper Ridge (JRA) #####

# Set data
data <- data.pPS.JRA

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.JRA.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.JRA.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.JRA, ylab = "pPS for JRA", xlab = "Temperature")
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)






##### 3i. Paso Robles (PAR) #####

# Set data
data <- data.pPS.PAR

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.PAR.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.PAR.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.PAR, ylab = "pPS for PAR", xlab = "Temperature")
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3j. Powder canyon (POW) #####
# Set data
data <- data.pPS.POW

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pPS.POW.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pPS.POW.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.POW, ylab = "pPS for POW", xlab = "Temperature")
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### For internal checking: 10 panel plot of fits using uniform priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.HOP, ylab = "pPS for Hop", xlab = "Temperature")
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.MAR1, ylab = "pPS for MAR1", xlab = "Temperature")
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.MAR2 , ylab = "pPS for MAR2 prior", xlab = "Temperature")
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.WAW , ylab = "pPS for WAW prior", xlab = "Temperature")
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.EUG , ylab = "pPS for EUG prior", xlab = "Temperature")
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.PLA , ylab = "pPS for PLA prior", xlab = "Temperature")
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.SB , ylab = "pPS for SB prior", xlab = "Temperature")
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.JRA , ylab = "pPS for JRA prior", xlab = "Temperature")
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.PAR , ylab = "pPS for PAR prior", xlab = "Temperature")
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pPS.POW , ylab = "pPS for POW prior", xlab = "Temperature")
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### 10 panel plots of fits with uniform priors, but raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$Pupal.Survival ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$Pupal.Survival ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$Pupal.Survival ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/pPS_UniformPriors.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

pPS.Hop.MeansSds = MeanSd(data.pPS.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.Hop.MeansSds, ylab = "pPS for Hop", xlab = "Temperature", pch= 16)
arrows(pPS.Hop.MeansSds$Temp.Treatment, pPS.Hop.MeansSds$mean-pPS.Hop.MeansSds$stderror, pPS.Hop.MeansSds$Temp.Treatment, pPS.Hop.MeansSds$mean + pPS.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.MAR1.MeansSds = MeanSd(data.pPS.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.MAR1.MeansSds, ylab = "pPS for MAR1", xlab = "Temperature", pch= 16)
arrows(pPS.MAR1.MeansSds$Temp.Treatment, pPS.MAR1.MeansSds$mean-pPS.MAR1.MeansSds$stderror, pPS.MAR1.MeansSds$Temp.Treatment, pPS.MAR1.MeansSds$mean + pPS.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.MAR2.MeansSds = MeanSd(data.pPS.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.MAR2.MeansSds, ylab = "pPS for MAR2", xlab = "Temperature", pch= 16)
arrows(pPS.MAR2.MeansSds$Temp.Treatment, pPS.MAR2.MeansSds$mean-pPS.MAR2.MeansSds$stderror, pPS.MAR2.MeansSds$Temp.Treatment, pPS.MAR2.MeansSds$mean + pPS.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.WAW.MeansSds = MeanSd(data.pPS.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.WAW.MeansSds, ylab = "pPS for WAW", xlab = "Temperature", pch= 16)
arrows(pPS.WAW.MeansSds$Temp.Treatment, pPS.WAW.MeansSds$mean-pPS.WAW.MeansSds$stderror, pPS.WAW.MeansSds$Temp.Treatment, pPS.WAW.MeansSds$mean + pPS.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.EUG.MeansSds = MeanSd(data.pPS.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.EUG.MeansSds, ylab = "pPS for EUG", xlab = "Temperature", pch= 16)
arrows(pPS.EUG.MeansSds$Temp.Treatment, pPS.EUG.MeansSds$mean-pPS.EUG.MeansSds$stderror, pPS.EUG.MeansSds$Temp.Treatment, pPS.EUG.MeansSds$mean + pPS.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.PLA.MeansSds = MeanSd(data.pPS.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.PLA.MeansSds, ylab = "pPS for PLA", xlab = "Temperature", pch= 16)
arrows(pPS.PLA.MeansSds$Temp.Treatment, pPS.PLA.MeansSds$mean-pPS.PLA.MeansSds$stderror, pPS.PLA.MeansSds$Temp.Treatment, pPS.PLA.MeansSds$mean + pPS.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.SB.MeansSds = MeanSd(data.pPS.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.SB.MeansSds, ylab = "pPS for SB", xlab = "Temperature", pch= 16)
arrows(pPS.SB.MeansSds$Temp.Treatment, pPS.SB.MeansSds$mean-pPS.SB.MeansSds$stderror, pPS.SB.MeansSds$Temp.Treatment, pPS.SB.MeansSds$mean + pPS.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.JRA.MeansSds = MeanSd(data.pPS.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.JRA.MeansSds, ylab = "pPS for JRA", xlab = "Temperature", pch= 16)
arrows(pPS.JRA.MeansSds$Temp.Treatment, pPS.JRA.MeansSds$mean-pPS.JRA.MeansSds$stderror, pPS.JRA.MeansSds$Temp.Treatment, pPS.JRA.MeansSds$mean + pPS.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.PAR.MeansSds = MeanSd(data.pPS.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.PAR.MeansSds, ylab = "pPS for PAR", xlab = "Temperature", pch= 16)
arrows(pPS.PAR.MeansSds$Temp.Treatment, pPS.PAR.MeansSds$mean-pPS.PAR.MeansSds$stderror, pPS.PAR.MeansSds$Temp.Treatment, pPS.PAR.MeansSds$mean + pPS.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

pPS.POW.MeansSds = MeanSd(data.pPS.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.POW.MeansSds, ylab = "pPS for POW", xlab = "Temperature", pch= 16)
arrows(pPS.POW.MeansSds$Temp.Treatment, pPS.POW.MeansSds$mean-pPS.POW.MeansSds$stderror, pPS.POW.MeansSds$Temp.Treatment, pPS.POW.MeansSds$mean + pPS.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

dev.off()

###### 4. Get TPC parameters #####
###### 4a. Topt of pPS for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 4b. Pmax of pPS for each population ######

Pmax = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  PmaxVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {PmaxVec[i] = max(Matrix[i,6:86])}
  return(c(mean(PmaxVec), quantile(PmaxVec, c(0.025, 0.975))))
}

###### 4c. Tbreadth for each population #####

# Tbreadth defined as the temperature range where 
# performance remains >= 50% peak performance

Tbreadth = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  TbreadthVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {fiftypeak = (max(Matrix[i,6:86])/2)
  fiftypeakindexes = which(Matrix[i,6:86] >= fiftypeak)
  mintempindex = fiftypeakindexes[1]
  maxtempindex = fiftypeakindexes[length(fiftypeakindexes)]
  TbreadthVec[i] = Temp.xs[maxtempindex] - Temp.xs[mintempindex]
  }
  return(c(mean(TbreadthVec), quantile(TbreadthVec, c(0.025, 0.975))))
}


###### 4d. Compile prediction data for all populations #####
# T0, Tm, and Topt along with 95% credible intervals (2.5%, 97.5%)

# function to compile mean & 95% CrIs for all parameters for each population 
# for Tmin and Tmax, uses quantile method
# for Topt and Pmax - generates 95% CrIs using quantile method

library(HDInterval)

ParamCompile = function(x){
  DF = as.data.frame(matrix(,nrow = 5, ncol =4))
  colnames(DF) = c("mean", "lower95", "upper95", "param")
  DF[,4] = c("Tmin", "Tmax", "Topt", "Tbreadth", "Pmax")
  DF[1,1] = x$BUGSoutput$summary[1,1]
  DF[1,2:3] = hdi(x$BUGSoutput$sims.list$cf.T0, 0.95)[c(1,2)]
  DF[2,1] = x$BUGSoutput$summary[2,1]
  DF[2,2:3] = hdi(x$BUGSoutput$sims.list$cf.Tm, 0.95)[c(1,2)]
  DF[3,1:3] = as.numeric(Topt(x))
  DF[4,1:3] = as.numeric(Tbreadth(x))
  DF[5,1:3] = as.numeric(Pmax(x))
  return(DF)
}

pPSdf_uni <- rbind.data.frame(ParamCompile(pPS.HOP.out), ParamCompile(pPS.MAR1.out),
                              ParamCompile(pPS.MAR2.out), ParamCompile(pPS.WAW.out),
                              ParamCompile(pPS.EUG.out), ParamCompile(pPS.PLA.out),
                              ParamCompile(pPS.SB.out), ParamCompile(pPS.JRA.out),
                              ParamCompile(pPS.PAR.out), ParamCompile(pPS.POW.out))
pPSdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(pPSdf_uni, file = "pPS_meansd_uni.Rsave")


##### 5. Plotting  #####

library(ggplot2)
plot.new()
par(mfrow = c(1,1))
par(mar = c(2.5, 2.5, 2.5, 2.5))

load("pPS_meansd_uni.Rsave")

###### 5a. Plot of CTmax, CTmin, Topt for each population #####

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pPSdf_uni$Population = factor(pPSdf_uni$Population, levels = PopOrder)
pPSdf_uni = pPSdf_uni[order(pPSdf_uni$Population),]
pPSdf_uni = pPSdf_uni[pPSdf_uni$param != "Tbreadth",]
pPSdf_uni = pPSdf_uni[pPSdf_uni$param != "Pmax",]

# colors for each population 
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotpPS = ggplot(pPSdf_uni, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 35)) +
  geom_point(stat="identity", col=colors3, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors3,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("pPS") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5b. Plot of Pmax for each population #####

load("pPS_meansd_uni.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pPSdf_uni$Population = factor(pPSdf_uni$Population, levels = PopOrder)
pPSdf_uni = pPSdf_uni[order(pPSdf_uni$Population),]
pPSdf_uni = pPSdf_uni[pPSdf_uni$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotpPS_Pmax = ggplot(pPSdf_uni, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.7, 1)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("pPS") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5c. Plot pPS curves for all populations #####

load(file = "pPS_meansd_uni.Rsave")
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + pPS
plot(Pupal.Survival ~ Temp.Treatment, 
     xlim = c(5, 35), ylim = c(0,1), data = data.pPS.HOP, type = "n", bty = "n",
     ylab = "pPS", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "pPS", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(EUGdata ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = "#e0f3f8", lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = "#fee090", lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(POWdata ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(PARdata ~ Temp.xs, col = "#d73027", lwd = 1.5)
lines(SBdata ~ Temp.xs, col = "#a50026", lwd = 1.5)


###### 6. Generate informative priors: Leave-one-out with uniform priors ######
# Use same uniform priors settings as above
# Subset Data for leave-one-out
data.pPS.HOP.prior <- subset(data.pPS, Population != "HOP")
data.pPS.MAR1.prior <- subset(data.pPS, Population != "MARIN35")
data.pPS.MAR2.prior <- subset(data.pPS, Population != "MARIN29")
data.pPS.WAW.prior <- subset(data.pPS, Population != "WAW")
data.pPS.EUG.prior <- subset(data.pPS, Population != "EUG")
data.pPS.PLA.prior <- subset(data.pPS, Population != "PLA")
data.pPS.SB.prior <- subset(data.pPS, Population != "SB")
data.pPS.JRA.prior <- subset(data.pPS, Population != "JRA")
data.pPS.PAR.prior <- subset(data.pPS, Population != "PAR")
data.pPS.POW.prior <- subset(data.pPS, Population != "POW")


###### 6a. Hopland #####
# set data
data <- data.pPS.HOP.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.HOP.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.HOP.prior, ylab = "pPS for Hop prior", xlab = "Temperature")
lines(pPS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6b: Marin35 (MAR1) #####

# set data
data <- data.pPS.MAR1.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.MAR1.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.MAR1.prior, ylab = "pPS for MAR1 prior", xlab = "Temperature")
lines(pPS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.pPS.MAR2.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.MAR2.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.MAR2.prior, ylab = "pPS for MAR2 prior", xlab = "Temperature")
lines(pPS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.pPS.WAW.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.WAW.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.WAW.prior, ylab = "pPS for WAW prior", xlab = "Temperature")
lines(pPS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.pPS.EUG.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.EUG.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.EUG.prior, ylab = "pPS for EUG prior", xlab = "Temperature")
lines(pPS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.pPS.PLA.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.PLA.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.PLA.prior, ylab = "pPS for PLA prior", xlab = "Temperature")
lines(pPS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.pPS.SB.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.SB.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.SB.prior, ylab = "pPS for SB prior", xlab = "Temperature")
lines(pPS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.pPS.JRA.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.JRA.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.JRA.prior, ylab = "pPS for JRA prior", xlab = "Temperature")
lines(pPS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.pPS.PAR.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.PAR.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.PAR.prior, ylab = "pPS for PAR prior", xlab = "Temperature")
lines(pPS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyon (POW) ######

# set data
data <- data.pPS.POW.prior

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.POW.prior.out)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.POW.prior, ylab = "pPS for POW prior", xlab = "Temperature")
lines(pPS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### Save posterior distributions for informative priors #####

pPS.HOP.prior.cf.dists <- data.frame(q = as.vector(pPS.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.HOP.prior.gamma.fits = apply(pPS.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.MAR1.prior.cf.dists <- data.frame(q = as.vector(pPS.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pPS.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pPS.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.MAR1.prior.gamma.fits = apply(pPS.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.MAR2.prior.cf.dists <- data.frame(q = as.vector(pPS.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pPS.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pPS.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.MAR2.prior.gamma.fits = apply(pPS.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.WAW.prior.cf.dists <- data.frame(q = as.vector(pPS.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.WAW.prior.gamma.fits = apply(pPS.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.EUG.prior.cf.dists <- data.frame(q = as.vector(pPS.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.EUG.prior.gamma.fits = apply(pPS.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.PLA.prior.cf.dists <- data.frame(q = as.vector(pPS.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.PLA.prior.gamma.fits = apply(pPS.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.SB.prior.cf.dists <- data.frame(q = as.vector(pPS.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(pPS.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(pPS.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.SB.prior.gamma.fits = apply(pPS.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.JRA.prior.cf.dists <- data.frame(q = as.vector(pPS.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.JRA.prior.gamma.fits = apply(pPS.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.PAR.prior.cf.dists <- data.frame(q = as.vector(pPS.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.PAR.prior.gamma.fits = apply(pPS.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.POW.prior.cf.dists <- data.frame(q = as.vector(pPS.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.POW.prior.gamma.fits = apply(pPS.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pPS.hypers <- list(pPS.HOP.prior.gamma.fits, pPS.MAR1.prior.gamma.fits, pPS.MAR2.prior.gamma.fits,
                   pPS.WAW.prior.gamma.fits, pPS.EUG.prior.gamma.fits, pPS.PLA.prior.gamma.fits,
                   pPS.SB.prior.gamma.fits, pPS.JRA.prior.gamma.fits, pPS.PAR.prior.gamma.fits, pPS.POW.prior.gamma.fits)
save(pPS.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/pPShypers.Rsave")



###### Quadratic Model with gamma priors (except sigma) ######

sink("quad_inf.binom.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
      for(i in 1:N.obs)
      {trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
      trait[i] ~ dbin(ifelse(trait.mu[i] > 0.99, 0.99, trait.mu[i]), 1)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- ifelse(-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) > 1, 0.99, 
    -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]))
    }
    
    } # close model
    ",fill=T)
sink()



###### 7. Fit using informative priors #####
# Load hyperparameters to bypass above code of generating them 
load("pPShypers.Rsave")
pPS.HOP.prior.gamma.fits <- pPS.hypers[[1]]
pPS.MAR1.prior.gamma.fits <- pPS.hypers[[2]]
pPS.MAR2.prior.gamma.fits <- pPS.hypers[[3]]
pPS.WAW.prior.gamma.fits <- pPS.hypers[[4]]
pPS.EUG.prior.gamma.fits <- pPS.hypers[[5]]
pPS.PLA.prior.gamma.fits <- pPS.hypers[[6]]
pPS.SB.prior.gamma.fits <- pPS.hypers[[7]]
pPS.JRA.prior.gamma.fits <- pPS.hypers[[8]]
pPS.PAR.prior.gamma.fits <- pPS.hypers[[9]]
pPS.POW.prior.gamma.fits <- pPS.hypers[[10]]

###### 7a. Hopland ######

# Set data 
data <- data.pPS.HOP
hypers <- pPS.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.HOP.out.inf)

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.HOP, ylab = "pPS for Hopland", xlab = "Temperature", pch = 1)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7b. Marin35 ######

# Set data 
data <- data.pPS.MAR1
hypers <- pPS.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.MAR1, ylab = "pPS for Marin35", xlab = "Temperature", pch = 1)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.pPS.MAR2
hypers <- pPS.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.MAR2, ylab = "pPS for Marin29", xlab = "Temperature", pch = 1)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.pPS.WAW
hypers <- pPS.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.WAW, ylab = "pPS for Wawona", xlab = "Temperature", pch = 1)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.pPS.EUG
hypers <- pPS.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.EUG, ylab = "pPS for Eugene", xlab = "Temperature", pch = 1)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.pPS.PLA
hypers <- pPS.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.PLA, ylab = "pPS for Placer", xlab = "Temperature", pch = 1)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.pPS.SB
hypers <- pPS.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.SB, ylab = "pPS for Santa Barbara", xlab = "Temperature", pch = 1)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.pPS.JRA
hypers <- pPS.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.JRA, ylab = "pPS for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.pPS.PAR
hypers <- pPS.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.PAR, ylab = "pPS for Paso Robles", xlab = "Temperature", pch = 1)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon ######

# Set data 
data <- data.pPS.POW
hypers <- pPS.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pPS.POW, ylab = "pPS for Powder Canyon", xlab = "Temperature", pch = 1)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### For internal checking: 10 panel plot of fits using informative priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.HOP, ylab = "pPS for Hop", xlab = "Temperature")
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.MAR1, ylab = "pPS for MAR1", xlab = "Temperature")
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.MAR2 , ylab = "pPS for MAR2 prior", xlab = "Temperature")
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.WAW , ylab = "pPS for WAW prior", xlab = "Temperature")
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.EUG , ylab = "pPS for EUG prior", xlab = "Temperature")
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.PLA , ylab = "pPS for PLA prior", xlab = "Temperature")
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.SB , ylab = "pPS for SB prior", xlab = "Temperature")
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.JRA , ylab = "pPS for JRA prior", xlab = "Temperature")
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.PAR , ylab = "pPS for PAR prior", xlab = "Temperature")
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Pupal.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1), data = data.pPS.POW , ylab = "pPS for POW prior", xlab = "Temperature")
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### 10 panel plots of fits with raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$Pupal.Survival ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$Pupal.Survival ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$Pupal.Survival ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/pPS.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

pPS.Hop.MeansSds = MeanSd(data.pPS.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.Hop.MeansSds, ylab = "pPS for Hop", xlab = "Temperature", pch= 16)
arrows(pPS.Hop.MeansSds$Temp.Treatment, pPS.Hop.MeansSds$mean-pPS.Hop.MeansSds$stderror, pPS.Hop.MeansSds$Temp.Treatment, pPS.Hop.MeansSds$mean + pPS.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.MAR1.MeansSds = MeanSd(data.pPS.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.MAR1.MeansSds, ylab = "pPS for MAR1", xlab = "Temperature", pch= 16)
arrows(pPS.MAR1.MeansSds$Temp.Treatment, pPS.MAR1.MeansSds$mean-pPS.MAR1.MeansSds$stderror, pPS.MAR1.MeansSds$Temp.Treatment, pPS.MAR1.MeansSds$mean + pPS.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.MAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.MAR2.MeansSds = MeanSd(data.pPS.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.MAR2.MeansSds, ylab = "pPS for MAR2", xlab = "Temperature", pch= 16)
arrows(pPS.MAR2.MeansSds$Temp.Treatment, pPS.MAR2.MeansSds$mean-pPS.MAR2.MeansSds$stderror, pPS.MAR2.MeansSds$Temp.Treatment, pPS.MAR2.MeansSds$mean + pPS.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.WAW.MeansSds = MeanSd(data.pPS.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.WAW.MeansSds, ylab = "pPS for WAW", xlab = "Temperature", pch= 16)
arrows(pPS.WAW.MeansSds$Temp.Treatment, pPS.WAW.MeansSds$mean-pPS.WAW.MeansSds$stderror, pPS.WAW.MeansSds$Temp.Treatment, pPS.WAW.MeansSds$mean + pPS.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.EUG.MeansSds = MeanSd(data.pPS.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.EUG.MeansSds, ylab = "pPS for EUG", xlab = "Temperature", pch= 16)
arrows(pPS.EUG.MeansSds$Temp.Treatment, pPS.EUG.MeansSds$mean-pPS.EUG.MeansSds$stderror, pPS.EUG.MeansSds$Temp.Treatment, pPS.EUG.MeansSds$mean + pPS.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.PLA.MeansSds = MeanSd(data.pPS.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.PLA.MeansSds, ylab = "pPS for PLA", xlab = "Temperature", pch= 16)
arrows(pPS.PLA.MeansSds$Temp.Treatment, pPS.PLA.MeansSds$mean-pPS.PLA.MeansSds$stderror, pPS.PLA.MeansSds$Temp.Treatment, pPS.PLA.MeansSds$mean + pPS.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.SB.MeansSds = MeanSd(data.pPS.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.SB.MeansSds, ylab = "pPS for SB", xlab = "Temperature", pch= 16)
arrows(pPS.SB.MeansSds$Temp.Treatment, pPS.SB.MeansSds$mean-pPS.SB.MeansSds$stderror, pPS.SB.MeansSds$Temp.Treatment, pPS.SB.MeansSds$mean + pPS.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.JRA.MeansSds = MeanSd(data.pPS.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.JRA.MeansSds, ylab = "pPS for JRA", xlab = "Temperature", pch= 16)
arrows(pPS.JRA.MeansSds$Temp.Treatment, pPS.JRA.MeansSds$mean-pPS.JRA.MeansSds$stderror, pPS.JRA.MeansSds$Temp.Treatment, pPS.JRA.MeansSds$mean + pPS.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.PAR.MeansSds = MeanSd(data.pPS.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.PAR.MeansSds, ylab = "pPS for PAR", xlab = "Temperature", pch= 16)
arrows(pPS.PAR.MeansSds$Temp.Treatment, pPS.PAR.MeansSds$mean-pPS.PAR.MeansSds$stderror, pPS.PAR.MeansSds$Temp.Treatment, pPS.PAR.MeansSds$mean + pPS.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pPS.POW.MeansSds = MeanSd(data.pPS.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pPS.POW.MeansSds, ylab = "pPS for POW", xlab = "Temperature", pch= 16)
arrows(pPS.POW.MeansSds$Temp.Treatment, pPS.POW.MeansSds$mean-pPS.POW.MeansSds$stderror, pPS.POW.MeansSds$Temp.Treatment, pPS.POW.MeansSds$mean + pPS.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 36.5, y = 0.95, cex = 1.5)
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pPS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

dev.off()

###### 8. Get TPC parameters #####
###### 8a. Topt of pPS for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of pPS for each population ######

Pmax = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  PmaxVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {PmaxVec[i] = max(Matrix[i,6:86])}
  return(c(mean(PmaxVec), quantile(PmaxVec, c(0.025, 0.975))))
}

###### 8c. Tbreadth for each population #####

# Tbreadth defined as the temperature range where 
# performance remains >= 50% peak performance

Tbreadth = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  TbreadthVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {fiftypeak = (max(Matrix[i,6:86])/2)
  fiftypeakindexes = which(Matrix[i,6:86] >= fiftypeak)
  mintempindex = fiftypeakindexes[1]
  maxtempindex = fiftypeakindexes[length(fiftypeakindexes)]
  TbreadthVec[i] = Temp.xs[maxtempindex] - Temp.xs[mintempindex]
  }
  return(c(mean(TbreadthVec), quantile(TbreadthVec, c(0.025, 0.975))))
}


###### 8d. Compile prediction data for all populations #####
# T0, Tm, and Topt along with 95% credible intervals (2.5%, 97.5%)

# function to compile mean & 95% CrIs for all parameters for each population 
# for Tmin and Tmax, uses quantile method
# for Topt and Pmax - generates 95% CrIs using quantile method

ParamCompile = function(x){
  DF = as.data.frame(matrix(,nrow = 5, ncol =4))
  colnames(DF) = c("mean", "lower95", "upper95", "param")
  DF[,4] = c("Tmin", "Tmax", "Topt", "Tbreadth", "Pmax")
  DF[1,1] = x$BUGSoutput$summary[1,1]
  DF[1,2:3] = hdi(x$BUGSoutput$sims.list$cf.T0, 0.95)[c(1,2)]
  DF[2,1] = x$BUGSoutput$summary[2,1]
  DF[2,2:3] = hdi(x$BUGSoutput$sims.list$cf.Tm, 0.95)[c(1,2)]
  DF[3,1:3] = as.numeric(Topt(x))
  DF[4,1:3] = as.numeric(Tbreadth(x))
  DF[5,1:3] = as.numeric(Pmax(x))
  return(DF)
}

pPSdf_inf <- rbind.data.frame(ParamCompile(pPS.HOP.out.inf), ParamCompile(pPS.MAR1.out.inf),
                              ParamCompile(pPS.MAR2.out.inf), ParamCompile(pPS.WAW.out.inf),
                              ParamCompile(pPS.EUG.out.inf), ParamCompile(pPS.PLA.out.inf),
                              ParamCompile(pPS.SB.out.inf), ParamCompile(pPS.JRA.out.inf),
                              ParamCompile(pPS.PAR.out.inf), ParamCompile(pPS.POW.out.inf))
pPSdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(pPSdf_inf, file = "pPS_meansd_inf.Rsave")
###### 9. Plotting  #####

load("pPS_meansd_inf.Rsave")

library(ggplot2)
plot.new()
par(mfrow = c(1,1))
par(mar = c(2.5, 2.5, 2.5, 2.5))

###### 9a. Plot of CTmax, CTmin, Topt for each population #####
LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors = rep(rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fee090", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), each = 3)

pPSdf_inf$Population = factor(pPSdf_inf$Population, levels = LatOrder)
pPSdf_inf = pPSdf_inf[order(pPSdf_inf$Population),]
pPSdf_inf = pPSdf_inf[pPSdf_inf$param != "Tbreadth",]
pPSdf_inf = pPSdf_inf[pPSdf_inf$param != "Pmax",]

# colors for each population 
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotpPSinf = ggplot(pPSdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Pupal survival probability") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=12), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 

###### 9b. Plot of Pmax for each population #####

load("pPS_meansd_inf.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pPSdf_inf$Population = factor(pPSdf_inf$Population, levels = PopOrder)
pPSdf_inf = pPSdf_inf[order(pPSdf_inf$Population),]
pPSdf_inf = pPSdf_inf[pPSdf_inf$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotpPSinf_Pmax = ggplot(pPSdf_inf, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.7, 1)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("pPS (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 


###### 9c. Plot pPS curve for all populations #####

load(file = "pPS_meansd_inf.Rsave")
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(Pupal.Survival ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,1), data = data.pPS.HOP, type = "n", bty = "n",
     ylab = "probability", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Probability larval survival", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(EUGdata ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = "#e0f3f8", lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = "#fee090", lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(POWdata ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(PARdata ~ Temp.xs, col = "#d73027", lwd = 1.5)
lines(SBdata ~ Temp.xs, col = "#a50026", lwd = 1.5)


