# Trait Fits for Larval Development Rate #

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for larval development rate (LDR) for 
# 10 Aedes sierrensis populations with uniform prior 

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit larval dev rate thermal responses (Briere) using low info priors
#  4) Fit gamma distributions to LDR prior thermal responses
#  5) Fit larval dev rate thermal responses with data-informed priors
#  6) Plot parameters


###### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.LDR = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
data.LDR = data.LDR[!is.na(data.LDR$LarvalDevRate),]

# Plot trait data to see which function to fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.15), pch = 16,
     data = data.LDR, ylab = "Larval dev rate", xlab = "Temperature")

# only HOP, WAW, SB, POW had larvae survive the 32 C treatment 
# for all other populations, add a '0' for development rate at 32 
# (i.e., individuals can't develop if they can't survive)

ManualMar1 = c(NA,"MARIN35", 32.18, rep(NA,16),0, rep(NA,6))
ManualMar2 = c(NA,"MARIN29", 32.18, rep(NA,16),0, rep(NA,6))
ManualEug = c(NA,"EUG", 32.18, rep(NA,16),0, rep(NA,6))
ManualPLA = c(NA,"PLA", 32.18, rep(NA,16),0, rep(NA,6))
ManualJRA = c(NA,"JRA", 32.18, rep(NA,16),0, rep(NA,6))
ManualPAR = c(NA,"PAR", 32.18, rep(NA,16),0, rep(NA,6))

data.LDR = rbind(data.LDR, ManualMar1, ManualMar2, ManualEug, 
                 ManualPLA, ManualJRA, ManualPAR)

data.LDR$Temp.Treatment = as.numeric(data.LDR$Temp.Treatment)
data.LDR$LarvalDevRate = as.numeric(data.LDR$LarvalDevRate)

# Subset data by population
data.LDR.HOP <- subset(data.LDR, Population == "HOP")
data.LDR.MAR1 <- subset(data.LDR, Population == "MARIN35")
data.LDR.MAR2 <- subset(data.LDR, Population == "MARIN29")
data.LDR.WAW <- subset(data.LDR, Population == "WAW")
data.LDR.EUG <- subset(data.LDR, Population == "EUG")
data.LDR.PLA <- subset(data.LDR, Population == "PLA")
data.LDR.SB <- subset(data.LDR, Population == "SB")
data.LDR.JRA <- subset(data.LDR, Population == "JRA")
data.LDR.PAR <- subset(data.LDR, Population == "PAR")
data.LDR.POW <- subset(data.LDR, Population == "POW")

###### 2a. Set up JAGS model #####
# 2a. Briere Model with uniform priors #

sink("briere.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 35) # 35 set based on highest CTmax for larval survival
    cf.sigma ~ dunif(0, 1000) # std of data points
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau) # trait values drawn from normal dist with variance from ct tau
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()





###### 2b. Shared settings for all models ######

# Initial values function #
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

# Temp sequence for derived quantity calculations 
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(5, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


###### 3. Fit using uniform priors ####

# uniform priors, constrained over a biologically realistic range

###### 3a: Hopland (HOP) ######

# set data
data <- data.LDR.HOP

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.HOP.out$BUGSoutput$summary[1:5,]
mcmcplot(LDR.HOP.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

###### 3b: Marin35 (MAR1) ######

# set data
data <- data.LDR.MAR1 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.MAR1.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.MAR1.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.MAR1.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.MAR1 , ylab = "LDR for MAR1 prior", xlab = "Temperature")
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3c: Marin29 (MAR2) ######

# set data
data <- data.LDR.MAR2 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.MAR2.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.MAR2.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.MAR2.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.MAR2 , ylab = "LDR for MAR2 prior", xlab = "Temperature")
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3d: Wawona (WAW) ######

# set data
data <- data.LDR.WAW 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.WAW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.WAW.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.WAW.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.WAW , ylab = "LDR for WAW prior", xlab = "Temperature")
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.LDR.EUG 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.EUG.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.EUG.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.EUG.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.EUG , ylab = "LDR for EUG prior", xlab = "Temperature")
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.LDR.PLA 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.PLA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.PLA.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.PLA.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.PLA , ylab = "LDR for PLA prior", xlab = "Temperature")
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.LDR.SB 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.SB.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.SB.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.SB.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.SB , ylab = "LDR for SB prior", xlab = "Temperature")
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.LDR.JRA 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.JRA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.JRA.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.JRA.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.JRA , ylab = "LDR for JRA prior", xlab = "Temperature")
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.LDR.PAR 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.PAR.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.PAR.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.PAR.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.PAR , ylab = "LDR for PAR prior", xlab = "Temperature")
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

# set data
data <- data.LDR.POW 

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.POW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.POW.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.POW.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.POW , ylab = "LDR for POW prior", xlab = "Temperature")
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)






###### For internal checking: 10 panel plot of fits using uniform priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR1, ylab = "LDR for MAR1", xlab = "Temperature")
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR2 , ylab = "LDR for MAR2 prior", xlab = "Temperature")
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.WAW , ylab = "LDR for WAW prior", xlab = "Temperature")
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.EUG , ylab = "LDR for EUG prior", xlab = "Temperature")
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PLA , ylab = "LDR for PLA prior", xlab = "Temperature")
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.SB , ylab = "LDR for SB prior", xlab = "Temperature")
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.JRA , ylab = "LDR for JRA prior", xlab = "Temperature")
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PAR , ylab = "LDR for PAR prior", xlab = "Temperature")
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.POW , ylab = "LDR for POW prior", xlab = "Temperature")
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 10 panel plots of fits with uniform priors, but raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$LarvalDevRate ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$LarvalDevRate ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$LarvalDevRate ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/LDR_UniformPriors.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

LDR.Hop.MeansSds = MeanSd(data.LDR.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.Hop.MeansSds, ylab = "LDR for Hop", xlab = "Temperature", pch= 16)
arrows(LDR.Hop.MeansSds$Temp.Treatment, LDR.Hop.MeansSds$mean-LDR.Hop.MeansSds$stderror, LDR.Hop.MeansSds$Temp.Treatment, LDR.Hop.MeansSds$mean + LDR.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("HOP", x = 36.5, y = 0.95, cex = 1.5)

LDR.MAR1.MeansSds = MeanSd(data.LDR.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.MAR1.MeansSds, ylab = "LDR for MAR1", xlab = "Temperature", pch= 16)
arrows(LDR.MAR1.MeansSds$Temp.Treatment, LDR.MAR1.MeansSds$mean-LDR.MAR1.MeansSds$stderror, LDR.MAR1.MeansSds$Temp.Treatment, LDR.MAR1.MeansSds$mean + LDR.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("MAR1", x = 36.5, y = 0.95, cex = 1.5)

LDR.MAR2.MeansSds = MeanSd(data.LDR.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.MAR2.MeansSds, ylab = "LDR for MAR2", xlab = "Temperature", pch= 16)
arrows(LDR.MAR2.MeansSds$Temp.Treatment, LDR.MAR2.MeansSds$mean-LDR.MAR2.MeansSds$stderror, LDR.MAR2.MeansSds$Temp.Treatment, LDR.MAR2.MeansSds$mean + LDR.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("MAR2", x = 36.5, y = 0.95, cex = 1.5)


LDR.WAW.MeansSds = MeanSd(data.LDR.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.WAW.MeansSds, ylab = "LDR for WAW", xlab = "Temperature", pch= 16)
arrows(LDR.WAW.MeansSds$Temp.Treatment, LDR.WAW.MeansSds$mean-LDR.WAW.MeansSds$stderror, LDR.WAW.MeansSds$Temp.Treatment, LDR.WAW.MeansSds$mean + LDR.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("WAW", x = 36.5, y = 0.95, cex = 1.5)


LDR.EUG.MeansSds = MeanSd(data.LDR.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.EUG.MeansSds, ylab = "LDR for EUG", xlab = "Temperature", pch= 16)
arrows(LDR.EUG.MeansSds$Temp.Treatment, LDR.EUG.MeansSds$mean-LDR.EUG.MeansSds$stderror, LDR.EUG.MeansSds$Temp.Treatment, LDR.EUG.MeansSds$mean + LDR.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("EUG", x = 36.5, y = 0.95, cex = 1.5)


LDR.PLA.MeansSds = MeanSd(data.LDR.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.PLA.MeansSds, ylab = "LDR for PLA", xlab = "Temperature", pch= 16)
arrows(LDR.PLA.MeansSds$Temp.Treatment, LDR.PLA.MeansSds$mean-LDR.PLA.MeansSds$stderror, LDR.PLA.MeansSds$Temp.Treatment, LDR.PLA.MeansSds$mean + LDR.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("PLA", x = 36.5, y = 0.95, cex = 1.5)


LDR.SB.MeansSds = MeanSd(data.LDR.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.SB.MeansSds, ylab = "LDR for SB", xlab = "Temperature", pch= 16)
arrows(LDR.SB.MeansSds$Temp.Treatment, LDR.SB.MeansSds$mean-LDR.SB.MeansSds$stderror, LDR.SB.MeansSds$Temp.Treatment, LDR.SB.MeansSds$mean + LDR.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("SB", x = 36.5, y = 0.95, cex = 1.5)


LDR.JRA.MeansSds = MeanSd(data.LDR.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.JRA.MeansSds, ylab = "LDR for JRA", xlab = "Temperature", pch= 16)
arrows(LDR.JRA.MeansSds$Temp.Treatment, LDR.JRA.MeansSds$mean-LDR.JRA.MeansSds$stderror, LDR.JRA.MeansSds$Temp.Treatment, LDR.JRA.MeansSds$mean + LDR.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("JRA", x = 36.5, y = 0.95, cex = 1.5)


LDR.PAR.MeansSds = MeanSd(data.LDR.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.PAR.MeansSds, ylab = "LDR for PAR", xlab = "Temperature", pch= 16)
arrows(LDR.PAR.MeansSds$Temp.Treatment, LDR.PAR.MeansSds$mean-LDR.PAR.MeansSds$stderror, LDR.PAR.MeansSds$Temp.Treatment, LDR.PAR.MeansSds$mean + LDR.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("PAR", x = 36.5, y = 0.95, cex = 1.5)


LDR.POW.MeansSds = MeanSd(data.LDR.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.POW.MeansSds, ylab = "LDR for POW", xlab = "Temperature", pch= 16)
arrows(LDR.POW.MeansSds$Temp.Treatment, LDR.POW.MeansSds$mean-LDR.POW.MeansSds$stderror, LDR.POW.MeansSds$Temp.Treatment, LDR.POW.MeansSds$mean + LDR.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("POW", x = 36.5, y = 0.95, cex = 1.5)


dev.off()
###### 4. Get TPC parameters #####
###### 4a. Topt of fit for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 4b. Pmax of fit for each population ######

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

LDRdf_uni <- rbind.data.frame(ParamCompile(LDR.HOP.out), ParamCompile(LDR.MAR1.out),
                              ParamCompile(LDR.MAR2.out), ParamCompile(LDR.WAW.out),
                              ParamCompile(LDR.EUG.out), ParamCompile(LDR.PLA.out),
                              ParamCompile(LDR.SB.out), ParamCompile(LDR.JRA.out),
                              ParamCompile(LDR.PAR.out), ParamCompile(LDR.POW.out))
LDRdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(LDRdf_uni, file = "LDR_meansd_uni.Rsave")


###### 5. Plotting  #####
library(ggplot2)
load("LDR_meansd_uni.Rsave")
###### 5a. Plot of CTmax, CTmin, Topt for each population #####
# order by latitude
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

LDRdf_uni$Population = factor(LDRdf_uni$Population, levels = LatOrder)
LDRdf_uni = LDRdf_uni[order(LDRdf_uni$Population),]
LDRdf_uni = LDRdf_uni[LDRdf_uni$param != "Tbreadth",]
LDRdf_uni = LDRdf_uni[LDRdf_uni$param != "Pmax",]

plotLDRuni = ggplot(LDRdf_uni, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Larval development rate") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=14), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 

###### 5b. Plot of Pmax for each population #####

load("LDR_meansd_uni.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
LDRdf_uni$Population = factor(LDRdf_uni$Population, levels = PopOrder)
LDRdf_uni = LDRdf_uni[order(LDRdf_uni$Population),]
LDRdf_uni = LDRdf_uni[LDRdf_uni$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotLDR_Pmax = ggplot(LDRdf_uni, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.075, 0.125)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("LDR (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5c. Plot larval dev rates curve for all populations #####

load(file = "LDR_meansd_uni.Rsave")
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,0.11), data = data.LDR.HOP, type = "n", bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Larval development rates", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
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
data.LDR.HOP.prior <- subset(data.LDR, Population != "HOP")
data.LDR.MAR1.prior <- subset(data.LDR, Population != "MARIN35")
data.LDR.MAR2.prior <- subset(data.LDR, Population != "MARIN29")
data.LDR.WAW.prior <- subset(data.LDR, Population != "WAW")
data.LDR.EUG.prior <- subset(data.LDR, Population != "EUG")
data.LDR.PLA.prior <- subset(data.LDR, Population != "PLA")
data.LDR.SB.prior <- subset(data.LDR, Population != "SB")
data.LDR.JRA.prior <- subset(data.LDR, Population != "JRA")
data.LDR.PAR.prior <- subset(data.LDR, Population != "PAR")
data.LDR.POW.prior <- subset(data.LDR, Population != "POW")


###### 6a. Hopland #####
# set data
data <- data.LDR.HOP.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.HOP.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.HOP.prior, ylab = "LDR for Hop prior", xlab = "Temperature")
lines(LDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6b: Marin35 (MAR1) #####

# set data
data <- data.LDR.MAR1.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.MAR1.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.MAR1.prior, ylab = "LDR for MAR1 prior", xlab = "Temperature")
lines(LDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.LDR.MAR2.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.MAR2.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.MAR2.prior, ylab = "LDR for MAR2 prior", xlab = "Temperature")
lines(LDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.LDR.WAW.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.WAW.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.WAW.prior, ylab = "LDR for WAW prior", xlab = "Temperature")
lines(LDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.LDR.EUG.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.EUG.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.EUG.prior, ylab = "LDR for EUG prior", xlab = "Temperature")
lines(LDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.LDR.PLA.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.PLA.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.PLA.prior, ylab = "LDR for PLA prior", xlab = "Temperature")
lines(LDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.LDR.SB.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.SB.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.SB.prior, ylab = "LDR for SB prior", xlab = "Temperature")
lines(LDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.LDR.JRA.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.JRA.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.JRA.prior, ylab = "LDR for JRA prior", xlab = "Temperature")
lines(LDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.LDR.PAR.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.PAR.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.PAR.prior, ylab = "LDR for PAR prior", xlab = "Temperature")
lines(LDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyon (POW) ######

# set data
data <- data.LDR.POW.prior

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
LDR.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.POW.prior.out)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.POW.prior, ylab = "LDR for POW prior", xlab = "Temperature")
lines(LDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




###### Save posterior distributions for informative priors #####

LDR.HOP.prior.cf.dists <- data.frame(q = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.HOP.prior.gamma.fits = apply(LDR.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.MAR1.prior.cf.dists <- data.frame(q = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.MAR1.prior.gamma.fits = apply(LDR.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.MAR2.prior.cf.dists <- data.frame(q = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.MAR2.prior.gamma.fits = apply(LDR.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.WAW.prior.cf.dists <- data.frame(q = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.WAW.prior.gamma.fits = apply(LDR.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.EUG.prior.cf.dists <- data.frame(q = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.EUG.prior.gamma.fits = apply(LDR.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.PLA.prior.cf.dists <- data.frame(q = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.PLA.prior.gamma.fits = apply(LDR.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.SB.prior.cf.dists <- data.frame(q = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.SB.prior.gamma.fits = apply(LDR.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.JRA.prior.cf.dists <- data.frame(q = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.JRA.prior.gamma.fits = apply(LDR.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.PAR.prior.cf.dists <- data.frame(q = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.PAR.prior.gamma.fits = apply(LDR.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.POW.prior.cf.dists <- data.frame(q = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.POW.prior.gamma.fits = apply(LDR.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

LDR.hypers <- list(LDR.HOP.prior.gamma.fits, LDR.MAR1.prior.gamma.fits, LDR.MAR2.prior.gamma.fits,
                   LDR.WAW.prior.gamma.fits, LDR.EUG.prior.gamma.fits, LDR.PLA.prior.gamma.fits,
                   LDR.SB.prior.gamma.fits, LDR.JRA.prior.gamma.fits, LDR.PAR.prior.gamma.fits, LDR.POW.prior.gamma.fits)
save(LDR.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/LDRhypers.Rsave")


###### Briere Model with gamma priors (except sigma) ######

sink("briere_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

###### 7. Fit using informative priors #####
# Load hyperparameters to bypass above code of generating them 
load("LDRhypers.Rsave")
LDR.HOP.prior.gamma.fits <- LDR.hypers[[1]]
LDR.MAR1.prior.gamma.fits <- LDR.hypers[[2]]
LDR.MAR2.prior.gamma.fits <- LDR.hypers[[3]]
LDR.WAW.prior.gamma.fits <- LDR.hypers[[4]]
LDR.EUG.prior.gamma.fits <- LDR.hypers[[5]]
LDR.PLA.prior.gamma.fits <- LDR.hypers[[6]]
LDR.SB.prior.gamma.fits <- LDR.hypers[[7]]
LDR.JRA.prior.gamma.fits <- LDR.hypers[[8]]
LDR.PAR.prior.gamma.fits <- LDR.hypers[[9]]
LDR.POW.prior.gamma.fits <- LDR.hypers[[10]]




###### 7a. Hopland ######

# Set data 
data <- data.LDR.HOP
hypers <- LDR.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(LDR.HOP.out.inf)

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hopland", xlab = "Temperature", pch = 1)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 7b. Marin35 ######

# Set data 
data <- data.LDR.MAR1
hypers <- LDR.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.MAR1, ylab = "LDR for Marin35", xlab = "Temperature", pch = 1)
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.LDR.MAR2
hypers <- LDR.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.MAR2, ylab = "LDR for Marin29", xlab = "Temperature", pch = 1)
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.LDR.WAW
hypers <- LDR.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.WAW, ylab = "LDR for Wawona", xlab = "Temperature", pch = 1)
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.LDR.EUG
hypers <- LDR.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.EUG, ylab = "LDR for Eugene", xlab = "Temperature", pch = 1)
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.LDR.PLA
hypers <- LDR.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.PLA, ylab = "LDR for Placer", xlab = "Temperature", pch = 1)
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.LDR.SB
hypers <- LDR.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.SB, ylab = "LDR for Santa Barbara", xlab = "Temperature", pch = 1)
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.LDR.JRA
hypers <- LDR.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.JRA, ylab = "LDR for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.LDR.PAR
hypers <- LDR.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.PAR, ylab = "LDR for Paso Robles", xlab = "Temperature", pch = 1)
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon ######

# Set data 
data <- data.LDR.POW
hypers <- LDR.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$LarvalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
LDR.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
LDR.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.POW, ylab = "LDR for Powder Canyon", xlab = "Temperature", pch = 1)
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### For internal checking: 10 panel plot of fits using informative priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hop", xlab = "Temperature")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col ="red")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR1, ylab = "LDR for MAR1", xlab = "Temperature")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.MAR2 , ylab = "LDR for MAR2 prior", xlab = "Temperature")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")


plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.WAW , ylab = "LDR for WAW prior", xlab = "Temperature")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.EUG , ylab = "LDR for EUG prior", xlab = "Temperature")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PLA , ylab = "LDR for PLA prior", xlab = "Temperature")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.SB , ylab = "LDR for SB prior", xlab = "Temperature")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.JRA , ylab = "LDR for JRA prior", xlab = "Temperature")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.PAR , ylab = "LDR for PAR prior", xlab = "Temperature")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

plot(LarvalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = data.LDR.POW , ylab = "LDR for POW prior", xlab = "Temperature")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")


###### 10 panel plots of fits with raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$LarvalDevRate ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$LarvalDevRate ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$LarvalDevRate ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/LDR.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

LDR.Hop.MeansSds = MeanSd(data.LDR.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.Hop.MeansSds, ylab = "LDR for Hop", xlab = "Temperature", pch= 16)
arrows(LDR.Hop.MeansSds$Temp.Treatment, LDR.Hop.MeansSds$mean-LDR.Hop.MeansSds$stderror, LDR.Hop.MeansSds$Temp.Treatment, LDR.Hop.MeansSds$mean + LDR.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


LDR.MAR1.MeansSds = MeanSd(data.LDR.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.MAR1.MeansSds, ylab = "LDR for MAR1", xlab = "Temperature", pch= 16)
arrows(LDR.MAR1.MeansSds$Temp.Treatment, LDR.MAR1.MeansSds$mean-LDR.MAR1.MeansSds$stderror, LDR.MAR1.MeansSds$Temp.Treatment, LDR.MAR1.MeansSds$mean + LDR.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


LDR.MAR2.MeansSds = MeanSd(data.LDR.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.MAR2.MeansSds, ylab = "LDR for MAR2", xlab = "Temperature", pch= 16)
arrows(LDR.MAR2.MeansSds$Temp.Treatment, LDR.MAR2.MeansSds$mean-LDR.MAR2.MeansSds$stderror, LDR.MAR2.MeansSds$Temp.Treatment, LDR.MAR2.MeansSds$mean + LDR.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")


LDR.WAW.MeansSds = MeanSd(data.LDR.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.WAW.MeansSds, ylab = "LDR for WAW", xlab = "Temperature", pch= 16)
arrows(LDR.WAW.MeansSds$Temp.Treatment, LDR.WAW.MeansSds$mean-LDR.WAW.MeansSds$stderror, LDR.WAW.MeansSds$Temp.Treatment, LDR.WAW.MeansSds$mean + LDR.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")


LDR.EUG.MeansSds = MeanSd(data.LDR.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.EUG.MeansSds, ylab = "LDR for EUG", xlab = "Temperature", pch= 16)
arrows(LDR.EUG.MeansSds$Temp.Treatment, LDR.EUG.MeansSds$mean-LDR.EUG.MeansSds$stderror, LDR.EUG.MeansSds$Temp.Treatment, LDR.EUG.MeansSds$mean + LDR.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")


LDR.PLA.MeansSds = MeanSd(data.LDR.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.PLA.MeansSds, ylab = "LDR for PLA", xlab = "Temperature", pch= 16)
arrows(LDR.PLA.MeansSds$Temp.Treatment, LDR.PLA.MeansSds$mean-LDR.PLA.MeansSds$stderror, LDR.PLA.MeansSds$Temp.Treatment, LDR.PLA.MeansSds$mean + LDR.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")

LDR.SB.MeansSds = MeanSd(data.LDR.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.SB.MeansSds, ylab = "LDR for SB", xlab = "Temperature", pch= 16)
arrows(LDR.SB.MeansSds$Temp.Treatment, LDR.SB.MeansSds$mean-LDR.SB.MeansSds$stderror, LDR.SB.MeansSds$Temp.Treatment, LDR.SB.MeansSds$mean + LDR.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")


LDR.JRA.MeansSds = MeanSd(data.LDR.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.JRA.MeansSds, ylab = "LDR for JRA", xlab = "Temperature", pch= 16)
arrows(LDR.JRA.MeansSds$Temp.Treatment, LDR.JRA.MeansSds$mean-LDR.JRA.MeansSds$stderror, LDR.JRA.MeansSds$Temp.Treatment, LDR.JRA.MeansSds$mean + LDR.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")

LDR.PAR.MeansSds = MeanSd(data.LDR.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.PAR.MeansSds, ylab = "LDR for PAR", xlab = "Temperature", pch= 16)
arrows(LDR.PAR.MeansSds$Temp.Treatment, LDR.PAR.MeansSds$mean-LDR.PAR.MeansSds$stderror, LDR.PAR.MeansSds$Temp.Treatment, LDR.PAR.MeansSds$mean + LDR.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")


LDR.POW.MeansSds = MeanSd(data.LDR.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.12), data = LDR.POW.MeansSds, ylab = "LDR for POW", xlab = "Temperature", pch= 16)
arrows(LDR.POW.MeansSds$Temp.Treatment, LDR.POW.MeansSds$mean-LDR.POW.MeansSds$stderror, LDR.POW.MeansSds$Temp.Treatment, LDR.POW.MeansSds$mean + LDR.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 36.5, y = 0.105, cex = 1.5)
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(LDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")


dev.off()
###### 8. Get TPC parameters #####
###### 8a. Topt of LDR for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of LDR for each population ######

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

LDRdf_inf <- rbind.data.frame(ParamCompile(LDR.HOP.out.inf), ParamCompile(LDR.MAR1.out.inf),
                              ParamCompile(LDR.MAR2.out.inf), ParamCompile(LDR.WAW.out.inf),
                              ParamCompile(LDR.EUG.out.inf), ParamCompile(LDR.PLA.out.inf),
                              ParamCompile(LDR.SB.out.inf), ParamCompile(LDR.JRA.out.inf),
                              ParamCompile(LDR.PAR.out.inf), ParamCompile(LDR.POW.out.inf))
LDRdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(LDRdf_inf, file = "LDR_meansd_inf.Rsave")

###### 9. Plotting  #####

library(ggplot2)
load("LDR_meansd_inf.Rsave")

###### 9a. Plot of CTmax, CTmin, Topt for each population #####
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

LDRdf_inf$Population = factor(LDRdf_inf$Population, levels = LatOrder)
LDRdf_inf = LDRdf_inf[order(LDRdf_inf$Population),]
LDRdf_inf = LDRdf_inf[LDRdf_inf$param != "Tbreadth",]
LDRdf_inf = LDRdf_inf[LDRdf_inf$param != "Pmax",]

plotLDRinf = ggplot(LDRdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, size = 2, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors, size = 0.7,
                position=position_dodge(.9), width = 0.5) + 
  ggtitle("Larval development rate") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=14), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 

###### 9b. Plot of Pmax for each population #####

load("LDR_meansd_inf.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
LDRdf_inf$Population = factor(LDRdf_inf$Population, levels = PopOrder)
LDRdf_inf = LDRdf_inf[order(LDRdf_inf$Population),]
LDRdf_inf = LDRdf_inf[LDRdf_inf$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotLDRinf_Pmax = ggplot(LDRdf_inf, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.075, 0.125)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("LDR (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 


###### 9c. Plot larval dev rates curve for all populations #####

load(file = "LDR_meansd_inf.Rsave")

EUGdata = LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
par(mar = (c(4.5, 4.5, 2.5, 1)))

plot(LarvalDevRate ~ Temp.Treatment, 
     xlim = c(0, 45), ylim = c(0,0.11), data = data.LDR.HOP, type = "n", bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Larval development rates", cex.main = 1.7, cex.lab = 1.4, cex.axis = 1.2)
box(col = "black")
lines(EUGdata ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = "#7cae00", lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = "#fed439ff", lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(PARdata ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(SBdata ~ Temp.xs, col = "#ec3c30", lwd = 1.5)
lines(POWdata ~ Temp.xs, col = "#ab041b", lwd = 1.5)









