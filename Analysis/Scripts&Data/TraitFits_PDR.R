# Trait Fits for Pupal Development Rate ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for pupal development rate (PDR) for 
# 10 Aedes sierrensis populations with uniform priors 

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit pupal dev rate thermal responses (Briere) using low info priors
#  4) Fit gamma distributions to PDR prior thermal responses
#  5) Fit pupal dev rate thermal responses with data-informed priors
#  6) Plot parameters


###### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.PDR = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# Keep only rows with PDR data
data.PDR = data.PDR[!is.na(data.PDR$PupalDevRate),]

# only WAW, SB, POW had pupae survive the 32 C treatment 
# for all other populations, add a '0' for development rate at 32 
# (i.e., individuals can't develop if they can't survive)

ManualHop = c(NA,"HOP", 32.18, rep(NA,17),0, rep(NA,5))
ManualMar1 = c(NA,"MARIN35", 32.18, rep(NA,17),0, rep(NA,5))
ManualMar2 = c(NA,"MARIN29", 32.18, rep(NA,17),0, rep(NA,5))
ManualEug = c(NA,"EUG", 32.18, rep(NA,17),0, rep(NA,5))
ManualPLA = c(NA,"PLA", 32.18, rep(NA,17),0, rep(NA,5))
ManualJRA = c(NA,"JRA", 32.18, rep(NA,17),0, rep(NA,5))
ManualPAR = c(NA,"PAR", 32.18, rep(NA,17),0, rep(NA,5))

data.PDR = rbind(data.PDR, ManualHop, ManualMar1, ManualMar2, ManualEug, 
                 ManualPLA, ManualJRA, ManualPAR)

data.PDR$Temp.Treatment = as.numeric(data.PDR$Temp.Treatment)
data.PDR$PupalDevRate = as.numeric(data.PDR$PupalDevRate)

# Subset data by population
data.PDR.HOP <- subset(data.PDR, Population == "HOP")
data.PDR.MAR1 <- subset(data.PDR, Population == "MARIN35")
data.PDR.MAR2 <- subset(data.PDR, Population == "MARIN29")
data.PDR.WAW <- subset(data.PDR, Population == "WAW")
data.PDR.EUG <- subset(data.PDR, Population == "EUG")
data.PDR.PLA <- subset(data.PDR, Population == "PLA")
data.PDR.SB <- subset(data.PDR, Population == "SB")
data.PDR.JRA <- subset(data.PDR, Population == "JRA")
data.PDR.PAR <- subset(data.PDR, Population == "PAR")
data.PDR.POW <- subset(data.PDR, Population == "POW")

# Plot trait data
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), pch = 16,
     cex = 1.5, cex.axis = 1.5, cex.lab = 1.5,
     data = data.PDR, ylab = "Pupal dev rate", xlab = "Temperature")

###### 2a. Set up JAGS model #####
# 2a. Briere Model with uniform priors #

sink("briere.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 39) # 39 set based on highest CTmax for pupal survival 
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

############## Briere Model with gamma priors (except sigma)

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



###### 2b. Shared settings for all models ####

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


###### 3. Fit using uniform priors #####
# uniform priors constrained over biologically realistic range
###### 3a: Hopland (HOP) ######

# set data
data <- data.PDR.HOP 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.HOP.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.HOP.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.HOP , ylab = "PDR for Hop prior", xlab = "Temperature")
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b: Marin35 (MAR1) ######

# set data
data <- data.PDR.MAR1 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.MAR1.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.MAR1.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.MAR1.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.MAR1 , ylab = "PDR for MAR1 prior", xlab = "Temperature")
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3c: Marin29 (MAR2) ######

# set data
data <- data.PDR.MAR2 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.MAR2.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.MAR2.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.MAR2.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.MAR2 , ylab = "PDR for MAR2 prior", xlab = "Temperature")
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3d: Wawona (WAW) ######

# set data
data <- data.PDR.WAW 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.WAW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.WAW.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.WAW.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.WAW , ylab = "PDR for WAW prior", xlab = "Temperature")
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.PDR.EUG 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.EUG.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.EUG.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.EUG.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.EUG , ylab = "PDR for EUG prior", xlab = "Temperature")
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.PDR.PLA 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.PLA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.PLA.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.PLA.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.PLA , ylab = "PDR for PLA prior", xlab = "Temperature")
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.PDR.SB 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.SB.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.SB.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.SB.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.SB , ylab = "PDR for SB prior", xlab = "Temperature")
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.PDR.JRA 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.JRA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.JRA.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.JRA.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.JRA , ylab = "PDR for JRA prior", xlab = "Temperature")
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.PDR.PAR 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.PAR.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.PAR.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.PAR.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.PAR , ylab = "PDR for PAR prior", xlab = "Temperature")
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

# set data
data <- data.PDR.POW 

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.POW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.POW.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.POW.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.POW , ylab = "PDR for POW prior", xlab = "Temperature")
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### For internal checking: 10 panel plot of fits using uniform priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.HOP, ylab = "PDR for Hop", xlab = "Temperature")
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.MAR1, ylab = "PDR for MAR1", xlab = "Temperature")
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.MAR2 , ylab = "PDR for MAR2 prior", xlab = "Temperature")
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.WAW , ylab = "PDR for WAW prior", xlab = "Temperature")
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.EUG , ylab = "PDR for EUG prior", xlab = "Temperature")
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.PLA , ylab = "PDR for PLA prior", xlab = "Temperature")
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.SB , ylab = "PDR for SB prior", xlab = "Temperature")
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.JRA , ylab = "PDR for JRA prior", xlab = "Temperature")
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.PAR , ylab = "PDR for PAR prior", xlab = "Temperature")
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.POW , ylab = "PDR for POW prior", xlab = "Temperature")
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### 10 panel plots of fits with uniform priors, but raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$PupalDevRate ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$PupalDevRate ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$PupalDevRate ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/PDR_UniformPriors.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

PDR.Hop.MeansSds = MeanSd(data.PDR.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.Hop.MeansSds, ylab = "PDR for Hop", xlab = "Temperature", pch= 16)
arrows(PDR.Hop.MeansSds$Temp.Treatment, PDR.Hop.MeansSds$mean-PDR.Hop.MeansSds$stderror, PDR.Hop.MeansSds$Temp.Treatment, PDR.Hop.MeansSds$mean + PDR.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.MAR1.MeansSds = MeanSd(data.PDR.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.MAR1.MeansSds, ylab = "PDR for MAR1", xlab = "Temperature", pch= 16)
arrows(PDR.MAR1.MeansSds$Temp.Treatment, PDR.MAR1.MeansSds$mean-PDR.MAR1.MeansSds$stderror, PDR.MAR1.MeansSds$Temp.Treatment, PDR.MAR1.MeansSds$mean + PDR.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.MAR2.MeansSds = MeanSd(data.PDR.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.MAR2.MeansSds, ylab = "PDR for MAR2", xlab = "Temperature", pch= 16)
arrows(PDR.MAR2.MeansSds$Temp.Treatment, PDR.MAR2.MeansSds$mean-PDR.MAR2.MeansSds$stderror, PDR.MAR2.MeansSds$Temp.Treatment, PDR.MAR2.MeansSds$mean + PDR.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.WAW.MeansSds = MeanSd(data.PDR.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.WAW.MeansSds, ylab = "PDR for WAW", xlab = "Temperature", pch= 16)
arrows(PDR.WAW.MeansSds$Temp.Treatment, PDR.WAW.MeansSds$mean-PDR.WAW.MeansSds$stderror, PDR.WAW.MeansSds$Temp.Treatment, PDR.WAW.MeansSds$mean + PDR.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.EUG.MeansSds = MeanSd(data.PDR.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.EUG.MeansSds, ylab = "PDR for EUG", xlab = "Temperature", pch= 16)
arrows(PDR.EUG.MeansSds$Temp.Treatment, PDR.EUG.MeansSds$mean-PDR.EUG.MeansSds$stderror, PDR.EUG.MeansSds$Temp.Treatment, PDR.EUG.MeansSds$mean + PDR.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.PLA.MeansSds = MeanSd(data.PDR.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.PLA.MeansSds, ylab = "PDR for PLA", xlab = "Temperature", pch= 16)
arrows(PDR.PLA.MeansSds$Temp.Treatment, PDR.PLA.MeansSds$mean-PDR.PLA.MeansSds$stderror, PDR.PLA.MeansSds$Temp.Treatment, PDR.PLA.MeansSds$mean + PDR.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.SB.MeansSds = MeanSd(data.PDR.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.SB.MeansSds, ylab = "PDR for SB", xlab = "Temperature", pch= 16)
arrows(PDR.SB.MeansSds$Temp.Treatment, PDR.SB.MeansSds$mean-PDR.SB.MeansSds$stderror, PDR.SB.MeansSds$Temp.Treatment, PDR.SB.MeansSds$mean + PDR.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.JRA.MeansSds = MeanSd(data.PDR.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.JRA.MeansSds, ylab = "PDR for JRA", xlab = "Temperature", pch= 16)
arrows(PDR.JRA.MeansSds$Temp.Treatment, PDR.JRA.MeansSds$mean-PDR.JRA.MeansSds$stderror, PDR.JRA.MeansSds$Temp.Treatment, PDR.JRA.MeansSds$mean + PDR.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.PAR.MeansSds = MeanSd(data.PDR.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.PAR.MeansSds, ylab = "PDR for PAR", xlab = "Temperature", pch= 16)
arrows(PDR.PAR.MeansSds$Temp.Treatment, PDR.PAR.MeansSds$mean-PDR.PAR.MeansSds$stderror, PDR.PAR.MeansSds$Temp.Treatment, PDR.PAR.MeansSds$mean + PDR.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.POW.MeansSds = MeanSd(data.PDR.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = PDR.POW.MeansSds, ylab = "PDR for POW", xlab = "Temperature", pch= 16)
arrows(PDR.POW.MeansSds$Temp.Treatment, PDR.POW.MeansSds$mean-PDR.POW.MeansSds$stderror, PDR.POW.MeansSds$Temp.Treatment, PDR.POW.MeansSds$mean + PDR.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

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

PDRdf_uni <- rbind.data.frame(ParamCompile(PDR.HOP.out), ParamCompile(PDR.MAR1.out),
                              ParamCompile(PDR.MAR2.out), ParamCompile(PDR.WAW.out),
                              ParamCompile(PDR.EUG.out), ParamCompile(PDR.PLA.out),
                              ParamCompile(PDR.SB.out), ParamCompile(PDR.JRA.out),
                              ParamCompile(PDR.PAR.out), ParamCompile(PDR.POW.out))
PDRdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(PDRdf_uni, file = "PDR_meansd_uni.Rsave")

###### 5. Plotting #####
load("PDR_meansd_uni.Rsave")
library(ggplot2)

###### 5a. Plot of CTmax, CTmin, Topt for each population #####

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
PDRdf_uni$Population = factor(PDRdf_uni$Population, levels = PopOrder)
PDRdf_uni = PDRdf_uni[order(PDRdf_uni$Population),]
PDRdf_uni = PDRdf_uni[PDRdf_uni$param != "Tbreadth",]
PDRdf_uni = PDRdf_uni[PDRdf_uni$param != "Pmax",]

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotPDR = ggplot(PDRdf_uni, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors3,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Pupal development rate") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5b. Plot of Pmax for each population #####

load("PDR_meansd_uni.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
PDRdf_uni$Population = factor(PDRdf_uni$Population, levels = PopOrder)
PDRdf_uni = PDRdf_uni[order(PDRdf_uni$Population),]
PDRdf_uni = PDRdf_uni[PDRdf_uni$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotPDR_Pmax = ggplot(PDRdf_uni, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.25, 0.35)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("PDR (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5c. Plot pupal dev rates curve for all populations #####

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, 
     xlim = c(0, 50), ylim = c(0,0.4), data = data.PDR.HOP, 
     type = "n",  bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Pupal development rate", cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)
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
data.PDR.HOP.prior <- subset(data.PDR, Population != "HOP")
data.PDR.MAR1.prior <- subset(data.PDR, Population != "MARIN35")
data.PDR.MAR2.prior <- subset(data.PDR, Population != "MARIN29")
data.PDR.WAW.prior <- subset(data.PDR, Population != "WAW")
data.PDR.EUG.prior <- subset(data.PDR, Population != "EUG")
data.PDR.PLA.prior <- subset(data.PDR, Population != "PLA")
data.PDR.SB.prior <- subset(data.PDR, Population != "SB")
data.PDR.JRA.prior <- subset(data.PDR, Population != "JRA")
data.PDR.PAR.prior <- subset(data.PDR, Population != "PAR")
data.PDR.POW.prior <- subset(data.PDR, Population != "POW")

###### 6a. Hopland #####
# set data
data <- data.PDR.HOP.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.HOP.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.HOP.prior, ylab = "PDR for Hop prior", xlab = "Temperature")
lines(PDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6b: Marin35 (MAR1) #####

# set data
data <- data.PDR.MAR1.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.MAR1.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.MAR1.prior, ylab = "PDR for MAR1 prior", xlab = "Temperature")
lines(PDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.PDR.MAR2.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.MAR2.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.MAR2.prior, ylab = "PDR for MAR2 prior", xlab = "Temperature")
lines(PDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.PDR.WAW.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.WAW.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.WAW.prior, ylab = "PDR for WAW prior", xlab = "Temperature")
lines(PDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.PDR.EUG.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.EUG.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.EUG.prior, ylab = "PDR for EUG prior", xlab = "Temperature")
lines(PDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.PDR.PLA.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.PLA.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.PLA.prior, ylab = "PDR for PLA prior", xlab = "Temperature")
lines(PDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.PDR.SB.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.SB.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.SB.prior, ylab = "PDR for SB prior", xlab = "Temperature")
lines(PDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.PDR.JRA.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.JRA.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.JRA.prior, ylab = "PDR for JRA prior", xlab = "Temperature")
lines(PDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.PDR.PAR.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.PAR.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.PAR.prior, ylab = "PDR for PAR prior", xlab = "Temperature")
lines(PDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyono (POW) ######

# set data
data <- data.PDR.POW.prior

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
PDR.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.POW.prior.out)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.POW.prior, ylab = "PDR for POW prior", xlab = "Temperature")
lines(PDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### Save posterior distributions for informative priors #####

PDR.HOP.prior.cf.dists <- data.frame(q = as.vector(PDR.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.HOP.prior.gamma.fits = apply(PDR.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.MAR1.prior.cf.dists <- data.frame(q = as.vector(PDR.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(PDR.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(PDR.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.MAR1.prior.gamma.fits = apply(PDR.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.MAR2.prior.cf.dists <- data.frame(q = as.vector(PDR.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(PDR.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(PDR.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.MAR2.prior.gamma.fits = apply(PDR.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.WAW.prior.cf.dists <- data.frame(q = as.vector(PDR.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.WAW.prior.gamma.fits = apply(PDR.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.EUG.prior.cf.dists <- data.frame(q = as.vector(PDR.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.EUG.prior.gamma.fits = apply(PDR.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.PLA.prior.cf.dists <- data.frame(q = as.vector(PDR.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.PLA.prior.gamma.fits = apply(PDR.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.SB.prior.cf.dists <- data.frame(q = as.vector(PDR.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(PDR.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(PDR.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.SB.prior.gamma.fits = apply(PDR.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.JRA.prior.cf.dists <- data.frame(q = as.vector(PDR.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.JRA.prior.gamma.fits = apply(PDR.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.PAR.prior.cf.dists <- data.frame(q = as.vector(PDR.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.PAR.prior.gamma.fits = apply(PDR.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.POW.prior.cf.dists <- data.frame(q = as.vector(PDR.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(PDR.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(PDR.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
PDR.POW.prior.gamma.fits = apply(PDR.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.hypers <- list(PDR.HOP.prior.gamma.fits, PDR.MAR1.prior.gamma.fits, PDR.MAR2.prior.gamma.fits,
                   PDR.WAW.prior.gamma.fits, PDR.EUG.prior.gamma.fits, PDR.PLA.prior.gamma.fits,
                   PDR.SB.prior.gamma.fits, PDR.JRA.prior.gamma.fits, PDR.PAR.prior.gamma.fits, PDR.POW.prior.gamma.fits)
save(PDR.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/PDRhypers.Rsave")

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
load("PDRhypers.Rsave")
PDR.HOP.prior.gamma.fits <- PDR.hypers[[1]]
PDR.MAR1.prior.gamma.fits <- PDR.hypers[[2]]
PDR.MAR2.prior.gamma.fits <- PDR.hypers[[3]]
PDR.WAW.prior.gamma.fits <- PDR.hypers[[4]]
PDR.EUG.prior.gamma.fits <- PDR.hypers[[5]]
PDR.PLA.prior.gamma.fits <- PDR.hypers[[6]]
PDR.SB.prior.gamma.fits <- PDR.hypers[[7]]
PDR.JRA.prior.gamma.fits <- PDR.hypers[[8]]
PDR.PAR.prior.gamma.fits <- PDR.hypers[[9]]
PDR.POW.prior.gamma.fits <- PDR.hypers[[10]]
###### Note: priors downweighted for this trait (0.01 vs 0.1 for other traits) 
###### 7a. Hopland ######

# Set data 
data <- data.PDR.HOP
hypers <- PDR.HOP.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(PDR.HOP.out.inf)

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.HOP, ylab = "PDR for Hopland", xlab = "Temperature", pch = 1)
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 7b. Marin35 ######

# Set data 
data <- data.PDR.MAR1
hypers <- PDR.MAR1.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.MAR1, ylab = "PDR for Marin35", xlab = "Temperature", pch = 1)
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.PDR.MAR2
hypers <- PDR.MAR2.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.MAR2, ylab = "PDR for Marin29", xlab = "Temperature", pch = 1)
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.PDR.WAW
hypers <- PDR.WAW.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.WAW, ylab = "PDR for Wawona", xlab = "Temperature", pch = 1)
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.PDR.EUG
hypers <- PDR.EUG.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.EUG, ylab = "PDR for Eugene", xlab = "Temperature", pch = 1)
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.PDR.PLA
hypers <- PDR.PLA.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.PLA, ylab = "PDR for Placer", xlab = "Temperature", pch = 1)
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.PDR.SB
hypers <- PDR.SB.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.SB, ylab = "PDR for Santa Barbara", xlab = "Temperature", pch = 1)
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.PDR.JRA
hypers <- PDR.JRA.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.JRA, ylab = "PDR for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.PDR.PAR
hypers <- PDR.PAR.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.PAR, ylab = "PDR for Paso Robles", xlab = "Temperature", pch = 1)
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon ######

# Set data 
data <- data.PDR.POW
hypers <- PDR.POW.prior.gamma.fits * 0.01

# Organize Data for JAGS
trait <- data$PupalDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
PDR.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
PDR.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.POW, ylab = "PDR for Powder Canyon", xlab = "Temperature", pch = 1)
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### For internal checking: 10 panel plot of fits using informative priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.HOP, ylab = "PDR for Hop", xlab = "Temperature")
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.MAR1, ylab = "PDR for MAR1", xlab = "Temperature")
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.MAR2 , ylab = "PDR for MAR2 prior", xlab = "Temperature")
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.WAW , ylab = "PDR for WAW prior", xlab = "Temperature")
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.EUG , ylab = "PDR for EUG prior", xlab = "Temperature")
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.PLA , ylab = "PDR for PLA prior", xlab = "Temperature")
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.SB , ylab = "PDR for SB prior", xlab = "Temperature")
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.JRA , ylab = "PDR for JRA prior", xlab = "Temperature")
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.PAR , ylab = "PDR for PAR prior", xlab = "Temperature")
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(PupalDevRate ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.6), data = data.PDR.POW , ylab = "PDR for POW prior", xlab = "Temperature")
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### 10 panel plots of fits with raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$PupalDevRate ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$PupalDevRate ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$PupalDevRate ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/PDR.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

PDR.Hop.MeansSds = MeanSd(data.PDR.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.Hop.MeansSds, ylab = "PDR for Hop", xlab = "Temperature", pch= 16)
arrows(PDR.Hop.MeansSds$Temp.Treatment, PDR.Hop.MeansSds$mean-PDR.Hop.MeansSds$stderror, PDR.Hop.MeansSds$Temp.Treatment, PDR.Hop.MeansSds$mean + PDR.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(PDR.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.MAR1.MeansSds = MeanSd(data.PDR.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.MAR1.MeansSds, ylab = "PDR for MAR1", xlab = "Temperature", pch= 16)
arrows(PDR.MAR1.MeansSds$Temp.Treatment, PDR.MAR1.MeansSds$mean-PDR.MAR1.MeansSds$stderror, PDR.MAR1.MeansSds$Temp.Treatment, PDR.MAR1.MeansSds$mean + PDR.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

PDR.MAR2.MeansSds = MeanSd(data.PDR.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.MAR2.MeansSds, ylab = "PDR for MAR2", xlab = "Temperature", pch= 16)
arrows(PDR.MAR2.MeansSds$Temp.Treatment, PDR.MAR2.MeansSds$mean-PDR.MAR2.MeansSds$stderror, PDR.MAR2.MeansSds$Temp.Treatment, PDR.MAR2.MeansSds$mean + PDR.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.WAW.MeansSds = MeanSd(data.PDR.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.WAW.MeansSds, ylab = "PDR for WAW", xlab = "Temperature", pch= 16)
arrows(PDR.WAW.MeansSds$Temp.Treatment, PDR.WAW.MeansSds$mean-PDR.WAW.MeansSds$stderror, PDR.WAW.MeansSds$Temp.Treatment, PDR.WAW.MeansSds$mean + PDR.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.EUG.MeansSds = MeanSd(data.PDR.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.EUG.MeansSds, ylab = "PDR for EUG", xlab = "Temperature", pch= 16)
arrows(PDR.EUG.MeansSds$Temp.Treatment, PDR.EUG.MeansSds$mean-PDR.EUG.MeansSds$stderror, PDR.EUG.MeansSds$Temp.Treatment, PDR.EUG.MeansSds$mean + PDR.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.PLA.MeansSds = MeanSd(data.PDR.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.PLA.MeansSds, ylab = "PDR for PLA", xlab = "Temperature", pch= 16)
arrows(PDR.PLA.MeansSds$Temp.Treatment, PDR.PLA.MeansSds$mean-PDR.PLA.MeansSds$stderror, PDR.PLA.MeansSds$Temp.Treatment, PDR.PLA.MeansSds$mean + PDR.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.SB.MeansSds = MeanSd(data.PDR.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.SB.MeansSds, ylab = "PDR for SB", xlab = "Temperature", pch= 16)
arrows(PDR.SB.MeansSds$Temp.Treatment, PDR.SB.MeansSds$mean-PDR.SB.MeansSds$stderror, PDR.SB.MeansSds$Temp.Treatment, PDR.SB.MeansSds$mean + PDR.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.JRA.MeansSds = MeanSd(data.PDR.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.JRA.MeansSds, ylab = "PDR for JRA", xlab = "Temperature", pch= 16)
arrows(PDR.JRA.MeansSds$Temp.Treatment, PDR.JRA.MeansSds$mean-PDR.JRA.MeansSds$stderror, PDR.JRA.MeansSds$Temp.Treatment, PDR.JRA.MeansSds$mean + PDR.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.PAR.MeansSds = MeanSd(data.PDR.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.PAR.MeansSds, ylab = "PDR for PAR", xlab = "Temperature", pch= 16)
arrows(PDR.PAR.MeansSds$Temp.Treatment, PDR.PAR.MeansSds$mean-PDR.PAR.MeansSds$stderror, PDR.PAR.MeansSds$Temp.Treatment, PDR.PAR.MeansSds$mean + PDR.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

PDR.POW.MeansSds = MeanSd(data.PDR.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,0.4), data = PDR.POW.MeansSds, ylab = "PDR for POW", xlab = "Temperature", pch= 16)
arrows(PDR.POW.MeansSds$Temp.Treatment, PDR.POW.MeansSds$mean-PDR.POW.MeansSds$stderror, PDR.POW.MeansSds$Temp.Treatment, PDR.POW.MeansSds$mean + PDR.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 36.5, y = 0.35, cex = 1.5)
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(PDR.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

dev.off()

###### 8. Get TPC parameters #####
###### 8a. Topt of PDR for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of PDR for each population ######

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

PDRdf_inf <- rbind.data.frame(ParamCompile(PDR.HOP.out.inf), ParamCompile(PDR.MAR1.out.inf),
                              ParamCompile(PDR.MAR2.out.inf), ParamCompile(PDR.WAW.out.inf),
                              ParamCompile(PDR.EUG.out.inf), ParamCompile(PDR.PLA.out.inf),
                              ParamCompile(PDR.SB.out.inf), ParamCompile(PDR.JRA.out.inf),
                              ParamCompile(PDR.PAR.out.inf), ParamCompile(PDR.POW.out.inf))
PDRdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(PDRdf_inf, file = "PDR_meansd_inf.Rsave")


###### 9. Plotting  #####
library(ggplot2)
load("PDR_meansd_inf.Rsave")

###### 9a. Plot of CTmax, CTmin, Topt for each population #####

LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors = rep(rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fee090", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")), each = 3)

PDRdf_inf$Population = factor(PDRdf_inf$Population, levels = LatOrder)
PDRdf_inf = PDRdf_inf[order(PDRdf_inf$Population),]
PDRdf_inf = PDRdf_inf[PDRdf_inf$param != "Tbreadth",]
PDRdf_inf = PDRdf_inf[PDRdf_inf$param != "Pmax",]

plotPDRinf = ggplot(PDRdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Pupal development rate") + 
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

load("PDR_meansd_inf.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
PDRdf_inf$Population = factor(PDRdf_inf$Population, levels = PopOrder)
PDRdf_inf = PDRdf_inf[order(PDRdf_inf$Population),]
PDRdf_inf = PDRdf_inf[PDRdf_inf$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotPDRinf_Pmax = ggplot(PDRdf_inf, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.250, 0.325)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("PDR (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 9c. Plot Pupal dev rates curve for all populations #####

load(file = "PDR_meansd_inf.Rsave")
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = PDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = PDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = PDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = PDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = PDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = PDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = PDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = PDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = PDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = PDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(PupalDevRate ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.HOP, type = "n", bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Pupal development rates", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
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



