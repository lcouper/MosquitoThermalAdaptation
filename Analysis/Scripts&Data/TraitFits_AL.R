# Trait Fits for Adult Lifespan #

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for adult lifespan (AL) for 
# 10 Aedes sierrensis populations with uniform priors 

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit adult life span thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to adult lifespan prior thermal responses
#  5) Fit adult lifespan thermal responses with data-informed priors
#  6) Plot parameters

###### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.AL = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# if lifespan is blank, it's because the individual did not survive to adulthood
# remove 5C treatment since no individuals reached pupal/adult life stage

data.AL = data.AL[data.AL$Temp.Treatment != 5.49,]

# only WAW, SB, POW had adults survive the 32 C treatment 
# for all other populations, add a '0' for lifespan at 32 
# (i.e., individuals can't develop if they can't survive)

ManualHop = c(NA,"HOP", 32.18, rep(NA,20),0, rep(NA,2))
ManualMar1 = c(NA,"MARIN35", 32.18, rep(NA,20),0, rep(NA,2))
ManualMar2 = c(NA,"MARIN29", 32.18, rep(NA,20),0, rep(NA,2))
ManualEug = c(NA,"EUG", 32.18, rep(NA,20),0, rep(NA,2))
ManualPLA = c(NA,"PLA", 32.18, rep(NA,20),0, rep(NA,2))
ManualJRA = c(NA,"JRA", 32.18, rep(NA,20),0, rep(NA,2))
ManualPAR = c(NA,"PAR", 32.18, rep(NA,20),0, rep(NA,2))

data.AL= rbind(data.AL, ManualHop, ManualMar1, ManualMar2, ManualEug, 
               ManualPLA, ManualJRA, ManualPAR)

data.AL$Temp.Treatment = as.numeric(data.AL$Temp.Treatment)
data.AL$AdultLifespan = as.numeric(data.AL$AdultLifespan)

# Subset data by population
data.AL.HOP <- subset(data.AL, Population == "HOP")
data.AL.MAR1 <- subset(data.AL, Population == "MARIN35")
data.AL.MAR2 <- subset(data.AL, Population == "MARIN29")
data.AL.WAW <- subset(data.AL, Population == "WAW")
data.AL.EUG <- subset(data.AL, Population == "EUG")
data.AL.PLA <- subset(data.AL, Population == "PLA")
data.AL.SB <- subset(data.AL, Population == "SB")
data.AL.JRA <- subset(data.AL, Population == "JRA")
data.AL.PAR <- subset(data.AL, Population == "PAR")
data.AL.POW <- subset(data.AL, Population == "POW")

# Plot trait data
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), pch = 16,
     data = data.AL, ylab = "Adult lifespan", xlab = "Temperature")

###### 2a. Set up JAGS Models ####
# Quadratic Model with uniform priors

sink("quad.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


# Quadratic Model with gamma priors (except sigma) 

sink("quad_inf.txt")
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
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

###### 2b. Shared settings for all models ####

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



###### 3. Fit using uniform priors ######


###### 3a: Hopland (HOP) ######

# set data
data <- data.AL.HOP 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.HOP.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.HOP , ylab = "AL for Hop prior", xlab = "Temperature")
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b: Marin35 (MAR1) ######

# set data
data <- data.AL.MAR1 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.MAR1.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.MAR1.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.MAR1.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.MAR1 , ylab = "AL for MAR1", xlab = "Temperature")
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3c: Marin29 (MAR2) ######

# set data
data <- data.AL.MAR2 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.MAR2.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.MAR2.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.MAR2.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.MAR2 , ylab = "AL for MAR2 prior", xlab = "Temperature")
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3d: Wawona (WAW) ######

# set data
data <- data.AL.WAW 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.WAW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.WAW.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.WAW.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.WAW , ylab = "AL for WAW prior", xlab = "Temperature")
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.AL.EUG 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
jag.data1<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.EUG.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.EUG.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.EUG.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.EUG , ylab = "AL for EUG prior", xlab = "Temperature")
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3f: Placer (PLA) ######

# set data
data <- data.AL.PLA 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.PLA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.PLA.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.PLA.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.PLA , ylab = "AL for PLA prior", xlab = "Temperature")
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3g: Santa Barbara (SB) ######

# set data
data <- data.AL.SB 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.SB.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.SB.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.SB.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.SB , ylab = "AL for SB prior", xlab = "Temperature")
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.AL.JRA 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.JRA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.JRA.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.JRA.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.JRA , ylab = "AL for JRA prior", xlab = "Temperature")
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.AL.PAR 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.PAR.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine output
AL.PAR.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.PAR.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.PAR , ylab = "AL for PAR prior", xlab = "Temperature")
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3j: Powder Canyon (POW) ######

# set data
data <- data.AL.POW 

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.POW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.POW.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.POW.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.POW , ylab = "AL for POW prior", xlab = "Temperature")
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### For internal checking: 10 panel plot of fits using uniform priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.HOP, ylab = "AL for Hop", xlab = "Temperature")
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.MAR1, ylab = "AL for MAR1", xlab = "Temperature")
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.MAR2 , ylab = "AL for MAR2 prior", xlab = "Temperature")
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.WAW , ylab = "AL for WAW prior", xlab = "Temperature")
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.EUG , ylab = "AL for EUG prior", xlab = "Temperature")
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.PLA , ylab = "AL for PLA prior", xlab = "Temperature")
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.SB , ylab = "AL for SB prior", xlab = "Temperature")
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.JRA , ylab = "AL for JRA prior", xlab = "Temperature")
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.PAR , ylab = "AL for PAR prior", xlab = "Temperature")
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.POW , ylab = "AL for POW prior", xlab = "Temperature")
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### 10 panel plots of fits with uniform priors, but raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$AdultLifespan ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$AdultLifespan ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$AdultLifespan ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/AL_UniformPriors.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

AL.Hop.MeansSds = MeanSd(data.AL.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.Hop.MeansSds, ylab = "AL for Hop", xlab = "Temperature", pch= 16)
arrows(AL.Hop.MeansSds$Temp.Treatment, AL.Hop.MeansSds$mean-AL.Hop.MeansSds$stderror, AL.Hop.MeansSds$Temp.Treatment, AL.Hop.MeansSds$mean + AL.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.MAR1.MeansSds = MeanSd(data.AL.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.MAR1.MeansSds, ylab = "AL for MAR1", xlab = "Temperature", pch= 16)
arrows(AL.MAR1.MeansSds$Temp.Treatment, AL.MAR1.MeansSds$mean-AL.MAR1.MeansSds$stderror, AL.MAR1.MeansSds$Temp.Treatment, AL.MAR1.MeansSds$mean + AL.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.MAR2.MeansSds = MeanSd(data.AL.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.MAR2.MeansSds, ylab = "AL for MAR2", xlab = "Temperature", pch= 16)
arrows(AL.MAR2.MeansSds$Temp.Treatment, AL.MAR2.MeansSds$mean-AL.MAR2.MeansSds$stderror, AL.MAR2.MeansSds$Temp.Treatment, AL.MAR2.MeansSds$mean + AL.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.WAW.MeansSds = MeanSd(data.AL.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.WAW.MeansSds, ylab = "AL for WAW", xlab = "Temperature", pch= 16)
arrows(AL.WAW.MeansSds$Temp.Treatment, AL.WAW.MeansSds$mean-AL.WAW.MeansSds$stderror, AL.WAW.MeansSds$Temp.Treatment, AL.WAW.MeansSds$mean + AL.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.EUG.MeansSds = MeanSd(data.AL.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.EUG.MeansSds, ylab = "AL for EUG", xlab = "Temperature", pch= 16)
arrows(AL.EUG.MeansSds$Temp.Treatment, AL.EUG.MeansSds$mean-AL.EUG.MeansSds$stderror, AL.EUG.MeansSds$Temp.Treatment, AL.EUG.MeansSds$mean + AL.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.PLA.MeansSds = MeanSd(data.AL.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.PLA.MeansSds, ylab = "AL for PLA", xlab = "Temperature", pch= 16)
arrows(AL.PLA.MeansSds$Temp.Treatment, AL.PLA.MeansSds$mean-AL.PLA.MeansSds$stderror, AL.PLA.MeansSds$Temp.Treatment, AL.PLA.MeansSds$mean + AL.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.SB.MeansSds = MeanSd(data.AL.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.SB.MeansSds, ylab = "AL for SB", xlab = "Temperature", pch= 16)
arrows(AL.SB.MeansSds$Temp.Treatment, AL.SB.MeansSds$mean-AL.SB.MeansSds$stderror, AL.SB.MeansSds$Temp.Treatment, AL.SB.MeansSds$mean + AL.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.JRA.MeansSds = MeanSd(data.AL.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.JRA.MeansSds, ylab = "AL for JRA", xlab = "Temperature", pch= 16)
arrows(AL.JRA.MeansSds$Temp.Treatment, AL.JRA.MeansSds$mean-AL.JRA.MeansSds$stderror, AL.JRA.MeansSds$Temp.Treatment, AL.JRA.MeansSds$mean + AL.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.PAR.MeansSds = MeanSd(data.AL.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.PAR.MeansSds, ylab = "AL for PAR", xlab = "Temperature", pch= 16)
arrows(AL.PAR.MeansSds$Temp.Treatment, AL.PAR.MeansSds$mean-AL.PAR.MeansSds$stderror, AL.PAR.MeansSds$Temp.Treatment, AL.PAR.MeansSds$mean + AL.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

AL.POW.MeansSds = MeanSd(data.AL.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.POW.MeansSds, ylab = "AL for POW", xlab = "Temperature", pch= 16)
arrows(AL.POW.MeansSds$Temp.Treatment, AL.POW.MeansSds$mean-AL.POW.MeansSds$stderror, AL.POW.MeansSds$Temp.Treatment, AL.POW.MeansSds$mean + AL.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

dev.off()




####### 4. Get TPC parameters #####
###### 4a. Topt of AL for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 4b. Pmax of AL for each population ######

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

ALdf_uni <- rbind.data.frame(ParamCompile(AL.HOP.out), ParamCompile(AL.MAR1.out),
                             ParamCompile(AL.MAR2.out), ParamCompile(AL.WAW.out),
                             ParamCompile(AL.EUG.out), ParamCompile(AL.PLA.out),
                             ParamCompile(AL.SB.out), ParamCompile(AL.JRA.out),
                             ParamCompile(AL.PAR.out), ParamCompile(AL.POW.out))
ALdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(ALdf_uni, file = "AL_meansd_uni.Rsave")

###### 5. Plotting ######
load("AL_meansd_uni.Rsave")
library(ggplot2)
###### 5a. Plot of CTmax, CTmin, Topt for each population ######
# order by latitude
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

ALdf_uni$Population = factor(ALdf_uni$Population, levels = LatOrder)
ALdf_uni = ALdf_uni[order(ALdf_uni$Population),]
ALdf_uni = ALdf_uni[ALdf_uni$param != "Tbreadth",]
ALdf_uni = ALdf_uni[ALdf_uni$param != "Pmax",]

plotALuni = ggplot(ALdf_uni, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Adult lifespan") + 
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

load("AL_meansd_uni.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
ALdf_uni$Population = factor(ALdf_uni$Population, levels = PopOrder)
ALdf_uni = ALdf_uni[order(ALdf_uni$Population),]
ALdf_uni = ALdf_uni[ALdf_uni$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotAL_Pmax = ggplot(ALdf_uni, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0, 20)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Adult lifespan (days)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5c. Plot adult lifespan curve for all populations #####

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]


# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, 
     xlim = c(0, 50), ylim = c(0,17), data = data.AL.HOP, type = "n", bty = "n",
     ylab = "lifespan (days)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Adult lifespan", cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)
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
data.AL.HOP.prior <- subset(data.AL, Population != "HOP")
data.AL.MAR1.prior <- subset(data.AL, Population != "MARIN35")
data.AL.MAR2.prior <- subset(data.AL, Population != "MARIN29")
data.AL.WAW.prior <- subset(data.AL, Population != "WAW")
data.AL.EUG.prior <- subset(data.AL, Population != "EUG")
data.AL.PLA.prior <- subset(data.AL, Population != "PLA")
data.AL.SB.prior <- subset(data.AL, Population != "SB")
data.AL.JRA.prior <- subset(data.AL, Population != "JRA")
data.AL.PAR.prior <- subset(data.AL, Population != "PAR")
data.AL.POW.prior <- subset(data.AL, Population != "POW")


###### 6a. Hopland #####
# set data
data <- data.AL.HOP.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.HOP.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.HOP.prior, ylab = "AL for Hop prior", xlab = "Temperature")
lines(AL.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6b: Marin35 (MAR1) #####

# set data
data <- data.AL.MAR1.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.MAR1.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.MAR1.prior, ylab = "AL for MAR1 prior", xlab = "Temperature")
lines(AL.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.AL.MAR2.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.MAR2.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.MAR2.prior, ylab = "AL for MAR2 prior", xlab = "Temperature")
lines(AL.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.AL.WAW.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.WAW.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.WAW.prior, ylab = "AL for WAW prior", xlab = "Temperature")
lines(AL.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.AL.EUG.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.EUG.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.EUG.prior, ylab = "AL for EUG prior", xlab = "Temperature")
lines(AL.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.AL.PLA.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.PLA.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.PLA.prior, ylab = "AL for PLA prior", xlab = "Temperature")
lines(AL.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.AL.SB.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.SB.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.SB.prior, ylab = "AL for SB prior", xlab = "Temperature")
lines(AL.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.AL.JRA.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.JRA.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.JRA.prior, ylab = "AL for JRA prior", xlab = "Temperature")
lines(AL.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.AL.PAR.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.PAR.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.PAR.prior, ylab = "AL for PAR prior", xlab = "Temperature")
lines(AL.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyon (POW) ######

# set data
data <- data.AL.POW.prior

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
AL.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.POW.prior.out)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.POW.prior, ylab = "AL for POW prior", xlab = "Temperature")
lines(AL.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### Save posterior distributions for informative priors #####

AL.HOP.prior.cf.dists <- data.frame(q = as.vector(AL.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.HOP.prior.gamma.fits = apply(AL.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.MAR1.prior.cf.dists <- data.frame(q = as.vector(AL.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(AL.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(AL.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.MAR1.prior.gamma.fits = apply(AL.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.MAR2.prior.cf.dists <- data.frame(q = as.vector(AL.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(AL.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(AL.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.MAR2.prior.gamma.fits = apply(AL.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.WAW.prior.cf.dists <- data.frame(q = as.vector(AL.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.WAW.prior.gamma.fits = apply(AL.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.EUG.prior.cf.dists <- data.frame(q = as.vector(AL.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.EUG.prior.gamma.fits = apply(AL.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.PLA.prior.cf.dists <- data.frame(q = as.vector(AL.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.PLA.prior.gamma.fits = apply(AL.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.SB.prior.cf.dists <- data.frame(q = as.vector(AL.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(AL.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(AL.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.SB.prior.gamma.fits = apply(AL.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.JRA.prior.cf.dists <- data.frame(q = as.vector(AL.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.JRA.prior.gamma.fits = apply(AL.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.PAR.prior.cf.dists <- data.frame(q = as.vector(AL.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.PAR.prior.gamma.fits = apply(AL.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.POW.prior.cf.dists <- data.frame(q = as.vector(AL.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.POW.prior.gamma.fits = apply(AL.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

AL.hypers <- list(AL.HOP.prior.gamma.fits, AL.MAR1.prior.gamma.fits, AL.MAR2.prior.gamma.fits,
                  AL.WAW.prior.gamma.fits, AL.EUG.prior.gamma.fits, AL.PLA.prior.gamma.fits,
                  AL.SB.prior.gamma.fits, AL.JRA.prior.gamma.fits, AL.PAR.prior.gamma.fits, AL.POW.prior.gamma.fits)
save(AL.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/ALhypers.Rsave")
###### Quadratic Model with gamma priors (except sigma) ######

sink("quad_inf.txt")
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
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

###### 7. Fit using informative priors #####
# Load hyperparameters to bypass above code of generating them 
load("ALhypers.Rsave")
AL.HOP.prior.gamma.fits <- AL.hypers[[1]]
AL.MAR1.prior.gamma.fits <- AL.hypers[[2]]
AL.MAR2.prior.gamma.fits <- AL.hypers[[3]]
AL.WAW.prior.gamma.fits <- AL.hypers[[4]]
AL.EUG.prior.gamma.fits <- AL.hypers[[5]]
AL.PLA.prior.gamma.fits <- AL.hypers[[6]]
AL.SB.prior.gamma.fits <- AL.hypers[[7]]
AL.JRA.prior.gamma.fits <- AL.hypers[[8]]
AL.PAR.prior.gamma.fits <- AL.hypers[[9]]
AL.POW.prior.gamma.fits <- AL.hypers[[10]]


###### 7a. Hopland ######

# Set data 
data <- data.AL.HOP
hypers <- AL.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(AL.HOP.out.inf)

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.HOP, ylab = "AL for Hopland", xlab = "Temperature", pch = 1)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 7b. Marin35 ######

# Set data 
data <- data.AL.MAR1
hypers <- AL.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.MAR1, ylab = "AL for Marin35", xlab = "Temperature", pch = 1)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.AL.MAR2
hypers <- AL.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.MAR2, ylab = "AL for Marin29", xlab = "Temperature", pch = 1)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.AL.WAW
hypers <- AL.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.WAW, ylab = "AL for Wawona", xlab = "Temperature", pch = 1)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.AL.EUG
hypers <- AL.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.EUG, ylab = "AL for Eugene", xlab = "Temperature", pch = 1)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.AL.PLA
hypers <- AL.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.PLA, ylab = "AL for Placer", xlab = "Temperature", pch = 1)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.AL.SB
hypers <- AL.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.SB, ylab = "AL for Santa Barbara", xlab = "Temperature", pch = 1)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.AL.JRA
hypers <- AL.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.JRA, ylab = "AL for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.AL.PAR
hypers <- AL.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.PAR, ylab = "AL for Paso Robles", xlab = "Temperature", pch = 1)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon ######

# Set data 
data <- data.AL.POW
hypers <- AL.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$AdultLifespan
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
AL.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
AL.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,35), data = data.AL.POW, ylab = "AL for Powder Canyon", xlab = "Temperature", pch = 1)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 10 panel plots of fits with raw data as means & sderrors  #####
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$AdultLifespan ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$AdultLifespan ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$AdultLifespan ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/AL.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

AL.Hop.MeansSds = MeanSd(data.AL.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.Hop.MeansSds, ylab = "AL for Hop", xlab = "Temperature", pch= 16)
arrows(AL.Hop.MeansSds$Temp.Treatment, AL.Hop.MeansSds$mean-AL.Hop.MeansSds$stderror, AL.Hop.MeansSds$Temp.Treatment, AL.Hop.MeansSds$mean + AL.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 18, cex = 1.5)
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.MAR1.MeansSds = MeanSd(data.AL.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.MAR1.MeansSds, ylab = "AL for MAR1", xlab = "Temperature", pch= 16)
arrows(AL.MAR1.MeansSds$Temp.Treatment, AL.MAR1.MeansSds$mean-AL.MAR1.MeansSds$stderror, AL.MAR1.MeansSds$Temp.Treatment, AL.MAR1.MeansSds$mean + AL.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 18, cex = 1.5)
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.MAR2.MeansSds = MeanSd(data.AL.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.MAR2.MeansSds, ylab = "AL for MAR2", xlab = "Temperature", pch= 16)
arrows(AL.MAR2.MeansSds$Temp.Treatment, AL.MAR2.MeansSds$mean-AL.MAR2.MeansSds$stderror, AL.MAR2.MeansSds$Temp.Treatment, AL.MAR2.MeansSds$mean + AL.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 36.5, y = 18, cex = 1.5)
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.WAW.MeansSds = MeanSd(data.AL.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.WAW.MeansSds, ylab = "AL for WAW", xlab = "Temperature", pch= 16)
arrows(AL.WAW.MeansSds$Temp.Treatment, AL.WAW.MeansSds$mean-AL.WAW.MeansSds$stderror, AL.WAW.MeansSds$Temp.Treatment, AL.WAW.MeansSds$mean + AL.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 36.5, y = 18, cex = 1.5)
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.EUG.MeansSds = MeanSd(data.AL.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.EUG.MeansSds, ylab = "AL for EUG", xlab = "Temperature", pch= 16)
arrows(AL.EUG.MeansSds$Temp.Treatment, AL.EUG.MeansSds$mean-AL.EUG.MeansSds$stderror, AL.EUG.MeansSds$Temp.Treatment, AL.EUG.MeansSds$mean + AL.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG", x = 36.5, y = 18, cex = 1.5)
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.PLA.MeansSds = MeanSd(data.AL.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.PLA.MeansSds, ylab = "AL for PLA", xlab = "Temperature", pch= 16)
arrows(AL.PLA.MeansSds$Temp.Treatment, AL.PLA.MeansSds$mean-AL.PLA.MeansSds$stderror, AL.PLA.MeansSds$Temp.Treatment, AL.PLA.MeansSds$mean + AL.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 36.5, y = 18, cex = 1.5)
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.SB.MeansSds = MeanSd(data.AL.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.SB.MeansSds, ylab = "AL for SB", xlab = "Temperature", pch= 16)
arrows(AL.SB.MeansSds$Temp.Treatment, AL.SB.MeansSds$mean-AL.SB.MeansSds$stderror, AL.SB.MeansSds$Temp.Treatment, AL.SB.MeansSds$mean + AL.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 36.5, y = 18, cex = 1.5)
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.JRA.MeansSds = MeanSd(data.AL.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.JRA.MeansSds, ylab = "AL for JRA", xlab = "Temperature", pch= 16)
arrows(AL.JRA.MeansSds$Temp.Treatment, AL.JRA.MeansSds$mean-AL.JRA.MeansSds$stderror, AL.JRA.MeansSds$Temp.Treatment, AL.JRA.MeansSds$mean + AL.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 36.5, y = 18, cex = 1.5)
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.PAR.MeansSds = MeanSd(data.AL.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.PAR.MeansSds, ylab = "AL for PAR", xlab = "Temperature", pch= 16)
arrows(AL.PAR.MeansSds$Temp.Treatment, AL.PAR.MeansSds$mean-AL.PAR.MeansSds$stderror, AL.PAR.MeansSds$Temp.Treatment, AL.PAR.MeansSds$mean + AL.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 36.5, y = 18, cex = 1.5)
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

AL.POW.MeansSds = MeanSd(data.AL.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,20), data = AL.POW.MeansSds, ylab = "AL for POW", xlab = "Temperature", pch= 16)
arrows(AL.POW.MeansSds$Temp.Treatment, AL.POW.MeansSds$mean-AL.POW.MeansSds$stderror, AL.POW.MeansSds$Temp.Treatment, AL.POW.MeansSds$mean + AL.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 36.5, y = 18, cex = 1.5)
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(AL.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

dev.off()

###### 8. Get TPC parameters #####
###### 8a. Topt of AL for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of AL for each population ######

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

ALdf_inf <- rbind.data.frame(ParamCompile(AL.HOP.out.inf), ParamCompile(AL.MAR1.out.inf),
                             ParamCompile(AL.MAR2.out.inf), ParamCompile(AL.WAW.out.inf),
                             ParamCompile(AL.EUG.out.inf), ParamCompile(AL.PLA.out.inf),
                             ParamCompile(AL.SB.out.inf), ParamCompile(AL.JRA.out.inf),
                             ParamCompile(AL.PAR.out.inf), ParamCompile(AL.POW.out.inf))
ALdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(ALdf_inf, file = "AL_meansd_inf.Rsave")

###### 9. Plotting  #####
library(ggplot2)
load("AL_meansd_inf.Rsave")
###### 9a. Plot of CTmax, CTmin, Topt for each population #####
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

ALdf_inf$Population = factor(ALdf_inf$Population, levels = LatOrder)
ALdf_inf = ALdf_inf[order(ALdf_inf$Population),]
ALdf_inf = ALdf_inf[ALdf_inf$param != "Tbreadth",]
ALdf_inf = ALdf_inf[ALdf_inf$param != "Pmax",]

plotALinf = ggplot(ALdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, size = 2, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors, size = 0.7,
                position=position_dodge(.9), width = 0.5) + 
  ggtitle("Adult lifespan") + 
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

load("AL_meansd_inf.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
ALdf_inf$Population = factor(ALdf_inf$Population, levels = PopOrder)
ALdf_inf = ALdf_inf[order(ALdf_inf$Population),]
ALdf_inf = ALdf_inf[ALdf_inf$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotALinf_Pmax = ggplot(ALdf_inf, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0, 20)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Adult lifespan (days)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 9c. Plot adult lifespans curve for all populations #####
load(file = "AL_meansd_inf.Rsave")

EUGdata = AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
par(mar = (c(4.5, 4.5, 2.5, 1)))
plot(AdultLifespan ~ Temp.Treatment, 
     xlim = c(5, 40), ylim = c(0,15), data = data.AL.HOP, type = "n", bty = "n",
     ylab = "lifespan (days)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Adult lifespan", cex.main = 1.7, cex.lab = 1.4, cex.axis = 1.2)
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





