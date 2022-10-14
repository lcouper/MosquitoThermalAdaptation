# Trait Fits for Larval Development Rate ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for larval development rate (LDR) for 
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit larval dev rate thermal responses (Briere) using low info priors
#  4) Fit gamma distributions to LDR prior thermal responses
#  5) Fit larval dev rate thermal responses with data-informed priors
#  6) Plot parameters


#### 1. Set up workspace, load packages, get data, etc. ####

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

# Subset data by population
data.LDR.HOP <- subset(data.LDR, Population == "HOP")
data.LDR.MAR1 <- subset(data.LDR, Population == "MAR35")
data.LDR.MAR2 <- subset(data.LDR, Population == "MAR29")
data.LDR.WAW <- subset(data.LDR, Population == "WAW")
data.LDR.EUG <- subset(data.LDR, Population == "EUG")
data.LDR.PLA <- subset(data.LDR, Population == "PLA")
data.LDR.SB <- subset(data.LDR, Population == "SB")
data.LDR.JRA <- subset(data.LDR, Population == "JRA")
data.LDR.PAR <- subset(data.LDR, Population == "PAR")
data.LDR.POW <- subset(data.LDR, Population == "POW")


##### 2. Set up JAGS model #####
# 2a. Briere Model with uniform priors #

sink("briere.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
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



###### 3. Shared settings for all models

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


#### 3. Generate data-informed priors using leave-one-out approach ####
# # For each population, use data for all *other* populations to obtain
# distributions of each TPC parameter
# For this initial fitting, we use uniform priors, 
# constrained over a biologically realistic range

# Subset Data - leave-one-out
data.LDR.HOP.prior <- subset(data.LDR, Population != "HOP")
data.LDR.MAR1.prior <- subset(data.LDR, Population != "MAR1")
data.LDR.MAR2.prior <- subset(data.LDR, Population != "MAR2")
data.LDR.WAW.prior <- subset(data.LDR, Population != "WAW")
data.LDR.EUG.prior <- subset(data.LDR, Population != "EUG")
data.LDR.PLA.prior <- subset(data.LDR, Population != "PLA")
data.LDR.SB.prior <- subset(data.LDR, Population != "SB")
data.LDR.JRA.prior <- subset(data.LDR, Population != "JRA")
data.LDR.PAR.prior <- subset(data.LDR, Population != "PAR")
data.LDR.POW.prior <- subset(data.LDR, Population != "POW")

###### 3a: Hopland (HOP) ######

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

###### 3b: Marin35 (MAR1) ######

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


###### 3c: Marin29 (MAR2) ######

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


###### 3d: Wawona (WAW) ######

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

###### 3e: Eugene (EUG) ######

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

###### 3f: Placer (PLA) ######

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


###### 3g: Santa Barbara (SB) ######

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


###### 3h: Jasper Ridge (JRA) ######

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


###### 3i: Paso Robles (PAR) ######

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


###### 3j: Powder Canyono (POW) ######

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











#### 4. Fit gamma distributions to LDR prior thermal responses ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists

###### 4a. Hopland ######

LDR.HOP.prior.cf.dists <- data.frame(q = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                T0 = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                Tm = as.vector(LDR.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.HOP.prior.gamma.fits = apply(LDR.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4b. Marin35 (MAR1) ######

LDR.MAR1.prior.cf.dists <- data.frame(q = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.MAR1.prior.gamma.fits = apply(LDR.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4c. Marin29 (MAR2) ######

LDR.MAR2.prior.cf.dists <- data.frame(q = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(LDR.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.MAR2.prior.gamma.fits = apply(LDR.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

LDR.WAW.prior.cf.dists <- data.frame(q = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(LDR.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.WAW.prior.gamma.fits = apply(LDR.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

LDR.EUG.prior.cf.dists <- data.frame(q = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.EUG.prior.gamma.fits = apply(LDR.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

LDR.PLA.prior.cf.dists <- data.frame(q = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.PLA.prior.gamma.fits = apply(LDR.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

LDR.SB.prior.cf.dists <- data.frame(q = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.SB.prior.gamma.fits = apply(LDR.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

LDR.JRA.prior.cf.dists <- data.frame(q = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(LDR.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.JRA.prior.gamma.fits = apply(LDR.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

LDR.PAR.prior.cf.dists <- data.frame(q = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(LDR.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.PAR.prior.gamma.fits = apply(LDR.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

LDR.POW.prior.cf.dists <- data.frame(q = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(LDR.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
LDR.POW.prior.gamma.fits = apply(LDR.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
LDR.hypers <- list(LDR.HOP.prior.gamma.fits, LDR.MAR1.prior.gamma.fits, LDR.MAR2.prior.gamma.fits,
                   LDR.WAW.prior.gamma.fits, LDR.EUG.prior.gamma.fits, LDR.PLA.prior.gamma.fits,
                   LDR.SB.prior.gamma.fits, LDR.JRA.prior.gamma.fits, LDR.PAR.prior.gamma.fits, LDR.POW.prior.gamma.fits)
save(LDR.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/LDRhypers.Rsave")



#### 5. Fit LDR thermal responses with data-informed priors ####

# to bypass above code of generating hyperparameters
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


###### 5a. Hopland ######

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
#save(LDR.HOP.out.inf, file = "jagsout_LDR.HOPp_inf.Rdata")

# Plot data + fit
plot(LarvalDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.12), data = data.LDR.HOP, ylab = "LDR for Hopland", xlab = "Temperature", pch = 1)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 5b. Marin35 ######

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

###### 5c. Marin29 ######

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

###### 5d. Wawona ######

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


###### 5e. Eugene ######

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

###### 5f. Placer ######

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


###### 5g. Santa Barbara ######

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

###### 5g. Jasper Ridge ######

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


###### 5i. Paso Robles ######

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


###### 5j. Powder canyon ######

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

###### Topt of LDR for each population ######
ToptHOP =Temp.xs[which.max(as.vector(LDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptMAR1 =Temp.xs[which.max(as.vector(LDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5
ToptMAR2 =Temp.xs[which.max(as.vector(LDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26.5
ToptWAW =Temp.xs[which.max(as.vector(LDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28
ToptEUG =Temp.xs[which.max(as.vector(LDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5
ToptPLA =Temp.xs[which.max(as.vector(LDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptSB =Temp.xs[which.max(as.vector(LDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5
ToptJRA =Temp.xs[which.max(as.vector(LDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26.5
ToptPAR =Temp.xs[which.max(as.vector(LDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptPOW =Temp.xs[which.max(as.vector(LDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5

##### Calculate Tbreadth for each population #####

# Tbreadth defined as the temperature range where performance remains >= 50% peak performance

# 1. Write Tbreadth function
Tbreadth = function(x)
{index = which.max(as.vector(x$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))
peak = x$BUGSoutput$summary[index + 5, "50%"]
fiftypeak = peak/2
fiftypeakindexes = which(as.vector(x$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]) >= fiftypeak)
mintempindex = fiftypeakindexes[1]
maxtempindex = fiftypeakindexes[length(fiftypeakindexes)]
Tbreadth = Temp.xs[maxtempindex] - Temp.xs[mintempindex]
return(Tbreadth)}

# 2. Apply to each population
TbreadthHOP = Tbreadth(LDR.HOP.out.inf) # 18
TbreadthMAR1 = Tbreadth(LDR.MAR1.out.inf) # 18
TbreadthMAR2 = Tbreadth(LDR.MAR2.out.inf) # 17.5
TbreadthWAW = Tbreadth(LDR.WAW.out.inf) # 18.5
TbreadthEUG = Tbreadth(LDR.EUG.out.inf) # 18
TbreadthPLA = Tbreadth(LDR.PLA.out.inf) # 18
TbreadthSB = Tbreadth(LDR.SB.out.inf) # 18
TbreadthJRA = Tbreadth(LDR.JRA.out.inf) # 17
TbreadthPAR = Tbreadth(LDR.PAR.out.inf) # 17.5
TbreadthPOW = Tbreadth(LDR.POW.out.inf) # 18


#### Compile prediction data for all populations
# T0, Tm, and Topt along with 95% credible intervals (2.5%, 97.5%)

LDRmeansd <- list(LDR.HOP.out.inf$BUGSoutput$summary[1:2,1:2], LDR.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                LDR.MAR2.out.inf$BUGSoutput$summary[1:2,1:2], LDR.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                LDR.EUG.out.inf$BUGSoutput$summary[1:2,1:2], LDR.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                LDR.SB.out.inf$BUGSoutput$summary[1:2,1:2], LDR.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                LDR.POW.out.inf$BUGSoutput$summary[1:2,1:2], LDR.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
LDRtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptMAR2, ToptWAW, ToptEUG, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 10)), nrow = 10, ncol = 2))
LDRtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthMAR2, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                      TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 10)), nrow = 10, ncol = 2))
colnames(LDRtops) = c("mean", "sd")
colnames(LDRtbreadths) = c("mean", "sd")
LDRdf2 = do.call(rbind.data.frame, LDRmeansd)
LDRdf_1 = rbind(LDRdf2, LDRtops)
LDRdf = rbind(LDRdf_1, LDRtbreadths)
LDRdf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                     rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
LDRdf$Param = c(rep(c("Tmin", "Tmax"), 10), rep("Topt", 10), rep("Tbreadth", 10))
save(LDRdf, file = "LDR_meansd.Rsave")


#### Compile prediction data for all populations ###### 

LDRmeansd <- list(LDR.HOP.out.inf$BUGSoutput$summary[1:2,1:2], LDR.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                  LDR.MAR2.out.inf$BUGSoutput$summary[1:2,1:2], LDR.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                  LDR.EUG.out.inf$BUGSoutput$summary[1:2,1:2], LDR.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                  LDR.SB.out.inf$BUGSoutput$summary[1:2,1:2], LDR.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                  LDR.POW.out.inf$BUGSoutput$summary[1:2,1:2], LDR.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
LDRtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptMAR2, ToptWAW, ToptEUG, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 10)), nrow = 10, ncol = 2))
LDRtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthMAR2, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                      TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 10)), nrow = 10, ncol = 2))
colnames(LDRtops) = c("mean", "sd")
colnames(LDRtbreadths) = c("mean", "sd")
LDRdf2 = do.call(rbind.data.frame, LDRmeansd)
LDRdf_1 = rbind(LDRdf2, LDRtops)
LDRdf = rbind(LDRdf_1, LDRtbreadths)
LDRdf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                     rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
LDRdf$Param = c(rep(c("Tmin", "Tmax"), 10), rep("Topt", 10), rep("Tbreadth", 10))
save(LDRdf, file = "LDR_meansd.Rsave")

#### 6a. Plot TPC parameters for each population #####

library(ggplot2)
load("LDR_meansd.Rsave")

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
LDRdf$Population = factor(LDRdf$Population, levels = PopOrder)
LDRdf = LDRdf[order(LDRdf$Population),]
LDRdf = LDRdf[LDRdf$Param != "Tbreadth",]


# colors for each population 
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotLDR = ggplot(LDRdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Larval development rate") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 


#### 6b. Plot larval dev rates curve for all populaions #####

load(file = "LDR_meansd.Rsave")
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

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
plot(LarvalDevRate ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,0.11), data = data.LDR.HOP, type = "n", bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Larval development rates", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(EUGdata ~ Temp.xs, col = "#313695")
lines(WAWdata ~ Temp.xs, col = "#4575b4")
lines(PLAdata ~ Temp.xs, col = "#74add1")
lines(HOPdata ~ Temp.xs, col = "#abd9e9")
lines(JRAdata ~ Temp.xs, col = "#e0f3f8")
lines(MAR2data ~ Temp.xs, col = "#fee090")
lines(MAR1data ~ Temp.xs, col = "#fdae61")
lines(POWdata ~ Temp.xs, col = "#f46d43")
lines(PARdata ~ Temp.xs, col = "#d73027")
lines(SBdata ~ Temp.xs, col = "#a50026")


