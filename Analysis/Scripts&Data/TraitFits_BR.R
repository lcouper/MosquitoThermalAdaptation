# Trait Fits for Biting Rate ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for biting rate for
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit biting rate thermal responses (Briere) using low info priors
#  4) Fit gamma distributions to biting rate prior thermal responses
#  5) Fit biting rate thermal responses with data-informed priors
#  6) Plot parameters

#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.BR = read.csv("ReproExp_Data.csv")[,-1]
# Calculate biting rate (1 / gonotrophic cycle length)
data.BR$BitingRate = 1/data.BR$CycleLength
# remove rows with no data
data.BR = data.BR[!is.na(data.BR$BitingRate),]

# Plot trait data to see which function to fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), pch = 16,
     data = data.BR, ylab = "Biting rate", xlab = "Temperature")

# Subset data by population
data.BR.HOP <- subset(data.BR, Population == "HOP")
data.BR.MAR1 <- subset(data.BR, Population == "MAR")
data.BR.WAW <- subset(data.BR, Population == "WAW")
data.BR.EUG <- subset(data.BR, Population == "EUG")
data.BR.PLA <- subset(data.BR, Population == "PLA")
data.BR.SB <- subset(data.BR, Population == "SB")
data.BR.JRA <- subset(data.BR, Population == "JRA")
data.BR.PAR <- subset(data.BR, Population == "PAR")
data.BR.POW <- subset(data.BR, Population == "POW")

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
data.BR.HOP.prior <- subset(data.BR, Population != "HOP")
data.BR.MAR1.prior <- subset(data.BR, Population != "MAR")
data.BR.WAW.prior <- subset(data.BR, Population != "WAW")
data.BR.EUG.prior <- subset(data.BR, Population != "EUG")
data.BR.PLA.prior <- subset(data.BR, Population != "PLA")
data.BR.SB.prior <- subset(data.BR, Population != "SB")
data.BR.JRA.prior <- subset(data.BR, Population != "JRA")
data.BR.PAR.prior <- subset(data.BR, Population != "PAR")
data.BR.POW.prior <- subset(data.BR, Population != "POW")


###### 3a: Hopland (HOP) ######

# set data
data <- data.BR.HOP.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.HOP.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.17), data = data.BR.HOP.prior, ylab = "BR for Hop prior", xlab = "Temperature")
lines(BR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b: Marin35 (MAR1) ######

# set data
data <- data.BR.MAR1.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.MAR1.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.MAR1.prior, ylab = "BR for MAR1 prior", xlab = "Temperature")
lines(BR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3c. Wawona ######

# set data
data <- data.BR.WAW.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.WAW.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.WAW.prior, ylab = "BR for WAW prior", xlab = "Temperature")
lines(BR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3d. Eugene #####

# set data
data <- data.BR.EUG.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.EUG.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.EUG.prior, ylab = "BR for EUG prior", xlab = "Temperature")
lines(BR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3e. 



##### 3e. Placer #####

# set data
data <- data.BR.PLA.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.PLA.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.PLA.prior, ylab = "BR for PLA prior", xlab = "Temperature")
lines(BR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3f. Santa Barbara #####

# set data
data <- data.BR.SB.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.SB.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.SB.prior, ylab = "BR for SB prior", xlab = "Temperature")
lines(BR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3g. Jasper Ridge ######

# set data
data <- data.BR.JRA.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.JRA.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.JRA.prior, ylab = "BR for JRA prior", xlab = "Temperature")
lines(BR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
##### 3h. Paso Robles ######

# set data
data <- data.BR.PAR.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.PAR.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.PAR.prior, ylab = "BR for PAR prior", xlab = "Temperature")
lines(BR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3i. Powder canyon #####
# set data
data <- data.BR.POW.prior

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
BR.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(BR.POW.prior.out)

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.POW.prior, ylab = "BR for POW prior", xlab = "Temperature")
lines(BR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




#### 4. Fit gamma distributions to BR prior thermal responses ####

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists

###### 4a. Hopland ######

BR.HOP.prior.cf.dists <- data.frame(q = as.vector(BR.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.HOP.prior.gamma.fits = apply(BR.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4b. Marin35 (MAR1) ######

BR.MAR1.prior.cf.dists <- data.frame(q = as.vector(BR.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(BR.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(BR.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.MAR1.prior.gamma.fits = apply(BR.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

BR.WAW.prior.cf.dists <- data.frame(q = as.vector(BR.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.WAW.prior.gamma.fits = apply(BR.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

BR.EUG.prior.cf.dists <- data.frame(q = as.vector(BR.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.EUG.prior.gamma.fits = apply(BR.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

BR.PLA.prior.cf.dists <- data.frame(q = as.vector(BR.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.PLA.prior.gamma.fits = apply(BR.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

BR.SB.prior.cf.dists <- data.frame(q = as.vector(BR.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(BR.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(BR.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.SB.prior.gamma.fits = apply(BR.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

BR.JRA.prior.cf.dists <- data.frame(q = as.vector(BR.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.JRA.prior.gamma.fits = apply(BR.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

BR.PAR.prior.cf.dists <- data.frame(q = as.vector(BR.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.PAR.prior.gamma.fits = apply(BR.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

BR.POW.prior.cf.dists <- data.frame(q = as.vector(BR.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(BR.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(BR.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
BR.POW.prior.gamma.fits = apply(BR.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
BR.hypers <- list(BR.HOP.prior.gamma.fits, BR.MAR1.prior.gamma.fits, 
                  BR.WAW.prior.gamma.fits, BR.EUG.prior.gamma.fits, BR.PLA.prior.gamma.fits,
                  BR.SB.prior.gamma.fits, BR.JRA.prior.gamma.fits, BR.PAR.prior.gamma.fits, BR.POW.prior.gamma.fits)
save(BR.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BRhypers.Rsave")



#### 5. Fit BR thermal responses with data-informed priors ####

# to bypass above code of generating hyperparameters
load("BRhypers.Rsave")
BR.HOP.prior.gamma.fits <- BR.hypers[[1]]
BR.MAR1.prior.gamma.fits <- BR.hypers[[2]]
BR.WAW.prior.gamma.fits <- BR.hypers[[3]]
BR.EUG.prior.gamma.fits <- BR.hypers[[4]]
BR.PLA.prior.gamma.fits <- BR.hypers[[5]]
BR.SB.prior.gamma.fits <- BR.hypers[[6]]
BR.JRA.prior.gamma.fits <- BR.hypers[[7]]
BR.PAR.prior.gamma.fits <- BR.hypers[[8]]
BR.POW.prior.gamma.fits <- BR.hypers[[9]]

###### 5a. Hopland ######

# Set data 
data <- data.BR.HOP
hypers <- BR.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.HOP.out.inf$BUGSoutput$summary[1:5,]
#save(BR.HOP.out.inf, file = "jagsout_BR.HOPp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.HOP, ylab = "BR for Hopland", xlab = "Temperature", pch = 1)
lines(BR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)



###### 5b. Marin35 #####
# Set data 
data <- data.BR.MAR1
hypers <- BR.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.MAR1.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.MAR1.out.inf)
#save(BR.MAR1.out.inf, file = "jagsout_BR.MAR1p_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.12), data = data.BR.MAR1, ylab = "BR for MAR1land", xlab = "Temperature", pch = 1)
lines(BR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


###### 5c. Wawona #####

# Set data 
data <- data.BR.WAW
hypers <- BR.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.WAW.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.WAW.out.inf)
#save(BR.WAW.out.inf, file = "jagsout_BR.WAWp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.WAW, ylab = "BR for WAWland", xlab = "Temperature", pch = 1)
lines(BR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 5d. Eugene ######

# Set data 
data <- data.BR.EUG
hypers <- BR.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.EUG.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.EUG.out.inf)
#save(BR.EUG.out.inf, file = "jagsout_BR.EUGp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.EUG, ylab = "BR for EUGland", xlab = "Temperature", pch = 1)
lines(BR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


###### 5e. Placer #######
# Set data 
data <- data.BR.PLA
hypers <- BR.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.PLA.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.PLA.out.inf)
#save(BR.PLA.out.inf, file = "jagsout_BR.PLAp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.PLA, ylab = "BR for PLAland", xlab = "Temperature", pch = 1)
lines(BR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 5f. Santa Barbara ######

# Set data 
data <- data.BR.SB
hypers <- BR.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.SB, ylab = "BR for SBland", xlab = "Temperature", pch = 1)
lines(BR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 5g. Jasper Ridge ######
# Set data 
data <- data.BR.JRA
hypers <- BR.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.JRA.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.JRA.out.inf)
#save(BR.JRA.out.inf, file = "jagsout_BR.JRAp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.JRA, ylab = "BR for JRAland", xlab = "Temperature", pch = 1)
lines(BR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 5h. Paso Robles #####

# Set data 
data <- data.BR.PAR
hypers <- BR.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.PAR.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.PAR.out.inf)
#save(BR.PAR.out.inf, file = "jagsout_BR.PARp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.PAR, ylab = "BR for PARland", xlab = "Temperature", pch = 1)
lines(BR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


###### 5i. Powder canyon ######

# Set data 
data <- data.BR.POW
hypers <- BR.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$BitingRate
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
BR.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
BR.POW.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(BR.POW.out.inf)
#save(BR.POW.out.inf, file = "jagsout_BR.POWp_inf.Rdata")

# Plot data + fit
plot(BitingRate ~ Temperature, xlim = c(5, 45), ylim = c(0,0.15), data = data.BR.POW, ylab = "BR for POWland", xlab = "Temperature", pch = 1)
lines(BR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(BR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(BR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


###### Topt of BR for each population ######
ToptHOP =Temp.xs[which.max(as.vector(BR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28.5
ToptMAR1 =Temp.xs[which.max(as.vector(BR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26.5
ToptWAW =Temp.xs[which.max(as.vector(BR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26.5
ToptEUG =Temp.xs[which.max(as.vector(BR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptPLA =Temp.xs[which.max(as.vector(BR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptSB =Temp.xs[which.max(as.vector(BR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 23
ToptJRA =Temp.xs[which.max(as.vector(BR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptPAR =Temp.xs[which.max(as.vector(BR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26
ToptPOW =Temp.xs[which.max(as.vector(BR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28

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
TbreadthHOP = Tbreadth(BR.HOP.out.inf) # 19
TbreadthMAR1 = Tbreadth(BR.MAR1.out.inf) # 17.5
TbreadthWAW = Tbreadth(BR.WAW.out.inf) # 17
TbreadthEUG = Tbreadth(BR.EUG.out.inf) # 17
TbreadthPLA = Tbreadth(BR.PLA.out.inf) # 18
TbreadthSB = Tbreadth(BR.SB.out.inf) # 15.5
TbreadthJRA = Tbreadth(BR.JRA.out.inf) # 17
TbreadthPAR = Tbreadth(BR.PAR.out.inf) # 18
TbreadthPOW = Tbreadth(BR.POW.out.inf) # 17.5

#### Compile prediction data for all populations
# T0, Tm, and Topt along with 95% credible intervals (2.5%, 97.5%)

BRmeansd <- list(BR.HOP.out.inf$BUGSoutput$summary[1:2,1:2], BR.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                 BR.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                 BR.EUG.out.inf$BUGSoutput$summary[1:2,1:2], BR.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                 BR.SB.out.inf$BUGSoutput$summary[1:2,1:2], BR.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                 BR.POW.out.inf$BUGSoutput$summary[1:2,1:2], BR.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
BRtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptWAW, ToptEUG, ToptPLA, 
                                ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 9)), nrow = 9, ncol = 2))
BRtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 9)), nrow = 9, ncol = 2))
colnames(BRtops) = c("mean", "sd")
colnames(BRtbreadths) = c("mean", "sd")
BRdf2 = do.call(rbind.data.frame, BRmeansd)
BRdf_1 = rbind(BRdf2, BRtops)
BRdf = rbind(BRdf_1, BRtbreadths)
BRdf$Population = c(rep(c("HOP", "MAR1",  "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                   rep(c("HOP", "MAR1", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
BRdf$Param = c(rep(c("Tmin", "Tmax"), 9), rep("Topt", 9), rep("Tbreadth",9))
save(BRdf, file = "BR_meansd.Rsave")

#### 6a. Plot TPC parameters for each population #####

library(ggplot2)
load("BR_meansd.Rsave")

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB",  "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA",  "MAR1", "POW", "PAR", "SB")
BRdf$Population = factor(BRdf$Population, levels = PopOrder)
BRdf = BRdf[order(BRdf$Population),]
BRdf = BRdf[BRdf$Param != "Tbreadth",]

# colors for each population 
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotBR = ggplot(BRdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Biting rate") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

#### 6b. Plot biting dev rates curve for all populations #####

colors3 = rep(rev(c("#a50026", "#d73027","#f46d43","#fdae61","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB",  "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA",  "MAR1", "POW", "PAR", "SB")

EUGdata = BR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = BR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = BR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = BR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = BR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = BR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = BR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = BR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = BR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(BitingRate ~ Temperature,
     xlim = c(5, 45), ylim = c(0,0.11), data = data.BR.HOP, type = "n", bty = "n",
     ylab = "biting rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Biting rate", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(EUGdata ~ Temp.xs, col = "#313695")
lines(WAWdata ~ Temp.xs, col = "#4575b4")
lines(PLAdata ~ Temp.xs, col = "#74add1")
lines(HOPdata ~ Temp.xs, col = "#abd9e9")
lines(JRAdata ~ Temp.xs, col = "#e0f3f8")
lines(MAR1data ~ Temp.xs, col = "#fdae61")
lines(POWdata ~ Temp.xs, col = "#f46d43")
lines(PARdata ~ Temp.xs, col = "#d73027")
lines(SBdata ~ Temp.xs, col = "#a50026")



