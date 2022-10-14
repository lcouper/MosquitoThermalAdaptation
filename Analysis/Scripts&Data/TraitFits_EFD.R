# Trait Fits for Eggs / Female / Day ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for eggs/female/day (EFD)
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit EFD thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to EFD prior thermal responses
#  5) Fit EFD thermal responses with data-informed priors
#  6) Plot parameters

#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.EFD_int = read.csv("ReproExp_Data.csv")[,-1]
# Calculate EFD for each population and temperature treatment
# Equal to individual EFD (in data sheet) x proportion that oviposited

# step 1: calculate proportion that oviposited for each population
# remove NAs
data.EFD_int = data.EFD_int[!is.na(data.EFD_int$Oviposited),]
data.EFD_int$Oviposited = as.numeric(recode(data.EFD_int$Oviposited, 'y' = 1, 'n' = 0))
data.EFD_oviposited = aggregate(Oviposited ~ Population + Temperature, data.EFD_int, FUN = sum)
data.EFD_total = aggregate(Oviposited ~ Population + Temperature, data.EFD_int, FUN = length)
data.EFD_propovi = cbind.data.frame(data.EFD_oviposited[,1:2], data.EFD_oviposited$Oviposited/ data.EFD_total$Oviposited)
colnames(data.EFD_propovi)[3] = "PropOviposited"

# step 2: calculate average EFD by population and temperature treatment
data.EFD_mean = aggregate(EFD ~ Population + Temperature, data = data.EFD_int, FUN = mean)
data.EFD = merge(data.EFD_propovi, data.EFD_mean, all.x = T)
# convert NAs to 0s
data.EFD$EFD[is.na(data.EFD$EFD)] <-0

# step 3. Multiply average EVG by proportion oviposited
data.EFD$EFD_trait = data.EFD$PropOviposited * data.EFD$EFD

# Plot
plot(EFD_trait ~ Temperature, xlim = c(5, 45), data = data.EFD, ylab = "EFD", xlab = "Temperature")
# Appears left skewed, so will use Briere function (as in Mordecai et al. 2017)

# Subset data by population
data.EFD.HOP <- subset(data.EFD, Population == "HOP")
data.EFD.MAR1 <- subset(data.EFD, Population == "MAR")
# note there is no MAR2 population from the reproduction experiment. 
data.EFD.WAW <- subset(data.EFD, Population == "WAW")
data.EFD.EUG <- subset(data.EFD, Population == "EUG")
data.EFD.PLA <- subset(data.EFD, Population == "PLA")
data.EFD.SB <- subset(data.EFD, Population == "SB")
data.EFD.JRA <- subset(data.EFD, Population == "JRA")
data.EFD.PAR <- subset(data.EFD, Population == "PAR")
data.EFD.POW <- subset(data.EFD, Population == "POW")

#### 2a. Set up JAGS model #####
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

# Briere Model with gamma priors (except sigma)

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

#### 2b. Shared settings for all models ####

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
data.EFD.HOP.prior <- subset(data.EFD, Population != "HOP")
data.EFD.MAR1.prior <- subset(data.EFD, Population != "MAR")
data.EFD.WAW.prior <- subset(data.EFD, Population != "WAW")
data.EFD.EUG.prior <- subset(data.EFD, Population != "EUG")
data.EFD.PLA.prior <- subset(data.EFD, Population != "PLA")
data.EFD.SB.prior <- subset(data.EFD, Population != "SB")
data.EFD.JRA.prior <- subset(data.EFD, Population != "JRA")
data.EFD.PAR.prior <- subset(data.EFD, Population != "PAR")
data.EFD.POW.prior <- subset(data.EFD, Population != "POW")


###### 3a: Hopland (HOP) ######

# set data
data <- data.EFD.HOP.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(EFD.HOP.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.HOP.prior, ylab = "EFD for Hop prior", xlab = "Temperature")
lines(EFD.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b: Marin35 (MAR1) ######

# set data
data <- data.EFD.MAR1.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(EFD.MAR.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.MAR1.prior, ylab = "EFD for MAR prior", xlab = "Temperature")
lines(EFD.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3c: Wawona (WAW) #####

# set data
data <- data.EFD.WAW.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.WAW.prior.out$BUGSoutput$sumWAWy[1:5,]
#mcmcplot(EFD.WAW.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.WAW.prior, ylab = "EFD for WAW prior", xlab = "Temperature")
lines(EFD.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)





###### 3d: Eugene (EUG) #####

# set data
data <- data.EFD.EUG.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.EUG.prior.out$BUGSoutput$sumEUGy[1:5,]
#mcmcplot(EFD.EUG.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.EUG.prior, ylab = "EFD for EUG prior", xlab = "Temperature")
lines(EFD.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e. Placer #####

# set data
data <- data.EFD.PLA.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.PLA.prior.out$BUGSoutput$sumPLAy[1:5,]
#mcmcplot(EFD.PLA.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.PLA.prior, ylab = "EFD for PLA prior", xlab = "Temperature")
lines(EFD.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f. SB #######
# set data
data <- data.EFD.SB.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.SB.prior.out$BUGSoutput$sumSBy[1:5,]
#mcmcplot(EFD.SB.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.SB.prior, ylab = "EFD for SB prior", xlab = "Temperature")
lines(EFD.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3g. Jasper Ridge #####

# set data
data <- data.EFD.JRA.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.JRA.prior.out$BUGSoutput$sumJRAy[1:5,]
#mcmcplot(EFD.JRA.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.JRA.prior, ylab = "EFD for JRA prior", xlab = "Temperature")
lines(EFD.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3h. P


###### 3h. Paso Robles (PAR) ######

# set data
data <- data.EFD.PAR.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.PAR.prior.out$BUGSoutput$sumPARy[1:5,]
#mcmcplot(EFD.PAR.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.PAR.prior, ylab = "EFD for PAR prior", xlab = "Temperature")
lines(EFD.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i. Powder Canyon (POW) ######

# set data
data <- data.EFD.POW.prior

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
EFD.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.POW.prior.out$BUGSoutput$sumPOWy[1:5,]
#mcmcplot(EFD.POW.prior.out)

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 40), ylim = c(0,15), data = data.EFD.POW.prior, ylab = "EFD for POW prior", xlab = "Temperature")
lines(EFD.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


#### 4. Fit gamma distributions to EFD prior thermal responses ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists

###### 4a. Hopland ######

EFD.HOP.prior.cf.dists <- data.frame(q = as.vector(EFD.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.HOP.prior.gamma.fits = apply(EFD.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4b. Marin35 (MAR1) ######

EFD.MAR1.prior.cf.dists <- data.frame(q = as.vector(EFD.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(EFD.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(EFD.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.MAR1.prior.gamma.fits = apply(EFD.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

EFD.WAW.prior.cf.dists <- data.frame(q = as.vector(EFD.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.WAW.prior.gamma.fits = apply(EFD.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

EFD.EUG.prior.cf.dists <- data.frame(q = as.vector(EFD.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.EUG.prior.gamma.fits = apply(EFD.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

EFD.PLA.prior.cf.dists <- data.frame(q = as.vector(EFD.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.PLA.prior.gamma.fits = apply(EFD.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

EFD.SB.prior.cf.dists <- data.frame(q = as.vector(EFD.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(EFD.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(EFD.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.SB.prior.gamma.fits = apply(EFD.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

EFD.JRA.prior.cf.dists <- data.frame(q = as.vector(EFD.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.JRA.prior.gamma.fits = apply(EFD.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

EFD.PAR.prior.cf.dists <- data.frame(q = as.vector(EFD.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.PAR.prior.gamma.fits = apply(EFD.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

EFD.POW.prior.cf.dists <- data.frame(q = as.vector(EFD.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EFD.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EFD.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
EFD.POW.prior.gamma.fits = apply(EFD.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
EFD.hypers <- list(EFD.HOP.prior.gamma.fits, EFD.MAR1.prior.gamma.fits, 
                   EFD.WAW.prior.gamma.fits, EFD.EUG.prior.gamma.fits, EFD.PLA.prior.gamma.fits,
                   EFD.SB.prior.gamma.fits, EFD.JRA.prior.gamma.fits, EFD.PAR.prior.gamma.fits, EFD.POW.prior.gamma.fits)
save(EFD.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/EFDhypers.Rsave")

#### 5. Fit LDR thermal responses with data-informed priors ####

# to bypass above code of generating hyperparameters
load("EFDhypers.Rsave")
EFD.HOP.prior.gamma.fits <- EFD.hypers[[1]]
EFD.MAR1.prior.gamma.fits <- EFD.hypers[[2]]
EFD.WAW.prior.gamma.fits <- EFD.hypers[[3]]
EFD.EUG.prior.gamma.fits <- EFD.hypers[[4]]
EFD.PLA.prior.gamma.fits <- EFD.hypers[[5]]
EFD.SB.prior.gamma.fits <- EFD.hypers[[6]]
EFD.JRA.prior.gamma.fits <- EFD.hypers[[7]]
EFD.PAR.prior.gamma.fits <- EFD.hypers[[8]]
EFD.POW.prior.gamma.fits <- EFD.hypers[[9]]


#### 5. Fit EFD thermal responses with data-informed priors ####

# to bypass above code of generating hyperparameters
load("EFDhypers.Rsave")
EFD.HOP.prior.gamma.fits <- EFD.hypers[[1]]
EFD.MAR1.prior.gamma.fits <- EFD.hypers[[2]]
EFD.MAR2.prior.gamma.fits <- EFD.hypers[[3]]
EFD.WAW.prior.gamma.fits <- EFD.hypers[[4]]
EFD.EUG.prior.gamma.fits <- EFD.hypers[[5]]
EFD.PLA.prior.gamma.fits <- EFD.hypers[[6]]
EFD.SB.prior.gamma.fits <- EFD.hypers[[7]]
EFD.JRA.prior.gamma.fits <- EFD.hypers[[8]]
EFD.PAR.prior.gamma.fits <- EFD.hypers[[9]]
EFD.POW.prior.gamma.fits <- EFD.hypers[[10]]


###### 5a. Hopland ######

# Set data 
data <- data.EFD.HOP
hypers <- EFD.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(EFD.HOP.out.inf)
#save(EFD.HOP.out.inf, file = "jagsout_EFD.HOPp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.HOP, ylab = "EFD for Hopland", xlab = "Temperature", pch = 1)
lines(EFD.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


##### 5b. Marin35 (MAR1) #####

# Set data 
data <- data.EFD.MAR1
hypers <- EFD.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.MAR1.out.inf$BUGSoutput$summary[1:5,]
#save(EFD.MAR1.out.inf, file = "jagsout_EFD.MAR1p_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.MAR1, ylab = "EFD for MAR1", xlab = "Temperature", pch = 1)
lines(EFD.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


##### 5c. Wawona #####

# Set data 
data <- data.EFD.WAW
hypers <- EFD.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.WAW.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(EFD.WAW.out.inf)
#save(EFD.WAW.out.inf, file = "jagsout_EFD.WAWp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.WAW, ylab = "EFD for Wawona", xlab = "Temperature", pch = 1)
lines(EFD.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


##### 5d. Eugene: REMOVE - TWO FEW POINTS #####
# oviposition only occurred at 2 temperatures (18 and 24), so can't fit model

# Set data 
data <- data.EFD.EUG
hypers <- EFD.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.EUG.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(EFD.EUG.out.inf)
#save(EFD.EUG.out.inf, file = "jagsout_EFD.EUG_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.EUG, ylab = "EFD for EUGland", xlab = "Temperature", pch = 1)
lines(EFD.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

##### 5e. Placer #####

# Set data 
data <- data.EFD.PLA
hypers <- EFD.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.PLA.out.inf$BUGSoutput$summary[1:5,]
#save(EFD.PLA.out.inf, file = "jagsout_EFD.PLAp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.PLA, ylab = "EFD for PLA", xlab = "Temperature", pch = 1)
lines(EFD.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

##### 5f. SB #####
# Set data 
data <- data.EFD.SB
hypers <- EFD.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.SB.out.inf$BUGSoutput$summary[1:5,]
#save(EFD.SB.out.inf, file = "jagsout_EFD.SBp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.SB, ylab = "EFD for SB", xlab = "Temperature", pch = 1)
lines(EFD.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

##### 5g. Jasper Ridge #####

# Set data 
data <- data.EFD.JRA
hypers <- EFD.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.JRA.out.inf$BUGSoutput$summary[1:5,]
#save(EFD.JRA.out.inf, file = "jagsout_EFD.JRAp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.JRA, ylab = "EFD for JRA", xlab = "Temperature", pch = 1)
lines(EFD.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


##### 5h. Paso Robles #####

# Set data 
data <- data.EFD.PAR
hypers <- EFD.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.PAR.out.inf$BUGSoutput$summary[1:5,]
#save(EFD.PAR.out.inf, file = "jagsout_EFD.PARp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.PAR, ylab = "EFD for PAR", xlab = "Temperature", pch = 1)
lines(EFD.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

##### 5i. Powder canyon ######
# Set data 
data <- data.EFD.POW
hypers <- EFD.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EFD_trait
N.obs <- length(trait)
temp <- data$Temperature

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
EFD.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
EFD.POW.out.inf$BUGSoutput$summary[1:5,]
#save(EFD.POW.out.inf, file = "jagsout_EFD.POWp_inf.Rdata")

# Plot data + fit
plot(EFD_trait ~ Temperature, xlim = c(5, 45), ylim = c(0,15), data = data.EFD.POW, ylab = "EFD for POW", xlab = "Temperature", pch = 1)
lines(EFD.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### Topt of adult lifespan for each population ######
ToptHOP = Temp.xs[which.max(as.vector(EFD.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 23.5
ToptMAR1 = Temp.xs[which.max(as.vector(EFD.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 23
ToptWAW = Temp.xs[which.max(as.vector(EFD.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16
#ToptEUG = Temp.xs[which.max(as.vector(EFD.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16.5
ToptPLA = Temp.xs[which.max(as.vector(EFD.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 23
ToptSB = Temp.xs[which.max(as.vector(EFD.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 25.5
ToptJRA = Temp.xs[which.max(as.vector(EFD.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19.5
ToptPAR = Temp.xs[which.max(as.vector(EFD.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21.5
ToptPOW = Temp.xs[which.max(as.vector(EFD.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 25



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
TbreadthHOP = Tbreadth(EFD.HOP.out.inf) # 12.5
TbreadthMAR1 = Tbreadth(EFD.MAR1.out.inf) # 13
TbreadthWAW = Tbreadth(EFD.WAW.out.inf) # 9
TbreadthPLA = Tbreadth(EFD.PLA.out.inf) # 12
TbreadthSB = Tbreadth(EFD.SB.out.inf) # 14
TbreadthJRA = Tbreadth(EFD.JRA.out.inf) # 10
TbreadthPAR = Tbreadth(EFD.PAR.out.inf) # 10
TbreadthPOW = Tbreadth(EFD.POW.out.inf) # 7.5

##### Compile prediction data for EFD populations #####

EFDmeansd <- list(EFD.HOP.out.inf$BUGSoutput$summary[1:2,1:2], EFD.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                  EFD.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                  EFD.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                  EFD.SB.out.inf$BUGSoutput$summary[1:2,1:2], EFD.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                  EFD.POW.out.inf$BUGSoutput$summary[1:2,1:2], EFD.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
EFDtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptWAW, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 8)), nrow = 8, ncol = 2))
EFDtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthWAW, TbreadthPLA, 
                                      TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 8)), nrow = 8, ncol = 2))
colnames(EFDtops) = c("mean", "sd")
colnames(EFDtbreadths) = c("mean", "sd")
EFDdf2 = do.call(rbind.data.frame, EFDmeansd)
EFDdf_1 = rbind(EFDdf2, EFDtops)
EFDdf = rbind(EFDdf_1, EFDtbreadths)
EFDdf$Population = c(rep(c("HOP", "MAR1", "WAW",  "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                     rep(c("HOP", "MAR1", "WAW",  "PLA", "SB", "JRA", "PAR", "POW"), 2))
EFDdf$Param = c(rep(c("Tmin", "Tmax"), 8), rep("Topt", 8), rep("Tbreadth", 8))
save(EFDdf, file = "EFD_meansd.Rsave")


#### Plot TPC parameters for each population #####

load("EFD_meansd.Rsave")
library(ggplot2)

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("WAW","HOP", "PLA", "SB", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("WAW", "PLA", "HOP", "JRA", "MAR1", "POW", "PAR", "SB")
EFDdf$Population = factor(EFDdf$Population, levels = PopOrder)
EFDdf = EFDdf[order(EFDdf$Population),]
EFDdf = EFDdf[EFDdf$Param != "Tbreadth",]


# colors for each population
colors3 = rep(rev(c("#d73027","#f46d43","#fdae61","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)



plotEFD = ggplot(EFDdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, size = 2,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(), width = 0.3) + 
  ggtitle("EFD") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

#### 6b. Plot EFD curve for all populaions #####

colors3 = rep(rev(c("#d73027","#f46d43","#fdae61","#fee090","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("WAW","HOP", "PLA", "SB", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("WAW", "PLA", "HOP", "JRA", "MAR1", "POW", "PAR", "SB") 

WAWdata = EFD.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = EFD.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = EFD.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = EFD.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = EFD.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = EFD.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = EFD.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = EFD.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(EFD_trait ~ Temperature, 
     xlim = c(5, 45), ylim = c(0,8), data = data.EFD.HOP, type = "n", bty = "n",
     ylab = "eggs / female / day", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "EFD", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(WAWdata ~ Temp.xs, col = "#4575b4")
lines(PLAdata ~ Temp.xs, col = "#74add1")
lines(HOPdata ~ Temp.xs, col = "#abd9e9")
lines(JRAdata ~ Temp.xs, col = "#e0f3f8")
lines(MAR1data ~ Temp.xs, col = "#fdae61")
lines(POWdata ~ Temp.xs, col = "#f46d43")
lines(PARdata ~ Temp.xs, col = "#d73027")
lines(SBdata ~ Temp.xs, col = "#a50026")



