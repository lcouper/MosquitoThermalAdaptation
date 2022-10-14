# Trait Fits for Probability Egg to Adult Survival ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for eggs/female/day (EFD)
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit pEA thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to pEA prior thermal responses
#  5) Fit pEA thermal responses with data-informed priors
#  6) Plot parameters

#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.pEA = read.csv("MosquitoDensity_Data.csv", header = T)[,-1]
# remove 5C treatment since no pupae developed

# Plot
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), data = data.pEA, ylab = "probability", xlab = "Temperature")

# Subset data by population
data.pEA.HOP <- subset(data.pEA, Population == "HOP")
data.pEA.MAR <- subset(data.pEA, Population == "MAR")
data.pEA.WAW <- subset(data.pEA, Population == "WAW")
data.pEA.EUG <- subset(data.pEA, Population == "EUG")
data.pEA.PLA <- subset(data.pEA, Population == "PLA")
data.pEA.SB <- subset(data.pEA, Population == "SB")
data.pEA.JRA <- subset(data.pEA, Population == "JRA")
data.pEA.PAR <- subset(data.pEA, Population == "PAR")
data.pEA.POW <- subset(data.pEA, Population == "POW")

#### 2a. JAGS Models ####

# Quadratic Model with uniform priors

sink("quad.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 50)
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

#### 2b. Shared settings for all models ####

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






#### 3. Generate data-informed priors using leave-one-out approach  ####

# Subset Data - leave-one-out
data.pEA.HOP.prior <- subset(data.pEA, Population != "HOP")
data.pEA.MAR1.prior <- subset(data.pEA, Population != "MAR")
data.pEA.WAW.prior <- subset(data.pEA, Population != "WAW")
data.pEA.EUG.prior <- subset(data.pEA, Population != "EUG")
data.pEA.PLA.prior <- subset(data.pEA, Population != "PLA")
data.pEA.SB.prior <- subset(data.pEA, Population != "SB")
data.pEA.JRA.prior <- subset(data.pEA, Population != "JRA")
data.pEA.PAR.prior <- subset(data.pEA, Population != "PAR")
data.pEA.POW.prior <- subset(data.pEA, Population != "POW")

###### 3a. Hopland (HOP) ####

# Set data
data <- data.pEA.HOP.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.HOP.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.HOP.prior, ylab = "probability for Hopland", xlab = "Temperature")
lines(pEA.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b. Marin35 (MAR) ####

# Set data
data <- data.pEA.MAR1.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.MAR1.prior, ylab = "probability for Marin35", xlab = "Temperature")
lines(pEA.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3d: Wawona (WAW) ######

# set data
data <- data.pEA.WAW.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.WAW.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.WAW.prior, ylab = "probability for WAW prior", xlab = "Temperature")
lines(pEA.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.pEA.EUG.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pEA.EUG.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.EUG.prior, ylab = "probability for EUG prior", xlab = "Temperature")
lines(pEA.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.pEA.PLA.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pEA.PLA.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.PLA.prior, ylab = "probability for PLA prior", xlab = "Temperature")
lines(pEA.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.pEA.SB.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pEA.SB.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.SB.prior, ylab = "probability for SB prior", xlab = "Temperature")
lines(pEA.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.pEA.JRA.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pEA.JRA.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.JRA.prior, ylab = "probability for JRA prior", xlab = "Temperature")
lines(pEA.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.pEA.PAR.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pEA.PAR.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.PAR.prior, ylab = "probability for PAR prior", xlab = "Temperature")
lines(pEA.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

# set data
data <- data.pEA.POW.prior

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pEA.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pEA.POW.prior.out)

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.pEA.POW.prior, ylab = "probability for POW prior", xlab = "Temperature")
lines(pEA.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




#### 4. Fit gamma distributions to pEA prior thermal responses ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists
###### 4a. Hopland ######
pEA.HOP.prior.cf.dists <- data.frame(q = as.vector(pEA.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))

pEA.HOP.prior.gamma.fits = apply(pEA.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4b. Marin35 (MAR1) ######

pEA.MAR1.prior.cf.dists <- data.frame(q = as.vector(pEA.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pEA.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pEA.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.MAR1.prior.gamma.fits = apply(pEA.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

pEA.WAW.prior.cf.dists <- data.frame(q = as.vector(pEA.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.WAW.prior.gamma.fits = apply(pEA.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

pEA.EUG.prior.cf.dists <- data.frame(q = as.vector(pEA.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.EUG.prior.gamma.fits = apply(pEA.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

pEA.PLA.prior.cf.dists <- data.frame(q = as.vector(pEA.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.PLA.prior.gamma.fits = apply(pEA.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

pEA.SB.prior.cf.dists <- data.frame(q = as.vector(pEA.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(pEA.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(pEA.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.SB.prior.gamma.fits = apply(pEA.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

pEA.JRA.prior.cf.dists <- data.frame(q = as.vector(pEA.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.JRA.prior.gamma.fits = apply(pEA.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

pEA.PAR.prior.cf.dists <- data.frame(q = as.vector(pEA.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.PAR.prior.gamma.fits = apply(pEA.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

pEA.POW.prior.cf.dists <- data.frame(q = as.vector(pEA.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pEA.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pEA.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
pEA.POW.prior.gamma.fits = apply(pEA.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list #####
pEA.hypers <- list(pEA.HOP.prior.gamma.fits, pEA.MAR1.prior.gamma.fits,
                   pEA.WAW.prior.gamma.fits, pEA.EUG.prior.gamma.fits, pEA.PLA.prior.gamma.fits,
                   pEA.SB.prior.gamma.fits, pEA.JRA.prior.gamma.fits, pEA.PAR.prior.gamma.fits, pEA.POW.prior.gamma.fits)
save(pEA.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/pEAhypers.Rsave")


#### 5. Fit pEA thermal responses with data-informed priors #####

# to bypass above code of generating hyperparameters
#load("pEAhypers.Rsave")
pEA.HOP.prior.gamma.fits <- pEA.hypers[[1]]
pEA.MAR1.prior.gamma.fits <- pEA.hypers[[2]]
pEA.WAW.prior.gamma.fits <- pEA.hypers[[3]]
pEA.EUG.prior.gamma.fits <- pEA.hypers[[4]]
pEA.PLA.prior.gamma.fits <- pEA.hypers[[5]]
pEA.SB.prior.gamma.fits <- pEA.hypers[[6]]
pEA.JRA.prior.gamma.fits <- pEA.hypers[[7]]
pEA.PAR.prior.gamma.fits <- pEA.hypers[[8]]
pEA.POW.prior.gamma.fits <- pEA.hypers[[9]]

###### 5a. Hopland ######

# Set data 
data <- data.pEA.HOP
hypers <- pEA.HOP.prior.gamma.fits * .1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.HOP.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2.0), data = data.pEA.HOP, ylab = "probability for Hopland", xlab = "Temperature", pch = 1)
lines(pEA.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### 5b. Marin35 ######

# Set data 
data <- data.pEA.MAR
hypers <- pEA.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.MAR, ylab = "probability for Marin35", xlab = "Temperature", pch = 1)
lines(pEA.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5d. Wawona ######

# Set data 
data <- data.pEA.WAW
hypers <- pEA.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.WAW, ylab = "probability for Wawona", xlab = "Temperature", pch = 1)
lines(pEA.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5e. Eugene ######

# Set data 
data <- data.pEA.EUG
hypers <- pEA.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.EUG, ylab = "probability for Eugene", xlab = "Temperature", pch = 1)
lines(pEA.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5f. Placer ######

# Set data 
data <- data.pEA.PLA
hypers <- pEA.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.PLA, ylab = "probability for Placer", xlab = "Temperature", pch = 1)
lines(pEA.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5g. Santa Barbara ######

# Set data 
data <- data.pEA.SB
hypers <- pEA.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.SB, ylab = "probability for Santa Barbara", xlab = "Temperature", pch = 1)
lines(pEA.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5g. Jasper Ridge ######

# Set data 
data <- data.pEA.JRA
hypers <- pEA.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.JRA, ylab = "probability for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(pEA.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5i. Paso Robles ######

# Set data 
data <- data.pEA.PAR
hypers <- pEA.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.PAR, ylab = "probability for Paso Robles", xlab = "Temperature", pch = 1)
lines(pEA.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5j. Powder canyon ######

# Set data 
data <- data.pEA.POW
hypers <- pEA.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.EggAdult.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pEA.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pEA.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.EggAdult.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pEA.POW, ylab = "probability for Powder Canyon", xlab = "Temperature", pch = 1)
lines(pEA.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)





###### Topt of pEA for each population ######
ToptHOP = Temp.xs[which.max(as.vector(pEA.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 23.5
ToptMAR1 = Temp.xs[which.max(as.vector(pEA.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 22
ToptWAW = Temp.xs[which.max(as.vector(pEA.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21.5
ToptEUG = Temp.xs[which.max(as.vector(pEA.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptPLA = Temp.xs[which.max(as.vector(pEA.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21
ToptSB = Temp.xs[which.max(as.vector(pEA.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21.5
ToptJRA = Temp.xs[which.max(as.vector(pEA.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 17.5
ToptPAR = Temp.xs[which.max(as.vector(pEA.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19
ToptPOW = Temp.xs[which.max(as.vector(pEA.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 17

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
TbreadthHOP = Tbreadth(pEA.HOP.out.inf) # 19.5
TbreadthMAR1 = Tbreadth(pEA.MAR1.out.inf) # 19.5
TbreadthWAW = Tbreadth(pEA.WAW.out.inf) # 19.5
TbreadthEUG = Tbreadth(pEA.EUG.out.inf) # 17.5
TbreadthPLA = Tbreadth(pEA.PLA.out.inf) # 19
TbreadthSB = Tbreadth(pEA.SB.out.inf) # 20
TbreadthJRA = Tbreadth(pEA.JRA.out.inf) # 20
TbreadthPAR = Tbreadth(pEA.PAR.out.inf) # 19.5
TbreadthPOW = Tbreadth(pEA.POW.out.inf) # 19.5

#### Compile prediction data for all populations ###### 

pEAmeansd <- list(pEA.HOP.out.inf$BUGSoutput$summary[1:2,1:2], pEA.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                  pEA.WAW.out.inf$BUGSoutput$summary[1:2,1:2], pEA.EUG.out.inf$BUGSoutput$summary[1:2,1:2], 
                  pEA.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                  pEA.SB.out.inf$BUGSoutput$summary[1:2,1:2], pEA.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                  pEA.POW.out.inf$BUGSoutput$summary[1:2,1:2], pEA.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
pEAtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptWAW, ToptEUG, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 9)), nrow = 9, ncol = 2))
pEAtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                    TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 9)), nrow = 9, ncol = 2))
colnames(pEAtops) = c("mean", "sd")
colnames(pEAtbreadths) = c("mean", "sd")
pEAdf2 = do.call(rbind.data.frame, pEAmeansd)
pEAdf_1 = rbind(pEAdf2, pEAtops)
pEAdf = rbind(pEAdf_1, pEAtbreadths)
pEAdf$Population = c(rep(c("HOP", "MAR1",  "WAW",  "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                   rep(c("HOP", "MAR1", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
pEAdf$Param = c(rep(c("Tmin", "Tmax"), 9), rep("Topt", 9), rep("Tbreadth", 9))
save(pEAdf, file = "pEA_meansd.Rsave")


#### Plot TPC parameters for each population #####

load("pEA_meansd.Rsave")
library(ggplot2)

# order from coldest to hottest temp in wettest quarter
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR1", "POW", "PAR", "SB")
pEAdf$Population = factor(pEAdf$Population, levels = PopOrder)
pEAdf = pEAdf[order(pEAdf$Population),]
pEAdf = pEAdf[pEAdf$Param != "Tbreadth",]

# colors for each population
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotpEA = ggplot(pEAdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, size = 2,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(), width = 0.3) + 
  ggtitle("Probability egg to adult survival") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

#### Plot Pupal survival curve for all populations #####

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#e0f3f8", "#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA",  "MAR1", "POW", "PAR", "SB")

plot(prob.EggAdult.Survival ~ Temp.Treatment, 
     xlim = c(0, 45), ylim = c(0,1.1), data = data.pEA.HOP, type = "n", bty = "n",
     ylab = "probability", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Probability egg to adult survival", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(pEA.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#313695", lwd = 2.5)
lines(pEA.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#4575b4", lwd = 2.5) 
lines(pEA.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#74add1", lwd = 2.5)
lines(pEA.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#abd9e9", lwd = 2.5)
lines(pEA.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#e0f3f8", lwd = 2.5)
lines(pEA.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fdae61", lwd = 2.5)
lines(pEA.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#f46d43", lwd = 2.5)
lines(pEA.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#d73027", lwd = 2.5)
lines(pEA.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#a50026", lwd = 2.5)

