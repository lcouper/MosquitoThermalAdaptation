# Trait Fits for Mosquito Density ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for eggs/female/day (EFD)
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit M thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to M prior thermal responses
#  5) Fit M thermal responses with data-informed priors
#  6) Plot parameters

#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.M = read.csv("MosquitoDensity_Data.csv", header = T)[,-1]

# convert NAs to 0s
data.M$M[is.na(data.M$M)] <-0 


# Plot
plot(M ~ Temp.Treatment, xlim = c(5, 45), data = data.M, ylab = "Mosquito Density", xlab = "Temperature")
# Appears quadratic

# Subset data by population
data.M.HOP <- subset(data.M, Population == "HOP")
data.M.MAR <- subset(data.M, Population == "MAR")
data.M.WAW <- subset(data.M, Population == "WAW")
data.M.EUG <- subset(data.M, Population == "EUG")
data.M.PLA <- subset(data.M, Population == "PLA")
data.M.SB <- subset(data.M, Population == "SB")
data.M.JRA <- subset(data.M, Population == "JRA")
data.M.PAR <- subset(data.M, Population == "PAR")
data.M.POW <- subset(data.M, Population == "POW")

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
data.M.HOP.prior <- subset(data.M, Population != "HOP")
data.M.MAR.prior <- subset(data.M, Population != "MAR")
data.M.WAW.prior <- subset(data.M, Population != "WAW")
data.M.EUG.prior <- subset(data.M, Population != "EUG")
data.M.PLA.prior <- subset(data.M, Population != "PLA")
data.M.SB.prior <- subset(data.M, Population != "SB")
data.M.JRA.prior <- subset(data.M, Population != "JRA")
data.M.PAR.prior <- subset(data.M, Population != "PAR")
data.M.POW.prior <- subset(data.M, Population != "POW")

###### 3a. Hopland (HOP) ####

# Set data
data <- data.M.HOP.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.HOP.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.HOP.prior, ylab = "M for Hopland", xlab = "Temperature")
lines(M.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b. Marin35 (MAR) ####

# Set data
data <- data.M.MAR.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.MAR.prior, ylab = "M for Marin35", xlab = "Temperature")
lines(M.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3d: Wawona (WAW) ######

# set data
data <- data.M.WAW.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.WAW.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.WAW.prior, ylab = "M for WAW prior", xlab = "Temperature")
lines(M.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.M.EUG.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(M.EUG.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.EUG.prior, ylab = "M for EUG prior", xlab = "Temperature")
lines(M.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.M.PLA.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(M.PLA.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.PLA.prior, ylab = "M for PLA prior", xlab = "Temperature")
lines(M.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.M.SB.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(M.SB.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.SB.prior, ylab = "M for SB prior", xlab = "Temperature")
lines(M.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.M.JRA.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(M.JRA.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.JRA.prior, ylab = "M for JRA prior", xlab = "Temperature")
lines(M.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.M.PAR.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(M.PAR.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.PAR.prior, ylab = "M for PAR prior", xlab = "Temperature")
lines(M.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

# set data
data <- data.M.POW.prior

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
M.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(M.POW.prior.out)

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2), data = data.M.POW.prior, ylab = "M for POW prior", xlab = "Temperature")
lines(M.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




#### 4. Fit gamma distributions to M prior thermal responses ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists

###### 4a. Hopland ######
M.HOP.prior.cf.dists <- data.frame(q = as.vector(M.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))

M.HOP.prior.gamma.fits = apply(M.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


###### 4b. Marin35 (MAR1) ######

M.MAR1.prior.cf.dists <- data.frame(q = as.vector(M.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(M.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(M.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
M.MAR1.prior.gamma.fits = apply(M.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

M.WAW.prior.cf.dists <- data.frame(q = as.vector(M.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
M.WAW.prior.gamma.fits = apply(M.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

M.EUG.prior.cf.dists <- data.frame(q = as.vector(M.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
M.EUG.prior.gamma.fits = apply(M.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

M.PLA.prior.cf.dists <- data.frame(q = as.vector(M.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
M.PLA.prior.gamma.fits = apply(M.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

M.SB.prior.cf.dists <- data.frame(q = as.vector(M.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                  T0 = as.vector(M.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                  Tm = as.vector(M.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
M.SB.prior.gamma.fits = apply(M.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

M.JRA.prior.cf.dists <- data.frame(q = as.vector(M.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
M.JRA.prior.gamma.fits = apply(M.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

M.PAR.prior.cf.dists <- data.frame(q = as.vector(M.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
M.PAR.prior.gamma.fits = apply(M.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

M.POW.prior.cf.dists <- data.frame(q = as.vector(M.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(M.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(M.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
M.POW.prior.gamma.fits = apply(M.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
M.hypers <- list(M.HOP.prior.gamma.fits, M.MAR1.prior.gamma.fits,
                 M.WAW.prior.gamma.fits, M.EUG.prior.gamma.fits, M.PLA.prior.gamma.fits,
                 M.SB.prior.gamma.fits, M.JRA.prior.gamma.fits, M.PAR.prior.gamma.fits, M.POW.prior.gamma.fits)
save(M.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/Mhypers.Rsave")


#### 5. Fit M thermal responses with data-informed priors #####

# to bypass above code of generating hyperparameters
#load("Mhypers.Rsave")
M.HOP.prior.gamma.fits <- M.hypers[[1]]
M.MAR1.prior.gamma.fits <- M.hypers[[2]]
M.WAW.prior.gamma.fits <- M.hypers[[3]]
M.EUG.prior.gamma.fits <- M.hypers[[4]]
M.PLA.prior.gamma.fits <- M.hypers[[5]]
M.SB.prior.gamma.fits <- M.hypers[[6]]
M.JRA.prior.gamma.fits <- M.hypers[[7]]
M.PAR.prior.gamma.fits <- M.hypers[[8]]
M.POW.prior.gamma.fits <- M.hypers[[9]]

###### 5a. Hopland ######

# Set data 
data <- data.M.HOP
hypers <- M.HOP.prior.gamma.fits * .1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.HOP.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,2.0), data = data.M.HOP, ylab = "mosquito density for Hopland", xlab = "Temperature", pch = 1)
lines(M.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5b. Marin35 ######

# Set data 
data <- data.M.MAR
hypers <- M.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.MAR, ylab = "M for Marin35", xlab = "Temperature", pch = 1)
lines(M.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5d. Wawona ######

# Set data 
data <- data.M.WAW
hypers <- M.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.WAW, ylab = "M for Wawona", xlab = "Temperature", pch = 1)
lines(M.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5e. Eugene (not enough points) ######

# Set data 
data <- data.M.EUG
hypers <- M.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.EUG, ylab = "M for Eugene", xlab = "Temperature", pch = 1)
lines(M.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5f. Placer ######

# Set data 
data <- data.M.PLA
hypers <- M.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.PLA, ylab = "M for Placer", xlab = "Temperature", pch = 1)
lines(M.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5g. Santa Barbara ######

# Set data 
data <- data.M.SB
hypers <- M.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.SB, ylab = "M for Santa Barbara", xlab = "Temperature", pch = 1)
lines(M.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5g. Jasper Ridge ######

# Set data 
data <- data.M.JRA
hypers <- M.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.JRA, ylab = "M for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(M.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5i. Paso Robles ######

# Set data 
data <- data.M.PAR
hypers <- M.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.PAR, ylab = "M for Paso Robles", xlab = "Temperature", pch = 1)
lines(M.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5j. Powder canyon ######

# Set data 
data <- data.M.POW
hypers <- M.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$M
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
M.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
M.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(M ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.M.POW, ylab = "M for Powder Canyon", xlab = "Temperature", pch = 1)
lines(M.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(M.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(M.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### Topt of M for each population ######
ToptHOP = Temp.xs[which.max(as.vector(M.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 23.5
ToptMAR1 = Temp.xs[which.max(as.vector(M.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 22
ToptWAW = Temp.xs[which.max(as.vector(M.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21.5
#ToptEUG = Temp.xs[which.max(as.vector(M.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptPLA = Temp.xs[which.max(as.vector(M.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21
ToptSB = Temp.xs[which.max(as.vector(M.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 21.5
ToptJRA = Temp.xs[which.max(as.vector(M.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 17.5
ToptPAR = Temp.xs[which.max(as.vector(M.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19
ToptPOW = Temp.xs[which.max(as.vector(M.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 17

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
TbreadthHOP = Tbreadth(M.HOP.out.inf) # 1.5
TbreadthMAR1 = Tbreadth(M.MAR1.out.inf) # 7
TbreadthWAW = Tbreadth(M.WAW.out.inf) # 0.5
#TbreadthEUG = Tbreadth(M.EUG.out.inf) # 
TbreadthPLA = Tbreadth(M.PLA.out.inf) # 3
TbreadthSB = Tbreadth(M.SB.out.inf) # 1
TbreadthJRA = Tbreadth(M.JRA.out.inf) # 2.5
TbreadthPAR = Tbreadth(M.PAR.out.inf) # 8
TbreadthPOW = Tbreadth(M.POW.out.inf) # 3


#### Compile prediction data for all populations ###### 

Mmeansd <- list(M.HOP.out.inf$BUGSoutput$summary[1:2,1:2], M.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                M.WAW.out.inf$BUGSoutput$summary[1:2,1:2],  M.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                M.SB.out.inf$BUGSoutput$summary[1:2,1:2], M.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                M.POW.out.inf$BUGSoutput$summary[1:2,1:2], M.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
Mtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptWAW, ToptPLA, 
                               ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 8)), nrow = 8, ncol = 2))
Mtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1, TbreadthWAW, TbreadthPLA, 
                                    TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 8)), nrow = 8, ncol = 2))
colnames(Mtops) = c("mean", "sd")
colnames(Mtbreadths) = c("mean", "sd")
Mdf2 = do.call(rbind.data.frame, Mmeansd)
Mdf_1 = rbind(Mdf2, Mtops)
Mdf = rbind(Mdf_1, Mtbreadths)
Mdf$Population = c(rep(c("HOP", "MAR1",  "WAW",  "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                   rep(c("HOP", "MAR1", "WAW",  "PLA", "SB", "JRA", "PAR", "POW"), 2))
Mdf$Param = c(rep(c("Tmin", "Tmax"), 8), rep("Topt", 8), rep("Tbreadth", 8))
save(Mdf, file = "M_meansd.Rsave")

#### Plot TPC parameters for each population #####

load("M_meansd.Rsave")
library(ggplot2)

# order from coldest to hottest temp in wettest quarter
PopOrder = c("WAW", "PLA", "HOP", "JRA", "MAR1", "POW", "PAR", "SB")
Mdf$Population = factor(Mdf$Population, levels = PopOrder)
Mdf = Mdf[order(Mdf$Population),]
Mdf = Mdf[Mdf$Param != "Tbreadth",]

# colors for each population
colors3 = rep(rev(c("#d73027","#f46d43","#fdae61","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)


plotM = ggplot(Mdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 40)) +
  geom_point(stat="identity", col=colors3, size = 2,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(), width = 0.3) + 
  ggtitle("Mosquito density") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

#### Plot mosquito density curve for all populations #####

colors3 = rep(rev(c("#d73027","#f46d43","#fdae61","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
PopOrder = c("WAW", "PLA", "HOP", "JRA", "MAR1", "POW", "PAR", "SB")

#### Unscaled plot #####
plot(M ~ Temp.Treatment, 
     xlim = c(5, 35), ylim = c(0,1.2), data = data.M.HOP, type = "n", bty = "n",
     ylab = " ", xlab = "Temperature (\u00B0C)", pch = 1, cex = 1.5,
     main = "Mosquito density", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(M.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#4575b4", lwd = 2.5)
lines(M.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#74add1",lwd = 2.5)
lines(M.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#abd9e9", lwd = 2.5)
lines(M.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#e0f3f8", lwd = 2.5)
lines(M.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fdae61", lwd = 2.5)
lines(M.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#f46d43", lwd = 2.5)
lines(M.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#d73027", lwd = 2.5)
lines(M.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#a50026", lwd = 2.5)


