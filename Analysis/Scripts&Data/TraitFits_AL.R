# Trait Fits for Adult Lifespan #

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for adult lifespan (AL) for 
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit adult life span thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to adult lifespan prior thermal responses
#  5) Fit adult lifespan thermal responses with data-informed priors
#  6) Plot parameters


#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.AL = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# Keep only rows with AL data
data.AL = data.AL[!is.na(data.AL$AdultLifespan),]

# Subset data by population
data.AL.HOP <- subset(data.AL, Population == "HOP")
data.AL.MAR1 <- subset(data.AL, Population == "MAR35")
data.AL.MAR2 <- subset(data.AL, Population == "MAR29")
data.AL.WAW <- subset(data.AL, Population == "WAW")
data.AL.EUG <- subset(data.AL, Population == "EUG")
data.AL.PLA <- subset(data.AL, Population == "PLA")
data.AL.SB <- subset(data.AL, Population == "SB")
data.AL.JRA <- subset(data.AL, Population == "JRA")
data.AL.PAR <- subset(data.AL, Population == "PAR")
data.AL.POW <- subset(data.AL, Population == "POW")

# Plot trait data
plot(AdultLifespan ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,40), pch = 16,
     data = data.AL, ylab = "Adult lifespan", xlab = "Temperature")

#### 2a. Set up JAGS Models ####
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
data.AL.HOP.prior <- subset(data.AL, Population != "HOP")
data.AL.MAR1.prior <- subset(data.AL, Population != "MAR35")
data.AL.MAR2.prior <- subset(data.AL, Population != "MAR29")
data.AL.WAW.prior <- subset(data.AL, Population != "WAW")
data.AL.EUG.prior <- subset(data.AL, Population != "EUG")
data.AL.PLA.prior <- subset(data.AL, Population != "PLA")
data.AL.SB.prior <- subset(data.AL, Population != "SB")
data.AL.JRA.prior <- subset(data.AL, Population != "JRA")
data.AL.PAR.prior <- subset(data.AL, Population != "PAR")
data.AL.POW.prior <- subset(data.AL, Population != "POW")

###### 3a: Hopland (HOP) ######

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
lines(AL.HOP.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b: Marin35 (MAR1) ######

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
lines(AL.MAR1.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3c: Marin29 (MAR2) ######

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
lines(AL.MAR2.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3d: Wawona (WAW) ######

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
lines(AL.WAW.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

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
lines(AL.EUG.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

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
lines(AL.PLA.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

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
lines(AL.SB.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

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
lines(AL.JRA.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

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
lines(AL.PAR.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

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
lines(AL.POW.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)













#### 4. Fit gamma distributions to adult lifespan prior thermal responses ####
###### 4a. Hopland ######

AL.HOP.prior.cf.dists <- data.frame(q = as.vector(AL.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(AL.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))

pPS.HOP.prior.gamma.fits = apply(pPS.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

##### 4b. Marin35 (MAR1) ######

AL.MAR1.prior.cf.dists <- data.frame(q = as.vector(AL.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(AL.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(AL.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.MAR1.prior.gamma.fits = apply(AL.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4c. Marin29 (MAR2) ######

AL.MAR2.prior.cf.dists <- data.frame(q = as.vector(AL.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(AL.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(AL.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.MAR2.prior.gamma.fits = apply(AL.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

AL.WAW.prior.cf.dists <- data.frame(q = as.vector(AL.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.WAW.prior.gamma.fits = apply(AL.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

AL.EUG.prior.cf.dists <- data.frame(q = as.vector(AL.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.EUG.prior.gamma.fits = apply(AL.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

AL.PLA.prior.cf.dists <- data.frame(q = as.vector(AL.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.PLA.prior.gamma.fits = apply(AL.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

AL.SB.prior.cf.dists <- data.frame(q = as.vector(AL.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                   T0 = as.vector(AL.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                   Tm = as.vector(AL.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.SB.prior.gamma.fits = apply(AL.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

AL.JRA.prior.cf.dists <- data.frame(q = as.vector(AL.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.JRA.prior.gamma.fits = apply(AL.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

AL.PAR.prior.cf.dists <- data.frame(q = as.vector(AL.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.PAR.prior.gamma.fits = apply(AL.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

AL.POW.prior.cf.dists <- data.frame(q = as.vector(AL.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(AL.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(AL.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
AL.POW.prior.gamma.fits = apply(AL.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
AL.hypers <- list(AL.HOP.prior.gamma.fits, AL.MAR1.prior.gamma.fits, AL.MAR2.prior.gamma.fits,
                  AL.WAW.prior.gamma.fits, AL.EUG.prior.gamma.fits, AL.PLA.prior.gamma.fits,
                  AL.SB.prior.gamma.fits, AL.JRA.prior.gamma.fits, AL.PAR.prior.gamma.fits, AL.POW.prior.gamma.fits)
save(AL.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/ALhypers.Rsave")


#### 5. Fit adult life span thermal responses with data-informed priors ####

# to bypass above code of generating hyperparameters
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
###### 5a. Hopland ######

# Set data 
data <- data.AL.HOP
hypers <- AL.HOP.prior.gamma.fits * .1

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

#save(pPS.HOP.out.inf, file = "jagsout_pPS_HOP_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.HOP, ylab = "prob Pupal survival for Hopland", xlab = "Temperature", pch = 1)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 5b. Marin35 #####

# Set data 
data <- data.AL.MAR1
hypers <- AL.MAR1.prior.gamma.fits * .1

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
#save(pPS.MAR1.out.inf, file = "jagsout_pPS_MAR1_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.MAR1, ylab = "prob Pupal survival for MAR1land", xlab = "Temperature", pch = 1)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 5c. Marin29 #####

# Set data 
data <- data.AL.MAR2
hypers <- AL.MAR2.prior.gamma.fits * .1

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
#save(pPS.MAR2.out.inf, file = "jagsout_pPS_MAR2_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.MAR2, ylab = "prob Pupal survival for MAR2land", xlab = "Temperature", pch = 1)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 5d. Wawona (WAW) ####

# Set data 
data <- data.AL.WAW
hypers <- AL.WAW.prior.gamma.fits * .1

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

#save(pPS.WAW.out.inf, file = "jagsout_pPS_WAW_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.WAW, ylab = "prob Pupal survival for WAWland", xlab = "Temperature", pch = 1)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##### 5e. Eugene (EUG) #####

# Set data 
data <- data.AL.EUG
hypers <- AL.EUG.prior.gamma.fits * .1

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

#save(pPS.EUG.out.inf, file = "jagsout_pPS_EUG_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.EUG, ylab = "prob Pupal survival for EUGland", xlab = "Temperature", pch = 1)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##### 5f. Placer ######

# Set data 
data <- data.AL.PLA
hypers <- AL.PLA.prior.gamma.fits * .1

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

#save(pPS.PLA.out.inf, file = "jagsout_pPS_PLA_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.PLA, ylab = "prob Pupal survival for PLAland", xlab = "Temperature", pch = 1)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




##### 5g. Santa Barbara ######
# Set data 
data <- data.AL.SB
hypers <- AL.SB.prior.gamma.fits * .1

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

#save(pPS.SB.out.inf, file = "jagsout_pPS_SB_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.SB, ylab = "prob Pupal survival for SBland", xlab = "Temperature", pch = 1)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 5h. Jasper Ridge #####

# Set data 
data <- data.AL.JRA
hypers <- AL.JRA.prior.gamma.fits * .1

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

#save(pPS.JRA.out.inf, file = "jagsout_pPS_JRA_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.JRA, ylab = "prob Pupal survival for JRAland", xlab = "Temperature", pch = 1)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 5i. Paso Robles #####

# Set data 
data <- data.AL.PAR
hypers <- AL.PAR.prior.gamma.fits * .1

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

#save(pPS.PAR.out.inf, file = "jagsout_pPS_PAR_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.PAR, ylab = "prob Pupal survival for PARland", xlab = "Temperature", pch = 1)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##### 5j. Powder canyon #####

# Set data 
data <- data.AL.POW
hypers <- AL.POW.prior.gamma.fits * .1

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
#save(pPS.POW.out.inf, file = "jagsout_pPS_POW_inf.Rdata")

# Plot data + fit
plot(AdultLifespan ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,35), data = data.AL.POW, ylab = "prob Pupal survival for POWland", xlab = "Temperature", pch = 1)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### Topt of adult lifespan for each population ######
ToptHOP = Temp.xs[which.max(as.vector(AL.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16
ToptMAR1 = Temp.xs[which.max(as.vector(AL.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 15.5
ToptMAR2 = Temp.xs[which.max(as.vector(AL.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16
ToptWAW = Temp.xs[which.max(as.vector(AL.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16.5
ToptEUG = Temp.xs[which.max(as.vector(AL.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16.5
ToptPLA = Temp.xs[which.max(as.vector(AL.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 15.5
ToptSB = Temp.xs[which.max(as.vector(AL.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16
ToptJRA = Temp.xs[which.max(as.vector(AL.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16.5
ToptPAR = Temp.xs[which.max(as.vector(AL.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16
ToptPOW = Temp.xs[which.max(as.vector(AL.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 16.5

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
TbreadthHOP = Tbreadth(AL.HOP.out.inf) # 22.5
TbreadthMAR1 = Tbreadth(AL.MAR1.out.inf) # 21
TbreadthMAR2 = Tbreadth(AL.MAR2.out.inf) # 22
TbreadthWAW = Tbreadth(AL.WAW.out.inf) # 22
TbreadthEUG = Tbreadth(AL.EUG.out.inf) # 22.5
TbreadthPLA = Tbreadth(AL.PLA.out.inf) # 21
TbreadthSB = Tbreadth(AL.SB.out.inf) # 22
TbreadthJRA = Tbreadth(AL.JRA.out.inf) # 22.5
TbreadthPAR = Tbreadth(AL.PAR.out.inf) # 21.5
TbreadthPOW = Tbreadth(AL.POW.out.inf) # 22.5

##### Compile prediction data for all populations #####

ALmeansd <- list(AL.HOP.out.inf$BUGSoutput$summary[1:2,1:2], AL.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                 AL.MAR2.out.inf$BUGSoutput$summary[1:2,1:2], AL.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                 AL.EUG.out.inf$BUGSoutput$summary[1:2,1:2], AL.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                 AL.SB.out.inf$BUGSoutput$summary[1:2,1:2], AL.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                 AL.POW.out.inf$BUGSoutput$summary[1:2,1:2], AL.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
ALtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptMAR2, ToptWAW, ToptEUG, ToptPLA, 
                                ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 10)), nrow = 10, ncol = 2))
ALtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthMAR2, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                     TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 10)), nrow = 10, ncol = 2))
colnames(ALtops) = c("mean", "sd")
colnames(ALtbreadths) = c("mean", "sd")
ALdf2 = do.call(rbind.data.frame, ALmeansd)
ALdf_1 = rbind(ALdf2, ALtops)
ALdf = rbind(ALdf_1, ALtbreadths)
ALdf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                    rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
ALdf$Param = c(rep(c("Tmin", "Tmax"), 10), rep("Topt", 10), rep("Tbreadth", 10))
save(ALdf, file = "AL_meansd.Rsave")

#### Plot TPC parameters for each population #####

load("AL_meansd.Rsave")
library(ggplot2)

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
ALdf$Population = factor(ALdf$Population, levels = PopOrder)
ALdf = ALdf[order(ALdf$Population),]
ALdf = ALdf[ALdf$Param != "Tbreadth",]

# colors for each population
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotAL = ggplot(ALdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, size = 2,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(), width = 0.3) + 
  ggtitle("Adult Lifespan") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

#### 6b. Plot Adult lifespan curve for all populaions #####

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

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
plot(AdultLifespan ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,15), data = data.AL.HOP, type = "n", bty = "n",
     ylab = "lifespan (days)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Adult Lifespan", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
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


