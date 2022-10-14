# Trait Fits for Mosquito Development Rate ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for pupal development rate (PDR) for 
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit mosquito dev rate thermal responses (Briere) using low info priors
#  4) Fit gamma distributions to MDR prior thermal responses
#  5) Fit mosquito dev rate thermal responses with data-informed priors
#  6) Plot parameters


#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.MDR = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# Keep only rows with PDR data
data.MDR = data.MDR[!is.na(data.MDR$JuvenileDevRate),]

# Subset data by population
# Subset data by population
data.MDR.HOP <- subset(data.MDR, Population == "HOP")
data.MDR.MAR1 <- subset(data.MDR, Population == "MAR35")
data.MDR.MAR2 <- subset(data.MDR, Population == "MAR29")
data.MDR.WAW <- subset(data.MDR, Population == "WAW")
data.MDR.EUG <- subset(data.MDR, Population == "EUG")
data.MDR.PLA <- subset(data.MDR, Population == "PLA")
data.MDR.SB <- subset(data.MDR, Population == "SB")
data.MDR.JRA <- subset(data.MDR, Population == "JRA")
data.MDR.PAR <- subset(data.MDR, Population == "PAR")
data.MDR.POW <- subset(data.MDR, Population == "POW")

# Plot trait data
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), pch = 16,
     data = data.MDR, ylab = "Juvenile dev rate", xlab = "Temperature")

##### 2a. Set up JAGS model #####
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



##### 2b. Shared settings for all models ####

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
# goal = use info from other populations to decrease uncertainty in parameter estimates

# Subset Data - leave-one-out
data.MDR.HOP.prior <- subset(data.MDR, Population != "HOP")
data.MDR.MAR1.prior <- subset(data.MDR, Population != "MAR35")
data.MDR.MAR2.prior <- subset(data.MDR, Population != "MAR29")
data.MDR.WAW.prior <- subset(data.MDR, Population != "WAW")
data.MDR.EUG.prior <- subset(data.MDR, Population != "EUG")
data.MDR.PLA.prior <- subset(data.MDR, Population != "PLA")
data.MDR.SB.prior <- subset(data.MDR, Population != "SB")
data.MDR.JRA.prior <- subset(data.MDR, Population != "JRA")
data.MDR.PAR.prior <- subset(data.MDR, Population != "PAR")
data.MDR.POW.prior <- subset(data.MDR, Population != "POW")

###### 3a: Hopland (HOP) ######

# set data
data <- data.MDR.HOP.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.HOP.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.HOP.prior, ylab = "MDR for Hop prior", xlab = "Temperature")
lines(MDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3b: Marin35 (MAR1) ######

# set data
data <- data.MDR.MAR1.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.MAR1.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.MAR1.prior, ylab = "MDR for MAR1 prior", xlab = "Temperature")
lines(MDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)





###### 3c: Marin29 (MAR2) ######

# set data
data <- data.MDR.MAR2.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.MAR2.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.MAR2.prior, ylab = "MDR for MAR2 prior", xlab = "Temperature")
lines(MDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3d: Wawona (WAW) ######

# set data
data <- data.MDR.WAW.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.WAW.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.WAW.prior, ylab = "MDR for WAW prior", xlab = "Temperature")
lines(MDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.MDR.EUG.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.EUG.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.EUG.prior, ylab = "MDR for EUG prior", xlab = "Temperature")
lines(MDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.MDR.PLA.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.PLA.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.PLA.prior, ylab = "MDR for PLA prior", xlab = "Temperature")
lines(MDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.MDR.SB.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.SB.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.SB.prior, ylab = "MDR for SB prior", xlab = "Temperature")
lines(MDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.MDR.JRA.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.JRA.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.JRA.prior, ylab = "MDR for JRA prior", xlab = "Temperature")
lines(MDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.MDR.PAR.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.PAR.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.PAR.prior, ylab = "MDR for PAR prior", xlab = "Temperature")
lines(MDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyono (POW) ######

# set data
data <- data.MDR.POW.prior

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
MDR.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.POW.prior.out)

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.POW.prior, ylab = "MDR for POW prior", xlab = "Temperature")
lines(MDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



#### 4. Fit gamma distributions to MDR prior thermal responses ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists

###### 4a. Hopland ######

MDR.HOP.prior.cf.dists <- data.frame(q = as.vector(MDR.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(MDR.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(MDR.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.HOP.prior.gamma.fits = apply(MDR.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)





###### 4b. Marin35 (MAR1) ######

MDR.MAR1.prior.cf.dists <- data.frame(q = as.vector(MDR.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(MDR.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(MDR.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.MAR1.prior.gamma.fits = apply(MDR.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4c. Marin29 (MAR2) ######

MDR.MAR2.prior.cf.dists <- data.frame(q = as.vector(MDR.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(MDR.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(MDR.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.MAR2.prior.gamma.fits = apply(MDR.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

MDR.WAW.prior.cf.dists <- data.frame(q = as.vector(MDR.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(MDR.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.WAW.prior.gamma.fits = apply(MDR.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

MDR.EUG.prior.cf.dists <- data.frame(q = as.vector(MDR.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(MDR.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.EUG.prior.gamma.fits = apply(MDR.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

MDR.PLA.prior.cf.dists <- data.frame(q = as.vector(MDR.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(MDR.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.PLA.prior.gamma.fits = apply(MDR.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

MDR.SB.prior.cf.dists <- data.frame(q = as.vector(MDR.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(MDR.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(MDR.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.SB.prior.gamma.fits = apply(MDR.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

MDR.JRA.prior.cf.dists <- data.frame(q = as.vector(MDR.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(MDR.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.JRA.prior.gamma.fits = apply(MDR.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

MDR.PAR.prior.cf.dists <- data.frame(q = as.vector(MDR.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(MDR.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.PAR.prior.gamma.fits = apply(MDR.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

MDR.POW.prior.cf.dists <- data.frame(q = as.vector(MDR.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(MDR.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(MDR.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
MDR.POW.prior.gamma.fits = apply(MDR.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


###### Hyperparameter list #####
MDR.hypers <- list(MDR.HOP.prior.gamma.fits, MDR.MAR1.prior.gamma.fits, MDR.MAR2.prior.gamma.fits,
                  MDR.WAW.prior.gamma.fits, MDR.EUG.prior.gamma.fits, MDR.PLA.prior.gamma.fits,
                  MDR.SB.prior.gamma.fits, MDR.JRA.prior.gamma.fits, MDR.PAR.prior.gamma.fits, MDR.POW.prior.gamma.fits)
save(MDR.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/MDRhypers.Rsave")


#### 5. Fit MDR thermal responses with data-informed priors ####

# to bypass above code of generating hyperparameters
load("MDRhypers.Rsave")
MDR.HOP.prior.gamma.fits <- MDR.hypers[[1]]
MDR.MAR1.prior.gamma.fits <- MDR.hypers[[2]]
MDR.MAR2.prior.gamma.fits <- MDR.hypers[[3]]
MDR.WAW.prior.gamma.fits <- MDR.hypers[[4]]
MDR.EUG.prior.gamma.fits <- MDR.hypers[[5]]
MDR.PLA.prior.gamma.fits <- MDR.hypers[[6]]
MDR.SB.prior.gamma.fits <- MDR.hypers[[7]]
MDR.JRA.prior.gamma.fits <- MDR.hypers[[8]]
MDR.PAR.prior.gamma.fits <- MDR.hypers[[9]]
MDR.POW.prior.gamma.fits <- MDR.hypers[[10]]

###### 5a. Hopland ######

# Set data 
data <- data.MDR.HOP
hypers <- MDR.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(MDR.HOP.out.inf)
#save(MDR.HOP.out.inf, file = "jagsout_MDR.HOPp_inf.Rdata")

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.HOP, ylab = "MDR for Hopland", xlab = "Temperature", pch = 1)
lines(MDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 5b. Marin35 ######

# Set data 
data <- data.MDR.MAR1
hypers <- MDR.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.MAR1, ylab = "MDR for Marin35", xlab = "Temperature", pch = 1)
lines(MDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5c. Marin29 ######

# Set data 
data <- data.MDR.MAR2
hypers <- MDR.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.MAR2, ylab = "MDR for Marin29", xlab = "Temperature", pch = 1)
lines(MDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5d. Wawona ######

# Set data 
data <- data.MDR.WAW
hypers <- MDR.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.WAW, ylab = "MDR for Wawona", xlab = "Temperature", pch = 1)
lines(MDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5e. Eugene ######

# Set data 
data <- data.MDR.EUG
hypers <- MDR.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.EUG, ylab = "MDR for Eugene", xlab = "Temperature", pch = 1)
lines(MDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5f. Placer ######

# Set data 
data <- data.MDR.PLA
hypers <- MDR.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.PLA, ylab = "MDR for Placer", xlab = "Temperature", pch = 1)
lines(MDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5g. Santa Barbara ######

# Set data 
data <- data.MDR.SB
hypers <- MDR.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.SB, ylab = "MDR for Santa Barbara", xlab = "Temperature", pch = 1)
lines(MDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5g. Jasper Ridge ######

# Set data 
data <- data.MDR.JRA
hypers <- MDR.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.JRA, ylab = "MDR for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(MDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5i. Paso Robles ######

# Set data 
data <- data.MDR.PAR
hypers <- MDR.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.PAR, ylab = "MDR for Paso Robles", xlab = "Temperature", pch = 1)
lines(MDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5j. Powder canyon ######

# Set data 
data <- data.MDR.POW
hypers <- MDR.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$JuvenileDevRate
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
MDR.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
MDR.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,0.1), data = data.MDR.POW, ylab = "MDR for Powder Canyon", xlab = "Temperature", pch = 1)
lines(MDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




###### Topt of MDR for each population ######
ToptHOP =Temp.xs[which.max(as.vector(MDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28
ToptMAR1 =Temp.xs[which.max(as.vector(MDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28
ToptMAR2 =Temp.xs[which.max(as.vector(MDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26.5
ToptWAW =Temp.xs[which.max(as.vector(MDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28
ToptEUG =Temp.xs[which.max(as.vector(MDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5
ToptPLA =Temp.xs[which.max(as.vector(MDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5
ToptSB =Temp.xs[which.max(as.vector(MDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27.5
ToptJRA =Temp.xs[which.max(as.vector(MDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 26
ToptPAR =Temp.xs[which.max(as.vector(MDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 27
ToptPOW =Temp.xs[which.max(as.vector(MDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 28

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
TbreadthHOP = Tbreadth(MDR.HOP.out.inf) # 18
TbreadthMAR1 = Tbreadth(MDR.MAR1.out.inf) # 17.5
TbreadthMAR2 = Tbreadth(MDR.MAR2.out.inf) # 17
TbreadthWAW = Tbreadth(MDR.WAW.out.inf) # 18
TbreadthEUG = Tbreadth(MDR.EUG.out.inf) # 18
TbreadthPLA = Tbreadth(MDR.PLA.out.inf) # 18
TbreadthSB = Tbreadth(MDR.SB.out.inf) # 18
TbreadthJRA = Tbreadth(MDR.JRA.out.inf) # 17
TbreadthPAR = Tbreadth(MDR.PAR.out.inf) # 17.5
TbreadthPOW = Tbreadth(MDR.POW.out.inf) # 17.5

#### Compile prediction data for all populations ###### 

MDRmeansd <- list(MDR.HOP.out.inf$BUGSoutput$summary[1:2,1:2], MDR.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                  MDR.MAR2.out.inf$BUGSoutput$summary[1:2,1:2], MDR.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                  MDR.EUG.out.inf$BUGSoutput$summary[1:2,1:2], MDR.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                  MDR.SB.out.inf$BUGSoutput$summary[1:2,1:2], MDR.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                  MDR.POW.out.inf$BUGSoutput$summary[1:2,1:2], MDR.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
MDRtops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptMAR2, ToptWAW, ToptEUG, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 10)), nrow = 10, ncol = 2))
MDRtbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthMAR2, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                      TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 10)), nrow = 10, ncol = 2))
colnames(MDRtops) = c("mean", "sd")
colnames(MDRtbreadths) = c("mean", "sd")
MDRdf2 = do.call(rbind.data.frame, MDRmeansd)
MDRdf_1 = rbind(MDRdf2, MDRtops)
MDRdf = rbind(MDRdf_1, MDRtbreadths)
MDRdf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                     rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
MDRdf$Param = c(rep(c("Tmin", "Tmax"), 10), rep("Topt", 10), rep("Tbreadth", 10))
save(MDRdf, file = "MDR_meansd.Rsave")





#### 6a. Plot TPC parameters for each population #####

load("MDR_meansd.Rsave")
library(ggplot2)

# order from coldest to hottest temp in wettest quarter
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
MDRdf$Population = factor(MDRdf$Population, levels = PopOrder)
MDRdf = MDRdf[order(MDRdf$Population),]
MDRdf = MDRdf[MDRdf$Param != "Tbreadth",]

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotMDR = ggplot(MDRdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Mosquito development rate") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

#### 6b. Plot Mosquito dev rates curve for all populations #####

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = MDR.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = MDR.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = MDR.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = MDR.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = MDR.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = MDR.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = MDR.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = MDR.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = MDR.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = MDR.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
plot(JuvenileDevRate ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,0.08), data = data.MDR.HOP, type = "n", bty = "n",
     ylab = "rate (1/day)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Mosquito development rates", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
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

