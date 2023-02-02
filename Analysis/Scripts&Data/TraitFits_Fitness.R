# Trait Fits for Fitness (for individuals with non-zero fitness) #

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for fitness 
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit fitness thermal responses (Briere) using low info priors
#  4) Fit gamma distributions to fitness prior thermal responses
#  5) Fit fitness thermal responses with data-informed priors
#  6) Plot parameters

#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.fit = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]


# Subset data by population
data.fit.HOP <- subset(data.fit, Population == "HOP")
data.fit.MAR1 <- subset(data.fit, Population == "MARIN35")
data.fit.MAR2 <- subset(data.fit, Population == "MARIN29")
data.fit.WAW <- subset(data.fit, Population == "WAW")
data.fit.EUG <- subset(data.fit, Population == "EUG")
data.fit.PLA <- subset(data.fit, Population == "PLA")
data.fit.SB <- subset(data.fit, Population == "SB")
data.fit.JRA <- subset(data.fit, Population == "JRA")
data.fit.PAR <- subset(data.fit, Population == "PAR")
data.fit.POW <- subset(data.fit, Population == "POW")

# Plot trait data
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), pch = 16,
     data = data.fit, ylab = "Fitness", xlab = "Temperature")


plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), pch = 16,
     data = data.fit, ylab = "Fitness", xlab = "Temperature")


##### 2a. Set up JAGS model #####
# 2a. Briere Model with uniform priors #

sink("briere.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 40)
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






##### 3. Fit using uniform priors ####
###### 3a: Hopland (HOP) ######

# set data
data <- data.fit.HOP

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.HOP.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.HOP.prior.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.HOP, ylab = "fit for Hopr", xlab = "Temperature")
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3b: Marin35 (MAR1) ######

# set data
data <- data.fit.MAR1 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.MAR1.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.MAR1.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.MAR1.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.MAR1 , ylab = "fit for MAR1 prior", xlab = "Temperature")
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3c: Marin29 (MAR2) ######

# set data
data <- data.fit.MAR2 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.MAR2.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.MAR2.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.MAR2.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.MAR2 , ylab = "fit for MAR2 prior", xlab = "Temperature")
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3d: Wawona (WAW) ######

# set data
data <- data.fit.WAW 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.WAW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.WAW.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.WAW.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.WAW , ylab = "fit for WAW prior", xlab = "Temperature")
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.fit.EUG 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.EUG.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.EUG.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.EUG.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.EUG , ylab = "fit for EUG prior", xlab = "Temperature")
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.fit.PLA 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.PLA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.PLA.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.PLA.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.PLA , ylab = "fit for PLA prior", xlab = "Temperature")
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.fit.SB 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.SB.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.SB.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.SB.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.SB , ylab = "fit for SB prior", xlab = "Temperature")
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.fit.JRA 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.JRA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.JRA.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.JRA.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.JRA , ylab = "fit for JRA prior", xlab = "Temperature")
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.fit.PAR 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.PAR.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.PAR.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.PAR.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.PAR , ylab = "fit for PAR prior", xlab = "Temperature")
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

# set data
data <- data.fit.POW 

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.POW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.POW.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.POW.out)

# Plot data + fit
plot(fitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.POW , ylab = "fit for POW prior", xlab = "Temperature")
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### For internal checking: 10 panel plot of fits using uniform priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.HOP, ylab = "fit for Hop", xlab = "Temperature")
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.MAR1, ylab = "fit for MAR1", xlab = "Temperature")
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.MAR2 , ylab = "fit for MAR2 prior", xlab = "Temperature")
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.WAW , ylab = "fit for WAW prior", xlab = "Temperature")
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.EUG , ylab = "fit for EUG prior", xlab = "Temperature")
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.PLA , ylab = "fit for PLA prior", xlab = "Temperature")
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.SB , ylab = "fit for SB prior", xlab = "Temperature")
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.JRA , ylab = "fit for JRA prior", xlab = "Temperature")
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.PAR , ylab = "fit for PAR prior", xlab = "Temperature")
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,50), data = data.fit.POW , ylab = "fit for POW prior", xlab = "Temperature")
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 10 panel plots of fits with uniform priors, but raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$SexWeightedFitness ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$SexWeightedFitness ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$SexWeightedFitness ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/fitness_UniformPriors.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

fit.HOP.MeansSds = MeanSd(data.fit.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.HOP.MeansSds, ylab = "fitness for Hop", xlab = "Temperature", pch= 16)
arrows(fit.HOP.MeansSds$Temp.Treatment, fit.HOP.MeansSds$mean-fit.HOP.MeansSds$stderror, fit.HOP.MeansSds$Temp.Treatment, fit.HOP.MeansSds$mean + fit.HOP.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.MAR1.MeansSds = MeanSd(data.fit.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.MAR1.MeansSds, ylab = "fitness for MAR1", xlab = "Temperature", pch= 16)
arrows(fit.MAR1.MeansSds$Temp.Treatment, fit.MAR1.MeansSds$mean-fit.MAR1.MeansSds$stderror, fit.MAR1.MeansSds$Temp.Treatment, fit.MAR1.MeansSds$mean + fit.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.MAR2.MeansSds = MeanSd(data.fit.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.MAR2.MeansSds, ylab = "fitness for MAR2", xlab = "Temperature", pch= 16)
arrows(fit.MAR2.MeansSds$Temp.Treatment, fit.MAR2.MeansSds$mean-fit.MAR2.MeansSds$stderror, fit.MAR2.MeansSds$Temp.Treatment, fit.MAR2.MeansSds$mean + fit.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.WAW.MeansSds = MeanSd(data.fit.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.WAW.MeansSds, ylab = "fitness for WAW", xlab = "Temperature", pch= 16)
arrows(fit.WAW.MeansSds$Temp.Treatment, fit.WAW.MeansSds$mean-fit.WAW.MeansSds$stderror, fit.WAW.MeansSds$Temp.Treatment, fit.WAW.MeansSds$mean + fit.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.EUG.MeansSds = MeanSd(data.fit.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.EUG.MeansSds, ylab = "fitness for EUG", xlab = "Temperature", pch= 16)
arrows(fit.EUG.MeansSds$Temp.Treatment, fit.EUG.MeansSds$mean-fit.EUG.MeansSds$stderror, fit.EUG.MeansSds$Temp.Treatment, fit.EUG.MeansSds$mean + fit.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.PLA.MeansSds = MeanSd(data.fit.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.PLA.MeansSds, ylab = "fitness for PLA", xlab = "Temperature", pch= 16)
arrows(fit.PLA.MeansSds$Temp.Treatment, fit.PLA.MeansSds$mean-fit.PLA.MeansSds$stderror, fit.PLA.MeansSds$Temp.Treatment, fit.PLA.MeansSds$mean + fit.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.SB.MeansSds = MeanSd(data.fit.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.SB.MeansSds, ylab = "fitness for SB", xlab = "Temperature", pch= 16)
arrows(fit.SB.MeansSds$Temp.Treatment, fit.SB.MeansSds$mean-fit.SB.MeansSds$stderror, fit.SB.MeansSds$Temp.Treatment, fit.SB.MeansSds$mean + fit.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.JRA.MeansSds = MeanSd(data.fit.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.JRA.MeansSds, ylab = "fitness for JRA", xlab = "Temperature", pch= 16)
arrows(fit.JRA.MeansSds$Temp.Treatment, fit.JRA.MeansSds$mean-fit.JRA.MeansSds$stderror, fit.JRA.MeansSds$Temp.Treatment, fit.JRA.MeansSds$mean + fit.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.PAR.MeansSds = MeanSd(data.fit.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.PAR.MeansSds, ylab = "fitness for PAR", xlab = "Temperature", pch= 16)
arrows(fit.PAR.MeansSds$Temp.Treatment, fit.PAR.MeansSds$mean-fit.PAR.MeansSds$stderror, fit.PAR.MeansSds$Temp.Treatment, fit.PAR.MeansSds$mean + fit.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

fit.POW.MeansSds = MeanSd(data.fit.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,12), data = fit.POW.MeansSds, ylab = "fitness for POW", xlab = "Temperature", pch= 16)
arrows(fit.POW.MeansSds$Temp.Treatment, fit.POW.MeansSds$mean-fit.POW.MeansSds$stderror, fit.POW.MeansSds$Temp.Treatment, fit.POW.MeansSds$mean + fit.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

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

fitdf_uni <- rbind.data.frame(ParamCompile(fit.HOP.out), ParamCompile(fit.MAR1.out),
                              ParamCompile(fit.MAR2.out), ParamCompile(fit.WAW.out),
                              ParamCompile(fit.EUG.out), ParamCompile(fit.PLA.out),
                              ParamCompile(fit.SB.out), ParamCompile(fit.JRA.out),
                              ParamCompile(fit.PAR.out), ParamCompile(fit.POW.out))
fitdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(fitdf_uni, file = "fit_meansd_uni.Rsave")



###### 5. Plotting  #####

library(ggplot2)
load("fit_meansd_uni.Rsave")

###### 5a. Plot of CTmax, CTmin, Topt for each population #####
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

fitdf_uni$Population = factor(fitdf_uni$Population, levels = LatOrder)
fitdf_uni = fitdf_uni[order(fitdf_uni$Population),]
fitdf_uni = fitdf_uni[fitdf_uni$param != "Tbreadth",]
fitdf_uni = fitdf_uni[fitdf_uni$param != "Pmax",]

ggplot(fitdf_uni, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 35)) +
  geom_point(stat="identity", col=LatColors, size = 2,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors, 
                position=position_dodge(.9), width = 0.3, size = 0.7) +
  ggtitle("Thermal limits of fitness (uniform priors)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (\u00B0C)") +
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA))

###### 5b. Plot of Pmax for each population #####

load("fit_meansd_uni.Rsave")

LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors2 = c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")

fitdf_uni$Population = factor(fitdf_uni$Population, levels = LatOrder)
fitdf_uni = fitdf_uni[order(fitdf_uni$Population),]
fitdf_uni = fitdf_uni[fitdf_uni$param == "Pmax",]

ggplot(fitdf_uni, aes(x=Population, y=mean, col = LatColors2)) + 
  scale_y_continuous(limits = c(0, 13)) +
  geom_point(stat="identity", col=LatColors2, size = 2,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors2,
                position=position_dodge(.9), width = 0.3, size = 0.7) + 
  ggtitle("Max fitness (uniform priors)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Max fitness") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 

###### 5c. Plot fitness curves for all populations #####
load(file = "fit_meansd_uni.Rsave")

LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors2 = rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"))

EUGdata = fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit 
plot(SexWeightedFitness ~ Temp.Treatment, 
     xlim = c(5, 35), ylim = c(0,12), data = data.fit.HOP, type = "n", bty = "n",
     ylab = "fitness (offspring per individual)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Fitness (uniform priors)", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
legend("topright", legend = LatOrder, bty = "n", cex  = 1.0,
       col = LatColors2, lty = 1, lwd = 2)
lines(POWdata ~ Temp.xs, col = LatColors2[10], lwd = 1.5)
lines(SBdata ~ Temp.xs, col = LatColors2[9], lwd = 1.5)
lines(PARdata ~ Temp.xs, col = LatColors2[8], lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = LatColors2[7], lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = LatColors2[6], lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = LatColors2[5], lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = LatColors2[4], lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = LatColors2[3], lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = LatColors2[2], lwd = 1.5)
lines(EUGdata ~ Temp.xs, col = LatColors2[1], lwd = 1.5)
box(col = "black")



###### 5d. Plot Tbreadth for all populations #####

load("fit_meansd_uni.Rsave")

LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors2 = c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")

fitdf_uni$Population = factor(fitdf_uni$Population, levels = LatOrder)
fitdf_uni = fitdf_uni[order(fitdf_uni$Population),]
fitdf_uni = fitdf_uni[fitdf_uni$param == "Tbreadth",]

ggplot(fitdf_uni, aes(x=Population, y=mean, col = LatColors2)) + 
  scale_y_continuous(limits = c(6, 16)) +
  geom_point(stat="identity", col=LatColors2, size = 2,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors2,
                position=position_dodge(.9), width = 0.3, size = 0.7) + 
  ggtitle("Thermal breadth (uniform priors)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 


###### 6. Generate informative priors: Leave-one-out with uniform priors ######
# Use same uniform priors settings as above
# Subset Data for leave-one-out
data.fit.HOP.prior <- subset(data.fit, Population != "HOP")
data.fit.MAR1.prior <- subset(data.fit, Population != "MARIN35")
data.fit.MAR2.prior <- subset(data.fit, Population != "MARIN29")
data.fit.WAW.prior <- subset(data.fit, Population != "WAW")
data.fit.EUG.prior <- subset(data.fit, Population != "EUG")
data.fit.PLA.prior <- subset(data.fit, Population != "PLA")
data.fit.SB.prior <- subset(data.fit, Population != "SB")
data.fit.JRA.prior <- subset(data.fit, Population != "JRA")
data.fit.PAR.prior <- subset(data.fit, Population != "PAR")
data.fit.POW.prior <- subset(data.fit, Population != "POW")

###### 6a. Hopland #####
# set data
data <- data.fit.HOP.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.HOP.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.HOP.prior, ylab = "fit for Hop prior", xlab = "Temperature")
lines(fit.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### 6b: Marin35 (MAR1) #####

# set data
data <- data.fit.MAR1.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.MAR1.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.MAR1.prior, ylab = "fit for MAR1 prior", xlab = "Temperature")
lines(fit.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.fit.MAR2.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.MAR2.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.MAR2.prior, ylab = "fit for MAR2 prior", xlab = "Temperature")
lines(fit.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.fit.WAW.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.WAW.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.WAW.prior, ylab = "fit for WAW prior", xlab = "Temperature")
lines(fit.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.fit.EUG.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.EUG.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.EUG.prior, ylab = "fit for EUG prior", xlab = "Temperature")
lines(fit.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.fit.PLA.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.PLA.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.PLA.prior, ylab = "fit for PLA prior", xlab = "Temperature")
lines(fit.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.fit.SB.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.SB.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.SB.prior, ylab = "fit for SB prior", xlab = "Temperature")
lines(fit.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.fit.JRA.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.JRA.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.JRA.prior, ylab = "fit for JRA prior", xlab = "Temperature")
lines(fit.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.fit.PAR.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.PAR.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.PAR.prior, ylab = "fit for PAR prior", xlab = "Temperature")
lines(fit.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyon (POW) ######

# set data
data <- data.fit.POW.prior

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fit.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fit.POW.prior.out)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.POW.prior, ylab = "fit for POW prior", xlab = "Temperature")
lines(fit.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### Save posterior distributions for informative priors #####

fit.HOP.prior.cf.dists <- data.frame(q = as.vector(fit.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.HOP.prior.gamma.fits = apply(fit.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.MAR1.prior.cf.dists <- data.frame(q = as.vector(fit.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(fit.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(fit.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.MAR1.prior.gamma.fits = apply(fit.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.MAR2.prior.cf.dists <- data.frame(q = as.vector(fit.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(fit.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(fit.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.MAR2.prior.gamma.fits = apply(fit.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.WAW.prior.cf.dists <- data.frame(q = as.vector(fit.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.WAW.prior.gamma.fits = apply(fit.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.EUG.prior.cf.dists <- data.frame(q = as.vector(fit.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.EUG.prior.gamma.fits = apply(fit.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.PLA.prior.cf.dists <- data.frame(q = as.vector(fit.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.PLA.prior.gamma.fits = apply(fit.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.SB.prior.cf.dists <- data.frame(q = as.vector(fit.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(fit.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(fit.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.SB.prior.gamma.fits = apply(fit.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.JRA.prior.cf.dists <- data.frame(q = as.vector(fit.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.JRA.prior.gamma.fits = apply(fit.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.PAR.prior.cf.dists <- data.frame(q = as.vector(fit.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.PAR.prior.gamma.fits = apply(fit.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.POW.prior.cf.dists <- data.frame(q = as.vector(fit.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fit.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fit.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
fit.POW.prior.gamma.fits = apply(fit.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fit.hypers <- list(fit.HOP.prior.gamma.fits, fit.MAR1.prior.gamma.fits, fit.MAR2.prior.gamma.fits,
                   fit.WAW.prior.gamma.fits, fit.EUG.prior.gamma.fits, fit.PLA.prior.gamma.fits,
                   fit.SB.prior.gamma.fits, fit.JRA.prior.gamma.fits, fit.PAR.prior.gamma.fits, fit.POW.prior.gamma.fits)
save(fit.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/fithypers.Rsave")


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
load("fithypers.Rsave")
fit.HOP.prior.gamma.fits <- fit.hypers[[1]]
fit.MAR1.prior.gamma.fits <- fit.hypers[[2]]
fit.MAR2.prior.gamma.fits <- fit.hypers[[3]]
fit.WAW.prior.gamma.fits <- fit.hypers[[4]]
fit.EUG.prior.gamma.fits <- fit.hypers[[5]]
fit.PLA.prior.gamma.fits <- fit.hypers[[6]]
fit.SB.prior.gamma.fits <- fit.hypers[[7]]
fit.JRA.prior.gamma.fits <- fit.hypers[[8]]
fit.PAR.prior.gamma.fits <- fit.hypers[[9]]
fit.POW.prior.gamma.fits <- fit.hypers[[10]]




###### 7a. Hopland ######

# Set data 
data <- data.fit.HOP
hypers <- fit.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(fit.HOP.out.inf)

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.HOP, ylab = "fit for Hopland", xlab = "Temperature", pch = 1)
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 7b. Marin35 ######

# Set data 
data <- data.fit.MAR1
hypers <- fit.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.MAR1, ylab = "fit for Marin35", xlab = "Temperature", pch = 1)
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.fit.MAR2
hypers <- fit.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.MAR2, ylab = "fit for Marin29", xlab = "Temperature", pch = 1)
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.fit.WAW
hypers <- fit.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.WAW, ylab = "fit for Wawona", xlab = "Temperature", pch = 1)
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.fit.EUG
hypers <- fit.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.EUG, ylab = "fit for Eugene", xlab = "Temperature", pch = 1)
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.fit.PLA
hypers <- fit.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.PLA, ylab = "fit for Placer", xlab = "Temperature", pch = 1)
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.fit.SB
hypers <- fit.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.SB, ylab = "fit for Santa Barbara", xlab = "Temperature", pch = 1)
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.fit.JRA
hypers <- fit.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.JRA, ylab = "fit for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.fit.PAR
hypers <- fit.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.PAR, ylab = "fit for Paso Robles", xlab = "Temperature", pch = 1)
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon (POW) ######

# Set data 
data <- data.fit.POW
hypers <- fit.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$SexWeightedFitness
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fit.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fit.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,50), data = data.fit.POW, ylab = "fit for Powder Canyon", xlab = "Temperature", pch = 1)
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### For internal checking: 10 panel plot of fits using informative priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.HOP, ylab = "fit for Hop", xlab = "Temperature")
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.MAR1, ylab = "fit for MAR1", xlab = "Temperature")
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.MAR2 , ylab = "fit for MAR2 prior", xlab = "Temperature")
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.WAW , ylab = "fit for WAW prior", xlab = "Temperature")
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.EUG , ylab = "fit for EUG prior", xlab = "Temperature")
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.PLA , ylab = "fit for PLA prior", xlab = "Temperature")
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.SB , ylab = "fit for SB prior", xlab = "Temperature")
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.JRA , ylab = "fit for JRA prior", xlab = "Temperature")
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.PAR , ylab = "fit for PAR prior", xlab = "Temperature")
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(SexWeightedFitness ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,50), data = data.fit.POW , ylab = "fit for POW prior", xlab = "Temperature")
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 10 panel plots of fits with raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$SexWeightedFitness ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$SexWeightedFitness ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$SexWeightedFitness ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/fitness_HigherWeight.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

fit.Hop.MeansSds = MeanSd(data.fit.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.Hop.MeansSds, ylab = "SexWeightedFitness for Hop", xlab = "Temperature", pch= 16)
arrows(fit.Hop.MeansSds$Temp.Treatment, fit.Hop.MeansSds$mean-fit.Hop.MeansSds$stderror, fit.Hop.MeansSds$Temp.Treatment, fit.Hop.MeansSds$mean + fit.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 12, cex = 1.5)
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.MAR1.MeansSds = MeanSd(data.fit.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.MAR1.MeansSds, ylab = "SexWeightedFitness for MAR1", xlab = "Temperature", pch= 16)
arrows(fit.MAR1.MeansSds$Temp.Treatment, fit.MAR1.MeansSds$mean-fit.MAR1.MeansSds$stderror, fit.MAR1.MeansSds$Temp.Treatment, fit.MAR1.MeansSds$mean + fit.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 12, cex = 1.5)
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.MAR2.MeansSds = MeanSd(data.fit.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.MAR2.MeansSds, ylab = "SexWeightedFitness for MAR2", xlab = "Temperature", pch= 16)
arrows(fit.MAR2.MeansSds$Temp.Treatment, fit.MAR2.MeansSds$mean-fit.MAR2.MeansSds$stderror, fit.MAR2.MeansSds$Temp.Treatment, fit.MAR2.MeansSds$mean + fit.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 36.5, y = 12, cex = 1.5)
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.WAW.MeansSds = MeanSd(data.fit.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.WAW.MeansSds, ylab = "SexWeightedFitness for WAW", xlab = "Temperature", pch= 16)
arrows(fit.WAW.MeansSds$Temp.Treatment, fit.WAW.MeansSds$mean-fit.WAW.MeansSds$stderror, fit.WAW.MeansSds$Temp.Treatment, fit.WAW.MeansSds$mean + fit.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 36.5, y = 12, cex = 1.5)
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.EUG.MeansSds = MeanSd(data.fit.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.EUG.MeansSds, ylab = "SexWeightedFitness for EUG", xlab = "Temperature", pch= 16)
arrows(fit.EUG.MeansSds$Temp.Treatment, fit.EUG.MeansSds$mean-fit.EUG.MeansSds$stderror, fit.EUG.MeansSds$Temp.Treatment, fit.EUG.MeansSds$mean + fit.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG", x = 36.5, y = 12, cex = 1.5)
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.PLA.MeansSds = MeanSd(data.fit.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.PLA.MeansSds, ylab = "SexWeightedFitness for PLA", xlab = "Temperature", pch= 16)
arrows(fit.PLA.MeansSds$Temp.Treatment, fit.PLA.MeansSds$mean-fit.PLA.MeansSds$stderror, fit.PLA.MeansSds$Temp.Treatment, fit.PLA.MeansSds$mean + fit.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 36.5, y = 12, cex = 1.5)
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.SB.MeansSds = MeanSd(data.fit.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.SB.MeansSds, ylab = "SexWeightedFitness for SB", xlab = "Temperature", pch= 16)
arrows(fit.SB.MeansSds$Temp.Treatment, fit.SB.MeansSds$mean-fit.SB.MeansSds$stderror, fit.SB.MeansSds$Temp.Treatment, fit.SB.MeansSds$mean + fit.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 36.5, y = 12, cex = 1.5)
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.JRA.MeansSds = MeanSd(data.fit.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.JRA.MeansSds, ylab = "SexWeightedFitness for JRA", xlab = "Temperature", pch= 16)
arrows(fit.JRA.MeansSds$Temp.Treatment, fit.JRA.MeansSds$mean-fit.JRA.MeansSds$stderror, fit.JRA.MeansSds$Temp.Treatment, fit.JRA.MeansSds$mean + fit.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 36.5, y = 12, cex = 1.5)
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.PAR.MeansSds = MeanSd(data.fit.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.PAR.MeansSds, ylab = "SexWeightedFitness for PAR", xlab = "Temperature", pch= 16)
arrows(fit.PAR.MeansSds$Temp.Treatment, fit.PAR.MeansSds$mean-fit.PAR.MeansSds$stderror, fit.PAR.MeansSds$Temp.Treatment, fit.PAR.MeansSds$mean + fit.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 36.5, y = 12, cex = 1.5)
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fit.POW.MeansSds = MeanSd(data.fit.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,13), data = fit.POW.MeansSds, ylab = "SexWeightedFitness for POW", xlab = "Temperature", pch= 16)
arrows(fit.POW.MeansSds$Temp.Treatment, fit.POW.MeansSds$mean-fit.POW.MeansSds$stderror, fit.POW.MeansSds$Temp.Treatment, fit.POW.MeansSds$mean + fit.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 36.5, y = 12, cex = 1.5)
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(fit.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

dev.off()





###### 8. Get TPC parameters #####
###### 8a. Topt of fit for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of fit for each population ######

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
  
fitdf_inf <- rbind.data.frame(ParamCompile(fit.HOP.out.inf), ParamCompile(fit.MAR1.out.inf),
                                  ParamCompile(fit.MAR2.out.inf), ParamCompile(fit.WAW.out.inf),
                                  ParamCompile(fit.EUG.out.inf), ParamCompile(fit.PLA.out.inf),
                                  ParamCompile(fit.SB.out.inf), ParamCompile(fit.JRA.out.inf),
                                  ParamCompile(fit.PAR.out.inf), ParamCompile(fit.POW.out.inf))
fitdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(fitdf_inf, file = "fit_meansd_inf.Rsave")


###### 9. Plotting  #####
library(ggplot2)
load("fit_meansd_inf.Rsave")

plot.new()
par(mfrow = c(1,1))
par(mar = c(5, 5, 3, 2))


###### 9a. Plot of CTmax, CTmin, Topt for each population #####
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

fitdf_inf$Population = factor(fitdf_inf$Population, levels = LatOrder)
fitdf_inf = fitdf_inf[order(fitdf_inf$Population),]
fitdf_inf = fitdf_inf[fitdf_inf$param != "Tbreadth",]
fitdf_inf = fitdf_inf[fitdf_inf$param != "Pmax",]

ggplot(fitdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 32)) +
  geom_point(stat="identity", col=LatColors, size = 2,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors, 
                position=position_dodge(.9), width = 0.3, size = 0.7) + 
  ggtitle("Thermal limits of fitness") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (\u00B0C)") +
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA))

###### 9b. Plot Pmax for each population #####
load("fit_meansd_inf.Rsave")

LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors2 = c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")

fitdf_inf$Population = factor(fitdf_inf$Population, levels = LatOrder)
fitdf_inf = fitdf_inf[order(fitdf_inf$Population),]
fitdf_inf = fitdf_inf[fitdf_inf$param == "Pmax",]

ggplot(fitdf_inf, aes(x=Population, y=mean, col = LatColors2)) + 
  scale_y_continuous(limits = c(0, 13)) +
  geom_point(stat="identity", col=LatColors2, size = 2,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors2,
                position=position_dodge(.9), width = 0.3, size = 0.7) + 
  ggtitle("fitness") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Max fitness") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 


###### 9c. Plot fitness curves for all populations #####

load(file = "fit_meansd_inf.Rsave")

LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors2 = rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"))

EUGdata = fit.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = fit.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = fit.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = fit.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = fit.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = fit.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = fit.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = fit.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = fit.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = fit.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit 
par(mar = (c(4.5, 4.5, 2.5, 1)))
plot(SexWeightedFitness ~ Temp.Treatment, 
     xlim = c(5, 35), ylim = c(0,12), data = data.fit.HOP, type = "n", bty = "n",
     ylab = "fitness (offspring per individual)", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Fitness", cex.main = 1.7, cex.lab = 1.4, cex.axis = 1.2)
box(col = "black")
lines(POWdata ~ Temp.xs, col = LatColors2[10], lwd = 1.5)
lines(SBdata ~ Temp.xs, col = LatColors2[9], lwd = 1.5)
lines(PARdata ~ Temp.xs, col = LatColors2[8], lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = LatColors2[7], lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = LatColors2[6], lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = LatColors2[5], lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = LatColors2[4], lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = LatColors2[3], lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = LatColors2[2], lwd = 1.5)
lines(EUGdata ~ Temp.xs, col = LatColors2[1], lwd = 1.5)
box(col = "black")
legend("topright", LatOrder, col = LatColors2, lty = 1, cex = 0.8, lwd = 2, bty = "n")


##### 9d. Plot Tbreadth for all populations #####
load("fit_meansd_inf.Rsave")

LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors2 = c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")

fitdf_inf$Population = factor(fitdf_inf$Population, levels = LatOrder)
fitdf_inf = fitdf_inf[order(fitdf_inf$Population),]
fitdf_inf = fitdf_inf[fitdf_inf$param == "Tbreadth",]

ggplot(fitdf_inf, aes(x=Population, y=mean, col = LatColors2)) + 
  scale_y_continuous(limits = c(6, 16)) +
  geom_point(stat="identity", col=LatColors2, size = 2,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors2,
                position=position_dodge(.9), width = 0.3, size = 0.7) + 
  ggtitle("Thermal breadth") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA)) 










