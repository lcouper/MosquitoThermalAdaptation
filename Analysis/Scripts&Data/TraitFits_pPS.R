# Trait Fits for Pupal Survival ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for Pupal survival pPS 
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit pupal survival thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to pupal survival prior thermal responses
#  5) Fit oupal survival thermal responses with data-informed priors
#  6) Plot parameters

#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.pPS_int = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# Calculate survival for each population and temp treatment
data.pPS_surv = aggregate(Pupal.Survival ~ Population + Temp.Treatment, data = data.pPS_int, FUN = sum)
data.pPS_all = aggregate(Pupal.Survival ~ Population + Temp.Treatment, data = data.pPS_int, FUN = length)
data.pPS_per = data.pPS_surv$Pupal.Survival/data.pPS_all$Pupal.Survival 

data.pPS = cbind.data.frame(data.pPS_surv[,1:2], data.pPS_per)
colnames(data.pPS)[3] = "prob.Pupal.Survival"

# Plot
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), data = data.pPS, ylab = "prob. pupal survival", xlab = "Temperature")

# Subset data by population
data.pPS.HOP <- subset(data.pPS, Population == "HOP")
data.pPS.MAR1 <- subset(data.pPS, Population == "MAR35")
data.pPS.MAR2 <- subset(data.pPS, Population == "MAR29")
data.pPS.WAW <- subset(data.pPS, Population == "WAW")
data.pPS.EUG <- subset(data.pPS, Population == "EUG")
data.pPS.PLA <- subset(data.pPS, Population == "PLA")
data.pPS.SB <- subset(data.pPS, Population == "SB")
data.pPS.JRA <- subset(data.pPS, Population == "JRA")
data.pPS.PAR <- subset(data.pPS, Population == "PAR")
data.pPS.POW <- subset(data.pPS, Population == "POW")

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
data.pPS.HOP.prior <- subset(data.pPS, Population != "HOP")
data.pPS.MAR1.prior <- subset(data.pPS, Population != "MAR35")
data.pPS.MAR2.prior <- subset(data.pPS, Population != "MAR29")
data.pPS.WAW.prior <- subset(data.pPS, Population != "WAW")
data.pPS.EUG.prior <- subset(data.pPS, Population != "EUG")
data.pPS.PLA.prior <- subset(data.pPS, Population != "PLA")
data.pPS.SB.prior <- subset(data.pPS, Population != "SB")
data.pPS.JRA.prior <- subset(data.pPS, Population != "JRA")
data.pPS.PAR.prior <- subset(data.pPS, Population != "PAR")
data.pPS.POW.prior <- subset(data.pPS, Population != "POW")

###### 3a. Hopland (HOP) ####

# Set data
data <- data.pPS.HOP.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.HOP.prior, ylab = "pPS for Hopland", xlab = "Temperature")
lines(pPS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b. Marin35 (MAR1) ####

# Set data
data <- data.pPS.MAR1.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.MAR1.prior, ylab = "pPS for Marin35", xlab = "Temperature")
lines(pPS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




###### 3c. Marin29 (MAR2) ####

# Set data
data <- data.pPS.MAR2.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR2.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.HOP.prior, ylab = "pPS for Marin29", xlab = "Temperature")
lines(pPS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3d: Wawona (WAW) ######

# set data
data <- data.pPS.WAW.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.WAW.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.WAW.prior, ylab = "pPS for WAW prior", xlab = "Temperature")
lines(pPS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.pPS.EUG.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.EUG.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.EUG.prior, ylab = "pPS for EUG prior", xlab = "Temperature")
lines(pPS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.pPS.PLA.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.PLA.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.PLA.prior, ylab = "pPS for PLA prior", xlab = "Temperature")
lines(pPS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.pPS.SB.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.SB.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.SB.prior, ylab = "pPS for SB prior", xlab = "Temperature")
lines(pPS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.pPS.JRA.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.JRA.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.JRA.prior, ylab = "pPS for JRA prior", xlab = "Temperature")
lines(pPS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.pPS.PAR.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.PAR.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.PAR.prior, ylab = "pPS for PAR prior", xlab = "Temperature")
lines(pPS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyon (POW) ######

# set data
data <- data.pPS.POW.prior

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pPS.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pPS.POW.prior.out)

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.POW.prior, ylab = "pPS for POW prior", xlab = "Temperature")
lines(pPS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)














#### 4. Fit gamma distributions to pPS prior thermal responses ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists

###### 4a. Hopland ######
pPS.HOP.prior.cf.dists <- data.frame(q = as.vector(pPS.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))

pPS.HOP.prior.gamma.fits = apply(pPS.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


###### 4b. Marin35 (MAR1) ######

pPS.MAR1.prior.cf.dists <- data.frame(q = as.vector(pPS.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pPS.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pPS.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.MAR1.prior.gamma.fits = apply(pPS.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4c. Marin29 (MAR2) ######

pPS.MAR2.prior.cf.dists <- data.frame(q = as.vector(pPS.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pPS.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pPS.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.MAR2.prior.gamma.fits = apply(pPS.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

pPS.WAW.prior.cf.dists <- data.frame(q = as.vector(pPS.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.WAW.prior.gamma.fits = apply(pPS.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

pPS.EUG.prior.cf.dists <- data.frame(q = as.vector(pPS.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.EUG.prior.gamma.fits = apply(pPS.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

pPS.PLA.prior.cf.dists <- data.frame(q = as.vector(pPS.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.PLA.prior.gamma.fits = apply(pPS.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

pPS.SB.prior.cf.dists <- data.frame(q = as.vector(pPS.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(pPS.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(pPS.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.SB.prior.gamma.fits = apply(pPS.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

pPS.JRA.prior.cf.dists <- data.frame(q = as.vector(pPS.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.JRA.prior.gamma.fits = apply(pPS.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

pPS.PAR.prior.cf.dists <- data.frame(q = as.vector(pPS.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.PAR.prior.gamma.fits = apply(pPS.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

pPS.POW.prior.cf.dists <- data.frame(q = as.vector(pPS.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pPS.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pPS.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
pPS.POW.prior.gamma.fits = apply(pPS.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
pPS.hypers <- list(pPS.HOP.prior.gamma.fits, pPS.MAR1.prior.gamma.fits, pPS.MAR2.prior.gamma.fits,
                   pPS.WAW.prior.gamma.fits, pPS.EUG.prior.gamma.fits, pPS.PLA.prior.gamma.fits,
                   pPS.SB.prior.gamma.fits, pPS.JRA.prior.gamma.fits, pPS.PAR.prior.gamma.fits, pPS.POW.prior.gamma.fits)
save(pPS.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/pPShypers.Rsave")


#### 6. Fit pPS thermal responses with data-informed priors

# to bypass above code of generating hyperparameters
#load("pPShypers.Rsave")
pPS.HOP.prior.gamma.fits <- pPS.hypers[[1]]
pPS.MAR1.prior.gamma.fits <- pPS.hypers[[2]]
pPS.MAR2.prior.gamma.fits <- pPS.hypers[[3]]
pPS.WAW.prior.gamma.fits <- pPS.hypers[[4]]
pPS.EUG.prior.gamma.fits <- pPS.hypers[[5]]
pPS.PLA.prior.gamma.fits <- pPS.hypers[[6]]
pPS.SB.prior.gamma.fits <- pPS.hypers[[7]]
pPS.JRA.prior.gamma.fits <- pPS.hypers[[8]]
pPS.PAR.prior.gamma.fits <- pPS.hypers[[9]]
pPS.POW.prior.gamma.fits <- pPS.hypers[[10]]

#### 5. Fit pPS thermal responses with data-informed priors ####
###### 5a. Hopland ######

# Set data 
data <- data.pPS.HOP
hypers <- pPS.HOP.prior.gamma.fits * .1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.HOP.out.inf$BUGSoutput$summary[1:5,]

#save(pPS.HOP.out.inf, file = "jagsout_pPS_HOP_inf.Rdata")

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.HOP, ylab = "prob Pupal survival for Hopland", xlab = "Temperature", pch = 1)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5b. Marin35 ######

# Set data 
data <- data.pPS.MAR1
hypers <- pPS.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.MAR1, ylab = "pPS for Marin35", xlab = "Temperature", pch = 1)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5c. Marin29 ######

# Set data 
data <- data.pPS.MAR2
hypers <- pPS.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.MAR2, ylab = "pPS for Marin29", xlab = "Temperature", pch = 1)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5d. Wawona ######

# Set data 
data <- data.pPS.WAW
hypers <- pPS.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.WAW, ylab = "pPS for Wawona", xlab = "Temperature", pch = 1)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5e. Eugene ######

# Set data 
data <- data.pPS.EUG
hypers <- pPS.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.EUG, ylab = "pPS for Eugene", xlab = "Temperature", pch = 1)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5f. Placer ######

# Set data 
data <- data.pPS.PLA
hypers <- pPS.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.PLA, ylab = "pPS for Placer", xlab = "Temperature", pch = 1)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5g. Santa Barbara ######

# Set data 
data <- data.pPS.SB
hypers <- pPS.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.SB, ylab = "pPS for Santa Barbara", xlab = "Temperature", pch = 1)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5g. Jasper Ridge ######

# Set data 
data <- data.pPS.JRA
hypers <- pPS.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.JRA, ylab = "pPS for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5i. Paso Robles ######

# Set data 
data <- data.pPS.PAR
hypers <- pPS.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.PAR, ylab = "pPS for Paso Robles", xlab = "Temperature", pch = 1)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5j. Powder canyon ######

# Set data 
data <- data.pPS.POW
hypers <- pPS.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Pupal.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pPS.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pPS.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Pupal.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.4), data = data.pPS.POW, ylab = "pPS for Powder Canyon", xlab = "Temperature", pch = 1)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




###### Topt of pPS for each population ######
ToptHOP = Temp.xs[which.max(as.vector(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptMAR1 = Temp.xs[which.max(as.vector(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptMAR2 = Temp.xs[which.max(as.vector(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptWAW = Temp.xs[which.max(as.vector(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptEUG = Temp.xs[which.max(as.vector(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptPLA = Temp.xs[which.max(as.vector(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19.5
ToptSB = Temp.xs[which.max(as.vector(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptJRA = Temp.xs[which.max(as.vector(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19
ToptPAR = Temp.xs[which.max(as.vector(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 20
ToptPOW = Temp.xs[which.max(as.vector(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19.5

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
TbreadthHOP = Tbreadth(pPS.HOP.out.inf) # 24
TbreadthMAR1 = Tbreadth(pPS.MAR1.out.inf) # 24.5
TbreadthMAR2 = Tbreadth(pPS.MAR2.out.inf) # 22.5
TbreadthWAW = Tbreadth(pPS.WAW.out.inf) # 23.5
TbreadthEUG = Tbreadth(pPS.EUG.out.inf) # 23.5
TbreadthPLA = Tbreadth(pPS.PLA.out.inf) # 25
TbreadthSB = Tbreadth(pPS.SB.out.inf) # 24
TbreadthJRA = Tbreadth(pPS.JRA.out.inf) # 25
TbreadthPAR = Tbreadth(pPS.PAR.out.inf) # 24
TbreadthPOW = Tbreadth(pPS.POW.out.inf) # 23.5

#### Compile prediction data for all populations ###### 

pPSmeansd <- list(pPS.HOP.out.inf$BUGSoutput$summary[1:2,1:2], pPS.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                  pPS.MAR2.out.inf$BUGSoutput$summary[1:2,1:2], pPS.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                  pPS.EUG.out.inf$BUGSoutput$summary[1:2,1:2], pPS.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                  pPS.SB.out.inf$BUGSoutput$summary[1:2,1:2], pPS.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                  pPS.POW.out.inf$BUGSoutput$summary[1:2,1:2], pPS.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
pPStops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptMAR2, ToptWAW, ToptEUG, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 10)), nrow = 10, ncol = 2))
pPStbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthMAR2, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                      TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 10)), nrow = 10, ncol = 2))
colnames(pPStops) = c("mean", "sd")
colnames(pPStbreadths) = c("mean", "sd")
pPSdf2 = do.call(rbind.data.frame, pPSmeansd)
pPSdf_1 = rbind(pPSdf2, pPStops)
pPSdf = rbind(pPSdf_1, pPStbreadths)
pPSdf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                     rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
pPSdf$Param = c(rep(c("Tmin", "Tmax"), 10), rep("Topt", 10), rep("Tbreadth", 10))
save(pPSdf, file = "pPS_meansd.Rsave")



#### Plot TPC parameters for each population #####

load("pPS_meansd.Rsave")
library(ggplot2)

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pPSdf$Population = factor(pPSdf$Population, levels = PopOrder)
pPSdf = pPSdf[order(pPSdf$Population),]
pPSdf = pPSdf[pPSdf$Param != "Tbreadth",]

# colors for each population
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

plotpPS = ggplot(pPSdf, aes(x=Population, y=mean, col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, size = 2,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(), width = 0.3) + 
  ggtitle("Pupal survival") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 


#### Plot Pupal survival curve for all populations #####

colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")


#### Unscaled plot #####

plot(prob.Pupal.Survival ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,1.1), data = data.pPS.HOP, type = "n", bty = "n",
     ylab = "probability", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Probability pupal survival", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(pPS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#313695")
lines(pPS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#4575b4")
lines(pPS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#74add1")
lines(pPS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#abd9e9")
lines(pPS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#e0f3f8")
lines(pPS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fee090")
lines(pPS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fdae61")
lines(pPS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#f46d43")
lines(pPS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#d73027")
lines(pPS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#a50026")

