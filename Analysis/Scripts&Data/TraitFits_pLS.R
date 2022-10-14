# Trait Fits for Larval Survival ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for larval survival pLS 
# 10 Aedes sierrensis populations - 1) with uniform priors and 2) data-informed priors

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit Larval survival thermal responses using low info priors
#  4) Fit gamma distributions to larval survival prior thermal responses
#  5) Fit larval survival thermal responses with data-informed priors


#### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraries for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

data.pLS_int = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# Calculate survival for each population and temp treatment
data.pLS_surv = aggregate(Larval.Survival ~ Population + Temp.Treatment, data = data.pLS_int, FUN = sum)
data.pLS_all = aggregate(Larval.Survival ~ Population + Temp.Treatment, data = data.pLS_int, FUN = length)
data.pLS_per = data.pLS_surv$Larval.Survival/data.pLS_all$Larval.Survival 

data.pLS = cbind.data.frame(data.pLS_surv[,1:2], data.pLS_per)
colnames(data.pLS)[3] = "prob.Larval.Survival"

# excluding the coldest temp treatment since larvae did survive, they just didn't pupate
data.pLS = data.pLS[data.pLS$Temp.Treatment != "5.49",]

# Plot
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), data = data.pLS, ylab = "prob. larval survival", xlab = "Temperature")

# Subset data by population
data.pLS.HOP <- subset(data.pLS, Population == "HOP")
data.pLS.MAR1 <- subset(data.pLS, Population == "MAR35")
data.pLS.MAR2 <- subset(data.pLS, Population == "MAR29")
data.pLS.WAW <- subset(data.pLS, Population == "WAW")
data.pLS.EUG <- subset(data.pLS, Population == "EUG")
data.pLS.PLA <- subset(data.pLS, Population == "PLA")
data.pLS.SB <- subset(data.pLS, Population == "SB")
data.pLS.JRA <- subset(data.pLS, Population == "JRA")
data.pLS.PAR <- subset(data.pLS, Population == "PAR")
data.pLS.POW <- subset(data.pLS, Population == "POW")


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


#### 3. Generate data-informed priors using leave-one-out approach #####
# For each population, use data for all *other* populations to obtain
# distributions of each TPC parameter
# For this initial fitting, we use uniform priors, 
# constrained over a biologically realistic range

# Subset Data - leave-one-out
data.pLS.HOP.prior <- subset(data.pLS, Population != "HOP")
data.pLS.MAR1.prior <- subset(data.pLS, Population != "MAR1")
data.pLS.MAR2.prior <- subset(data.pLS, Population != "MAR2")
data.pLS.WAW.prior <- subset(data.pLS, Population != "WAW")
data.pLS.EUG.prior <- subset(data.pLS, Population != "EUG")
data.pLS.PLA.prior <- subset(data.pLS, Population != "PLA")
data.pLS.SB.prior <- subset(data.pLS, Population != "SB")
data.pLS.JRA.prior <- subset(data.pLS, Population != "JRA")
data.pLS.PAR.prior <- subset(data.pLS, Population != "PAR")
data.pLS.POW.prior <- subset(data.pLS, Population != "POW")

###### 3a. Hopland (HOP) ####

# Set data
data <- data.pLS.HOP.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.HOP.prior, ylab = "pLS for Hopland", xlab = "Temperature")
lines(pLS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b. Marin35 (MAR1) ####

# Set data
data <- data.pLS.MAR1.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.MAR1.prior, ylab = "pLS for Marin35", xlab = "Temperature")
lines(pLS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




###### 3c. Marin29 (MAR2) ####

# Set data
data <- data.pLS.MAR2.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR2.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.HOP.prior, ylab = "pLS for Marin29", xlab = "Temperature")
lines(pLS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3d: Wawona (WAW) ######

# set data
data <- data.pLS.WAW.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.WAW.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.WAW.prior, ylab = "pLS for WAW prior", xlab = "Temperature")
lines(pLS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.pLS.EUG.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.EUG.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.EUG.prior, ylab = "pLS for EUG prior", xlab = "Temperature")
lines(pLS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3f: Placer (PLA) ######

# set data
data <- data.pLS.PLA.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.PLA.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.PLA.prior, ylab = "pLS for PLA prior", xlab = "Temperature")
lines(pLS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3g: Santa Barbara (SB) ######

# set data
data <- data.pLS.SB.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.SB.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.SB.prior, ylab = "pLS for SB prior", xlab = "Temperature")
lines(pLS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.pLS.JRA.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.JRA.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.JRA.prior, ylab = "pLS for JRA prior", xlab = "Temperature")
lines(pLS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.pLS.PAR.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.PAR.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.PAR.prior, ylab = "pLS for PAR prior", xlab = "Temperature")
lines(pLS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3j: Powder Canyono (POW) ######

# set data
data <- data.pLS.POW.prior

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.POW.prior.out)

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.POW.prior, ylab = "pLS for POW prior", xlab = "Temperature")
lines(pLS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)













#### 4. Fit gamma distributions to these TPC parameters ####
# Get the posterior dists for 3 main parameters (not sigma) into a data frame
# Fit gamma distributions for each parameter posterior dists
# q = scaling coefficient 

###### 4a. Hopland ######
pLS.HOP.prior.cf.dists <- data.frame(q = as.vector(pLS.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLS.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLS.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))

pLS.HOP.prior.gamma.fits = apply(pLS.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4b. Marin35 (MAR1) ######

pLS.MAR1.prior.cf.dists <- data.frame(q = as.vector(pLS.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLS.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLS.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.MAR1.prior.gamma.fits = apply(pLS.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4c. Marin29 (MAR2) ######

pLS.MAR2.prior.cf.dists <- data.frame(q = as.vector(pLS.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLS.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLS.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.MAR2.prior.gamma.fits = apply(pLS.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4d. Wawona ######

pLS.WAW.prior.cf.dists <- data.frame(q = as.vector(pLS.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.WAW.prior.gamma.fits = apply(pLS.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4e. Eugene ######

pLS.EUG.prior.cf.dists <- data.frame(q = as.vector(pLS.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.EUG.prior.gamma.fits = apply(pLS.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4f. Placer ######

pLS.PLA.prior.cf.dists <- data.frame(q = as.vector(pLS.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.PLA.prior.gamma.fits = apply(pLS.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4g. Santa Barbara ######

pLS.SB.prior.cf.dists <- data.frame(q = as.vector(pLS.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(pLS.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(pLS.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.SB.prior.gamma.fits = apply(pLS.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4h. Jasper Ridge ######

pLS.JRA.prior.cf.dists <- data.frame(q = as.vector(pLS.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.JRA.prior.gamma.fits = apply(pLS.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4i. Paso Robles ######

pLS.PAR.prior.cf.dists <- data.frame(q = as.vector(pLS.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.PAR.prior.gamma.fits = apply(pLS.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### 4j. Powder Canyon ######

pLS.POW.prior.cf.dists <- data.frame(q = as.vector(pLS.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.POW.prior.gamma.fits = apply(pLS.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###### Hyperparameter list
pLS.hypers <- list(pLS.HOP.prior.gamma.fits, pLS.MAR1.prior.gamma.fits, pLS.MAR2.prior.gamma.fits,
                   pLS.WAW.prior.gamma.fits, pLS.EUG.prior.gamma.fits, pLS.PLA.prior.gamma.fits,
                   pLS.SB.prior.gamma.fits, pLS.JRA.prior.gamma.fits, pLS.PAR.prior.gamma.fits, pLS.POW.prior.gamma.fits)
save(pLS.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/pLShypers.Rsave")





#### 5. Fit pLS thermal responses with data-informed priors #####

# to bypass above code of generating hyperparameters
load("pLShypers.Rsave")
pLS.HOP.prior.gamma.fits <- pLS.hypers[[1]]
pLS.MAR1.prior.gamma.fits <- pLS.hypers[[2]]
pLS.MAR2.prior.gamma.fits <- pLS.hypers[[3]]
pLS.WAW.prior.gamma.fits <- pLS.hypers[[4]]
pLS.EUG.prior.gamma.fits <- pLS.hypers[[5]]
pLS.PLA.prior.gamma.fits <- pLS.hypers[[6]]
pLS.SB.prior.gamma.fits <- pLS.hypers[[7]]
pLS.JRA.prior.gamma.fits <- pLS.hypers[[8]]
pLS.PAR.prior.gamma.fits <- pLS.hypers[[9]]
pLS.POW.prior.gamma.fits <- pLS.hypers[[10]]



###### 5a. Hopland ######

# Set data 
data <- data.pLS.HOP
hypers <- pLS.HOP.prior.gamma.fits * .1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.HOP.out.inf$BUGSoutput$summary[1:5,]

#save(pLS.HOP.out.inf, file = "jagsout_pLS_HOP_inf.Rdata")

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.HOP, ylab = "prob larval survival for Hopland", xlab = "Temperature", pch = 1)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5b. Marin35 ######

# Set data 
data <- data.pLS.MAR1
hypers <- pLS.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.MAR1, ylab = "pLS for Marin35", xlab = "Temperature", pch = 1)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5c. Marin29 ######

# Set data 
data <- data.pLS.MAR2
hypers <- pLS.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.MAR2, ylab = "pLS for Marin29", xlab = "Temperature", pch = 1)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5d. Wawona ######

# Set data 
data <- data.pLS.WAW
hypers <- pLS.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.WAW, ylab = "pLS for Wawona", xlab = "Temperature", pch = 1)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5e. Eugene ######

# Set data 
data <- data.pLS.EUG
hypers <- pLS.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.EUG, ylab = "pLS for Eugene", xlab = "Temperature", pch = 1)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5f. Placer ######

# Set data 
data <- data.pLS.PLA
hypers <- pLS.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.PLA, ylab = "pLS for Placer", xlab = "Temperature", pch = 1)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5g. Santa Barbara ######

# Set data 
data <- data.pLS.SB
hypers <- pLS.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.SB, ylab = "pLS for Santa Barbara", xlab = "Temperature", pch = 1)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 5g. Jasper Ridge ######

# Set data 
data <- data.pLS.JRA
hypers <- pLS.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.JRA, ylab = "pLS for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5i. Paso Robles ######

# Set data 
data <- data.pLS.PAR
hypers <- pLS.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.PAR, ylab = "pLS for Paso Robles", xlab = "Temperature", pch = 1)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 5j. Powder canyon ######

# Set data 
data <- data.pLS.POW
hypers <- pLS.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$prob.Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(prob.Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLS.POW, ylab = "pLS for Powder Canyon", xlab = "Temperature", pch = 1)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### Topt of pLS for each population ######
ToptHOP = Temp.xs[which.max(as.vector(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18
ToptMAR1 = Temp.xs[which.max(as.vector(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18
ToptMAR2 = Temp.xs[which.max(as.vector(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptWAW = Temp.xs[which.max(as.vector(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19
ToptEUG = Temp.xs[which.max(as.vector(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 19
ToptPLA = Temp.xs[which.max(as.vector(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptSB = Temp.xs[which.max(as.vector(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptJRA = Temp.xs[which.max(as.vector(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptPAR = Temp.xs[which.max(as.vector(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18.5
ToptPOW = Temp.xs[which.max(as.vector(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))] # 18

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
TbreadthHOP = Tbreadth(pLS.HOP.out.inf) # 20.5
TbreadthMAR1 = Tbreadth(pLS.MAR1.out.inf) # 21
TbreadthMAR2 = Tbreadth(pLS.MAR2.out.inf) # 20.5
TbreadthWAW = Tbreadth(pLS.WAW.out.inf) # 20.5
TbreadthEUG = Tbreadth(pLS.EUG.out.inf) # 21
TbreadthPLA = Tbreadth(pLS.PLA.out.inf) # 21
TbreadthSB = Tbreadth(pLS.SB.out.inf) # 21
TbreadthJRA = Tbreadth(pLS.JRA.out.inf) # 21
TbreadthPAR = Tbreadth(pLS.PAR.out.inf) # 21.5
TbreadthPOW = Tbreadth(pLS.POW.out.inf) # 21.5

#### Compile prediction data for all populations ###### 

pLSmeansd <- list(pLS.HOP.out.inf$BUGSoutput$summary[1:2,1:2], pLS.MAR1.out.inf$BUGSoutput$summary[1:2,1:2], 
                  pLS.MAR2.out.inf$BUGSoutput$summary[1:2,1:2], pLS.WAW.out.inf$BUGSoutput$summary[1:2,1:2], 
                  pLS.EUG.out.inf$BUGSoutput$summary[1:2,1:2], pLS.PLA.out.inf$BUGSoutput$summary[1:2,1:2],
                  pLS.SB.out.inf$BUGSoutput$summary[1:2,1:2], pLS.JRA.out.inf$BUGSoutput$summary[1:2,1:2],
                  pLS.POW.out.inf$BUGSoutput$summary[1:2,1:2], pLS.PAR.out.inf$BUGSoutput$summary[1:2,1:2])
pLStops = as.data.frame(matrix(c(ToptHOP, ToptMAR1, ToptMAR2, ToptWAW, ToptEUG, ToptPLA, 
                                 ToptSB, ToptJRA, ToptPAR, ToptPOW, rep(NA, 10)), nrow = 10, ncol = 2))
pLStbreadths = as.data.frame(matrix(c(TbreadthHOP, TbreadthMAR1,TbreadthMAR2, TbreadthWAW, TbreadthEUG, TbreadthPLA, 
                                      TbreadthSB, TbreadthJRA, TbreadthPAR, TbreadthPOW, rep(NA, 10)), nrow = 10, ncol = 2))
colnames(pLStops) = c("mean", "sd")
colnames(pLStbreadths) = c("mean", "sd")
pLSdf2 = do.call(rbind.data.frame, pLSmeansd)
pLSdf_1 = rbind(pLSdf2, pLStops)
pLSdf = rbind(pLSdf_1, pLStbreadths)
pLSdf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 2),
                     rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), 2))
pLSdf$Param = c(rep(c("Tmin", "Tmax"), 10), rep("Topt", 10), rep("Tbreadth", 10))
save(pLSdf, file = "pLS_meansd.Rsave")






#### Plot TPC parameters for each population #####

load("pLS_meansd.Rsave")
library(ggplot2)
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)

# order from coldest to hottest temp in wettest quarter
#OldPopOrder = c("EUG", "WAW","HOP", "PLA", "SB", "MAR2", "MAR1", "PAR", "JRA", "POW") 
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pLSdf$Population = factor(pLSdf$Population, levels = PopOrder)
pLSdf = pLSdf[order(pLSdf$Population),]
pLSdf = pLSdf[pLSdf$Param != "Tbreadth",]

# Parameter plot
plotpLS = ggplot(pLSdf, aes(x=Population, y=as.numeric(mean), col = colors3)) + 
  scale_y_continuous(limits = c(0, 45)) +
  geom_point(stat="identity", col=colors3, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), col = colors3,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Larval survival") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Temperature (C)") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none")

#### Plot larval survival curve for all populations #####

###### Unscaled plot ######
plot(prob.Larval.Survival ~ Temp.Treatment, 
     xlim = c(5, 45), ylim = c(0,1.1), data = data.pLS.HOP, type = "n", bty = "n",
     ylab = "probability", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Larval survival", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#313695")
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#4575b4")
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#74add1")
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#abd9e9")
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#e0f3f8")
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fee090")
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#fdae61")
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#f46d43")
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#d73027")
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "#a50026")

