# Trait Fits for Fecundity #

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for estimated fecundity (fec)
# 10 Aedes sierrensis populations with uniform priors 

#  1) Set-up,load packages, get data, etc.
#  2) Set up JAGS model and settings
#  3) Fit adult life span thermal responses (quadratic) using low info priors
#  4) Fit gamma distributions to adult lifespan prior thermal responses
#  5) Fit adult lifespan thermal responses with data-informed priors
#  6) Plot parameters

###### 1. Set up workspace, load packages, get data, etc. ####

# Set working directory
setwd("~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')
library(HDInterval)

data.fec = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# remove rows with missing data (e.g,, did not survive to adulthood )

data.fec = data.fec[!is.na(data.fec$EggCount),]

# only WAW, SB, POW had adults survive the 32 C treatment 
# for all other populations, add a '0' for egg count at 32 
# (i.e., individuals can't develop if they can't survive)

ManualHop = c(NA,"HOP", 32.18, rep(NA,20),0, rep(NA,2))
ManualMar1 = c(NA,"MARIN35", 32.18, rep(NA,20),0, rep(NA,2))
ManualMar2 = c(NA,"MARIN29", 32.18, rep(NA,20),0, rep(NA,2))
ManualEug = c(NA,"EUG", 32.18, rep(NA,20),0, rep(NA,2))
ManualPLA = c(NA,"PLA", 32.18, rep(NA,20),0, rep(NA,2))
ManualJRA = c(NA,"JRA", 32.18, rep(NA,20),0, rep(NA,2))
ManualPAR = c(NA,"PAR", 32.18, rep(NA,20),0, rep(NA,2))

data.fec= rbind(data.fec, ManualHop, ManualMar1, ManualMar2, ManualEug, 
               ManualPLA, ManualJRA, ManualPAR)

data.fec$Temp.Treatment = as.numeric(data.fec$Temp.Treatment)
data.fec$EggCount = as.numeric(data.fec$EggCount)

# Subset data by population
data.fec.HOP <- subset(data.fec, Population == "HOP")
data.fec.MAR1 <- subset(data.fec, Population == "MARIN35")
data.fec.MAR2 <- subset(data.fec, Population == "MARIN29")
data.fec.WAW <- subset(data.fec, Population == "WAW")
data.fec.EUG <- subset(data.fec, Population == "EUG")
data.fec.PLA <- subset(data.fec, Population == "PLA")
data.fec.SB <- subset(data.fec, Population == "SB")
data.fec.JRA <- subset(data.fec, Population == "JRA")
data.fec.PAR <- subset(data.fec, Population == "PAR")
data.fec.POW <- subset(data.fec, Population == "POW")

# Plot trait data
plot(EggCount ~ Temp.Treatment, xlim = c(5, 40), ylim = c(0,110), pch = 16,
     data = data.fec, ylab = "Fecundity (estimated)", xlab = "Temperature")
# appears roughly quadratic

###### 2a. Set up JAGS Models ####
# Quadratic Model with uniform priors

sink("quad.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
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

###### 2b. Shared settings for all models ####

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




###### 3. Fit using uniform priors ######
###### 3a: Hopland (HOP) ######

# set data
data <- data.fec.HOP 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.HOP.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.HOP.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.HOP , ylab = "fecundity for Hop prior", xlab = "Temperature")
lines(fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3b: Marin35 (MAR1) ######

# set data
data <- data.fec.MAR1 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.MAR1.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.MAR1.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.MAR1.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.MAR1 , ylab = "fec for MAR1", xlab = "Temperature")
lines(fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3c: Marin29 (MAR2) ######

# set data
data <- data.fec.MAR2 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.MAR2.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.MAR2.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.MAR2.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.MAR2 , ylab = "fec for MAR2 prior", xlab = "Temperature")
lines(fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3d: Wawona (WAW) ######

# set data
data <- data.fec.WAW 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.WAW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.WAW.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.WAW.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.WAW , ylab = "fec for WAW prior", xlab = "Temperature")
lines(fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3e: Eugene (EUG) ######

# set data
data <- data.fec.EUG 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
jag.data1<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.EUG.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.EUG.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.EUG.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.EUG , ylab = "fec for EUG prior", xlab = "Temperature")
lines(fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3f: Placer (PLA) ######

# set data
data <- data.fec.PLA 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.PLA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.PLA.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.PLA.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.PLA , ylab = "fec for PLA prior", xlab = "Temperature")
lines(fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3g: Santa Barbara (SB) ######

# set data
data <- data.fec.SB 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.SB.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.SB.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.SB.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.SB , ylab = "fec for SB prior", xlab = "Temperature")
lines(fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3h: Jasper Ridge (JRA) ######

# set data
data <- data.fec.JRA 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.JRA.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.JRA.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.JRA.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.JRA , ylab = "fec for JRA prior", xlab = "Temperature")
lines(fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 3i: Paso Robles (PAR) ######

# set data
data <- data.fec.PAR 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.PAR.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine output
fec.PAR.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.PAR.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.PAR , ylab = "fec for PAR prior", xlab = "Temperature")
lines(fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 3j: Powder Canyon (POW) ######

# set data
data <- data.fec.POW 

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.POW.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.POW.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.POW.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.POW , ylab = "fec for POW prior", xlab = "Temperature")
lines(fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

####### 4. Get TPC parameters #####
###### 4a. Topt of AL for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 4b. Pmax of AL for each population ######

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
# T0, Tm, and Topt along with 95% credible intervfecs (2.5%, 97.5%)

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

fecdf_uni <- rbind.data.frame(ParamCompile(fec.HOP.out), ParamCompile(fec.MAR1.out),
                              ParamCompile(fec.MAR2.out), ParamCompile(fec.WAW.out),
                              ParamCompile(fec.EUG.out), ParamCompile(fec.PLA.out),
                              ParamCompile(fec.SB.out), ParamCompile(fec.JRA.out),
                              ParamCompile(fec.PAR.out), ParamCompile(fec.POW.out))
fecdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(fecdf_uni, file = "fec_meansd_uni.Rsave")

###### 5. Plotting ######
load("fec_meansd_uni.Rsave")
library(ggplot2)
###### 5a. Plot of CTmax, CTmin, Topt for each population ######
# order by latitude
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

fecdf_uni$Population = factor(fecdf_uni$Population, levels = LatOrder)
fecdf_uni = fecdf_uni[order(fecdf_uni$Population),]
fecdf_uni = fecdf_uni[fecdf_uni$param != "Tbreadth",]
fecdf_uni = fecdf_uni[fecdf_uni$param != "Pmax",]

plotfecuni = ggplot(fecdf_uni, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Fecundity (estimated)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=14), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA))

###### 5b. Plot of Pmax for each population #####

load("fec_meansd_uni.Rsave")
LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors = rep(rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")), each = 3)

fecdf_uni$Population = factor(fecdf_uni$Population, levels = LatOrder)
fecdf_uni = fecdf_uni[order(fecdf_uni$Population),]
fecdf_uni = fecdf_uni[fecdf_uni$param == "Pmax",]

plotfec_Pmax = ggplot(ALdf_uni, aes(x=Population, y=mean, col = LatColors2)) + 
  scale_y_continuous(limits = c(0, 20)) +
  geom_point(stat="identity", col=LatColors2, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col =LatColors2,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Fecundity (estimated)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5c. Plot fecundity curves for all populations #####

EUGdata = fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
par(mar = (c(4.5, 4.5, 2.5, 1)))

plot(EggCount ~ Temp.Treatment, 
     xlim = c(0, 50), ylim = c(0,80), data = data.fec.HOP, type = "n", bty = "n",
     ylab = "# eggs laid", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Fecundity (estimate)", cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)
lines(EUGdata ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = "#7cae00", lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = "#fed439ff", lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(PARdata ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(SBdata ~ Temp.xs, col = "#ec3c30", lwd = 1.5)
lines(POWdata ~ Temp.xs, col = "#ab041b", lwd = 1.5)






###### 6. Generate informative priors: Leave-one-out with uniform priors ######
# Use same uniform priors settings as above
# Subset Data for leave-one-out
data.fec.HOP.prior <- subset(data.fec, Population != "HOP")
data.fec.MAR1.prior <- subset(data.fec, Population != "MARIN35")
data.fec.MAR2.prior <- subset(data.fec, Population != "MARIN29")
data.fec.WAW.prior <- subset(data.fec, Population != "WAW")
data.fec.EUG.prior <- subset(data.fec, Population != "EUG")
data.fec.PLA.prior <- subset(data.fec, Population != "PLA")
data.fec.SB.prior <- subset(data.fec, Population != "SB")
data.fec.JRA.prior <- subset(data.fec, Population != "JRA")
data.fec.PAR.prior <- subset(data.fec, Population != "PAR")
data.fec.POW.prior <- subset(data.fec, Population != "POW")


###### 6a. Hopland #####
# set data
data <- data.fec.HOP.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.HOP.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(AL.HOP.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.HOP.prior, ylab = "Fecundity for Hop prior", xlab = "Temperature")
lines(fec.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6b: Marin35 (MAR1) #####

# set data
data <- data.fec.MAR1.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.MAR1.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.MAR1.prior, ylab = "fec for MAR1 prior", xlab = "Temperature")
lines(fec.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.fec.MAR2.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.MAR2.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.MAR2.prior, ylab = "fec for MAR2 prior", xlab = "Temperature")
lines(fec.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.fec.WAW.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.WAW.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.WAW.prior, ylab = "fec for WAW prior", xlab = "Temperature")
lines(fec.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.fec.EUG.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.EUG.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.EUG.prior, ylab = "fec for EUG prior", xlab = "Temperature")
lines(fec.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.fec.PLA.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.PLA.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.PLA.prior, ylab = "fec for PLA prior", xlab = "Temperature")
lines(fec.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.fec.SB.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.SB.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.SB.prior, ylab = "fec for SB prior", xlab = "Temperature")
lines(fec.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.fec.JRA.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.JRA.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.JRA.prior, ylab = "fec for JRA prior", xlab = "Temperature")
lines(fec.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.fec.PAR.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.PAR.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.PAR.prior, ylab = "fec for PAR prior", xlab = "Temperature")
lines(fec.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyon (POW) ######

# set data
data <- data.fec.POW.prior

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
fec.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(fec.POW.prior.out)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.POW.prior, ylab = "fec for POW prior", xlab = "Temperature")
lines(fec.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### Save posterior distributions for informative priors #####

fec.HOP.prior.cf.dists <- data.frame(q = as.vector(fec.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.HOP.prior.gamma.fits = apply(fec.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.MAR1.prior.cf.dists <- data.frame(q = as.vector(fec.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(fec.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(fec.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.MAR1.prior.gamma.fits = apply(fec.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.MAR2.prior.cf.dists <- data.frame(q = as.vector(fec.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(fec.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(fec.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.MAR2.prior.gamma.fits = apply(fec.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.WAW.prior.cf.dists <- data.frame(q = as.vector(fec.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.WAW.prior.gamma.fits = apply(fec.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.EUG.prior.cf.dists <- data.frame(q = as.vector(fec.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.EUG.prior.gamma.fits = apply(fec.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.PLA.prior.cf.dists <- data.frame(q = as.vector(fec.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.PLA.prior.gamma.fits = apply(fec.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.SB.prior.cf.dists <- data.frame(q = as.vector(fec.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(fec.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(fec.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.SB.prior.gamma.fits = apply(fec.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.JRA.prior.cf.dists <- data.frame(q = as.vector(fec.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.JRA.prior.gamma.fits = apply(fec.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.PAR.prior.cf.dists <- data.frame(q = as.vector(fec.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.PAR.prior.gamma.fits = apply(fec.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.POW.prior.cf.dists <- data.frame(q = as.vector(fec.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(fec.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(fec.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
fec.POW.prior.gamma.fits = apply(fec.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

fec.hypers <- list(fec.HOP.prior.gamma.fits, fec.MAR1.prior.gamma.fits, fec.MAR2.prior.gamma.fits,
                   fec.WAW.prior.gamma.fits, fec.EUG.prior.gamma.fits, fec.PLA.prior.gamma.fits,
                   fec.SB.prior.gamma.fits, fec.JRA.prior.gamma.fits, fec.PAR.prior.gamma.fits, fec.POW.prior.gamma.fits)
save(fec.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/fechypers.Rsave")
###### Quadratic Model with gamma priors (except sigma) ######

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

###### 7. Fit using informative priors #####
# Load hyperparameters to bypass above code of generating them 
load("fechypers.Rsave")
fec.HOP.prior.gamma.fits <- fec.hypers[[1]]
fec.MAR1.prior.gamma.fits <- fec.hypers[[2]]
fec.MAR2.prior.gamma.fits <- fec.hypers[[3]]
fec.WAW.prior.gamma.fits <- fec.hypers[[4]]
fec.EUG.prior.gamma.fits <- fec.hypers[[5]]
fec.PLA.prior.gamma.fits <- fec.hypers[[6]]
fec.SB.prior.gamma.fits <- fec.hypers[[7]]
fec.JRA.prior.gamma.fits <- fec.hypers[[8]]
fec.PAR.prior.gamma.fits <- fec.hypers[[9]]
fec.POW.prior.gamma.fits <- fec.hypers[[10]]

###### 7a. Hopland ######

# Set data 
data <- data.fec.HOP
hypers <- fec.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(fec.HOP.out.inf)

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.HOP, ylab = "fec for Hopland", xlab = "Temperature", pch = 1)
lines(fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

###### 7b. Marin35 ######

# Set data 
data <- data.fec.MAR1
hypers <- fec.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.MAR1, ylab = "fec for Marin35", xlab = "Temperature", pch = 1)
lines(fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.fec.MAR2
hypers <- fec.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.MAR2, ylab = "fec for Marin29", xlab = "Temperature", pch = 1)
lines(fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.fec.WAW
hypers <- fec.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.WAW, ylab = "fec for Wawona", xlab = "Temperature", pch = 1)
lines(fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.fec.EUG
hypers <- fec.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.EUG, ylab = "fec for Eugene", xlab = "Temperature", pch = 1)
lines(fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.fec.PLA
hypers <- fec.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.PLA, ylab = "fec for Placer", xlab = "Temperature", pch = 1)
lines(fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.fec.SB
hypers <- fec.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.SB, ylab = "fec for Santa Barbara", xlab = "Temperature", pch = 1)
lines(fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.fec.JRA
hypers <- fec.JRA.prior.gamma.fits * 0.05

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.JRA, ylab = "fec for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.fec.PAR
hypers <- fec.PAR.prior.gamma.fits * 0.05

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.PAR, ylab = "fec for Paso Robles", xlab = "Temperature", pch = 1)
lines(fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon ######

# Set data 
data <- data.fec.POW
hypers <- fec.POW.prior.gamma.fits * 0.05

# Organize Data for JAGS
trait <- data$EggCount
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
fec.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
fec.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(EggCount ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,110), data = data.fec.POW, ylab = "fec for Powder Canyon", xlab = "Temperature", pch = 1)
lines(fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###### 10 panel plots of fits with raw data as means & sderrors  #####
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$EggCount ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$EggCount ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$EggCount ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/fec.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

fec.Hop.MeansSds = MeanSd(data.fec.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.Hop.MeansSds, ylab = "fec for Hop", xlab = "Temperature", pch= 16)
arrows(fec.Hop.MeansSds$Temp.Treatment, fec.Hop.MeansSds$mean-fec.Hop.MeansSds$stderror, fec.Hop.MeansSds$Temp.Treatment, fec.Hop.MeansSds$mean + fec.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 37, y = 74, cex = 1.5)
lines(fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.MAR1.MeansSds = MeanSd(data.fec.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.MAR1.MeansSds, ylab = "fec for MAR1", xlab = "Temperature", pch= 16)
arrows(fec.MAR1.MeansSds$Temp.Treatment, fec.MAR1.MeansSds$mean-fec.MAR1.MeansSds$stderror, fec.MAR1.MeansSds$Temp.Treatment, fec.MAR1.MeansSds$mean + fec.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 37, y = 74, cex = 1.5)
lines(fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.MAR2.MeansSds = MeanSd(data.fec.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.MAR2.MeansSds, ylab = "fec for MAR2", xlab = "Temperature", pch= 16)
arrows(fec.MAR2.MeansSds$Temp.Treatment, fec.MAR2.MeansSds$mean-fec.MAR2.MeansSds$stderror, fec.MAR2.MeansSds$Temp.Treatment, fec.MAR2.MeansSds$mean + fec.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 37, y = 74, cex = 1.5)
lines(fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.WAW.MeansSds = MeanSd(data.fec.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.WAW.MeansSds, ylab = "fec for WAW", xlab = "Temperature", pch= 16)
arrows(fec.WAW.MeansSds$Temp.Treatment, fec.WAW.MeansSds$mean-fec.WAW.MeansSds$stderror, fec.WAW.MeansSds$Temp.Treatment, fec.WAW.MeansSds$mean + fec.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 37, y = 74, cex = 1.5)
lines(fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.EUG.MeansSds = MeanSd(data.fec.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.EUG.MeansSds, ylab = "fec for EUG", xlab = "Temperature", pch= 16)
arrows(fec.EUG.MeansSds$Temp.Treatment, fec.EUG.MeansSds$mean-fec.EUG.MeansSds$stderror, fec.EUG.MeansSds$Temp.Treatment, fec.EUG.MeansSds$mean + fec.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG",x = 37, y = 74, cex = 1.5)
lines(fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.PLA.MeansSds = MeanSd(data.fec.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.PLA.MeansSds, ylab = "fec for PLA", xlab = "Temperature", pch= 16)
arrows(fec.PLA.MeansSds$Temp.Treatment, fec.PLA.MeansSds$mean-fec.PLA.MeansSds$stderror, fec.PLA.MeansSds$Temp.Treatment, fec.PLA.MeansSds$mean + fec.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 37, y = 74, cex = 1.5)
lines(fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.SB.MeansSds = MeanSd(data.fec.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.SB.MeansSds, ylab = "fec for SB", xlab = "Temperature", pch= 16)
arrows(fec.SB.MeansSds$Temp.Treatment, fec.SB.MeansSds$mean-fec.SB.MeansSds$stderror, fec.SB.MeansSds$Temp.Treatment, fec.SB.MeansSds$mean + fec.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 37, y = 74, cex = 1.5)
lines(fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.JRA.MeansSds = MeanSd(data.fec.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.JRA.MeansSds, ylab = "fec for JRA", xlab = "Temperature", pch= 16)
arrows(fec.JRA.MeansSds$Temp.Treatment, fec.JRA.MeansSds$mean-fec.JRA.MeansSds$stderror, fec.JRA.MeansSds$Temp.Treatment, fec.JRA.MeansSds$mean + fec.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 37, y = 74, cex = 1.5)
lines(fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.PAR.MeansSds = MeanSd(data.fec.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.PAR.MeansSds, ylab = "fec for PAR", xlab = "Temperature", pch= 16)
arrows(fec.PAR.MeansSds$Temp.Treatment, fec.PAR.MeansSds$mean-fec.PAR.MeansSds$stderror, fec.PAR.MeansSds$Temp.Treatment, fec.PAR.MeansSds$mean + fec.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 37, y = 74, cex = 1.5)
lines(fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

fec.POW.MeansSds = MeanSd(data.fec.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,80), data = fec.POW.MeansSds, ylab = "fec for POW", xlab = "Temperature", pch= 16)
arrows(fec.POW.MeansSds$Temp.Treatment, fec.POW.MeansSds$mean-fec.POW.MeansSds$stderror, fec.POW.MeansSds$Temp.Treatment, fec.POW.MeansSds$mean + fec.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 37, y = 74, cex = 1.5)
lines(fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col =  "grey35")
lines(fec.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

dev.off()



###### 8. Get TPC parameters #####
###### 8a. Topt of AL for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of AL for each population ######

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

# function to compile mean & 95% CrIs for alll parameters for each population 
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

fecdf_inf <- rbind.data.frame(ParamCompile(fec.HOP.out.inf), ParamCompile(fec.MAR1.out.inf),
                              ParamCompile(fec.MAR2.out.inf), ParamCompile(fec.WAW.out.inf),
                              ParamCompile(fec.EUG.out.inf), ParamCompile(fec.PLA.out.inf),
                              ParamCompile(fec.SB.out.inf), ParamCompile(fec.JRA.out.inf),
                              ParamCompile(fec.PAR.out.inf), ParamCompile(fec.POW.out.inf))
fecdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(fecdf_inf, file = "fec_meansd_inf.Rsave")



###### 9. Plotting  #####
library(ggplot2)
load("fec_meansd_inf.Rsave")

###### 9a. Plot of CTmax, CTmin, Topt for each population #####

LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors = rep(rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695")), each = 3)

fecdf_inf$Population = factor(fecdf_inf$Population, levels = LatOrder)
fecdf_inf = fecdf_inf[order(fecdf_inf$Population),]
fecdf_inf = fecdf_inf[fecdf_inf$param != "Tbreadth",]
fecdf_inf = fecdf_inf[fecdf_inf$param != "Pmax",]

plotfecinf = ggplot(fecdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, size = 2, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors, size = 0.7,
                position=position_dodge(.9), width = 0.5) + 
  ggtitle("Fecundity") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=14), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA))

###### 9c. Plot adult lifespans curve for all populations #####
load(file = "fec_meansd_inf.Rsave")

EUGdata = fec.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = fec.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = fec.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = fec.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = fec.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = fec.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = fec.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = fec.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = fec.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = fec.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
par(mar = (c(4.5, 4.5, 2.5, 1)))
plot(EggCount ~ Temp.Treatment, 
     xlim = c(5, 40), ylim = c(0,80), data = data.fec.HOP, type = "n", bty = "n",
     ylab = "# eggs laid", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Fecundity", cex.main = 1.7, cex.lab = 1.4, cex.axis = 1.2)
box(col = "black")
lines(EUGdata ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = "#7cae00", lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = "#fed439ff", lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(PARdata ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(SBdata ~ Temp.xs, col = "#ec3c30", lwd = 1.5)
lines(POWdata ~ Temp.xs, col = "#ab041b", lwd = 1.5)




 

