# Trait Fits for Larval Survival ##

# Lisa Couper, Stanford University. Modified from Shocket et al. 2020 
# Code summary:
# Use Bayesian Inference (JAGS) to fit temperature-dependent functions for larval survival pLS 
# 10 Aedes sierrensis populations with uniform priors 

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

data.pLS = read.csv("LifeHistoryTraitExp_Data.csv")[,-1]
# remove rows with missing larval survival data
data.pLS = data.pLS[!is.na(data.pLS$Larval.Survival),]

# Subset data by population
data.pLS.HOP <- subset(data.pLS, Population == "HOP")
data.pLS.MAR1 <- subset(data.pLS, Population == "MARIN35")
data.pLS.MAR2 <- subset(data.pLS, Population == "MARIN29")
data.pLS.WAW <- subset(data.pLS, Population == "WAW")
data.pLS.EUG <- subset(data.pLS, Population == "EUG")
data.pLS.PLA <- subset(data.pLS, Population == "PLA")
data.pLS.SB <- subset(data.pLS, Population == "SB")
data.pLS.JRA <- subset(data.pLS, Population == "JRA")
data.pLS.PAR <- subset(data.pLS, Population == "PAR")
data.pLS.POW <- subset(data.pLS, Population == "POW")

# Plot
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), data = data.pLS, ylab = "Larval survival", xlab = "Temperature")

##### 2a. JAGS Models ####

# Quadratic Model with uniform priors and binomial likelihood model :

sink("quad.binom.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 40)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
      for(i in 1:N.obs)
      {trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
      trait[i] ~ dbin(ifelse(trait.mu[i] > 0.99, 0.99, trait.mu[i]), 1)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- ifelse(-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) > 1, 0.99, 
    -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]))
    }
    
    } # close model 
    ",fill=TRUE)
sink()


##### 2b. Shared settings for all models ####

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


##### 3. Fit using uniform priors #####
##### 3a. Hopland #####

# Set data
data <- data.pLS.HOP

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom1.txt
pLS.HOP.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.HOP.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.HOP, ylab = "pLS for Hopland", xlab = "Temperature")
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3b. Marin35 (MAR1) ######

# Set data
data <- data.pLS.MAR1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom1.txt
pLS.MAR1.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.MAR1.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.MAR1, ylab = "pLS for MAR1land", xlab = "Temperature")
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3c. Marin29 (MAR2) #####

# Set data
data <- data.pLS.MAR2

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.MAR2.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.MAR2.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.MAR2, ylab = "pLS for MAR2land", xlab = "Temperature")
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3d. Wawona (WAW) #####

# Set data
data <- data.pLS.WAW

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.WAW.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.WAW.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.WAW, ylab = "pLS for WAW", xlab = "Temperature")
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3e. Eugene (EUG) #####

# Set data
data <- data.pLS.EUG

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.EUG.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.EUG.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.EUG, ylab = "pLS for EUG", xlab = "Temperature")
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3f. Placer (PLA) ######
# Set data
data <- data.pLS.PLA

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.PLA.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.PLA.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.PLA, ylab = "pLS for PLA", xlab = "Temperature")
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3g. Santa Barbara (SB) #####
# Set data
data <- data.pLS.SB

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.SB.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.SB.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.SB, ylab = "pLS for SB", xlab = "Temperature")
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

##### 3h. Jasper Ridge (JRA) #####

# Set data
data <- data.pLS.JRA

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.JRA.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.JRA.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.JRA, ylab = "pLS for JRA", xlab = "Temperature")
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)






##### 3i. Paso Robles (PAR) #####

# Set data
data <- data.pLS.PAR

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.PAR.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.PAR.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.PAR, ylab = "pLS for PAR", xlab = "Temperature")
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##### 3j. Powder canyon (POW) #####
# Set data
data <- data.pLS.POW

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS using quad.binom.txt
pLS.POW.out <- jags(data=jag.data, inits = inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd(), quiet = FALSE)

# Examine Output
pLS.POW.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.POW, ylab = "pLS for POW", xlab = "Temperature")
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### For checking fits: 10 panel plot of fits using uniform priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.HOP, ylab = "pLS for Hop", xlab = "Temperature")
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.MAR1, ylab = "pLS for MAR1", xlab = "Temperature")
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.MAR2 , ylab = "pLS for MAR2 prior", xlab = "Temperature")
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.WAW , ylab = "pLS for WAW prior", xlab = "Temperature")
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.EUG , ylab = "pLS for EUG prior", xlab = "Temperature")
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.PLA , ylab = "pLS for PLA prior", xlab = "Temperature")
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.SB , ylab = "pLS for SB prior", xlab = "Temperature")
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.JRA , ylab = "pLS for JRA prior", xlab = "Temperature")
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.PAR , ylab = "pLS for PAR prior", xlab = "Temperature")
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = data.pLS.POW , ylab = "pLS for POW prior", xlab = "Temperature")
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 10 panel plots of fits with uniform priors, but raw data as means & sderrorss ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$Larval.Survival ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$Larval.Survival ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$Larval.Survival ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/pLS_UniformPriors.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

pLS.Hop.MeansSds = MeanSd(data.pLS.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.Hop.MeansSds, ylab = "pLS for Hop", xlab = "Temperature", pch= 16)
arrows(pLS.Hop.MeansSds$Temp.Treatment, pLS.Hop.MeansSds$mean-pLS.Hop.MeansSds$stderror, pLS.Hop.MeansSds$Temp.Treatment, pLS.Hop.MeansSds$mean + pLS.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("HOP", x = 36.5, y = 0.95, cex = 1.5)

pLS.MAR1.MeansSds = MeanSd(data.pLS.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.MAR1.MeansSds, ylab = "pLS for MAR1", xlab = "Temperature", pch= 16)
arrows(pLS.MAR1.MeansSds$Temp.Treatment, pLS.MAR1.MeansSds$mean-pLS.MAR1.MeansSds$stderror, pLS.MAR1.MeansSds$Temp.Treatment, pLS.MAR1.MeansSds$mean + pLS.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("MAR1", x = 36.5, y = 0.95, cex = 1.5)

pLS.MAR2.MeansSds = MeanSd(data.pLS.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.MAR2.MeansSds, ylab = "pLS for MAR2", xlab = "Temperature", pch= 16)
arrows(pLS.MAR2.MeansSds$Temp.Treatment, pLS.MAR2.MeansSds$mean-pLS.MAR2.MeansSds$stderror, pLS.MAR2.MeansSds$Temp.Treatment, pLS.MAR2.MeansSds$mean + pLS.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("MAR2", x = 36.5, y = 0.95, cex = 1.5)

pLS.WAW.MeansSds = MeanSd(data.pLS.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.WAW.MeansSds, ylab = "pLS for WAW", xlab = "Temperature", pch= 16)
arrows(pLS.WAW.MeansSds$Temp.Treatment, pLS.WAW.MeansSds$mean-pLS.WAW.MeansSds$stderror, pLS.WAW.MeansSds$Temp.Treatment, pLS.WAW.MeansSds$mean + pLS.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("WAW", x = 36.5, y = 0.95, cex = 1.5)

pLS.EUG.MeansSds = MeanSd(data.pLS.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.EUG.MeansSds, ylab = "pLS for EUG", xlab = "Temperature", pch= 16)
arrows(pLS.EUG.MeansSds$Temp.Treatment, pLS.EUG.MeansSds$mean-pLS.EUG.MeansSds$stderror, pLS.EUG.MeansSds$Temp.Treatment, pLS.EUG.MeansSds$mean + pLS.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("EUG", x = 36.5, y = 0.95, cex = 1.5)

pLS.PLA.MeansSds = MeanSd(data.pLS.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.PLA.MeansSds, ylab = "pLS for PLA", xlab = "Temperature", pch= 16)
arrows(pLS.PLA.MeansSds$Temp.Treatment, pLS.PLA.MeansSds$mean-pLS.PLA.MeansSds$stderror, pLS.PLA.MeansSds$Temp.Treatment, pLS.PLA.MeansSds$mean + pLS.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("PLA", x = 36.5, y = 0.95, cex = 1.5)

pLS.SB.MeansSds = MeanSd(data.pLS.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.SB.MeansSds, ylab = "pLS for SB", xlab = "Temperature", pch= 16)
arrows(pLS.SB.MeansSds$Temp.Treatment, pLS.SB.MeansSds$mean-pLS.SB.MeansSds$stderror, pLS.SB.MeansSds$Temp.Treatment, pLS.SB.MeansSds$mean + pLS.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("SB", x = 36.5, y = 0.95, cex = 1.5)

pLS.JRA.MeansSds = MeanSd(data.pLS.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.JRA.MeansSds, ylab = "pLS for JRA", xlab = "Temperature", pch= 16)
arrows(pLS.JRA.MeansSds$Temp.Treatment, pLS.JRA.MeansSds$mean-pLS.JRA.MeansSds$stderror, pLS.JRA.MeansSds$Temp.Treatment, pLS.JRA.MeansSds$mean + pLS.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("JRA", x = 36.5, y = 0.95, cex = 1.5)

pLS.PAR.MeansSds = MeanSd(data.pLS.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.PAR.MeansSds, ylab = "pLS for PAR", xlab = "Temperature", pch= 16)
arrows(pLS.PAR.MeansSds$Temp.Treatment, pLS.PAR.MeansSds$mean-pLS.PAR.MeansSds$stderror, pLS.PAR.MeansSds$Temp.Treatment, pLS.PAR.MeansSds$mean + pLS.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("PAR", x = 36.5, y = 0.95, cex = 1.5)

pLS.POW.MeansSds = MeanSd(data.pLS.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.POW.MeansSds, ylab = "pLS for POW", xlab = "Temperature", pch= 16)
arrows(pLS.POW.MeansSds$Temp.Treatment, pLS.POW.MeansSds$mean-pLS.POW.MeansSds$stderror, pLS.POW.MeansSds$Temp.Treatment, pLS.POW.MeansSds$mean + pLS.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
text("POW", x = 36.5, y = 0.95, cex = 1.5)

dev.off()

###### 4. Get TPC parameters #####
###### 4a. Topt of pLS for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 4b. Pmax of pLS for each population ######

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

pLSdf_uni <- rbind.data.frame(ParamCompile(pLS.HOP.out), ParamCompile(pLS.MAR1.out),
                              ParamCompile(pLS.MAR2.out), ParamCompile(pLS.WAW.out),
                              ParamCompile(pLS.EUG.out), ParamCompile(pLS.PLA.out),
                              ParamCompile(pLS.SB.out), ParamCompile(pLS.JRA.out),
                              ParamCompile(pLS.PAR.out), ParamCompile(pLS.POW.out))
pLSdf_uni$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(pLSdf_uni, file = "pLS_meansd_uni.Rsave")

###### 5. Plotting  #####

library(ggplot2)
plot.new()
par(mfrow = c(1,1))
par(mar = c(2.5, 2.5, 2.5, 2.5))

load("pLS_meansd_uni.Rsave")

###### 5a. Plot of CTmax, CTmin, Topt for each population #####

# order by latitude of collection
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)
pLSdf_uni$Population = factor(pLSdf_uni$Population, levels = LatOrder)
pLSdf_uni = pLSdf_uni[order(pLSdf_uni$Population),]
pLSdf_uni = pLSdf_uni[pLSdf_uni$param != "Tbreadth",]
pLSdf_uni = pLSdf_uni[pLSdf_uni$param != "Pmax",]

plotpLSuni = ggplot(pLSdf_uni, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("Larval survival probability") + 
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

load("pLS_meansd_uni.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pLSdf_uni$Population = factor(pLSdf_uni$Population, levels = PopOrder)
pLSdf_uni = pLSdf_uni[order(pLSdf_uni$Population),]
pLSdf_uni = pLSdf_uni[pLSdf_uni$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotpLS_Pmax = ggplot(pLSdf_uni, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.7, 1)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("pLS") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 

###### 5c. Plot pLS curves for all populations #####

load(file = "pLS_meansd_uni.Rsave")
colors3 = rep(rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695")), each = 3)
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")

EUGdata = pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + pLS
plot(Larval.Survival ~ Temp.Treatment, 
     xlim = c(5, 35), ylim = c(0,1), data = data.pLS.HOP, type = "n", bty = "n",
     ylab = "pLS", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "pLS", cex.main = 2, cex.lab = 1.4, cex.axis = 1.2)
lines(EUGdata ~ Temp.xs, col = "#313695", lwd = 1.5)
lines(WAWdata ~ Temp.xs, col = "#4575b4", lwd = 1.5)
lines(PLAdata ~ Temp.xs, col = "#74add1", lwd = 1.5)
lines(HOPdata ~ Temp.xs, col = "#abd9e9", lwd = 1.5)
lines(JRAdata ~ Temp.xs, col = "#e0f3f8", lwd = 1.5)
lines(MAR2data ~ Temp.xs, col = "#fee090", lwd = 1.5)
lines(MAR1data ~ Temp.xs, col = "#fdae61", lwd = 1.5)
lines(POWdata ~ Temp.xs, col = "#f46d43", lwd = 1.5)
lines(PARdata ~ Temp.xs, col = "#d73027", lwd = 1.5)
lines(SBdata ~ Temp.xs, col = "#a50026", lwd = 1.5)

###### 6. Generate informative priors: Leave-one-out with uniform priors ######
# Use same uniform priors settings as above
# Subset Data for leave-one-out
data.pLS.HOP.prior <- subset(data.pLS, Population != "HOP")
data.pLS.MAR1.prior <- subset(data.pLS, Population != "MARIN35")
data.pLS.MAR2.prior <- subset(data.pLS, Population != "MARIN29")
data.pLS.WAW.prior <- subset(data.pLS, Population != "WAW")
data.pLS.EUG.prior <- subset(data.pLS, Population != "EUG")
data.pLS.PLA.prior <- subset(data.pLS, Population != "PLA")
data.pLS.SB.prior <- subset(data.pLS, Population != "SB")
data.pLS.JRA.prior <- subset(data.pLS, Population != "JRA")
data.pLS.PAR.prior <- subset(data.pLS, Population != "PAR")
data.pLS.POW.prior <- subset(data.pLS, Population != "POW")

###### 6a. Hopland #####
# set data
data <- data.pLS.HOP.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.HOP.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.HOP.prior.out$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.HOP.prior, ylab = "pLS for Hop prior", xlab = "Temperature")
lines(pLS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6b: Marin35 (MAR1) #####

# set data
data <- data.pLS.MAR1.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.MAR1.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR1.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.MAR1.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.MAR1.prior, ylab = "pLS for MAR1 prior", xlab = "Temperature")
lines(pLS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6c: Marin29 (MAR2) ######

# set data
data <- data.pLS.MAR2.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.MAR2.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR2.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.MAR2.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.MAR2.prior, ylab = "pLS for MAR2 prior", xlab = "Temperature")
lines(pLS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6d: Wawona (WAW) ######

# set data
data <- data.pLS.WAW.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.WAW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.WAW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.WAW.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.WAW.prior, ylab = "pLS for WAW prior", xlab = "Temperature")
lines(pLS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6e: Eugene (EUG) ######

# set data
data <- data.pLS.EUG.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.EUG.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.EUG.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.EUG.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.EUG.prior, ylab = "pLS for EUG prior", xlab = "Temperature")
lines(pLS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 6f: Placer (PLA) ######

# set data
data <- data.pLS.PLA.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.PLA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PLA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.PLA.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.PLA.prior, ylab = "pLS for PLA prior", xlab = "Temperature")
lines(pLS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6g: Santa Barbara (SB) ######

# set data
data <- data.pLS.SB.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.SB.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.SB.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.SB.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.SB.prior, ylab = "pLS for SB prior", xlab = "Temperature")
lines(pLS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6h: Jasper Ridge (JRA) ######

# set data
data <- data.pLS.JRA.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.JRA.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.JRA.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.JRA.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.JRA.prior, ylab = "pLS for JRA prior", xlab = "Temperature")
lines(pLS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6i: Paso Robles (PAR) ######

# set data
data <- data.pLS.PAR.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.PAR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PAR.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.PAR.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.PAR.prior, ylab = "pLS for PAR prior", xlab = "Temperature")
lines(pLS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 6j: Powder Canyon (POW) ######

# set data
data <- data.pLS.POW.prior

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Run JAGS
pLS.POW.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.binom.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.POW.prior.out$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.POW.prior.out)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.POW.prior, ylab = "pLS for POW prior", xlab = "Temperature")
lines(pLS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### Save posterior distributions for informative priors #####

pLS.HOP.prior.cf.dists <- data.frame(q = as.vector(pLS.HOP.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.HOP.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.HOP.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.HOP.prior.gamma.fits = apply(pLS.HOP.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.MAR1.prior.cf.dists <- data.frame(q = as.vector(pLS.MAR1.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLS.MAR1.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLS.MAR1.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.MAR1.prior.gamma.fits = apply(pLS.MAR1.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.MAR2.prior.cf.dists <- data.frame(q = as.vector(pLS.MAR2.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLS.MAR2.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLS.MAR2.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.MAR2.prior.gamma.fits = apply(pLS.MAR2.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.WAW.prior.cf.dists <- data.frame(q = as.vector(pLS.WAW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.WAW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.WAW.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.WAW.prior.gamma.fits = apply(pLS.WAW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.EUG.prior.cf.dists <- data.frame(q = as.vector(pLS.EUG.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.EUG.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.EUG.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.EUG.prior.gamma.fits = apply(pLS.EUG.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.PLA.prior.cf.dists <- data.frame(q = as.vector(pLS.PLA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.PLA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.PLA.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.PLA.prior.gamma.fits = apply(pLS.PLA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.SB.prior.cf.dists <- data.frame(q = as.vector(pLS.SB.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(pLS.SB.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(pLS.SB.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.SB.prior.gamma.fits = apply(pLS.SB.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.JRA.prior.cf.dists <- data.frame(q = as.vector(pLS.JRA.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.JRA.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.JRA.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.JRA.prior.gamma.fits = apply(pLS.JRA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.PAR.prior.cf.dists <- data.frame(q = as.vector(pLS.PAR.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.PAR.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.PAR.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.PAR.prior.gamma.fits = apply(pLS.PAR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.POW.prior.cf.dists <- data.frame(q = as.vector(pLS.POW.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pLS.POW.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pLS.POW.prior.out$BUGSoutput$sims.list$cf.Tm))
pLS.POW.prior.gamma.fits = apply(pLS.POW.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

pLS.hypers <- list(pLS.HOP.prior.gamma.fits, pLS.MAR1.prior.gamma.fits, pLS.MAR2.prior.gamma.fits,
                   pLS.WAW.prior.gamma.fits, pLS.EUG.prior.gamma.fits, pLS.PLA.prior.gamma.fits,
                   pLS.SB.prior.gamma.fits, pLS.JRA.prior.gamma.fits, pLS.PAR.prior.gamma.fits, pLS.POW.prior.gamma.fits)
save(pLS.hypers, file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/pLShypers.Rsave")

###### Quadratic Model with gamma priors (except sigma) ######

sink("quad_inf.binom.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
      for(i in 1:N.obs)
      {trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
      trait[i] ~ dbin(ifelse(trait.mu[i] > 0.99, 0.99, trait.mu[i]), 1)
    }

    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- ifelse(-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) > 1, 0.99, 
    -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]))
    }
    
    } # close model
    ",fill=T)
sink()

###### 7. Fit using informative priors #####
# Load hyperparameters to bypass above code of generating them 
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

###### 7a. Hopland ######

# Set data 
data <- data.pLS.HOP
hypers <- pLS.HOP.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.HOP.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.HOP.out.inf$BUGSoutput$summary[1:5,]
#mcmcplot(pLS.HOP.out.inf)

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.HOP, ylab = "pLS for Hopland", xlab = "Temperature", pch = 1)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7b. Marin35 ######

# Set data 
data <- data.pLS.MAR1
hypers <- pLS.MAR1.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.MAR1.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR1.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.MAR1, ylab = "pLS for Marin35", xlab = "Temperature", pch = 1)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7c. Marin29 ######

# Set data 
data <- data.pLS.MAR2
hypers <- pLS.MAR2.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.MAR2.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.MAR2.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.MAR2, ylab = "pLS for Marin29", xlab = "Temperature", pch = 1)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7d. Wawona ######

# Set data 
data <- data.pLS.WAW
hypers <- pLS.WAW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.WAW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.WAW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.WAW, ylab = "pLS for Wawona", xlab = "Temperature", pch = 1)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7e. Eugene ######

# Set data 
data <- data.pLS.EUG
hypers <- pLS.EUG.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.EUG.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.EUG.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.EUG, ylab = "pLS for Eugene", xlab = "Temperature", pch = 1)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7f. Placer ######

# Set data 
data <- data.pLS.PLA
hypers <- pLS.PLA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.PLA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PLA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.PLA, ylab = "pLS for Placer", xlab = "Temperature", pch = 1)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7g. Santa Barbara ######

# Set data 
data <- data.pLS.SB
hypers <- pLS.SB.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.SB.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.SB.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.SB, ylab = "pLS for Santa Barbara", xlab = "Temperature", pch = 1)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###### 7h. Jasper Ridge ######

# Set data 
data <- data.pLS.JRA
hypers <- pLS.JRA.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.JRA.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.JRA.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.JRA, ylab = "pLS for Jasper Ridge", xlab = "Temperature", pch = 1)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7i. Paso Robles ######

# Set data 
data <- data.pLS.PAR
hypers <- pLS.PAR.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.PAR.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.PAR.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.PAR, ylab = "pLS for Paso Robles", xlab = "Temperature", pch = 1)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###### 7j. Powder canyon ######

# Set data 
data <- data.pLS.POW
hypers <- pLS.POW.prior.gamma.fits * 0.1

# Organize Data for JAGS
trait <- data$Larval.Survival
N.obs <- length(trait)
temp <- data$Temp.Treatment

# Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

# Run JAGS
pLS.POW.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.binom.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

# Examine Output
pLS.POW.out.inf$BUGSoutput$summary[1:5,]

# Plot data + fit
plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.POW, ylab = "pLS for Powder Canyon", xlab = "Temperature", pch = 1)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




###### For internal checking: 10 panel plot of fits using informative priors #####
plot.new()
par(mfrow = c(5,2))
par(mar = c(1, 2, 1, 1))

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.HOP, ylab = "pLS for Hop", xlab = "Temperature")
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.MAR1, ylab = "pLS for MAR1", xlab = "Temperature")
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.MAR2 , ylab = "pLS for MAR2 prior", xlab = "Temperature")
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.WAW , ylab = "pLS for WAW prior", xlab = "Temperature")
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.EUG , ylab = "pLS for EUG prior", xlab = "Temperature")
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.PLA , ylab = "pLS for PLA prior", xlab = "Temperature")
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.SB , ylab = "pLS for SB prior", xlab = "Temperature")
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.JRA , ylab = "pLS for JRA prior", xlab = "Temperature")
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.PAR , ylab = "pLS for PAR prior", xlab = "Temperature")
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(Larval.Survival ~ Temp.Treatment, xlim = c(5, 45), ylim = c(0,1), data = data.pLS.POW , ylab = "pLS for POW prior", xlab = "Temperature")
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
###### 10 panel plots of fits with raw data as means & sderrors ######
plot.new()

MeanSd = function(x) {
  aggxmean = aggregate(x$Larval.Survival ~ x$Temp.Treatment, FUN = mean, na.rm = T)
  aggxsd = aggregate(x$Larval.Survival ~ x$Temp.Treatment, FUN = sd, na.rm = T) 
  aggxss = aggregate(x$Larval.Survival ~ x$Temp.Treatment, FUN = length)
  aggxstderror = aggxsd[,2] / (sqrt(aggxss[,2]))
  df = cbind(aggxmean, aggxstderror)
  colnames(df) = c("Temp.Treatment", "mean", "stderror")
  return(df)
}

# Output to pdf
pdf(file = "~/Documents/Current_Projects/LifeHistoryTraitExp/Analysis_TraitFits/BayesianFitPlots/pLS.pdf",   
    width = 6, # The width of the plot in inches
    height = 8) # The height of the plot in inches

par(mfrow = c(5,2))
par(mar = c(2.5, 2, 1.5, 1))

pLS.Hop.MeansSds = MeanSd(data.pLS.HOP)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.Hop.MeansSds, ylab = "pLS for Hop", xlab = "Temperature", pch= 16)
arrows(pLS.Hop.MeansSds$Temp.Treatment, pLS.Hop.MeansSds$mean-pLS.Hop.MeansSds$stderror, pLS.Hop.MeansSds$Temp.Treatment, pLS.Hop.MeansSds$mean + pLS.Hop.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("HOP", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.HOP.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pLS.MAR1.MeansSds = MeanSd(data.pLS.MAR1)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.MAR1.MeansSds, ylab = "pLS for MAR1", xlab = "Temperature", pch= 16)
arrows(pLS.MAR1.MeansSds$Temp.Treatment, pLS.MAR1.MeansSds$mean-pLS.MAR1.MeansSds$stderror, pLS.MAR1.MeansSds$Temp.Treatment, pLS.MAR1.MeansSds$mean + pLS.MAR1.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR1", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.MAR1.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


pLS.MAR2.MeansSds = MeanSd(data.pLS.MAR2)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.MAR2.MeansSds, ylab = "pLS for MAR2", xlab = "Temperature", pch= 16)
arrows(pLS.MAR2.MeansSds$Temp.Treatment, pLS.MAR2.MeansSds$mean-pLS.MAR2.MeansSds$stderror, pLS.MAR2.MeansSds$Temp.Treatment, pLS.MAR2.MeansSds$mean + pLS.MAR2.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("MAR2", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.MAR2.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.MAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


pLS.WAW.MeansSds = MeanSd(data.pLS.WAW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.WAW.MeansSds, ylab = "pLS for WAW", xlab = "Temperature", pch= 16)
arrows(pLS.WAW.MeansSds$Temp.Treatment, pLS.WAW.MeansSds$mean-pLS.WAW.MeansSds$stderror, pLS.WAW.MeansSds$Temp.Treatment, pLS.WAW.MeansSds$mean + pLS.WAW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("WAW", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.WAW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


pLS.EUG.MeansSds = MeanSd(data.pLS.EUG)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.EUG.MeansSds, ylab = "pLS for EUG", xlab = "Temperature", pch= 16)
arrows(pLS.EUG.MeansSds$Temp.Treatment, pLS.EUG.MeansSds$mean-pLS.EUG.MeansSds$stderror, pLS.EUG.MeansSds$Temp.Treatment, pLS.EUG.MeansSds$mean + pLS.EUG.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("EUG", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.EUG.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


pLS.PLA.MeansSds = MeanSd(data.pLS.PLA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.PLA.MeansSds, ylab = "pLS for PLA", xlab = "Temperature", pch= 16)
arrows(pLS.PLA.MeansSds$Temp.Treatment, pLS.PLA.MeansSds$mean-pLS.PLA.MeansSds$stderror, pLS.PLA.MeansSds$Temp.Treatment, pLS.PLA.MeansSds$mean + pLS.PLA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PLA", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.PLA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


pLS.SB.MeansSds = MeanSd(data.pLS.SB)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.SB.MeansSds, ylab = "pLS for SB", xlab = "Temperature", pch= 16)
arrows(pLS.SB.MeansSds$Temp.Treatment, pLS.SB.MeansSds$mean-pLS.SB.MeansSds$stderror, pLS.SB.MeansSds$Temp.Treatment, pLS.SB.MeansSds$mean + pLS.SB.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("SB", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.SB.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


pLS.JRA.MeansSds = MeanSd(data.pLS.JRA)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.JRA.MeansSds, ylab = "pLS for JRA", xlab = "Temperature", pch= 16)
arrows(pLS.JRA.MeansSds$Temp.Treatment, pLS.JRA.MeansSds$mean-pLS.JRA.MeansSds$stderror, pLS.JRA.MeansSds$Temp.Treatment, pLS.JRA.MeansSds$mean + pLS.JRA.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("JRA", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.JRA.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.JRa.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pLS.PAR.MeansSds = MeanSd(data.pLS.PAR)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.PAR.MeansSds, ylab = "pLS for PAR", xlab = "Temperature", pch= 16)
arrows(pLS.PAR.MeansSds$Temp.Treatment, pLS.PAR.MeansSds$mean-pLS.PAR.MeansSds$stderror, pLS.PAR.MeansSds$Temp.Treatment, pLS.PAR.MeansSds$mean + pLS.PAR.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("PAR", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.PAR.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")

pLS.POW.MeansSds = MeanSd(data.pLS.POW)
plot(mean ~ Temp.Treatment, xlim = c(0, 40), ylim = c(0,1.1), data = pLS.POW.MeansSds, ylab = "pLS for POW", xlab = "Temperature", pch= 16)
arrows(pLS.POW.MeansSds$Temp.Treatment, pLS.POW.MeansSds$mean-pLS.POW.MeansSds$stderror, pLS.POW.MeansSds$Temp.Treatment, pLS.POW.MeansSds$mean + pLS.POW.MeansSds$stderror, length=0.05, angle=90, code=3)
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
text("POW", x = 36.5, y = 0.95, cex = 1.5)

lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "grey35")
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "grey35")
lines(pLS.POW.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "grey35")


dev.off()


###### 8. Get TPC parameters #####
###### 8a. Topt of pLS for each population ######

Topt = function(x) {
  Matrix = x[["BUGSoutput"]][["sims.matrix"]]
  ToptVec = rep(NA, nrow(Matrix))
  for (i in 1:nrow(Matrix))
  {ToptVec[i] = Temp.xs[which.max(Matrix[i,6:86])]}
  return(c(mean(ToptVec), quantile(ToptVec, c(0.025, 0.975))))
}

###### 8b. Pmax of pLS for each population ######

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

pLSdf_inf <- rbind.data.frame(ParamCompile(pLS.HOP.out.inf), ParamCompile(pLS.MAR1.out.inf),
                              ParamCompile(pLS.MAR2.out.inf), ParamCompile(pLS.WAW.out.inf),
                              ParamCompile(pLS.EUG.out.inf), ParamCompile(pLS.PLA.out.inf),
                              ParamCompile(pLS.SB.out.inf), ParamCompile(pLS.JRA.out.inf),
                              ParamCompile(pLS.PAR.out.inf), ParamCompile(pLS.POW.out.inf))
pLSdf_inf$Population = c(rep(c("HOP", "MAR1", "MAR2", "WAW", "EUG", "PLA", "SB", "JRA", "PAR", "POW"), each = 5))
save(pLSdf_inf, file = "pLS_meansd_inf.Rsave")

###### 9. Plotting  #####

load("pLS_meansd_inf.Rsave")

library(ggplot2)
plot.new()
par(mfrow = c(1,1))
par(mar = c(2.5, 2.5, 2.5, 2.5))

###### 9a. Plot of CTmax, CTmin, Topt for each population #####
# order by latitude of collection
LatOrder = rev(c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW"))
LatColors = rep(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#abd9e9", "#74add1", "#4575b4", "#313695"), each = 3)

pLSdf_inf$Population = factor(pLSdf_inf$Population, levels = LatOrder)
pLSdf_inf = pLSdf_inf[order(pLSdf_inf$Population),]
pLSdf_inf = pLSdf_inf[pLSdf_inf$param != "Tbreadth",]
pLSdf_inf = pLSdf_inf[pLSdf_inf$param != "Pmax",]

plotpLSinf = ggplot(pLSdf_inf, aes(x=Population, y=mean, col = LatColors)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(position = "top") +
  geom_point(stat="identity", col=LatColors, size = 2, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = LatColors, size = 0.7,
                position=position_dodge(.9), width = 0.5) + 
  ggtitle("Larval survival probability") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = "Temperature (\u00B0C)") + 
  theme(axis.text=element_text(size=14), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA))

###### 9b. Plot of Pmax for each population #####

load("pLS_meansd_inf.Rsave")
PopOrder = c("EUG", "WAW", "PLA", "HOP", "JRA", "MAR2", "MAR1", "POW", "PAR", "SB")
pLSdf_inf$Population = factor(pLSdf_inf$Population, levels = PopOrder)
pLSdf_inf = pLSdf_inf[order(pLSdf_inf$Population),]
pLSdf_inf = pLSdf_inf[pLSdf_inf$param == "Pmax",]

colors4 = rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1", "#4575b4","#313695"))

plotpLSinf_Pmax = ggplot(pLSdf_inf, aes(x=Population, y=mean, col = colors4)) + 
  scale_y_continuous(limits = c(0.7, 1)) +
  geom_point(stat="identity", col=colors4, 
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), col = colors4,
                position=position_dodge(.9), width = 0.3) + 
  ggtitle("pLS (1/day)") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = "Population", y = "Peak trait performance") + 
  theme(axis.text=element_text(size=15), 
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none") 


###### 9c. Plot pLS curve for all populations #####

load(file = "pLS_meansd_inf.Rsave")

LatOrder = c("EUG", "HOP", "PLA", "MAR2", "MAR1", "JRA", "WAW", "PAR", "SB", "POW")
LatColors = rep(rev(c("#ab041b", "#ec3c30", "#f46d43", "#fdae61", "#fed439ff", "#7cae00", 
                      "#abd9e9", "#74add1", "#4575b4", "#313695")), each = 3)

EUGdata = pLS.EUG.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
WAWdata = pLS.WAW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
HOPdata = pLS.HOP.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PLAdata = pLS.PLA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
SBdata = pLS.SB.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR2data = pLS.MAR2.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
MAR1data = pLS.MAR1.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
PARdata = pLS.PAR.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
JRAdata = pLS.JRA.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]
POWdata = pLS.POW.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]

# Plot data + fit
par(mar = (c(4.5, 4.5, 2.5, 1)))

plot(Larval.Survival ~ Temp.Treatment, 
     xlim = c(5, 40), ylim = c(0,1), data = data.pLS.HOP, type = "n", bty = "n",
     ylab = "probability", xlab = "Temperature (\u00B0C)", pch = 1,
     main = "Probability larval survival", cex.main = 1.7, cex.lab = 1.4, cex.axis = 1.2)
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
