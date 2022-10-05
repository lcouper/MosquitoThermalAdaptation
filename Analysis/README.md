
# Fitting thermal performance curves #

## Overview ##

Following pipeline from Shocket et al. 2020: https://github.com/mshocket/Six-Viruses-Temp
using R2jags package

Overview of trait-fitting procedure:


1. Generate data informed priors 
    1. For each population, use trait data from all other populations (‘leave-one-out’)
    2. Fit either quadratic or Briere curves (depending on trait) to the temperature-trait data using uninformative priors (i.e., uniform, but constrained to be over a realistic temperature range)
    3. Fit gamma distributions to the TPC parameters generated for each population with the leave-one-out approach
2. Use these data-informed priors to fit TPCs for each population 
3. Obtain estimates of CTmin, CTmax 


Previously:
* Briere: Biting rate (inverse of gonotrophic cycle length (i.e., time between bloodfeed and egg-laying), eggs per gonotrophic cycle, development rate
* Quadratic: Probability egg to adult survival, adult lifespan
* Equations (c = positive rate constant, T0 = Tmin, Tm = Tmax):
  * Briere: cT(T–T0)sqrt(Tm−T)        Quadratic: -c(T–T0)(T–Tm))  

## Currently working on ##



## Scripts & Data files ##
In LifeHistoryTraitExp > Analysis_TraitFits
* TraitFits_LDR.R
* TraitFits_PDR.R
* TraitFits_pLS.R 
* TraitFits_pPS.R
* TraitFits_AL.R 
* TraitFits_EFD.R
* TraitFits_BR.R
* LifeHistoryTraitExp_DataPreProcessing.R
* LifeHistoryTraitExp_Data.csv
# ReproExp_Data.csv

## References ##

- https://elifesciences.org/articles/58511
- https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568#pntd.0005568.ref052


