
# Fitting thermal performance curves #

## Overview ##

- Following pipeline from Shocket et al. 2020: https://github.com/mshocket/Six-Viruses-Temp
using R2jags package
- Using uniform priors. Then using leave-one-out approach (using all other populations)
- Fit either symmetric, quadratic or asymmetric, Brière functions to the temperature-trait data.
(Visually inspect data for that trait to determine which to use.)
Previously:

Biting rate (inverse of gonotrophic cycle length (i.e., time between bloodfeed and egg-laying): Briere
Eggs per gonotrophic cycle: Briere
Probability egg to adult survival: Quadratic
Development rate (1/days): Briere
Adult lifespan: Quadratic

Briere: cT(T–T0)sqrt(Tm−T)
Quadratic: -c(T–T0)(T–Tm))  
where c = positive rate constant, T0 = Tmin, Tm = Tmax

## Currently working on ##

Fitting TPCs for pupal dev rate for each population   

## Scripts & Data files ##
In LifeHistoryTraitExp > Analysis_TraitFits
- TraitFits_LDR.R
- TraitFits_PDR.R
- LifeHistoryTraitExp_DataPreProcessing.R
- LifeHistoryExp_TempData_PreProcessing.csv
- LifeHistoryExp_TempData.csv

## References ##

- https://elifesciences.org/articles/58511
- https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568#pntd.0005568.ref052


