
# Fitting thermal performance curves #

References: 
- https://elifesciences.org/articles/58511
- https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568#pntd.0005568.ref052

## Overview ##

Fit either symmetric, quadratic or asymmetric, Brière functions to the temperature-trait data.
Visually inspect data for that trait to determine which to use. Previously:

Biting rate (inverse of gonotrophic cycle length (i.e., time between bloodfeed and egg-laying): Briere
Eggs per gonotrophic cycle: Briere
Probability egg to adult survival: Quadratic
Development rate (1/days): Briere
Adult lifespan: Quadratic

Briere: cT(T–T0)sqrt(Tm−T)
Quadratic: -c(T–T0)(T–Tm))  
where c = positive rate constant, T0 = Tmin, Tm = Tmax

Using R2jags package (R interface for JAGS package, which does MCMC sampling)

