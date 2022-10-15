# LifeHistoryTraitExp
Protocols, datasheets, and code for the *Aedes sierrensis* life history trait experiment 

## Links ##
Results slides: https://docs.google.com/presentation/d/1t0al5WpG7hgLubEsBmxENzvtfGpIB8K2Z7HGYpsTtd0/edit#slide=id.p   
Protocols & Descriptions: https://drive.google.com/drive/u/1/folders/1Mo3wjYjz4vL88z8FdgrOB-oj7w-4425v    
Datasheet: https://docs.google.com/spreadsheets/d/1oZZGo85oaEHsGGcL0MKQG68ncqCRL5oygw-cUxur2E0/edit#gid=336945607   
Background notes: https://docs.google.com/document/d/1AiJJt0yAJZJsXEDOziYk3WfPh45Odg9c9qZhwX3piQM/edit  
Wing Images: https://drive.google.com/drive/u/1/folders/1NtzC93sPJOMRPK9BPH_nkGAijmumqk0B  

## Projects ##

#### 1. Investigating evidence for mosquito thermal adaptation
Working doc: https://docs.google.com/document/d/1xFx6vmBDdw5-maTafi-_qo3OfKLJdVki4uI07BEJNdA/edit  

1. To what extent is there variation in mosquito thermal tolerance?   
1a.  Is the variation primarily in thermal *limits* or trait performance?  (evidence: for composite fitness measure, calculate variance in peak performance and in thermal limits)  
1b. Is the variation within-, between-population, or both?   (do above for both within (for dev rate, fecundity, and life span traits) and between-populations)   
1c.  Is the variation genetically-based? (evidence: results from sequencing on % genetic contribution to fitness)   
 
 
2. With what aspects of temperature/climate is this variation associated? (evidence: correlations between thermal limits and bioclim variables)  
3. Importance of each life history trait (struggling with the right question here)  
3a. What life history traits vary the most? (evidence: for each trait, calcualte between-population variance?)  
3b. What life history traits are most limiting to fitness?  (evidence: sensitivity analysis of fitness calculations to trait inputs as in Mordecai et al 2017 (see supplementals)? or just look at which traits have lowest TPC parameters)  
3c. What life history traits contribute most to local adaptation? (evidence: which show strongest correlation with climatic variable?)
3d. Potential add-on idea = measure correlational selection as in Gilbert 2017 (https://royalsocietypublishing.org/doi/full/10.1098/rspb.2017.0536#d3e1705)  
Model fitness ~ trait1+ trait2 + trait3 + trait1*trait2   
allows you to estimate correlations between traits

#### 2. Intra-specific variation in modeling climate responses 
Working doc: https://docs.google.com/document/d/1kU_lWkhXjSKcnZGYVIlGcciu3q279f0bcKXwr3HfMgQ/edit   
1. How do projections of mosquito climate responses vary when using species-wide vs population-level life history trait responses?
2. How do they vary when using iButton vs weather station data vs remotely sensed data?

#### 3. Genetic varaiation and potential for thermal adaptaiton ####
Working doc: https://docs.google.com/document/d/1E3hGLk8f2fVaF-348dDgzJ6As3sVeZl5RvhXSWXEKEI/edit
1. How much within- and between- genetic variation exists between populations? 
2. How much of the variaiton in fitness can be explained by genetics vs the environment?




## Currently working on ##
- finishing thermal performance experiments (4C treatment ongoing)
- fitting TPCs to prelim data using Bayesian methods
- DNA extraction

## Upcoming ##
- calculate fitness for each population & temperature treatment (i.e., probability of surviving to day [x] adult stage x fecundity (proxy based on wing length))
- calculate between-population variation in upper thermal limits for traits 
- pull other temperature data for sites (e.g., max temp in Spring, max temp in wettest quarter)
- wing length measurements & DNA extractions
- Estimate sequencing cost. Decide whether to sequence all individuals or just those with wing length measurements (i.e., those with available fecundity proxies)

### Approximate sample numbers to extract (9/7/22) ###
Refers to adults with complete data (i.e., will be sequenced)
725 in total

| Population | HOP | MAR1 | MAR2 | EUG | WAW | JRA | SB | PAR | POW | PLA | 
| :-----: | :---: | :---: |:---: | :---: |:---: | :---: |:---: | :---: |:---: | :---: |
| Approx. # | 66 | 67 | 60 | 47 | 289 | 91 | 73 | 78 | 83 | 76 | 84 |


 
