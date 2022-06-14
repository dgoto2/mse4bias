# Management Strategy Evaluation (MSE) for estimation bias

## Description
mse4bias is a management strategy evaluation (MSE) framework to evaluate management implications of persistent estimation bias in stock assessment. The framework was originally developed (2018) using the [Fisheries Library in R (FLR) mse package](https://github.com/flr/mse) for North Sea saithe (Pollachius virens) in Subareas 4, 6 and Division 3.a (North Sea, Rockall and West of Scotland, Skagerrak and Kattegat) as part of the [Workshop on North Sea stocks Management Strategy Evaluation (WKNSMSE)](https://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/Fisheries%20Resources%20Steering%20Group/2019/WKNSMSE/ICES%20WKNSMSE%20Report%202019.pdf).

This framework simulates population and harvest dynamics, surveys, assessments, and implementation of management strategies to explore trade-offs in achieving conservation-oriented (minimizing overexploitation risk) and harvest-oriented (maximizing yield) goals. The framework consists of submodels that simulate (1) true population and harvest dynamics at sea (operating model [OM]), from which observations through monitoring surveys and catch reporting (data generation) are made; and (2) management processes: assessments based on observations from the surveys and reported catch (using the State-space Assessment Model-SAM as estimation model (EM); https://github.com/fishfollower/SAM) and subsequent decisionmaking (management procedure, MP) based on the harvest control rule set for saithe ([ICES 2019]( https://www.ices.dk/sites/pub/Publication%20Reports/Advice/2019/2019/pok.27.3a46.pdf)).

<img src="https://github.com/dgoto2/mse4bias/blob/main/saithe.mse.png?raw=true" width="500"> 

###### (redrawn from https://github.com/ejardim; image credit: IAN Symbols, courtesy of the Integration and Application Network, University of Maryland Center for Environmental Science ([ian.umces.edu/symbols/](https://ian.umces.edu/media-library/symbols/)))

## Prerequisites
Install the following packages:
```r
install.packages(c("foreach", "doParallel", "dplyr", "tidyr", "data.table")) 

devtools::install_github(repo = "flr/FLCore", ref = "d55bc6570c0134c6bea6c3fc44be20378691e042")
devtools::install_github(repo = "flr/FLash", ref = "7c47560cf57627068259404bb553f2b644682726")
devtools::install_github(repo = "flr/FLBRP", ref = "5644cfccefb0ec3965b1d028090bbf75b1e59da2")
devtools::install_github(repo = "flr/FLAssess", ref = "f1e5acb98c106bcdfdc81034f1583f76bb485514")
devtools::install_github(repo = "flr/ggplotFL", ref = "e9e0d74e872815c1df3f172522da35ade5c70638")
devtools::install_github(repo = "flr/mse", ref = "e39ddd75cdb2bb693601e31428404d48ea810308")

# This MSE uses SAM as an assessment model and requires the biomassindex branch of the stockassessment R package:
install.packages("TMB") 
devtools::install_github("fishfollower/SAM/stockassessment", ref="biomassindex")

# To use State-space Assessment Model (SAM) within FLR, the following R package is required:
devtools::install_github("shfischer/FLfse/FLfse", ref = "c561f5bf28cbad0f711ef53a49bde7e9868dc257")

```

## Core scripts to run simulations
• OM.R creates the baseline operating model (OM)  
• flr_mse_WKNSMSE_funs.R contains a collection of functions and methods used for creating the OM and for running the MSE  
• run_mse.R is for running MSE scenarios (assessment bias scenarios can be set in prop_biasN and prop_biasF) and is called from a job submission script  
• run_mse.sh is a job submission script for a high performance computing cluster (HPC) to call run_mse.R  
• analyse_mse.R is for analyzing the MSE results  


<img src="https://github.com/dgoto2/mse4bias/blob/main/heatmap_optimHCR.png?raw=true" width="800"> 

###### Optimization of the harvest control rule parameters (Ftarget and Btrigger) under varying levels (10% to 50%) of estimation bias in stock assessment (overestimation of stock abundance and underestimation of fishing mortality rate).


## Reference
Goto, D., J.A. Devine, I. Umar, S.H. Fischer, J.A.A. De Oliveira, D. Howell, E. Jardim, I. Mosqueira, K. Ono. 2022. [Shaping sustainable harvest boundaries for marine populations despite estimation bias](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.3923). Ecosphere. 13(2): e3923. 
