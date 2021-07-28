# Management Strategy Evaluation (MSE) for estimation bias

## Description
mse4bias is a management strategy evaluation (MSE) framework to evaluate management implications of persistent estimation bias in stock assessment. The framework was originally developed using the [Fisheries Library in R (FLR) mse package](https://github.com/flr/mse) for North Sea saithe (Pollachius virens) in Subareas 4, 6 and Division 3.a (North Sea, Rockall and West of Scotland, Skagerrak and Kattegat) as part of the [Workshop on North Sea stocks Management Strategy Evaluation (WKNSMSE)](https://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/Fisheries%20Resources%20Steering%20Group/2019/WKNSMSE/ICES%20WKNSMSE%20Report%202019.pdf).

## Prerequisites
Install the following packages:
```r
install.packages(c("foreach", "DoParallel", "dplyr", "tidyr", "data.table")) 

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
• a4a_mse_WKNSMSE_funs.R contains a collection of functions and methods used for creating the OM and for running the MSE  
• run_mse.R is for running MSE scenarios and is called from a job submission script  
• run_google.sh is a job submission script used on a high performance computing cluster and call run_mse.R  
• run_mse_analyse.R is for analyzing the MSE results  
