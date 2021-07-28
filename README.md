# Management Strategy Evaluation (MSE) for assessment bias using FLR

## Description
mse4bias is a management strategy evaluation (MSE) framework using the Fisheries Library in R (FLR) mse package to evaluate management implications of persistent bias in stock assessment. The framework was originally developed for North Sea saithe (Pollachius virens) in Subareas 4, 6 and Division 3.a (North Sea, Rockall and West of Scotland, Skagerrak and Kattegat) as part of the [Workshop on North Sea stocks Management Strategy Evaluation (WKNSMSE)](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/viewer.html?pdfurl=https%3A%2F%2Fwww.ices.dk%2Fsites%2Fpub%2FPublication%2520Reports%2FExpert%2520Group%2520Report%2FFisheries%2520Resources%2520Steering%2520Group%2F2019%2FWKNSMSE%2FICES%2520WKNSMSE%2520Report%25202019.pdf&clen=48067636&chunk=true).

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

# To use SAM within FLR, the following R package is required:
devtools::install_github("shfischer/FLfse/FLfse", ref = "c561f5bf28cbad0f711ef53a49bde7e9868dc257")

```

## Core scripts to run simulations
• OM.R creates the baseline operating model (OM)  
• a4a_mse_WKNSMSE_funs.R contains a collection of functions and methods used for creating the OM and for running the MSE  
• run_mse.R is for running MSE scenarios and is called from a job submission script  
• run_google.sh is a job submission script used on a high performance computing cluster and call run_mse.R  
• run_mse_analyse.R is for analyzing the MSE results  
