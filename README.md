# Management Strategy Evaluation (MSE) for estimation bias

## Description
mse4bias is a management strategy evaluation (MSE) framework to evaluate management implications of persistent estimation bias in stock assessment and re-optimize harvest control rules (HCR). The framework was originally developed (2018) using the [Fisheries Library in R (FLR) mse package](https://github.com/flr/mse) for North Sea saithe (Pollachius virens) in Subareas 4, 6 and Division 3.a (North Sea, Rockall and West of Scotland, Skagerrak and Kattegat) as part of the [Workshop on North Sea stocks Management Strategy Evaluation (WKNSMSE)](https://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/Fisheries%20Resources%20Steering%20Group/2019/WKNSMSE/ICES%20WKNSMSE%20Report%202019.pdf).

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
• run_mse.R is for specifying and running MSE scenarios (HCR parameters, TAC contraint, and banking & borrowing, see [ICES 2019](https://ices-library.figshare.com/articles/report/Workshop_on_North_Sea_Stocks_Management_Strategy_Evaluation_WKNSMSE_/18621668)) and is also called from a job submission script to run on a high-performance computing system  
```r
### set HCR option: A, B, C
if (exists("HCRoption")) {
  input$ctrl.mp$ctrl.hcr@args$option <- switch(HCRoption, 
                                               "1" = "A", 
                                               "2" = "B", 
                                               "3" = "C",
                                               "4" = "A",
                                               "5" = "B",
                                               "6" = "C",
                                               "7" = "A")
  cat(paste0("\nSetting custom HCR option: HCRoption = ", HCRoption, 
             " => HCR ", input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
} else {
  cat(paste0("\nUsing default HCR option: HCR ", 
             input$ctrl.mp$ctrl.hcr@args$option, "\n\n"))
  HCRoption <- 0
}

### After combinations run  -- this section is turned on without running the above grid
### Btrigger and Ftrgt numbers must be updated with any changes 
if (HCRoption %in% 1:7) {
  comb_max <- switch(HCRoption, 
                     "1" = c(250000, 0.35), 
                     "2" = c(200000, 0.39), 
                     "3" = c(250000, 0.35), 
                     "4" = c(210000, 0.41),
                     "5" = c(220000, 0.39),
                     "6" = c(230000, 0.36),
                     "7" = c(230000, 0.36))
  hcr_vals <- expand.grid(Ftrgt = c(comb_max[2], 0.34, 0.33, 0.32, 0.31,
                                    0.36, 0.37, 0.38, 0.39, 0.40),    
                                Btrigger = comb_max[1])
}

### implement
if (exists("HCR_comb")) {
  
  ### set Btrigger
  Btrigger <- hcr_vals[HCR_comb, "Btrigger"]
  input$ctrl.mp$ctrl.phcr@args$Btrigger <- Btrigger
  input$ctrl.mp$ctrl.is@args$hcrpars$Btrigger <- Btrigger
  
  ### set Ftrgt
  Ftrgt <- hcr_vals[HCR_comb, "Ftrgt"]
  input$ctrl.mp$ctrl.phcr@args$Ftrgt <- Ftrgt
  input$ctrl.mp$ctrl.is@args$hcrpars$Ftrgt <- Ftrgt
  cat(paste0("\nSetting custom Btrigger/Ftrgt values.\n",
             "Using HCR_comb = ", HCR_comb, "\n",
             "Ftrgt = ", Ftrgt, "\n",
             "Btrigger = ", Btrigger, "\n\n"))
} else {
  cat(paste0("\nUsing default Btrigger/Ftrgt values.\n",
             "Ftrgt = ", input$ctrl.mp$ctrl.phcr@args$Ftrgt, "\n",
             "Btrigger = ", input$ctrl.mp$ctrl.phcr@args$Btrigger, "\n\n"))
}

### ------------------------------------------------------------------------ ###
### TAC constraint
input$ctrl.mp$ctrl.is@args$TAC_constraint <- FALSE

### check conditions
### either manually requested or as part of HCR options 4-7
if (exists("TAC_constraint")) {
  if (isTRUE(as.logical(TAC_constraint))) {
    input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
  }
}
if (HCRoption %in% 4:7) {
    input$ctrl.mp$ctrl.is@args$TAC_constraint <- TRUE
}
### implement
if (isTRUE(input$ctrl.mp$ctrl.is@args$TAC_constraint)) {
    if(HCRoption == 7){
      input$ctrl.mp$ctrl.is@args$lower <- 85
      input$ctrl.mp$ctrl.is@args$upper <- 115
      input$ctrl.mp$ctrl.is@args$Btrigger_cond <- TRUE
    } else {  
      input$ctrl.mp$ctrl.is@args$lower <- 80
      input$ctrl.mp$ctrl.is@args$upper <- 125
      input$ctrl.mp$ctrl.is@args$Btrigger_cond <- TRUE
    }
    cat(paste0("\nImplementing TAC constraint.\n\n"))
} else {
    cat(paste0("\nTAC constraint NOT implemented.\n\n"))
}

### ------------------------------------------------------------------------ ###
### banking & borrowing
input$ctrl.mp$ctrl.is@args$BB <- FALSE
input$iem <- NULL

### check conditions
### either manually requested or as part of HCR options 4-7
if (exists("BB")) {
  if (isTRUE(as.logical(BB))) {
    input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
    input$ctrl.mp$ctrl.is@args$BB <- TRUE
    input$ctrl.mp$ctrl.is@args$BB_check_hcr <- TRUE
    input$ctrl.mp$ctrl.is@args$BB_check_fc <- TRUE
    input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
  }
}
if (HCRoption %in% 4:7) {
  input$iem <- FLiem(method = iem_WKNSMSE, args = list(BB = TRUE))
  input$ctrl.mp$ctrl.is@args$BB <- TRUE
  input$ctrl.mp$ctrl.is@args$BB_rho <- c(-0.1, 0.1)
  input$ctrl.mp$ctrl.is@args$BB_check_hcr <- FALSE
  input$ctrl.mp$ctrl.is@args$BB_check_fc <- FALSE
  if (HCRoption %in% 4) {
    input$ctrl.mp$ctrl.is@args$BB_check_hcr <- TRUE
  } else if (HCRoption %in% 5:7) {
    input$ctrl.mp$ctrl.is@args$BB_check_fc <- TRUE
  }
}
if (!is.null(input$iem)) {
  cat(paste0("\nImplementing banking and borrowing.\n\n"))
} else {
  cat(paste0("\nBanking and borrowing NOT implemented.\n\n"))
}

```

• run_mse.sh is a job submission script for a high performance computing cluster (HPC) to call run_mse.R  
• analyse_mse.R is for analyzing the MSE results  


#### Example output of HCR (re)optimization for biased assessments
<img src="https://github.com/dgoto2/mse4bias/blob/main/heatmap_optimHCR.png?raw=true" width="800"> 

###### Optimization of the harvest control rule parameters (Ftarget and Btrigger) under varying levels (10% to 50%) of estimation bias in stock assessment (overestimation of stock abundance and underestimation of fishing mortality rate). Black boxes indicate maximum catches.


## References
Goto, D., J.A. Devine, I. Umar, S.H. Fischer, J.A.A. De Oliveira, D. Howell, E. Jardim, I. Mosqueira, K. Ono. 2022. [Shaping sustainable harvest boundaries for marine populations despite estimation bias](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.3923). Ecosphere. 13(2): e3923. 

ICES. 2020a. [The third Workshop on Guidelines for Management Strategy Evaluations (WKGMSE3)](https://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/Fisheries%20Resources%20Steering%20Group/2020/ICES%20WKGMSE3%202020.pdf). ICES Scientific Reports. 2:116. 112 pp. doi.org/10.17895/ices.pub.7627

ICES. 2020b. [Workshop on Catch Forecast from Biased Assessments (WKFORBIAS; outputs from 2019 meeting)](https://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/Fisheries%20Resources%20Steering%20Group/2020/WKFORBIAS_2019.pdf). ICES Scientific Reports. 2:28. 38 pp. http://doi.org/10.17895/ices.pub.5997 

ICES. 2019. [Workshop on North Sea stocks management strategy evaluation (WKNSMSE)](https://ices-library.figshare.com/articles/report/Workshop_on_North_Sea_Stocks_Management_Strategy_Evaluation_WKNSMSE_/18621668). ICES Scientific Reports. 1:12. 347 pp. doi.org/10.17895/ices.pub.5090


