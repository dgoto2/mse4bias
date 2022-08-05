### ------------------------------------------------------------------------ ###
### R script to run MSE on HPC ####
### ------------------------------------------------------------------------ ###
### This is designed to be called by a job submission script

### ------------------------------------------------------------------------ ###
### load arguments from job script ####
### ------------------------------------------------------------------------ ###
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))

  ## alternative OMs with different natural mortality rates (M)
  if( !exists("m_criteria") ) {
    # Assuming base
    m_criteria <- ""
  } else {
    
    # Remove zero
    m_criteria <- gsub("[.]", "", m_criteria)
    if(m_criteria %in% c("01", "02", "03")) {
      if(m_criteria == "01") { m_criteria <- "M01_" }
      else if(m_criteria == "02") { m_criteria <- "" }
      else if(m_criteria == "03") { m_criteria <- "M03_" }
    } else {
      stop("m_criteria must be 0.1, 0.2 or 0.3")
    }
  }
} else {
  stop("no argument passed to R")
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###
### load packages
library(FLfse)
library(stockassessment)
library(ggplotFL)
library(FLAssess)
library(mse)
library(FLash)
library(tidyr)
library(dplyr)

### load additional functions
source("flr_mse_WKNSMSE_funs.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###
if (exists("par_env") && exists("n_workers") && exists("nblocks") && nblocks > 1) {

  ### par_env=1 -> MPI (Rmpi, DoMPI)
  ### par_env=2 -> DoParallel
  if (par_env == 1) {
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)
    cl_length <- cl$workerCount
    pmethod <- "MPI"
    n_workers <- cl_length
  } else if (par_env == 2) {
    library(doParallel)
    cl <- makeCluster(n_workers, outfile="")
    registerDoParallel(cl)
    cl_length <- length(cl)
    pmethod <- "MultiCore"
  }
  print(paste0("Starting ", n_workers, " parallel workers with ", pmethod)) 

  ### load packages and functions into workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    #devtools::load_all("../mse/")
    library(mse)
    library(FLash)
    library(FLfse)
    library(stockassessment)
    library(foreach)
    library(doRNG)
    source("flr_mse_WKNSMSE_funs.R")
  }

  ### set random seed for reproducibility
  library(doRNG)
  registerDoRNG(123)
} else {
  nblocks <- 1
}

### ------------------------------------------------------------------------ ###
### load data for MSE ####
### ------------------------------------------------------------------------ ###
### data path
path_data <- paste0("input/pok/", iters, "_", years, "/", m_criteria)

### load input objects
input <- readRDS(paste0(path_data, "base_run.rds"))

### modify input for running in parallel
input$genArgs$nblocks <- nblocks

### ------------------------------------------------------------------------ ###
### set up HCR & options ####
### ------------------------------------------------------------------------ ###
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

### ------------------------------------------------------------------------ ###
### specify the levels of bias in assessment
input$ctrl.mp$ctrl.est@args$prop_biasN <- 1.0+prop_biasN ## bias in N-at-age from sam fit
input$ctrl.mp$ctrl.est@args$prop_biasF <- 1.0+prop_biasF ## bias in F-at-age from sam fit

### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###
### Use A1 for Option 7
if(HCRoption == 7) {
  HCRopt <- "A1"
} else {
  HCRopt <- input$ctrl.mp$ctrl.hcr@args$option
}

### Create output file name
path_out <- paste0("output/runs/pok/", iters, "_", years)
dir.create(path = path_out, recursive = TRUE)
file_out <- paste0("HCR-", HCRopt,
                   "_Ftrgt-", input$ctrl.mp$ctrl.phcr@args$Ftrgt,
                   "_Btrigger-", input$ctrl.mp$ctrl.phcr@args$Btrigger,
                   "_TACconstr-", input$ctrl.mp$ctrl.is@args$TAC_constraint,
                   "_BB-", input$ctrl.mp$ctrl.is@args$BB,
                   "_biasN-", input$ctrl.mp$ctrl.est@args$prop_biasN,
                   "_biasF-", input$ctrl.mp$ctrl.est@args$prop_biasF
            )
outFile <- paste0(path_out, "/", m_criteria, file_out, ".rds")

### Checking if output exists. If yes, quit
print(paste("Saving to: ", outFile))
if( file.exists(outFile) && !(exists("forceOverwrite") && forceOverwrite) ) {
  print("File exists. Not overwriting.")
  quit(save = "no")
}

### run MSE
res1 <- mp(om = input$om,
           oem = input$oem,
           iem = input$iem,
           ctrl.mp = input$ctrl.mp,
           genArgs = input$genArgs,
           tracking = input$tracking)

### save output
saveRDS(object = res1, file = outFile)

### ------------------------------------------------------------------------ ###
### terminate ####
### ------------------------------------------------------------------------ ###
### close R
### mpi.finalize() or mpi.quit() hang...
### -> kill R, the MPI processes stop afterwards
### clean up
if (exists("cl")) {
  if (par_env == 1 && exists("kill")) {
    system("bkill $LSB_JOBID")
  } else if (par_env == 2) {
    stopCluster(cl)
  }
}
quit(save = "no")
