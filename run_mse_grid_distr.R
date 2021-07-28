### ------------------------------------------------------------------------ ###
### R script to a grid search using MSEs ####
### ------------------------------------------------------------------------ ###

## REQUIRED ARGUMENTS
# iterations
iters <- 1000

# projection period
years <- 21

# M (natural mortality) selection
m_criterias <- c(0.2)

# HCR options
#hcrs <- c(1)

# HCR combs
#hcr_combs <- c(1:25)

# parallel workers
#n_workers <- 80

## ------------------
Rargs <- ""
extraArgsPar <- ""
if (exists("n_workers"))
  extraArgsPar <- paste0(" par_env=2 n_workers=", n_workers, " nblocks=", n_workers)
for( m_criteria in m_criterias ) {
  extraArgsAll <- ""
  extraArgsAll <- paste0(" m_criteria=", m_criteria)
  for(hcr in hcrs) {
    for(hcrcomb in hcr_combs)
      # Run scenario
      system(paste0("Rscript ", Rargs, " run_mse.R iters=", iters," years=", years, " HCRoption=", 
                    hcr, " HCR_comb=", hcrcomb," TAC_constraint=0 BB=0 ", extraArgsAll, extraArgsPar))
  }
}
