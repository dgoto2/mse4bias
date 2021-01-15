#### BASE RUN TEST

## REQUIRED ARGUMENTS

# Iterations
iters <- 1000

# Years
years <- 21

## OPTIONAL ARGUMENTS

# M selection
m_criteria <- 0.1

# parallel workers
n_workers <- 4

## ------------------

Rargs <- ""

extraArgsAll <- ""
extraArgsPar <- ""

if (exists("m_criteria"))
  extraArgsAll <- paste0(" m_criteria=", m_criteria)

if (exists("n_workers"))
  extraArgsPar <- paste0(" par_env=2 n_workers=", n_workers, " nblocks=", n_workers)

# Prepare objects
system(paste0("Rscript ", Rargs, " OM.R iters=", iters," years=", years, extraArgsAll))

# Switch to run mse with base_run after OM
if(exists("no_run") && !no_run) {
  # Run base MSE
  system(paste0("Rscript ", Rargs, " run_mse.R iters=", iters," years=", years, extraArgsAll, extraArgsPar))
}
