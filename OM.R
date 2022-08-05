### ------------------------------------------------------------------------ ###
### R script to specify an operating model (OM) for FLR MSEs ####
### ------------------------------------------------------------------------ ###

### install the stockassessment R package if needed
# source("http://flr-project.org/R/instFLR.R")
#devtools::install_github("fishfollower/SAM/stockassessment", ref="biomassindex") 

### load required packages
library(FLa4a)
library(doParallel)
library(FLfse)
library(stockassessment)
library(ggplotFL)
library(FLAssess)
library(mse)
library(FLash)
library(tidyr)
library(dplyr)

### create plots and print to screen?
verbose <- FALSE

###-------------------------------------------------
### load the stock objects and create the config file
###-------------------------------------------------
source("create_stock_object.R")
source("flr_mse_WKNSMSE_funs.R")

### grab stock objects from stockassessment.org (this works only on newer versions of SAM code, o/w build from scratch)
sei <- fitfromweb(NS_saithe_2018_rerun)
pok_conf <- loadConf(sei$data, "./conf/model.cfg")

### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###
### number of replicates
n_iter <- 1000

### years of projection period
n_years <- 21

### last data year
yr_data <- 2017

### ------------------------------------------------------------------------ ###
### fit SAM ####
### ------------------------------------------------------------------------ ###
### use input data provided in FLfse
### recreates the WGNSSK2018 saithe assessment
fit <- FLR_SAM(stk = pok_stk, idx = pok_idx, conf = pok_conf, conf_full = TRUE) 
if (isTRUE(verbose)) {
  is(fit)
  fit
  plot(fit)
}

### re-estimate Blim (to be used when calculating performance stats)
round(min(ssbtable(fit)[,1]))

### extract model parameters and use them in the simulation as starting values
sam_initial <- sam_getpar(fit)
sam_initial$logScale <- numeric(0)

### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration
stk <- SAM2FLStock(object = fit, stk = pok_stk)
if (isTRUE(verbose)) summary(stk)

### set units
units(stk)[1:17] <- as.list(c(rep(c("t", "1000", "kg"), 4), "", "", "f", "", ""))
if (isTRUE(verbose)) plot(stk)

### save for later comparison
stk_orig <- stk

### ------------------------------------------------------------------------ ###
### add uncertainty ####
### ------------------------------------------------------------------------ ###
### use variance-covariance matrix 
### add iteration dimension
stk <- FLCore::propagate(stk, n_iter)
dim(stk)

### add uncertainty estimated by SAM as iterations
set.seed(1)
uncertainty <- SAM_uncertainty(fit = fit, n = n_iter, print_screen = FALSE)
### SAM_uncertainty gets the covariance matrix of all model parameters and creates 
### random values based on estimation and covar sim states
### output = stock.n = stock.n, harvest = harvest,survey_catchability = catchability, 
### catch_sd = catch_sd,survey_sd = survey_sd, proc_error = SdLogN, sam_initial$logFpar

### add noise to stock
stock.n(stk)[] <- uncertainty$stock.n
stock(stk)[] <- computeStock(stk)

### add noise to F
harvest(stk)[] <- uncertainty$harvest

### catch noise added later
if (isTRUE(verbose)) plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
if (isTRUE(verbose)) plot(stk)

### maximum observed F
max(fbar(stk)); max(fbar(stk_orig))
max(harvest(stk)); max(harvest(stk_orig))

### get estimated catch numbers
catch_n <- uncertainty$catch_n

### ------------------------------------------------------------------------ ###
### extend stock projection for MSE simulations ####
### ------------------------------------------------------------------------ ###
stk_stf2017 <- stf(window(stk, end = 2017), n_years + 1)
stk_stf <- stk_stf2017

### ------------------------------------------------------------------------ ###
### biological data for OM ####
### ------------------------------------------------------------------------ ###
### the approach used in eqsim for NS saithe: 
### resample weights, maturity, M from the last 10 years (2008-2017)
### set up an array with one resampled year for each projection year
### (including intermediate year) and replicate
### use the same resampled year for all biological parameters
set.seed(2)

### use last 10 data years to sample biological parameters
bio_sample_yrs <- 2008:2017
sel_sample_yrs <- 2013:2017

### get year position of sample years
sample_yrs_pos_bio <- which(dimnames(stk_stf)$year %in% bio_sample_yrs)
sample_yrs_pos_sel <- which(dimnames(stk_stf)$year %in% sel_sample_yrs)

### create samples for biological data (weights, etc.)
### the historical biological parameters are identical for all iterations
### and consequently do not need to be treated individually (but keep age structure)
### create vector with resampled years
bio_samples <- sample(x = sample_yrs_pos_bio, size = (n_years + 1) * n_iter, replace = TRUE)

### do the same for selectivity
sel_samples <- sample(x = sample_yrs_pos_sel, size = (n_years + 1) * n_iter, replace = TRUE)

### populate year parameters
bio_yrs <- which(dimnames(stk_stf)$year %in% 2018:dims(stk_stf)$maxyear)

### insert values
catch.wt(stk_stf)[, bio_yrs] <- c(catch.wt(stk)[, bio_samples,,,, 1])
stock.wt(stk_stf)[, bio_yrs] <- c(stock.wt(stk)[, bio_samples,,,, 1])
landings.wt(stk_stf)[, bio_yrs] <- c(landings.wt(stk)[, bio_samples,,,, 1])
discards.wt(stk_stf)[, bio_yrs] <- c(discards.wt(stk)[, bio_samples,,,, 1])
m(stk_stf)[, bio_yrs] <- c(m(stk)[, bio_samples,,,, 1])
mat(stk_stf)[, bio_yrs] <- c(mat(stk)[, bio_samples,,,, 1])

### use different samples for selectivity
harvest(stk_stf)[, bio_yrs] <- c(harvest(stk)[, sel_samples,,,, 1])
if (isTRUE(verbose)) plot(stk_stf)

### ------------------------------------------------------------------------ ###
### recruitment ####
### ------------------------------------------------------------------------ ###
### fit the hockey-stick model
### get residuals from smoothed residuals
### use only recruitment data from 1998 and later (ssb from 1995)
sr <- as.FLSR(window(stk_stf, start = 1995), model = "segreg")

### fit model individually to each iteration and suppress output to screen
suppressWarnings(. <- capture.output(sr <- fmle(sr, fixed=list(b=min(ssb(stk_stf), na.rm=T)))))

if (isTRUE(verbose)) {
  plot(sr)  

	### check breakpoints
	summary(params(sr))

	### plot model and data
	as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb)) %>%
	  mutate(age = NULL, year = ifelse(qname == "SSB", year + 3, year)) %>%
	  tidyr::spread(key = qname, value = data) %>%
	  ggplot() +
	  geom_point(aes(x = SSB, y = rec, group = iter), 
				 alpha = 0.5, colour = "darkgrey", shape = 3) +
	  geom_line(aes(x = SSB, y = fitted, group = iter)) +
	  theme_bw() + xlim(0, NA) + ylim(0, NA)

	### Check extent of autocorrelation
	# cod code showed start year was 1997 and acf was +1 year for start
	acf(window(stock.n(stk_orig)[1], start = 1998))
	# Not significant, so no need to account for it in this OM

	### Check method proposed for generating recruitment compares with past recruitment estimates
	test <- as.data.frame(FLQuants(fitted = sr@fitted, rec = sr@rec, SSB = sr@ssb))
	test <- mutate(test, age = NULL, year = ifelse(qname == "SSB", year + 3, year))
	test <- tidyr::spread(test, key = qname, value = data)
	test <- test[complete.cases(test),]
	test$res <- rep(NA, nrow(test))

	# Generate residuals for future recruitments
	foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", .errorhandling = "pass") %do% {
	  set.seed(iter_i^2)
			  
    ### get residuals for current iteration
    res_i <- c(FLCore::iter(residuals(sr), iter_i))
    res_i <- res_i[!is.na(res_i)]
    
    ### calculate kernel density of residuals
    density <- density(x = res_i)
    
    ### sample residuals
    mu <- sample(x = res_i, size = length(res_i), replace = TRUE)
    
    ### "smooth", i.e. sample from density distribution
    test$res[test$iter==iter_i] <- rnorm(n = length(res_i), mean = mu, sd = density$bw)
	}
	
	# Generate future recruits from past SSBs and generated residuals
	test$future <- test$fitted * exp(test$res)

	# 10 randomly selected iters for plotting
	# Use less if n_iter is low
	n_samp <- 10
	if(n_iter < n_samp) n_samp <- round(n_iter/2)
	i_samp <- sample(seq(dim(sr)[6]), n_samp, replace=FALSE)

	# plot past and future stock recruit pairs for selected iters
	ggplot(test[is.element(test$iter, i_samp),]) +
	  geom_point(aes(x = SSB, y = rec), alpha = 0.5, colour = "red", shape = 19) +
	  geom_point(aes(x = SSB, y = future), alpha = 0.5, colour = "black", shape = 19) +
	  geom_line(aes(x = SSB, y = fitted)) +
	  facet_wrap(~iter) +
	  theme_bw() + xlim(0, NA) + ylim(0, NA)

	# empirical cumulative distributions for the same iters
	ggplot(test[is.element(test$iter, i_samp),]) +
	  stat_ecdf(aes(rec), geom = "step", colour = "red") +
	  stat_ecdf(aes(future), geom = "step", colour = "black") +
	  facet_wrap(~iter) +
	  theme_bw() + xlim(0, NA) + ylim(0, NA)

	# Combine previous two plots over all iters
	i_samp <- sample(seq(dim(sr)[6]), n_iter, replace=FALSE)

	# stock-recruit pairs - all
	ggplot(test[is.element(test$iter, i_samp),]) +
	  geom_point(aes(x = SSB, y = rec), alpha = 0.5, colour = "red", shape = 19) +
	  geom_point(aes(x = SSB, y = future), alpha = 0.5, colour = "black", shape = 19) +
	  theme_bw() + xlim(0, NA) + ylim(0, NA)
	
	# empirical cumulative distribution -all 
	ggplot(test[is.element(test$iter, i_samp),]) +
	  stat_ecdf(aes(rec), geom = "step", colour = "red") +
	  stat_ecdf(aes(future), geom = "step", colour = "black") +
	  theme_bw() + xlim(0, NA) + ylim(0, NA)
	rm(test, i_samp)
}

############################################################ 
### generate residuals
### years with missing residuals
yrs_res <- dimnames(sr)$year[which(is.na(iterMeans(rec(sr))))]
## check this as there was a note that NULL is produced sometimes - if find NULL, use the commented out code

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
res_new <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", .errorhandling = "pass") %dopar% { 
  set.seed(iter_i)
                     
 ### get residuals for current iteration
 res_i <- c(FLCore::iter(residuals(sr), iter_i))
 res_i <- res_i[!is.na(res_i)]
 
 ### calculate kernel density of residuals
 density <- density(x = res_i)
 
 ### sample residuals
 mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
 
 ### "smooth", i.e. sample from density distribution
 res_new <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
 return(res_new) 
}
summary(exp(unlist(res_new)))

### insert into the stock object
residuals(sr)[, yrs_res] <- unlist(res_new)

### exponeniate residuals
residuals(sr) <- exp(residuals(sr))
sr_res <- residuals(sr)
if (isTRUE(verbose)) plot(sr_res)

### ------------------------------------------------------------------------ ###
### process noise ####
### ------------------------------------------------------------------------ ###
### create FLQuant with process noise
### this will be added to the values obtained from fwd() in the MSE

### create noise for process error
set.seed(3)
proc_res <- stock.n(stk_stf) %=% 0 ### template FLQuant
proc_res[] <- stats::rnorm(n = length(proc_res), mean = 0, sd = uncertainty$proc_error)

### the proc_res values are on a normal scale,
### exponentiate to get log-normal 
proc_res <- exp(proc_res)
### proc_res is a factor by which the numbers at age are multiplied

### for historical period, numbers already include process error from SAM
### -> remove deviation
proc_res[, dimnames(proc_res)$year <= 2017] <- 1

### remove deviation for first age class (recruits)
proc_res[1, ] <- 1

### this gets passed on to the projection module
fitted(sr) <- proc_res
if (isTRUE(verbose)) plot(proc_res)

### ------------------------------------------------------------------------ ###
### stf for 2018: assume catch advice is taken ####
### ------------------------------------------------------------------------ ###
c2018 <- 116008
ctrl <- fwdControl(data.frame(year = 2018, quantity = "catch", val = c2018))

### project forward for intermediate year (2018)
stk_int <- stk_stf
stk_int[] <- fwd(stk_stf, ctrl = ctrl, sr = sr, sr.residuals = sr_res,
                 sr.residuals.mult = TRUE, maxF = 5)[]

### add process noise
stock.n(stk_int) <- stock.n(stk_int) * proc_res
stock(stk_int)[] <- computeStock(stk_int)

### create a stock object
stk_fwd <- stk_stf

### insert values for 2018
stk_fwd[, ac(2018)] <- stk_int[, ac(2018)]

### insert stock number for 2019 to calculate SSB at beginning of 2019
stock.n(stk_fwd)[, ac(2019)] <- stock.n(stk_int)[, ac(2019)]
stock(stk_fwd)[, ac(2019)] <- computeStock(stk_fwd[, ac(2019)])

### ------------------------------------------------------------------------ ###
### biological data for OEM ####
### ------------------------------------------------------------------------ ###
### based on the OM
stk_oem <- stk_fwd

### projection years
proj_yrs <- 2018:range(stk_oem)[["maxyear"]]

### use means of sampled values for projection period
catch.wt(stk_oem)[, ac(proj_yrs)] <- yearMeans(catch.wt(stk_oem)[, ac(bio_sample_yrs)])
landings.wt(stk_oem)[, ac(proj_yrs)] <- yearMeans(landings.wt(stk_oem)[, ac(bio_sample_yrs)])
discards.wt(stk_oem)[, ac(proj_yrs)] <- yearMeans(discards.wt(stk_oem)[, ac(bio_sample_yrs)])
stock.wt(stk_oem)[, ac(proj_yrs)] <- yearMeans(stock.wt(stk_oem)[, ac(bio_sample_yrs)])
m(stk_oem)[, ac(proj_yrs)] <- ifelse(m_criteria == "", yearMeans(m(stk_oem)[, ac(bio_sample_yrs)]), 0.2)
mat(stk_oem)[, ac(proj_yrs)] <- yearMeans(mat(stk_oem)[, ac(bio_sample_yrs)])

### remove stock assessment results
stock.n(stk_oem)[] <- stock(stk_oem)[] <- harvest(stk_oem)[] <- NA

### ------------------------------------------------------------------------ ###
### indices ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template (included in FLfse)
idx <- pok_idx

### second index (FSB) is a mid-year survey
range(idx[[2]])[c("startf", "endf")] <- 0.5

### extend for the simulation period
idx <- window(idx, end = yr_data + n_years)

### add iterations
idx <- lapply(idx, propagate, n_iter)

### insert catchability
for (idx_i in seq_along(idx)) {
  
  ### set catchability for projection
  index.q(idx[[idx_i]])[] <- uncertainty$survey_catchability[[idx_i]]
}
### create copy of index with original values
idx_raw <- lapply(idx ,index)

### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx)

### create deviances for indices; first, get template
idx_dev <- lapply(idx, index)

### create random noise based on covariance for IBTS-Q3
set.seed(4)
idx_dev[[1]][] <- unlist(lapply(lapply(uncertainty$survey_cov, "[[", 1), function(x){
  t(mvrnorm(n = dim(idx_dev[[1]])[2],  mu = rep(0, dim(idx_dev[[1]])[1]), Sigma = x))
}))

### exponentiate to convert to a log-normal scale
idx_dev[[1]] <- exp(idx_dev[[1]])

### create random noise based on sd for the ExIdx
set.seed(4)
### insert sd
idx_dev[[2]][] <- uncertainty$survey_sd[[2]]

### noise
idx_dev[[2]][] <- stats::rnorm(n = length(idx_dev[[2]]), mean = 0, sd = idx_dev[[2]])

### exponentiate toconvert to a log-normal scale
idx_dev[[2]] <- exp(idx_dev[[2]])

### modify residuals for historical period so that index values passed to 
### stock assessment are the ones observed in reality
### IBTS Q3, values up to 2017
idx_dev$DATRAS_Q3_3_8[, dimnames(idx_dev$DATRAS_Q3_3_8)$year <= 2017] <- 
  idx_raw$DATRAS_Q3_3_8[, dimnames(idx_raw$DATRAS_Q3_3_8)$year <= 2017] /
  index(idx$DATRAS_Q3_3_8)[, dimnames(idx$DATRAS_Q3_3_8@index)$year <= 2017]

### Exploitable_biomass_no_age_info, values up to 2017
idx_dev[[2]][, dimnames(idx_dev[[2]])$year <= 2017] <- 
  idx_raw[[2]][, dimnames(idx_raw[[2]])$year <= 2017] /
  index(idx[[2]])[, dimnames(idx[[2]]@index)$year <= 2017]

if (isTRUE(verbose)) {
	### compare simulated to original survey(s)
	as.data.frame(FLQuants(pok_q3 = index(pok_idx$DATRAS_Q3_3_8), 
						   pok_exidx = index(pok_idx$Exploitable_biomass_no_age_info),
						   sim_q3 = (index(idx$DATRAS_Q3_3_8)),
						   sim_exidx = (index(idx$Exploitable_biomass_no_age_info))
	)) %>%
	  mutate(survey = ifelse(grepl(x = qname, pattern = "*_q3$"), "IBTS Q3", "ExIdx"),
			 source = ifelse(grepl(x = qname, pattern = "^sim*"), "sim", "data")) %>%
	  filter(year <= 2019) %>%
	  ggplot(aes(x = year, y = data, colour = source)) +
	  facet_grid(paste("age", age) ~ survey, scales = "free_y") +
	  stat_summary(fun.y = quantile, fun.args = 0.25, geom = "line", alpha = 0.5) +
	  stat_summary(fun.y = quantile, fun.args = 0.75, geom = "line", alpha = 0.5) +
	  stat_summary(fun.y = median, geom = "line") +
	  theme_bw()
}

### ------------------------------------------------------------------------ ###
### catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch
### create noise for catch
set.seed(5)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, sd = uncertainty$catch_sd)

### catch_res values are on a normal scale,
### exponentiate to get log-normal 
catch_res <- exp(catch_res)
### catch_res is a factor by which the numbers at age are multiplied

### for historical period, pass on real observed catch
### -> remove deviation
catch_res[, dimnames(catch_res)$year <= 2017] <- 1
if (isTRUE(verbose)) plot(catch_res)

### ------------------------------------------------------------------------ ###
### save OM objects ####
### ------------------------------------------------------------------------ ###
### set a file path
input_path <- paste0("input/pok/", n_iter, "_", n_years, "/")
dir.create(input_path, recursive = TRUE)

# add stock M criteria
if(m_criteria != "")
  input_path <- paste0(input_path, m_criteria, "_")

### stock
saveRDS(stk_fwd, file = paste0(input_path, "stk.rds"))

### stock-recruitment model
saveRDS(sr, file = paste0(input_path, "sr.rds"))

### recruitment residuals
saveRDS(sr_res, file = paste0(input_path, "sr_res.rds"))

### surveys
saveRDS(idx, file = paste0(input_path, "idx.rds"))
saveRDS(idx_dev, file = paste0(input_path, "idx_dev.rds"))

### catch noise
saveRDS(catch_res, file = paste0(input_path, "catch_res.rds"))

### process error
saveRDS(proc_res, file = paste0(input_path, "proc_res.rds"))

### observed stock
saveRDS(stk_oem, file = paste0(input_path, "stk_oem.rds"))

### sam initial parameters
saveRDS(sam_initial, file = paste0(input_path, "sam_initial.rds"))

### sam configuration
saveRDS(pok_conf, file = paste0(input_path, "pok_conf"))

### catch numbers
saveRDS(catch_n, file = paste0(input_path, "catch_n.rds"))
save.image(file = paste0(input_path, "image.RData"))

### ------------------------------------------------------------------------ ###
### prepare objects for new a4a standard mse package ####
### ------------------------------------------------------------------------ ###
### https://github.com/flr/mse
### reference points  -- these do not change for the alt OMs 
refpts_mse <- list(Btrigger = 149098, Ftrgt = 0.363, Fpa = 0.446, Bpa = 149098, Blim = 107297)

### some specifications for short term forecast with SAM
### using last 3 years, as specified in advice sheet  
pok_stf_def <- list(fwd_yrs_average = -2:0,
                    fwd_yrs_rec_start = 1998,
                    fwd_yrs_sel = -2:0,
                    fwd_yrs_lf_remove = NULL,
                    fwd_splitLD = TRUE)

### arguments passed to mp()
genArgs <- list(fy =  dims(stk_fwd)$maxyear, ### final simulation year
                y0 = dims(stk_fwd)$minyear, ### first simulation year
                iy = yr_data+1, ### first simulation (intermediate) year
                nsqy = 3, ### not used, but has to provided
                nblocks = 90, ### block for parallel processing
                seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr, ### stock recruitment and precompiled residuals
           projection = mseCtrl(method = fwd_WKNSMSE, 
                                args = list(maxF = 2,
                                            ### process noise on stock.n (use fitted value from sr so that parallelization can work)
                                            proc_res = "fitted"
                                ))
)

### observation (error) model
oem <- FLoem(method = oem_WKNSMSE,
            observations = list(stk = stk_oem, idx = idx),
            deviances = list(stk = FLQuants(catch.dev = catch_res),
                             idx = idx_dev),
            args = list(idx_timing = c(-1, -1),		### index timing relative to assessment year: POK.idx end in 2017 both indices (same as assessment year)
                        catch_timing = -1,				### catch timing relative to ay
                        use_catch_residuals = TRUE, 	### use residuals for
                        use_idx_residuals = TRUE,		### observations
                        use_stk_oem = TRUE  ### biological parameters, wts etc
))

### default management procedure (MP)
ctrl_obj <- mpCtrl(list(
  ctrl.est = mseCtrl(method = SAM_wrapper,
                     args = c(### short term forecast specifications
                       forecast = TRUE, 
                       fwd_trgt = "fsq", fwd_yrs = 1, 
                       pok_stf_def,
                       prop_biasN = 1.0,## bias in N-at-age from SAM fit
                       prop_biasF = 1.0,## bias in F-at-age from SAM fit
                       
                       ### speeding SAM up
                       newtonsteps = 0, rel.tol = 0.001,
                       par_ini = list(sam_initial),
                       track_ini = TRUE, ### store ini for next year
                       
                       ### SAM model specifications
                       conf = list(pok_conf),
                       parallel = FALSE ### TESTING ONLY
                     )),
  ctrl.phcr = mseCtrl(method = phcr_WKNSMSE, args = refpts_mse),
  ctrl.hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  ctrl.is = mseCtrl(method = is_WKNSMSE, 
                    args = c(hcrpars = list(refpts_mse),
                             
                             ### for short term forecast
                             fwd_trgt = list(c("fsq", "fsq", "hcr")), fwd_yrs = 3,
                             pok_stf_def,
                             
                             ### TAC constraint
                             TAC_constraint = TRUE,
                             #lower = -Inf, upper = Inf,
                             #Btrigger_cond = FALSE,
                             
                             ### banking and borrowing 
                             BB = TRUE,
                             BB_conditional = TRUE,
                             BB_rho = list(c(-0.00001, 0.00001))
                    ))
))

### additional tracking metrics
tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")

### save mse objects
input <- list(om = om, oem = oem, ctrl.mp = ctrl_obj,
              genArgs = genArgs, tracking = tracking_add)
saveRDS(object = input, file = paste0(input_path, "base_run.rds"))
q("no")
