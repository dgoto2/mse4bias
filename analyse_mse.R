### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###

### load required packages
library(FLCore)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(doParallel)

### load additional functions
source("flr_mse_WKNSMSE_funs.R")

### ------------------------------------------------------------------------ ###
### get files ####
### ------------------------------------------------------------------------ ###
path_res <- "output/runs/pok/1000_21/"
files_res <- data.frame(file = list.files(path_res, pattern = "*.rds"), stringsAsFactors = FALSE)
files_res <- files_res[files_res$file != "stats.rds",, drop = FALSE]
files_res$OM <- sapply(strsplit(x = files_res$file, split = "\\HCR"), "[[", 1)
files_res$Ftrgt <- as.numeric(gsub(x = regmatches(x = files_res$file, m = regexpr(text = files_res$file, 
                                                                                  pattern = "Ftrgt-0.[0-9]{1,}")), 
                                   pattern = "Ftrgt-", replacement = ""))
files_res$Btrigger <- as.numeric(gsub(x = regmatches(x = files_res$file, m = regexpr(text = files_res$file,
                                                                                     pattern = "Btrigger-[0-9]{1,}[e+]{0,}[0-9]{0,}")),
                                      pattern = "Btrigger-", replacement = ""))
files_res$HCR <- gsub(x = regmatches(x = files_res$file, m = regexpr(text = files_res$file, pattern = "HCR-[A-Z]{1,}")),
                      pattern = "HCR-", replacement = "")
files_res$HCR[substr(files_res$file,6,6)==0] <- gsub(x = regmatches(x = files_res$file[substr(files_res$file,6,6)==0], 
                                                                    m = regexpr(text = files_res$file[substr(files_res$file,6,6)==0], 
                                                                                pattern = "HCR-[A-Z][0-9]{1,}")), pattern = "HCR-", replacement = "")
files_res$HCR[substr(files_res$file,6,6)==1] <- gsub(x = regmatches(x = files_res$file[substr(files_res$file,6,6)==1], 
                                                                    m = regexpr(text = files_res$file[substr(files_res$file,6,6)==1], 
                                                                                pattern = "HCR-[A-Z][0-9]{1,}")), pattern = "HCR-", replacement = "")
files_res$HCR[substr(files_res$file,10,10)==1] <- gsub(x = regmatches(x = files_res$file[substr(files_res$file,10,10)==1], 
                                                                      m = regexpr(text = files_res$file[substr(files_res$file,10,10)==1], 
                                                                                  pattern = "HCR-[A-Z][0-9]{1,}")), pattern = "HCR-", replacement = "")
files_res$TACconstr <- as.logical(gsub(x = regmatches(x = files_res$file, m = regexpr(text = files_res$file, 
                                                                                      pattern = "TACconstr-TRUE|TACconstr-FALSE")),
                                       pattern = "TACconstr-", replacement = ""))
files_res$BB <- as.logical(gsub(x = regmatches(x = files_res$file, m = regexpr(text = files_res$file, 
                                                                               pattern = "TACconstr-TRUE|TACconstr-FALSE")),
                                pattern = "TACconstr-", replacement = ""))
files_res$biasN <- as.numeric(gsub(x = regmatches(x = files_res$file,  m = regexpr(text = files_res$file, 
                                                                                   pattern = "biasN-[-|+|\\+]0.[0-9]{1,}")),
                                   pattern = "biasN-", replacement = ""))
files_res$biasF <- as.numeric(gsub(x = regmatches(x = files_res$file, m = regexpr(text = files_res$file,  
                                                                                  pattern = "biasF-[-|+|\\+]0.[0-9]{1,}")),
                                   pattern = "biasF-", replacement = ""))

if(file.exists(paste0(path_res, "stats.rds"))) {
  stats <- readRDS(paste0(path_res, "stats.rds"))
  stats_new <- merge(stats, files_res, all = TRUE)
  stats_new <- stats_new[!stats_new$file %in% stats$file, ]
  if( nrow(stats_new) < 1 ) stats_new <- files_res
} else {
  stats_new <- files_res
}

### set Blim depending on OM
stats_new$Blim <- sapply(stats_new$OM, function(x) {
  switch(x, "base" = 107297)})
res_list <- foreach(i = seq(nrow(stats_new)), .packages = "mse") %dopar% {
  tmp <- readRDS(paste0(path_res, stats_new$file[i]))
  tmp@oem <- FLoem()
  tmp
}

### ------------------------------------------------------------------------ ###
### calculate performance statistics ####
### ------------------------------------------------------------------------ ###
### calculate for short- (year 1-5), medium- (year 6-10) and long-term (year 11-20)
### risk: proportion of stock below Blim, maximum over iterations and period ('risk3' in ICES MSEs)
### catch: median catch in period (over years and iterations)
### iav: inter-annual variation of catch, average over years and iterations
### SSB
stats_new$ssb_median_long <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2029))  
}
stats_new$ssb_median_short <- foreach(x = res_list, .packages = "FLCore",  .combine = "c") %dopar% {
  median(window(ssb(x@stock), start = 2019, end = 2023))  
}
stats_new$ssb_median_medium <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {  
  median(window(ssb(x@stock), start = 2024, end = 2028)) 
}

### risk
stats_new$risk3_long <- foreach(x = res_list, Blim = stats_new$Blim[seq(nrow(stats_new))], .packages = "FLCore", .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2029) < Blim))
}
stats_new$risk3_short <- foreach(x = res_list, Blim = stats_new$Blim[seq(nrow(stats_new))], .packages = "FLCore", .combine = "c") %dopar% {
  max(iterMeans(window(ssb(x@stock), start = 2019, end = 2023) < Blim))
}
stats_new$risk3_medium <- foreach(x = res_list, Blim = stats_new$Blim[seq(nrow(stats_new))], .packages = "FLCore", .combine = "c") %dopar% {  
  max(iterMeans(window(ssb(x@stock), start = 2024, end = 2028) < Blim))  
}

### catch 
stats_new$catch_median_long <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2029))  
}
stats_new$catch_median_short <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2019, end = 2023))
}
stats_new$catch_median_medium <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(catch(x@stock), start = 2024, end = 2028)) 
}

### inter-annual variation of catch
stats_new$iav_long <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% { 
  iav(object = catch(window(stock(x), start = 2028)), summary_all = median)
}
stats_new$iav_short <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {  
  iav(object = catch(window(stock(x), start = 2018, end = 2023)), summary_all = median)
}
stats_new$iav_medium <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {     
  iav(object = catch(window(stock(x), start = 2023, end = 2028)),  summary_all = median)  
}

### Fbar
stats_new$fbar_median_long <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(fbar(x@stock), start = 2029))  
}
stats_new$fbar_median_short <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(fbar(x@stock), start = 2019, end = 2023)) 
}
stats_new$fbar_median_medium <- foreach(x = res_list, .packages = "FLCore", .combine = "c") %dopar% {
  median(window(fbar(x@stock), start = 2024, end = 2028))
}

### export output files as rds and csv
if(exists("stats")) {
  stats <- rbind(stats, stats_new)
  stats <- stats[order(stats$file), ]
} else {
  stats <- stats_new
}
saveRDS(object = stats, file = paste0(path_res, "stats.rds"))
write.csv(x = stats, file = paste0(path_res, "stats.csv"), row.names = FALSE)
