### ------------------------------------------------------------------------ ###
### process results ####
### ------------------------------------------------------------------------ ###

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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


### ------------------------------------------------------------------------ ###
### summary plots: compare HCR options ####
### ------------------------------------------------------------------------ ###

### select maximum yield combinations
combs <- data.frame(name = c("F0","A*", "base", "10", "20", "30", "40", "50"), 
  OM = c("base"),
  HCR = c("F0","A*", "A", "A", "A", "A", "A", "A"),
  BB = c(rep(TRUE, 8)),
  TACconstr = c(rep(TRUE, 8)),
  Btrigger = c(0, rep(Btgr, 7)), 
  Ftrgt = c(0, rep(Ftgt, 7)),
  biasN = c(0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
  biasF = c(0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5),
  scenario = 0) 					
combs<-combs[c(3:8),] # select only the scenarios read in from the output files

combs1 <- stats %>% left_join(combs)
combs1$name <- c("base","10", "20", "30", "40", "50")
combs2 <- gather(data = combs1, key = "key", value = "value",
                 catch_median_long, risk3_long, iav_long,
                 ssb_median_long, recovery_proportion, recovery_time)
combs2$name <- factor(combs2$name, levels = c("F0","A*", "base", "10", "20", "30", "40", "50"))
ggplot(data = combs2, 
       mapping = aes(x = name, y = value, group = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ key, scales = "free_y") +
  theme_bw()

### load entire distribution for stats (read in saved files)
stats_full <- function(data) {
  combs_full <- foreach(i = split(data, seq(nrow(data))), 
                        .packages = "FLCore", .combine = rbind) %dopar% {
                          
                          stk_i <- readRDS(paste0("output/runs/pok/1000_21/", i$file))
                          MSYBtrigger <- 149098
                          Blim <- ifelse(i$OM == "base", 107297, 
                                         ifelse(i$OM == "M01", 90094, 133650))
                          res <- rbind(
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "catch_long",
                                       value = c(window(catch(stk_i@stock), start = 2029))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "catch_medium",
                                       value = c(window(catch(stk_i@stock), start = 2024, end = 2028))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "catch_short",
                                       value = c(window(catch(stk_i@stock), start = 2019, end = 2023))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "risk1_long",
                                       value = mean(window(ssb(stk_i@stock), start = 2029) < Blim)),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "risk1_medium",
                                       value = mean(window(ssb(stk_i@stock), 
                                                           start = 2024, end = 2028) < Blim)),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "risk1_short",
                                       value = mean(window(ssb(stk_i@stock), 
                                                           start = 2019, end = 2023) < Blim)),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "risk3_long",
                                       value = max(iterMeans(window(ssb(stk_i@stock), 
                                                                    start = 2029) < Blim))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "risk3_medium",
                                       value = max(iterMeans(window(ssb(stk_i@stock), 
                                                                    start = 2024, end = 2028) < Blim))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "risk3_short",
                                       value = max(iterMeans(window(ssb(stk_i@stock), 
                                                                    start = 2019, end = 2023) < Blim))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iav_long",
                                       value = c(iav(object = catch(window(stock(stk_i), start = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iav_medium",
                                       value = c(iav(object = catch(window(stock(stk_i), 
                                                                           start = 2023, end = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iav_short",
                                       value = c(iav(object = catch(window(stock(stk_i), 
                                                                           start = 2018, end = 2023))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iavssb_long",
                                       value = c(iav(object = ssb(window(stock(stk_i), start = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iavssb_medium",
                                       value = c(iav(object = ssb(window(stock(stk_i), 
                                                                         start = 2023, end = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iavssb_short",
                                       value = c(iav(object = ssb(window(stock(stk_i), 
                                                                         start = 2018, end = 2023))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iavfbar_long",
                                       value = c(iav(object = fbar(window(stock(stk_i), start = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iavfbar_medium",
                                       value = c(iav(object = fbar(window(stock(stk_i), 
                                                                          start = 2023, end = 2028))))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "iavfbar_short",
                                       value = c(iav(object = fbar(window(stock(stk_i), 
                                                                          start = 2018, end = 2023))))),
                            
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "ssb_long",
                                       value = c(window(ssb(stk_i@stock), start = 2029))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "ssb_medium",
                                       value = c(window(ssb(stk_i@stock), start = 2024, end = 2028))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "ssb_short",
                                       value = c(window(ssb(stk_i@stock), start = 2019, end = 2023))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "F_median_long",
                                       value = c(window(fbar(stk_i@stock), start = 2029))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "F_median_medium",
                                       value = c(window(fbar(stk_i@stock), start = 2024, end = 2028))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "F_median_short",
                                       value = c(window(fbar(stk_i@stock), start = 2019, end = 2023))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "recovery_proportion",
                                       value = mean(apply(window(ssb(stk_i@stock), 
                                                                 start = 2019) >= MSYBtrigger, 6, max))),
                            data.frame(name = i$name, scenario = i$scenario, biasN = i$biasN, biasF = i$biasF,
                                       key = "recovery_time",
                                       value = c(apply(window(ssb(stk_i@stock), 
                                                              start = 2019)@.Data >= MSYBtrigger, 6,
                                                       function(x) {
                                                         if (any(x)) {which(x)[1]} else {Inf}})))
                          )
                          if (i$HCR == "F") {
                            res$value[res$key %in% c("catch_long", "catch_medium", "catch_short",
                                                     "iav_long", "iav_medium", "iav_short")] <- 0
                          }
                          res <- merge(res, i[, c("name", "OM", "HCR", "BB", "TACconstr", "Btrigger",
                                                  "Ftrgt", "biasN", "biasF")])
                          return(res)
                        }
  combs_full$name <- factor(combs_full$name,
                            levels = c("F0","A*", "base", "10", "20", "30", "40", "50"))
  return(combs_full)
}


### base OM
combs_base <- stats_full(data = combs1)
ggplot(data = combs_base, 
       mapping = aes(x = name, y = value, group = name)) +
  #geom_bar(stat = "identity") +
  geom_boxplot() + 
  facet_wrap(~ key, scales = "free_y") +
  theme_bw() +
  ylim(0, NA)
levels(combs_base$name) <- c("F0","A*","base", "10%", "20%", "30%", "40%", "50%")

### get median for option A* or base
combs_base_b <- left_join(combs_base_b, 
                          combs_base_b %>%
                            group_by(key, OM, name) %>%
                            summarise(value_median = median(value)) %>%
                            filter(name == "base") %>%
                            select(-name))
### perceived biased
p_iavssb_long_b <- ggplot(data = combs_base_b[combs_base_b$key == "iavssb_long", ], 
                          mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Interannual variability \nin SSB") +
  theme(axis.text.x = element_blank()) 
p_iavssb_short_b <- ggplot(data = combs_base_b[combs_base_b$key == "iavssb_short", ], 
                           mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.2, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Interannual variability \nin SSB") +
  theme(axis.text.x = element_blank()) 
p_iavfbar_long_b <- ggplot(data = combs_base_b[combs_base_b$key == "iavfbar_long", ], 
                           mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Variability in \nfishing mortality rate") +
  theme(axis.text.x = element_blank()) 
p_iavfbar_short_b <- ggplot(data = combs_base_b[combs_base_b$key == "iavfbar_short", ], 
                            mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.2, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Variability in \nfishing mortality rate") +
  theme(axis.text.x = element_blank())
p_ssb_long_b <- ggplot(data = combs_base_b[combs_base_b$key == "ssb_long", ], 
                       mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  coord_cartesian(ylim = c(0, 1200)) +
  labs(x = "", y = "") +
  labs(x = "", y = "SSB [1000 t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2))
p_ssb_short_b <- ggplot(data = combs_base_b[combs_base_b$key == "ssb_short", ], 
                        mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  coord_cartesian(ylim = c(0, 1000)) +
  labs(x = "", y = "") +
  labs(x = "", y = "SSB [1000 t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2),
                     limits = c(0, NA))
p_fbar_long_b <- ggplot(data = combs_base_b[combs_base_b$key == "F_median_long", ], 
                        mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  ylim(0, 2) +
  labs(x = "", y = "Fishing mortality rate") +
  theme(legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
p_fbar_short_b <- ggplot(data = combs_base_b[combs_base_b$key == "F_median_short", ], 
                         mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  ylim(0, 2) +
  labs(x = "", y = "Fishing mortality rate") +
  theme(legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 

### get median for option A* or base
combs_base <- left_join(combs_base, 
                        combs_base %>%
                          group_by(key, OM, name) %>%
                          summarise(value_median = median(value)) %>%
                          filter(name == "base") %>%
                          select(-name))
### om
p_iavssb_long <- ggplot(data = combs_base[combs_base$key == "iavssb_long", ], 
                        mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Interannual variability \nin SSB") +
  theme(axis.text.x = element_blank()) 
p_iavssb_short <- ggplot(data = combs_base[combs_base$key == "iavssb_short", ], 
                         mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.2, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 1)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Interannual variability \nin SSB") +
  theme(axis.text.x = element_blank()) 
p_iavfbar_long <- ggplot(data = combs_base[combs_base$key == "iavfbar_long", ], 
                         mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Variability in \nfishing mortality rate") +
  theme(axis.text.x = element_blank()) 
p_iavfbar_short <- ggplot(data = combs_base[combs_base$key == "iavfbar_short", ], 
                          mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.2, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Variability in \nfishing mortality rate") +
  theme(axis.text.x = element_blank()) 
p_ssb_long <- ggplot(data = combs_base[combs_base$key == "ssb_long", ], 
                     mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  coord_cartesian(ylim = c(0, 1200)) +
  labs(x = "", y = "") +
  labs(x = "", y = "SSB [1000 t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2))
p_ssb_medium <- ggplot(data = combs_base[combs_base$key == "ssb_medium", ], 
                       mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  coord_cartesian(ylim = c(0, 1100)) +
  labs(x = "", y = "") +
  labs(x = "", y = "SSB [t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2),
                     limits = c(0, NA))
p_ssb_short <- ggplot(data = combs_base[combs_base$key == "ssb_short", ], 
                      mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  coord_cartesian(ylim = c(0, 1000)) +
  labs(x = "", y = "") +
  labs(x = "", y = "SSB [1000 t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 2),
                     limits = c(0, NA))
p_fbar_long <- ggplot(data = combs_base[combs_base$key == "F_median_long", ], 
                      mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  ylim(0, 2) +
  labs(x = "", y = "Fishing mortality rate") +
  theme(legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
p_fbar_short <- ggplot(data = combs_base[combs_base$key == "F_median_short", ], 
                       mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  ylim(0, 2) + 
  labs(x = "", y = "Fishing mortality rate") +
  theme(legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
plot_grid(p_ssb_short, p_iavssb_short,p_fbar_short, p_iavfbar_short,
          ncol=1,align = "hv")
plot_grid( p_ssb_short_b, p_iavssb_short_b, p_fbar_short_b,p_iavfbar_short_b,
           ncol=1,align = "hv")
plot_grid(p_ssb_long, p_iavssb_short,p_fbar_long, p_iavfbar_short,
          ncol=1,align = "hv")
plot_grid(p_ssb_long_b, p_iavssb_short_b,p_fbar_long_b, p_iavfbar_short_b,
          ncol=1,align = "hv")

### management measures
p_catch_long <- ggplot(data = combs_base[combs_base$key == "catch_long", ], 
                       mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 550)) +
  labs(x = "", y = "") +
  labs(x = "", y = "Catch [1000 t]") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 1))
p_catch_medium <- ggplot(data = combs_base[combs_base$key == "catch_medium", ], 
                         mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  labs(x = "", y = "Catch [t]") +
  labs(x = "", y = "") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 1))
p_catch_short <- ggplot(data = combs_base[combs_base$key == "catch_short", ], 
                        mapping = aes(x = name, y = value/1000)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median/1000), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 450)) +
  theme(axis.text.x = element_blank()) +  
  labs(x = "", y = "") +
  labs(x = "", y = "Catch [1000 t]") +
  scale_y_continuous(labels = function(x) format(x, scientific = F, digits = 1))
p_risk3_long <- ggplot(data = combs_base[combs_base$key == "risk3_long", ], 
                       mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", lwd = 0.2, color = "black", aes(fill = name), show.legend = FALSE) + scale_fill_jco() +
  geom_blank(data = combs_base[combs_base$key == "risk1_long", ]) +
  geom_hline(aes(yintercept = 0.05), colour = "red", alpha = 0.5) +
  theme_classic() + ylim(0, 0.45) +
  labs(x = "", y = "") +
  labs(x = "", y = "Risk") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
p_risk3_medium <- ggplot(data = combs_base[combs_base$key == "risk3_medium", ], 
                         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", lwd = 0.2, color = "black", aes(fill = name), show.legend = FALSE) + scale_fill_jco() +
  geom_blank(data = combs_base[combs_base$key == "risk1_medium", ]) +
  geom_hline(aes(yintercept = 0.05), colour = "red", alpha = 0.5) +
  theme_classic() + ylim(0, 0.45) +
  labs(x = "", y = "Risk") +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
p_risk3_short <- ggplot(data = combs_base[combs_base$key == "risk3_short", ], 
                        mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", lwd = 0.2, color = "black", aes(fill = name), show.legend = FALSE) + scale_fill_jco() +
  geom_blank(data = combs_base[combs_base$key == "risk1_short", ]) +
  geom_hline(aes(yintercept = 0.05), colour = "red", alpha = 0.5) +
  theme_classic()+ ylim(0, 0.45) +
  labs(x = "", y = "") +
  labs(x = "", y = "Risk") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
p_iav_long <- ggplot(data = combs_base[combs_base$key == "iav_long", ], 
                     mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Interannual catch\n variability") +
  theme(axis.text.x = element_blank()) 
p_iav_medium <- ggplot(data = combs_base[combs_base$key == "iav_medium", ], 
                       mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_blank()) +  
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") 
p_iav_short <- ggplot(data = combs_base[combs_base$key == "iav_short", ], 
                      mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.2, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  coord_cartesian(ylim = c(0, 2)) + 
  labs(x = "", y = "") +
  labs(x = "", y = "Interannual catch\n variability") +
  theme(axis.text.x = element_blank()) 
p_recovery_proportion <- 
  ggplot(data = combs_base[combs_base$key == "recovery_proportion", ], 
         mapping = aes(x = name, y = value)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = name), lwd = 0.2,show.legend = FALSE) + scale_fill_jco() +
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme_bw() + ylim(0, NA) +
  labs(x = "", y = "") +
  labs(x = "", y = "recovery proportion")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) 
p_recovery_time <- 
  ggplot(data = combs_base[combs_base$key == "recovery_time", ], 
         mapping = aes(x = name, y = value)) +
  geom_violin(trim=FALSE, aes(fill = name), lwd = 0.2, show.legend = FALSE) + 
  geom_boxplot(width=0.25, outlier.size = 0.8, outlier.alpha = 0.5, lwd = 0.2, fill = "white") + 
  theme_classic() + ylim(0, NA) + scale_fill_jco() + 
  geom_hline(aes(yintercept = value_median), colour = "red", alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.2)) +
  labs(x = "", y = "recovery time [years]")

## grouping plots by period 
plot_grid(p_catch_short, p_risk1_short, p_iavssb_short,
          p_iav_short, p_risk3_short,
          ncol=1,align = "hv")
plot_grid(p_catch_long, p_risk1_long, p_iavssb_long,
          p_iav_long, p_risk3_long,
          ncol=1,align = "hv")
plot_grid(p_recovery_proportion, p_recovery_time,ncol=1,
          align = "hv")

