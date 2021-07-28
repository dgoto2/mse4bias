### ------------------------------------------------------------------------ ###
### R script to create a stock object for an operating model (OM) ####
### ------------------------------------------------------------------------ ###

### set a directory for input data
my.dir2 <- "./data/"

### read in the model files
pok <- FLStock(name="POK3a46 from May 2018 assessment", desc='SAM latest TMB version', 
               FLQuant(NA, dimnames=list(age=3:10, year=1967:2017)))
range(pok, c('minfbar', 'maxfbar')) <- c(4, 7)
discards.wt(pok) <- readVPAFile(paste(my.dir2, "dw.dat", sep=''))
catch.n(pok) <- readVPAFile(paste(my.dir2, "cn.dat", sep=''))
catch.wt(pok) <- readVPAFile(paste(my.dir2, "cw.dat", sep=''))
m.spwn(pok) <- readVPAFile(paste(my.dir2, 'pm.dat', sep=''))
harvest.spwn(pok) <- readVPAFile(paste(my.dir2, 'pf.dat', sep=''))
landings.wt(pok) <- readVPAFile(paste(my.dir2, "lw.dat", sep=''))
mat(pok) <- readVPAFile(paste(my.dir2, "mo.dat", sep=''))
m(pok) <- 0.2
stock.wt(pok) <- readVPAFile(paste(my.dir2, "cw.dat", sep=''))
lf <-readVPAFile(paste(my.dir2, 'lf.dat', sep=''))

## estimate landings and discards numbers
landings.n(pok) <- catch.n(pok)*lf
discards.n(pok) <- as.numeric(catch.n(pok)*(1-lf))
landings(pok) <- computeLandings(pok)      
discards(pok) <- computeDiscards(pok)     
catch(pok) <- computeCatch(pok, "catch")
stock(pok) <- computeStock(pok)
range(pok,c("minfbar", "maxfbar")) <- c(4, 7)
plusgroup <- 10
pok<-setPlusGroup(pok, plusgroup, na.rm=TRUE)
pok@m[1:8,,,,] <- 0.2
length(units(pok)); units(pok)
units(pok)[1:17] <- as.list(c(rep(c("tonnes", "thousands", "kg"), 4), "NA", "NA", "f", "NA", "NA"))
units(pok)

## read in input data files
indices <- readFLIndices("./data/survey.dat") #incl. DATRAS ages 3-8
lapply(indices, summary)
pok_stk <- pok
pok_idx <- indices

## alternative OMs with different Ms (natural mortality rates)
### M=0.1
pok@m[1:8,,,,] <- 0.1
pok_stkM01 <- pok

### M=0.3
pok@m[1:8,,,,] <- 0.3
pok_stkM03 <- pok

### export the stock objects
dir.create("./output")
save(pok_stk, pok_stkM03, pok_stkM01, pok_idx, file="./output/POK273a46.stock.object.Rdata")
rm(pok, plusgroup, indices, lf)
