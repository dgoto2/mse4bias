library(mse)

iters <- 1000
years <- 21

outDir <- "output/runs/pok/"
inDir <- "input/pok/"
plotDir <- "plots/"

prefix <- paste0(iters, "_", years)

m_crit <- c("", "M01", "M03")

suffix <- c("HCR-A_Ftrgt-0.363_Btrigger-149098_TACconstr-FALSE_BB-FALSE.rds")

dbs <- paste0(outDir,  prefix, "/", paste(m_crit, suffix, sep = "_"))
inputs <- paste0(inDir, prefix, "/", paste(m_crit, "base_run.rds", sep = "_"))

# cleanup texts
dbs <- gsub("/_", "/", dbs)
inputs <- gsub("/_", "/", inputs)

plotAll <- function(db) {

	if(!file.exists(db)) return(FALSE)

	x <- readRDS(db)

	dir.create(plotDir)

	pdf(paste0(plotDir, prefix, "_", tools::file_path_sans_ext(basename(db)), ".pdf"))
	print(plot(x) + ggtitle("MSE Result"))
	print(plot(index(x@oem@observations$idx[[1]])) + ggtitle("IBTS Index"))
	print(plot(index(x@oem@observations$idx[[2]])) + ggtitle("ExIdx Index"))
	print(plot(stock.n(x@stock)) + ggtitle("stock.n"))
	print(plot(x@tracking[c("F.est", "F.om")]) + ggtitle("F Performance"))
	print(plot(x@tracking[c("B.est", "B.om")]) + ggtitle("SSB Performance"))

	# Process input
	no <- (db == dbs)
	input <- readRDS(inputs[no])
	# ### get stock before simulation
	stk <- input$om@stock
	# ### add simulated data
	stk[, dimnames(x@stock)$year] <- x@stock
	# ### save
	saveRDS(object = stk, file = paste0(outDir, prefix, "_",
                                     m_crit[no], "_base_full_stk.rds"))
	# ### plot
	stkPlot <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) + 
	   xlab("year") + geom_vline(xintercept = 2018.5) +
	   geom_hline(data = data.frame(qname = "SSB", data = input$ctrl.mp$ctrl.phcr@args$Blim),
		      aes(yintercept = data), linetype = "dashed") +
	   geom_hline(data = data.frame(qname = "SSB", data = input$ctrl.mp$ctrl.phcr@args$Btrigger),
		      aes(yintercept = data), linetype = "solid") +
	   geom_hline(data = data.frame(qname = "F", data = input$ctrl.mp$ctrl.phcr@args$Fpa),
		      aes(yintercept = data), linetype = "dashed") +
	   geom_hline(data = data.frame(qname = "F", data = input$ctrl.mp$ctrl.phcr@args$Ftrgt),
		      aes(yintercept = data), linetype = "solid") +
	   theme_bw() + ggtitle("Full stock")

	print(stkPlot)

	dev.off()
	return(TRUE)
}

lapply(dbs, plotAll)

