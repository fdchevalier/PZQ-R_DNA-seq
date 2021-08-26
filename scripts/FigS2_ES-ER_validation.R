#!/usr/bin/env Rscript
# Title: FigS2_ER-ES_validation.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-08-26
# Modified in: 



#==========#
# Comments #
#==========#

# Module performing X-QTL analysis. The module must be called from the main script.



#==========#
# Versions #
#==========#

# v0.0 - 2021-08-26: creation



#======================#
# Packages and options #
#======================#

cat("Loading packages\n")

suppressMessages({
    library("magrittr")
})

options(stringsAsFactors = FALSE)



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Chromosome renaming
source("functions/rename_chr.R")

# Plotting function
source("functions/Sm.matplot.data.R")



#===========#
# Variables #
#===========#

# Folders
data_fd    <- "../data/"
results_fd <- "../results/2-QTL/2-ER-ES populations/"
graph_fd   <- "../graphs/"

# Load data
load("../results/2-QTL/2-ER-ES populations/myfreq.data.fltr.RData")



#=========#
# Figures #
#=========#

cat("Generating figures\n")

if (! dir.exists(graph_fd)) {dir.create(graph_fd, recursive=TRUE)}

# Columns to use
pv.vec2 <- 47
sf.vec2 <- 29

# Data reshaping
myfreq.data.fltr.r <- rename_chr_SmV7(myfreq.data.fltr, 1)
myfreq.data.fltr.r <- myfreq.data.fltr.r[ is.finite(myfreq.data.fltr.r[,pv.vec2]), ]

mymax <- -log10(myfreq.data.fltr.r[,pv.vec2]) %>% max() %>% floor()

# png(paste0(graph_fd, myfn, ".pv-sf.png"), width=72*12, height=72*10)
pdf(paste0(graph_fd, "Supp. Fig. 2 - ER-ES-pv-sf.pdf"), width = 12, height = 10)
layout(matrix(1:2, ncol = 1))
par(mar=c(5,4,1,1)+0.1) # For small size
matplot.data(myfreq.data.fltr.r, pv.vec2, "pvalue", ylim.min=0, ylim.max=mymax, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), by.pos=TRUE, type="p")
matplot.data(myfreq.data.fltr.r, sf.vec2, "freq", ylim.min=-1, ylim.max=1, abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), ylab="Difference in allele frequency", by.pos=TRUE, type="p")
dev.off()

