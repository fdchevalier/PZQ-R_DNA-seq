#!/usr/bin/env Rscript
# Title: Fig1_v7_v10_QTL.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2023-08-25
# Modified in:



#==========#
# Comments #
#==========#

# Draw the figure 1 of the manuscript on updated mapping using the S. mansoni v10 genome



#==========#
# Versions #
#==========#

# v0.0 - 2023-08-25: creation



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    library("magrittr")
    library("GenomicFeatures")
    library("rtracklayer")
    library("gplots")
    library("plotrix")
    library("Sushi")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Graphic
source("functions/rename_chr.R")
source("functions/line2user.R")

## Plotting function for the v7 genome
source("functions/Sm.matplot.data.R")
data.order_v7 <- data.order
matplot.data_v7 <- matplot.data

## Plotting function for the v10 genome
source("functions/Sm.matplot.data_v5.0.R")



#====================#
# Data and variables #
#====================#

# Folders
data_fd  <- "../data/"
graph_fd <- "../graphs/1-GWAS/"

mysfx    <- "v10"
graph_fd <- paste0(graph_fd, "/", mysfx, "/")

# Load data
mydata_v7  <- read.csv("../results/2-QTL/1-GWAS/peaks_report.myfreq.data.fltr.pv.sd-0.1.gq-0.rd-10-no_bf", header=TRUE, sep="\t", stringsAsFactors = FALSE)
mydata_v10 <- read.csv("../results/2-QTL/1-GWAS/v10/peaks_report.myfreq.data.fltr.pv.sd-0.1.gq-0.rd-10-no_bf", header=TRUE, sep="\t", stringsAsFactors = FALSE)

mylo <- read.delim("../results/2-QTL/liftover_v7_10_chr2.bed", header = FALSE)

mychr_v7  <- "SM_V7_3"
mychr_v10 <- "SM_V10_3"
exp.tag   <- c("Exp1","Exp2")

# Load GFF
mygff_fl <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS18.annotations.gff3")
mygff    <- readGFF(mygff_fl)
mygff.genes <- mygff[ mygff[,3] == "gene", ]

# Gene expression
myexpr_fl <- paste0(data_fd, "genome/TPM_isoforms_Sm_2020-08-07.tsv")
myexpr    <- read.csv(myexpr_fl, header = TRUE, sep = "\t")

# Chromosome of interest
coi_v7  <- c("SM_V7_2", "SM_V7_3")
coi_v10 <- c("SM_V10_3")

# Gene of interest
goi <- "Smp_246790"

myseed <- 553309



#=================#
# Data processing #
#=================#

#----------------#
# Table redesign #
#----------------#

# Get combined p.values and build table
## v7
myp.val.cln_v7 <- rev(grep(paste0(".*", paste0(exp.tag, sep = ".*", collapse = "")), colnames(mydata_v7), value = TRUE))[1]
myp.val.tb_v7  <- cbind(mydata_v7[, 1:2], mydata_v7[, myp.val.cln_v7])
myp.val.tb_v7  <- myp.val.tb_v7[! is.na(myp.val.tb_v7[, 3]), ]

## v10
myp.val.cln_v10 <- rev(grep(paste0(".*", paste0(exp.tag, sep = ".*", collapse = "")), colnames(mydata_v10), value = TRUE))[1]
myp.val.tb_v10  <- cbind(mydata_v10[, 1:2], mydata_v10[, myp.val.cln_v10])
myp.val.tb_v10  <- myp.val.tb_v10[! is.na(myp.val.tb_v10[, 3]), ]

# Bonferroni correction
## v7
bf.cor_v7 <- 0.05/nrow(myp.val.tb_v7)

bf.cor_v10 <- 0.05/nrow(myp.val.tb_v10)
my.qtl_v10 <- myp.val.tb_v10[myp.val.tb_v10[,3] <= bf.cor_v10 & myp.val.tb_v10[,1] == mychr_v10, 2]
bf.lim_v10 <- c(my.qtl_v10[1], my.qtl_v10[which(diff(my.qtl_v10) > 1e6)[1]])

# myarr <- sapply(coi, function(x) c(x, myp.val.tb[which.min(myp.val.tb[ myp.val.tb[, 1] == x, 3]), 2])) %>% t()
myarr_v7  <- sapply(coi_v7, function(x) c(x, myp.val.tb_v7[ myp.val.tb_v7[, 1] == x, 2][which.min(myp.val.tb_v7[ myp.val.tb_v7[, 1] == x, 3])])) %>% t()
myarr_v10 <- sapply(coi_v10, function(x) c(x, myp.val.tb_v10[ myp.val.tb_v10[, 1] == x, 2][which.min(myp.val.tb_v10[ myp.val.tb_v10[, 1] == x, 3])])) %>% t()



#=========#
# Figures #
#=========#

cat("\nDrawing graphs...\n")

if(! dir.exists(graph_fd)) { dir.create(graph_fd) }

my.x1 <- 1
my.x2 <- 7e6

mygff.genes <- mygff.genes[ mygff.genes[,1] == mychr_v10 & mygff.genes[,4] >= my.x1 & mygff.genes[,5] <= my.x2 , ]

# Gene expression
mygenes <- mygff.genes[, "Name"]
mycol   <- rep("gray", length(mygenes))
for (i in mygenes) {
    mygene.i <-  grepl(paste0(i, "\\."), myexpr[,1])
    if (any(mygene.i)) {
        exp.gene <- (myexpr[ grepl(paste0(i, "\\."), myexpr[,1]), 14:15] > 0) %>% sum(.) > 0
        if (exp.gene) {
            mycol[i] <- "black"
        }
    }
}

# Data restricted to QTL
myp.val.tb.chr_v10     <- myp.val.tb_v10[ myp.val.tb_v10[,1] == mychr_v10 & myp.val.tb_v10[,2] <= my.x2, ]
myp.val.tb.chr_v10[,3] <- -log10(myp.val.tb.chr_v10[,3])

myp.val.tb_v7 <- rename_chr_SmV7(myp.val.tb_v7, 1)
myp.val.tb_v7 <- data.order_v7(myp.val.tb_v7)
myarr_v7      <- rename_chr_SmV7(myarr_v7, 1)

# Graphic parameters
mycex.axis <- 1.2
mycex <- 1.5

pdf(paste0(graph_fd, "Fig. 1 v10.pdf"), width = 12, height = 12)

    layout(matrix(c(1, 2, 3), byrow = TRUE), height=c(1, 1))
    par(mar = c(5, 5, 3, 0) + 0.1)

    # Plot 1: version 7
    matplot.data_v7(myp.val.tb_v7, ncol(myp.val.tb_v7), "pvalue", xlab = "", ylim.min = 0, ylim.max = 22, abline.h = -log10(bf.cor_v7), abline.lwd = 2, arrow.head = myarr_v7, xlab.axis = c("1","2","3","4","5","6","7","Z","Unass. sc."), by.pos = TRUE, type = "p", cex.axis = mycex.axis, data.order = FALSE)

    ## Add title
    text(par("xaxp")[1], 22 * 0.98, labels = "version 7", adj = 0, cex = 1.2)

    ## Add letter
    mtext("A", side=3, line=1, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

    # Plot 2: the complete data
    matplot.data(myp.val.tb_v10, ncol(myp.val.tb_v10), "pvalue", ylim = c(0, 22), abline.h = -log10(bf.cor_v10), abline.lwd = 2, arrow.head = myarr_v10, by.pos = TRUE, type = "p", col = "red", cex.axis = mycex.axis, gap = 0)

    ## Compute shift to put box (length of the previous chr)
    myshift <- 0
    for (i in c("SM_V10_1", "SM_V10_2")) { myshift <- sum(myshift, myp.val.tb_v10[ tail(which(myp.val.tb_v10[,1] == i),1), 2 ]) }

    zoomsregion(c(my.x1, my.x2) + myshift)

    ## Add title
    text(par("xaxp")[1], 22 * 0.98, labels = "version 10", adj = 0, cex = 1.2)

    ## Add letter
    mtext("B", side=3, line=1, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

    # Plot 3: p-values
    plot(myp.val.tb.chr_v10[,2:3], xlim = c(my.x1, my.x2), type = "n", bty = "n", axes = FALSE, ann = FALSE)

    for (i in 1:nrow(mylo)) {
        rect(mylo[i, 2], 0, mylo[i, 3], par("usr")[4], col = adjustcolor("orange", alpha.f = 0.5), border = NA)
    }

    ## Add gene of interest
    g <- mygff.genes[, "Name"] == goi
    rect(mygff.genes[g, 4], 0, mygff.genes[g, 5], par("usr")[4], col = adjustcolor("lightgreen", alpha.f = 0.5), border=NA)

    points(myp.val.tb.chr_v10[,2:3], col = "red", pch = 19)

    ablineclip(h = -log10(bf.cor_v10), x1 = 0, x2 = my.x2,, col = "blue", lty = 3, lwd = 2)

    ## Annotations
    axis (2, at=format(seq(0,21,length=5),digits=2), cex.axis=mycex.axis)
    mtext(2, text=expression(-Log[10] ~ p-value), line=3) #, cex=mycex)
    mtext("C", side=3, line=1, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

    ## Add x axis
    axis (1, at=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1), labels=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1)/1e6, cex.axis=mycex.axis)
    mtext(1, text="Position (Mb)", line=3) #, cex=mycex)


invisible(dev.off())
