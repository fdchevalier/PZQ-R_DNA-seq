#!/usr/bin/env Rscript
# Title: PZQ-R_GWAS_CNV.R
# Version: 0.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-03-12
# Modified in: 2021-04-15



#==========#
# Comments #
#==========#

# v0.2 - 2021-04-15: adapt to the notebook environment
# v0.1 - 2019-10-17: adapt script to the whole genome sequencing data
# v0.0 - 2018-03-12: creation



#==========#
# Pacakges #
#==========#
cat("Loading packages...\n")

suppressMessages({
    library(cn.mops)
    library(Biobase)
    library(parallel)
})



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/1-QTL/"


exp.tag <- c("Exp1","Exp2")

BAMFolder <- paste0(data_fd, "libraries/1-GWAS/")

if(! dir.exists(graph_fd)) { dir.create(graph_fd) }

mychr <- paste0("SM_V7_", c(1:7, "ZW"))

# Number of core used
if (detectCores() > 120) { 
    n.cores <- 120
} else {
    n.cores <- detectCores()
}

# To store call table
call.tb.ls <- vector("list", length(exp.tag))
names(call.tb.ls) <- exp.tag

for (t in exp.tag) {
    
    # List of BAM files
    BAMFiles <- list.files(BAMFolder, pattern = paste0(t,".*bam$"), recursive = TRUE, full.name = TRUE)

    # Get read counts for the group
    X  <- getReadCountsFromBAM(BAMFiles, refSeqNames = mychr, parallel = n.cores)

    # Scan for CNVs
    resCNMOPS <- referencecn.mops(X[,4:6], X[,1:3], minWidth = 3, segAlgorithm = "fast", parallel = n.cores)

    # Number of regions with CNVs
    nb.regions <- length(resCNMOPS@cnvr)

    # Skip if no CNV detected
    if (nb.regions > 0) {
        resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

        print(resCNMOPS)

        regions <- as.vector(resCNMOPS@cnvr@seqnames)

        for (i in 1:nb.regions) {
            pdf(paste0(graph_fd, "CNV_region_", i, "_", regions[i], "_", t, ".pdf"))
                plot(resCNMOPS, which = i, toFile = TRUE)
            dev.off()
        }

        call.tb <- cbind(as.data.frame(resCNMOPS@gr)[,1:3], as.data.frame(resCNMOPS@individualCall))

        for (c in regions) {
            call.tb.chr <- call.tb[ grep(c, call.tb[,1]), ]

            for (i in 4:ncol(call.tb.chr)) {
                call.val <- call.tb.chr[,c(2:3,i)]
                myname   <- paste(strsplit(colnames(call.tb.chr)[i], "_")[[1]][1:4], collapse = " ")

                pdf(paste0(graph_fd, "CNV_", c, "_", myname, ".pdf"), width = 13, height = 8)
                    my.xlim <- max(call.val[,2])
                    my.ylim <- max(abs(call.val[,3]))
                    if (my.ylim < 1) { my.ylim <- 1 }
                    plot(0,xlim = c(1,my.xlim), ylim = c(-my.ylim, my.ylim), type = "n", xlab = "Position (bp)", ylab = "Log ratio of the call", main = myname)
                    segments(call.val[,1], call.val[,3], call.val[,2], lwd = 5, col = "red")
                dev.off()
            }
        }

        # Store table for use in downstream script
        call.tb.ls[[match(t,exp.tag)]] <- call.tb
    }
}
