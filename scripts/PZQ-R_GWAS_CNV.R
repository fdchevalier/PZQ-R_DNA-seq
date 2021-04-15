#!/usr/bin/env Rscript
# Title: PZQ-R_GWAS_CNV.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-03-12
# Modified in: 



#==========#
# Comments #
#==========#

# v0.0 - 2018-03-12: creation



#==========#
# Pacakges #
#==========#

library(cn.mops)
library(Biobase)
library(parallel)

exp.tag <- c("Exp1","Exp2")

BAMFolder <- "../1-Alignment/data/" 


segments  <- read.table("~/data/sm_exons/sma_agilent_baits.v7.0.chr_reorderd.bed", sep="\t", as.is=TRUE)

mychr <- c(paste0("Chr_",seq(1,7)), "Chr_W")

n.cores <- detectCores()

if (n.cores > 1) { n.cores <- round(n.cores*2/3) } else {n.cores <- 0}

# Test on the complete BAM
segments.chr  <- segments[ grep(paste0(mychr,"$",collapse="|"),segments[,1]) , ]

gr <- GRanges(segments.chr[,1],IRanges(segments.chr[,2],segments.chr[,3]))

# To store call table
call.tb.ls <- vector("list", length(exp.tag))
names(call.tb.ls) <- exp.tag

for (t in exp.tag) {
    
    # List of BAM files
    BAMFiles <- list.files(BAMFolder, pattern=paste0(t,".*bam$"), recursive=TRUE, full.name=TRUE)

    # Get read counts for the group
    X  <- getSegmentReadCountsFromBAM(BAMFiles, GR=gr, parallel=n.cores)

    # Remove redundant markers
    X <- X[isUnique(X),]

    # Scan for CNVs
    resCNMOPS <- exomecn.mops(X, parallel=n.cores)

    # Number of regions with CNVs
    nb.regions <- length(resCNMOPS@cnvr)

    # Skip if no CNV detected
    if (nb.regions > 0) {
        resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

        print(resCNMOPS)

        regions <- as.vector(resCNMOPS@cnvr@seqnames)

        for (i in 1:nb.regions) {
            pdf(paste0("graphs/CNV_region_",regions[i],"_",t,".pdf"))
                plot(resCNMOPS, which=i, toFile=TRUE)
            dev.off()
        }

        call.tb <- cbind(as.data.frame(resCNMOPS@gr)[,1:3], as.data.frame(resCNMOPS@individualCall))

        for (c in regions) {
            call.tb.chr <- call.tb[ grep(c, call.tb[,1]), ]

            for (i in 1:length(BAMFiles)) {
                call.val <- call.tb.chr[,c(2:3,3+i)]
                myname   <- paste(strsplit(colnames(call.tb.chr)[3+i], "_")[[1]][1:4], collapse=" ")

                png(paste0("graphs/CNV_",c,"_",myname,".png"), width=13*72, height=8*72)
                    my.xlim <- max(call.val[,2])
                    my.ylim <- max(abs(call.val[,3]))
                    if (my.ylim < 1) { my.ylim <- 1 }
                    plot(0,xlim=c(1,my.xlim), ylim=c(-my.ylim, my.ylim), type="n", xlab="Position (bp)", ylab="Log ratio of the call", main=myname)
                    segments(call.val[,1], call.val[,3], call.val[,2], lwd=5, col="red")
                dev.off()
            }
        }

        # Store table for use in downstream script
        call.tb.ls[[match(t,exp.tag)]] <- call.tb
    }
}
