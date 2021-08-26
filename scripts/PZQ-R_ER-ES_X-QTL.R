#!/usr/bin/env Rscript
# Title: PZQ-R_ER-ES_X-QTL.R
# Version: 0.3
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-14-13
# Modified in: 2021-08-25



#==========#
# Comments #
#==========#

# Analyze genomic data from ER and ES selected populations using X-QTL method.



#==========#
# Versions #
#==========#

# v0.3 - 2021-08-25: rename script / add X-QTL module / correct number of worms / clean code
# v0.2 - 2018-05-01: update source data to use VCF file directly
# v0.1 - 2018-03-15: update email address / remove unnecessary parts



#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/2-ER-ES populations/"
result_fd <- "../results/2-QTL/2-ER-ES populations/"

# Data file
myvcf_file <- paste0(data_fd, "calling/PZQ_ER-ES.vcf.gz")

# Annotation data
mygff_file <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
myann_file <- paste0(data_fd, "genome/Sm_transcript_table_gff-hhpred_2020-09-21.tsv")

# Expression data
myexpr     <- read.csv(paste0(data_fd, "genome/TPM_isoforms_Sm_2020-08-07.tsv"), header=TRUE, sep="\t")
myexpr.cln <- 9:12    # Column containing juvenile and adult expression data
myexpr.nm  <- paste0("TPM ", c("juv m", "juv f", "adt m", "adt f"))

# Baits
mybaits <- NULL

# Lib
mylib     <- c("SmLE-PZQ-ES_F1_m", "SmLE-PZQ-ER_F1_m", "SmLE-PZQ-ES_F1_f", "SmLE-PZQ-ER_F1_f")
mylib.ind <- c(100, 100, 25, 25)
mylib.rep <- mylib 

# Allele polarization (list of samples)
mypol.all <- matrix(c(
        "SmLE-PZQ-ES_F1_m", "SmLE-PZQ-ER_F1_m",
        "SmLE-PZQ-ER_F1_m", "SmLE-PZQ-ES_F1_m",
        "SmLE-PZQ-ES_F1_f", "SmLE-PZQ-ER_F1_f",
        "SmLE-PZQ-ER_F1_f", "SmLE-PZQ-ES_F1_f"
), ncol=2, byrow=TRUE)

# Z-score comparison (list of samples)
mycomp <- t(combn(mylib, 2))

# Z-score combination (list of computed z-scores)
mycomp.cb <- matrix(c(
        "SmLE-PZQ-ES_F1_m-SmLE-PZQ-ER_F1_m", "SmLE-PZQ-ES_F1_f-SmLE-PZQ-ER_F1_f"
            ), ncol=2, byrow=TRUE)

# Pattern for selecting column on which to apply filtering using Bonferroni correction
mypv.pat <- ".*ES.*ER.*ES"

# Pattern for selecting allele frequency column to report
myfrq.pat <- ".*ES.*ER.*"

# Run the X-QTL module
source("X-QTL_module.R")
