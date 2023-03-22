#!/usr/bin/env Rscript
# Title: PZQ-R_GWAS_X-QTL_v9.R
# Version: 0.4
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-14-13
# Modified in: 2023-03-22



#==========#
# Comments #
#==========#

# Analyze GWAS data using X-QTL method.



#==========#
# Versions #
#==========#

# v0.4 - 2023-03-22: adapt the PZQ-R_GWAS_X-QTL.R to S. mansoni v9 genome
# v0.3 - 2021-08-24: rename script / add X-QTL module / clean code
# v0.2 - 2018-05-01: update source data to use VCF file directly
# v0.1 - 2018-03-15: update email address / remove unnecessary parts



#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/1-GWAS/"
result_fd <- "../results/2-QTL/1-GWAS/"
mysfx     <- "v9"

# Data file
myvcf_file <- paste0(data_fd, "calling/PZQ_GWAS_v9_snpEff.vcf.gz")

# Annotation data
mygff_file <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS17.annotations.gff3")
myann_file <- paste0(data_fd, "genome/Sm_transcript_table_gff-hhpred_2020-09-21.tsv")

# Expression data
myexpr     <- read.csv(paste0(data_fd, "genome/TPM_isoforms_Sm_2020-08-07.tsv"), header=TRUE, sep="\t")
myexpr.cln <- 9:12    # Column containing juvenile and adult expression data
myexpr.nm  <- paste0("TPM ", c("juv m", "juv f", "adt m", "adt f"))

# Baits
mybaits <- NULL

# Library
mylib     <- c("SmLE-PZQ-R_Exp1_Contracted", "SmLE-PZQ-R_Exp1_Recovered", "SmLE-PZQ-R_Exp2_Contracted", "SmLE-PZQ-R_Exp2_Recovered")
mylib.ind <- c(116, 116, 137, 137)
mylib.rep <- mylib 

# Allele polarization (list of samples)
mypol.all <- matrix(c(
        "SmLE-PZQ-R_Exp1_Contracted", "SmLE-PZQ-R_Exp1_Recovered",
        "SmLE-PZQ-R_Exp1_Recovered", "SmLE-PZQ-R_Exp1_Recovered",
        "SmLE-PZQ-R_Exp2_Contracted", "SmLE-PZQ-R_Exp2_Recovered",
        "SmLE-PZQ-R_Exp2_Recovered", "SmLE-PZQ-R_Exp2_Recovered"
), ncol=2, byrow=TRUE)

# Z-score comparison (list of samples)
mycomp <- t(combn(mylib, 2))

# Z-score combination (list of computed z-scores)
mycomp.cb <- matrix(c(
        "SmLE-PZQ-R_Exp1_Contracted-SmLE-PZQ-R_Exp1_Recovered", "SmLE-PZQ-R_Exp2_Contracted-SmLE-PZQ-R_Exp2_Recovered"
            ), ncol=2, byrow=TRUE)

# Pattern for selecting column on which to apply filtering using Bonferroni correction
mypv.pat <- ".*Con.*Rec.*Con"

# Pattern for selecting allele frequency column to report
myfrq.pat <- ".*Con.*Rec.*"

# Run the X-QTL module
source("X-QTL_module.R")
