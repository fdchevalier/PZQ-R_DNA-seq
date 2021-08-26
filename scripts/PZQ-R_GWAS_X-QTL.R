#!/usr/bin/env Rscript
# Title: PZQ_analysis
# Version: 0.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-14-13
# Modified in: 2018-05-01



#==========#
# Comments #
#==========#

# v0.2 - 2018-05-01: update source data to use VCF file directly
# v0.1 - 2018-03-15: update email address / remove unnecessary parts



#===========#
# Variables #
#===========#

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/1-QTL GWAS"

# Data file
myvcf_file <- paste0(data_fd, "calling/PZQ_GWAS_snpEff.vcf.gz")

# Annotation data
mygff_file <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
myann_file <- paste0(data_fd, "genome/Sm_transcript_table_gff-hhpred_2020-09-21.tsv")

# Expression data
myexpr     <- read.csv(paste0(data_fd, "genome/TPM_isoforms_Sm_2020-08-07.tsv"), header=TRUE, sep="\t")
myexpr.cln <- 9:12    # Column containing juvenil and adult expression data
myexpr.nm  <- paste0("TPM ", c("juv m", "juv f", "adt m", "adt f"))

# Library
mylib     <- c("SmLE-PZQ-R_Exp1_Contracted", "SmLE-PZQ-R_Exp1_Recovered", "SmLE-PZQ-R_Exp2_Contracted", "SmLE-PZQ-R_Exp2_Recovered")
mylib.ind <- c(116, 116, 137, 137)
mylib.rep <- mylib 

#mycolnames <- colnames(mydata)

# Allele polarization (list of samples)
mypol.all <- matrix(c(
        "SmLE-PZQ-R_Exp1_Contracted", "SmLE-PZQ-R_Exp1_Recovered",
        "SmLE-PZQ-R_Exp1_Recovered", "SmLE-PZQ-R_Exp1_Recovered",
        "SmLE-PZQ-R_Exp2_Contracted", "SmLE-PZQ-R_Exp2_Recovered",
        "SmLE-PZQ-R_Exp2_Recovered", "SmLE-PZQ-R_Exp2_Recovered"
), ncol=2, byrow=TRUE)

# Z-score comparaison (list of samples)
mycomp <- t(combn(mylib, 2))

# Z-score combinaison (list of computed z-scores)
mycomp.cb <- matrix(c(
        "SmLE-PZQ-R_Exp1_Contracted-SmLE-PZQ-R_Exp1_Recovered", "SmLE-PZQ-R_Exp2_Contracted-SmLE-PZQ-R_Exp2_Recovered"
            ), ncol=2, byrow=TRUE)

# Pattern for selecting column on which to apply filtering using Bonferroni correction
mypv.pat <- ".*Con.*Rec.*Con"

# Pattern for selecting allele frequency column to report
myfrq.pat <- ".*Con.*Rec.*"
