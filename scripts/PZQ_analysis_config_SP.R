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
result_fd <- "../results/3-QTL Selected population/"

# Data file
myvcf_file <- paste0(data_fd, "calling/PZQ_ER-ES.vcf.gz")

# Annotation data
mygff_file <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
myann_file <- paste0(data_fd, "genome/Sm_transcript_table_gff-hhpred_2020-09-21.tsv")

# Expression data
myexpr     <- read.csv(paste0(data_fd, "genome/TPM_isoforms_Sm_2020-08-07.tsv"), header=TRUE, sep="\t")
myexpr.cln <- 9:12    # Column containing juvenile and adult expression data
myexpr.nm  <- paste0("TPM ", c("juv m", "juv f", "adt m", "adt f"))

mybaits <- read.csv("~/data/sm_exons/sma_agilent_baits.v7.0.chr_reorderd.bed", header=FALSE, sep="\t")

# Lib
mylib     <- c("SmLE-PZQ-ES_F1_m", "SmLE-PZQ-ER_F1_m", "SmLE-PZQ-ES_F1_f", "SmLE-PZQ-ER_F1_f")
mylib.ind <- c(116, 116, 137, 137)
mylib.rep <- mylib 

# Allele polarization (list of samples)
#mypol.all <- matrix(c(
#        "SmLE-PZQ-ES_F1_m", "SmLE-PZQ-ER_F1_m",
#        "SmLE-PZQ-ES_F1_f", "SmLE-PZQ-ER_F1_f",
#        "SmLE-PZQ-ES_F1_f", "SmLE-PZQ-ES_F1_m",
#        "SmLE-PZQ-ER_F1_f", "SmLE-PZQ-ER_F1_m"
#), ncol=2, byrow=TRUE)
mypol.all <- matrix(c(
        "SmLE-PZQ-ES_F1_m", "SmLE-PZQ-ER_F1_m",
        "SmLE-PZQ-ER_F1_m", "SmLE-PZQ-ES_F1_m",
        "SmLE-PZQ-ES_F1_f", "SmLE-PZQ-ER_F1_f",
        "SmLE-PZQ-ER_F1_f", "SmLE-PZQ-ES_F1_f"
), ncol=2, byrow=TRUE)

# Z-score comparaison (list of samples)
mycomp <- t(combn(mylib, 2))

# Z-score combinaison (list of computed z-scores)
mycomp.cb <- matrix(c(
        "SmLE-PZQ-ES_F1_m-SmLE-PZQ-ER_F1_m", "SmLE-PZQ-ES_F1_f-SmLE-PZQ-ER_F1_f"
            ), ncol=2, byrow=TRUE)

# Pattern for selecting column on which to apply filtering using Bonferroni correction
mypv.pat <- ".*ES.*ER.*ES"

# Pattern for selecting allele frequency column to report
myfrq.pat <- ".*ES.*ER.*"

# SmLE-PZQ-ER_F1_f.1
# SmLE-PZQ-ER_F1_f.2
# SmLE-PZQ-ER_F1_f.3
# SmLE-PZQ-ER_F1_m.1
# SmLE-PZQ-ER_F1_m.2
# SmLE-PZQ-ER_F1_m.3
# SmLE-PZQ-ES_F1_f.1
# SmLE-PZQ-ES_F1_f.2
# SmLE-PZQ-ES_F1_f.3
# SmLE-PZQ-ES_F1_m.1
# SmLE-PZQ-ES_F1_m.2
# SmLE-PZQ-ES_F1_m.3
# SmLE-PZQ-R_Exp1_Contracted_1_WGS
# SmLE-PZQ-R_Exp1_Contracted_2_WGS
# SmLE-PZQ-R_Exp1_Contracted_3_WGS
# SmLE-PZQ-R_Exp1_Recovered_1_WGS
# SmLE-PZQ-R_Exp1_Recovered_2_WGS
# SmLE-PZQ-R_Exp1_Recovered_3_WGS
# SmLE-PZQ-R_Exp2_Contracted_1_WGS
# SmLE-PZQ-R_Exp2_Contracted_2_WGS
# SmLE-PZQ-R_Exp2_Contracted_3_WGS
# SmLE-PZQ-R_Exp2_Recovered_1_WGS
# SmLE-PZQ-R_Exp2_Recovered_2_WGS
# SmLE-PZQ-R_Exp2_Recovered_3_WGS

