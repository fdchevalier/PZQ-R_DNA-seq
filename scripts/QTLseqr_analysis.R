#!/usr/bin/env Rscript
# Title: QTLseqr_analysis.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-05-13
# Modified in: 2021-08-25



#==========#
# Comments #
#==========#

# Perform QTL analysis using G' statistics.



#==========#
# Versions #
#==========#

# v0.1 - 2021-08-25: format script / correct minor bug
# v0.0 - 2021-05-13: creation



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    library("QTLseqr")
    library("vcfR")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("X-QTL_func.R")

source("functions/csvConnection.R")

#-------------------#
# Sm exome plotting #
#-------------------#

# Chromosome renaming
source("functions/rename_chr.R")

# Plotting function
source("functions/Sm.matplot.data.R")



#===========#
# Variables #
#===========#

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/1-GWAS/"
result_fd <- "../results/2-QTL/1-GWAS/"

# Data file
myvcf_file <- paste0(data_fd, "calling/PZQ_GWAS_snpEff.vcf.gz")

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

# Comparison (list of samples)
mycomp <- t(combn(mylib, 2))



#=================#
# Data processing #
#=================#

# Set core cluster
registerDoParallel(cores = detectCores())

# Load data
myvcf <- read.vcfR(myvcf_file)

# Input file name
myfn <- filename(myvcf_file, length=1)


#-----------------#
# First filtering #
#-----------------#

# Initial number of variants
nb.var <- nrow(myvcf)

# SNPs only
snp.rows <- ! is.indel(myvcf)
myvcf    <- myvcf[ snp.rows, ]
snp.rows <- sum(snp.rows)

# Biallelic sites only
bial.rows <- is.biallelic(myvcf)
myvcf     <- myvcf[ bial.rows, ]
bial.rows <- sum(bial.rows)


# Build read depth table
mydata.ad <- extract.gt(myvcf, "AD")
mydata.gq <- extract.gt(myvcf, "GQ")


mydata <- as.data.frame(cbind( getFIX(myvcf), getINFO(myvcf), matrix(NA, nrow=nrow(mydata.ad), ncol=ncol(mydata.ad)*4) ))

mymat <- foreach(i=1:ncol(mydata.ad), .combine='cbind') %dopar% {

    gq     <- as.numeric(mydata.gq[,i])
    ad     <- strsplit(mydata.ad[,i], ",")
    ad.ref <- lapply(ad, function(x) x[1]) %>% unlist() %>% as.numeric()
    ad.alt <- lapply(ad, function(x) x[2]) %>% unlist() %>% as.numeric()
    ad.tot <- ad.ref + ad.alt
    ad.frq <- ad.alt / ad.tot

    x <- cbind(gq, ad.alt, ad.tot, ad.frq)
    colnames(x) <- paste0(colnames(mydata.ad)[i], c("-GQ", "-nb.alt.allele.reads", "-nb.tot.allele.reads", "-freq.allele"))

    return(x)

}
mydata[,9:ncol(mydata)]  <- mymat 
colnames(mydata)[9:ncol(mydata)] <- colnames(mymat)

# Adjustments
mydata[,2] <- as.numeric(mydata[,2])
colnames(mydata)[8] <- "INFO"

mycolnames <- colnames(mydata)


#~~~~~~~~~~~~~~~~~~~~~~~~#
# Allele frequency stats #
#~~~~~~~~~~~~~~~~~~~~~~~~#

# Allele frequency replication (mean and sd)
cat("Generation of table from replicated libraries\n")
myaf.data <- data.frame(mydata[,1:8], matrix(nrow=nrow(mydata), ncol=length(mylib.rep)*2))
mycln <- 9:ncol(myaf.data)
myaf.data[, mycln] <- foreach (c=1:length(mycln), .combine='cbind') %dopar% {
    i <- mylib.rep[ ceiling(c/2) ]

    if (c %% 2 != 0) { # if c is odd
        x <- rowMeans(mydata[,grep(paste0(i,".*.freq"),mycolnames)])
    } else {
        x <- apply(mydata[,grep(paste0(i,".*.freq"),mycolnames)], 1, sd)
    }
    return(x)
}
colnames(myaf.data)[mycln] <- lapply(mylib.rep, function(x) paste(x, c("mean.freq","sd.freq"), sep=".")) %>% unlist()


# Allele frequency across all replicates

## Build primary table
myfreq.data <- mydata[,1:8]

### Run loop
mymat <- foreach(i=mylib.rep, .combine='cbind') %dopar% {
    
    # Method recalculating AF across all the samples
    mycln.nr <- grep(paste(i,".*.alt",sep=""), mycolnames)
    mycln.t  <- grep(paste(i,".*.tot",sep=""), mycolnames)

    my.nr <- mydata[, mycln.nr]
    my.nr <- suppressWarnings(as.data.frame(sapply(my.nr, as.numeric)))
    my.nr.sum <- rowSums(my.nr)

    my.t <- mydata[, mycln.t]
    my.t <- suppressWarnings(as.data.frame(sapply(my.t, as.numeric)))
    my.t.sum <- rowSums(my.t)
    

    x <- cbind(my.nr.sum, my.t.sum, my.nr.sum/my.t.sum)
    return(x)
}
colnames(mymat) <- lapply(mylib.rep, function(x) paste(x, c("alt.reads","tot.reads","freq"), sep=".")) %>% unlist()

myfreq.data <- cbind(myfreq.data,mymat)


# Genotype deduction
cat("Generation of genotypes\n")
mycln.REF <- grep("^REF$", colnames(myfreq.data))
mycln.ALT <- grep("^ALT$", colnames(myfreq.data))
if (length(mycln.REF) ==0 || (length(mycln.ALT) ==0)) {stop("REF and/or ALT columns missing.", call.=FALSE)}

Hm.REF <- 0.2
Hm.ALT <- 0.8
x <- paste(myfreq.data[,mycln.REF], myfreq.data[,mycln.ALT], sep="/")
mymat <- foreach(i=mylib, .combine='cbind') %dopar% {
    x[ which(myfreq.data[,paste0(i,".freq")] <= Hm.REF) ] <- paste(myfreq.data[ which(myfreq.data[,paste0(i,".freq")] <= Hm.REF) , mycln.REF], myfreq.data[ which(myfreq.data[,paste0(i,".freq")] <= Hm.REF) , mycln.REF], sep="/")
    x[ which(myfreq.data[,paste0(i,".freq")] >= Hm.ALT) ] <- paste(myfreq.data[ which(myfreq.data[,paste0(i,".freq")] >= Hm.ALT) , mycln.ALT], myfreq.data[ which(myfreq.data[,paste0(i,".freq")] >= Hm.ALT) , mycln.ALT], sep="/")
    return(x)
}
colnames(mymat) <- paste0("GT.", mylib)
myfreq.data <- cbind(myfreq.data, mymat)

# Allele frequency polarized toward PZQ-R origin
cat("Generation of polarized allele frequencies\n")
AF.trsh <- 0.5
mymat <- foreach(i=1:nrow(mypol.all), .combine='cbind') %dopar% {
    x <- myfreq.data[ ,paste0(mypol.all[i,1],".freq")]
    x[ which(myfreq.data[,paste0(mypol.all[i,2],".freq")] < AF.trsh) ] <- 1 - myfreq.data[ which(myfreq.data[,paste0(mypol.all[i,2],".freq")] < AF.trsh) , paste0(mypol.all[i,1],".freq") ]
    return(x)
}
colnames(mymat) <- paste0("pol.freq.", mypol.all[,1], "-", mypol.all[,2])
myfreq.data <- cbind(myfreq.data, mymat)

# Allele frequency subtraction
mymat <- foreach(i=1:nrow(mycomp), .combine='cbind') %dopar% {
    x <- myfreq.data[ ,grep(paste0("pol.freq.",mycomp[i,1]),colnames(myfreq.data))] - myfreq.data [ ,grep(paste0("pol.freq.",mycomp[i,2]),colnames(myfreq.data))]
}

colnames(mymat) <- paste0("sub.freq.", mycomp[,1], "-", mycomp[,2])
myfreq.data <- cbind(myfreq.data, mymat)


#~~~~~~~~~~~~~~~~#
# Data filtering #
#~~~~~~~~~~~~~~~~#

cat("Data filtering\n")


# Data filtering removing SNPs with high allele frequency variation within replicates
sd.trsh <- 0.10
sd.rows <- rowSums(myaf.data[, grep("*sd",colnames(myaf.data))] > sd.trsh) == 0
sd.rows[is.na(sd.rows)] <- FALSE

# Data filtering removing SNPs with low GQ
gq.trsh <- 0
gq.rows <- rowSums(mydata[, grep("*GQ",colnames(mydata))] < gq.trsh) == 0
gq.rows[is.na(gq.rows)] <- FALSE

# Data filtering removing SNPs with low read.depth
rd.trsh <- 10
rd.rows <- rowSums(mydata[, grep("*tot.allele.reads$",colnames(mydata))] < rd.trsh) == 0
rd.rows[is.na(rd.rows)] <- FALSE

myfreq.data.fltr <- myfreq.data[ apply( cbind(sd.rows, gq.rows, rd.rows), 1, all), ]


#~~~~~~~~#
# Report #
#~~~~~~~~#

# Summary
report <- c(
    paste("Nb of variable sites before filtering:",                         nb.var                                          ),
    paste("Nb of variable sites that are SNPs:",                            snp.rows                                        ),
    paste("Nb of variable sites that are biallelic:",                       bial.rows                                       ),
    paste("Nb of variable sites after filtering on sd between replicates:", sum(sd.rows)                                    ),
    paste("Nb of variable sites after filtering on read depth:",            (rowSums(cbind(rd.rows, sd.rows)) > 1) %>% sum()),
    paste("Nb of variable sites retained:",                                 nrow(myfreq.data.fltr)                          )
)

cat(" ", "Fitering report:", report, " ", sep="\n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Data processing for QTLseqr #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Recreate table compatible with the package
mydata2 <- as.data.frame(cbind( getFIX(myvcf)[, c(1, 2, 4, 5)],  matrix(NA, nrow=nrow(mydata.ad), ncol=ncol(mydata.ad)*3) ))

mymat2 <- foreach(i=1:ncol(mydata.ad), .combine='cbind') %dopar% {

    gq     <- as.numeric(mydata.gq[,i])
    ad     <- strsplit(mydata.ad[,i], ",")
    ad.ref <- lapply(ad, function(x) x[1]) %>% unlist() %>% as.numeric()
    ad.alt <- lapply(ad, function(x) x[2]) %>% unlist() %>% as.numeric()

    x <- cbind(ad.ref, ad.alt, gq)
    colnames(x) <- paste0(c("AD_REF.", "AD_ALT.", "GQ."), colnames(mydata.ad)[i])

    return(x)

}
mydata2[,5:ncol(mydata2)]  <- mymat2
colnames(mydata2)[5:ncol(mydata2)] <- colnames(mymat2)

mymat2 <- foreach(i=mylib.rep, .combine='cbind') %dopar% {
    
    # Method recalculating AF accross all the samples
    mycln.ref <- grep(paste0("AD_REF.", i), colnames(mydata2))
    mycln.alt <- grep(paste0("AD_ALT.", i), colnames(mydata2))
    mycln.gq <- grep(paste0("GQ.", i), colnames(mydata2))

    my.ref <- rowMeans(mydata2[, mycln.ref], na.rm = TRUE)
    my.alt <- rowMeans(mydata2[, mycln.alt], na.rm = TRUE)
    my.gq  <- rowMeans(mydata2[, mycln.gq ], na.rm = TRUE)

    x <- cbind(my.ref, my.alt, my.gq)
    return(x)
}
colnames(mymat2) <- lapply(mylib.rep, function(x) paste0(c("AD_REF.", "AD_ALT.", "GQ."), x)) %>% unlist()
mydata2 <- cbind(mydata2[, 1:4], mymat2)


# Combined
mygroup <- c("Recovered", "Contracted")
mymat3 <- foreach(i=mygroup, .combine='cbind') %dopar% {
    
    # Method recalculating AF across all the samples
    mycln.ref <- grep(paste0("AD_REF.*", i), colnames(mydata2))
    mycln.alt <- grep(paste0("AD_ALT.*", i), colnames(mydata2))
    mycln.gq <- grep(paste0("GQ.*", i), colnames(mydata2))

    my.ref <- rowMeans(mydata2[, mycln.ref], na.rm = TRUE)
    my.alt <- rowMeans(mydata2[, mycln.alt], na.rm = TRUE)
    my.gq  <- rowMeans(mydata2[, mycln.gq ], na.rm = TRUE)

    x <- cbind(my.ref, my.alt, my.gq)
    return(x)
}
colnames(mymat3) <- lapply(mygroup, function(x) paste0(c("AD_REF.", "AD_ALT.", "GQ."), x)) %>% unlist()
mydata3 <- cbind(mydata2[, 1:4], mymat3)


# Filter the data with the previous filter parameters
mydata2 <- mydata2[ apply( cbind(sd.rows, gq.rows, rd.rows), 1, all), ]
mydata3 <- mydata3[ apply( cbind(sd.rows, gq.rows, rd.rows), 1, all), ]

# Load the data with QTLseqr function using a connection (no disk I/O)
HighBulk <- grep("Recovered", colnames(mydata2), value = TRUE) %>% strsplit(., "\\.") %>% lapply(., function(x) x[2]) %>% unlist %>% unique()
LowBulk <- grep("Contracted", colnames(mydata2), value = TRUE) %>% strsplit(., "\\.") %>% lapply(., function(x) x[2]) %>% unlist %>% unique()
binary_csv <- csvConnection(mydata2)

for (i in 1:2) {
    x <- importFromTable(file = binary_csv, highBulk = HighBulk[i], lowBulk = LowBulk[i], chromList = paste0("SM_V7_", c(1:7, "ZW")))

    # Perform analysis
    ## G'
    y <- runGprimeAnalysis(x)

    # Save object
    assign(paste0("df", i), x)
    assign(paste0("df", i, ".q"), y)
}

# Combined
HighBulk <- grep("Recovered", colnames(mydata3), value = TRUE) %>% strsplit(., "\\.") %>% lapply(., function(x) x[2]) %>% unlist %>% unique()
LowBulk <- grep("Contracted", colnames(mydata3), value = TRUE) %>% strsplit(., "\\.") %>% lapply(., function(x) x[2]) %>% unlist %>% unique()
binary_csv <- csvConnection(mydata3)
df3 <- importFromTable(file = binary_csv, highBulk = HighBulk[1], lowBulk = LowBulk[1], chromList = paste0("SM_V7_", c(1:7, "ZW")))

# Perform analysis
## G'
df3.q <- runGprimeAnalysis(df3)



#=========#
# Figures #
#=========#

if (! dir.exists(graph_fd)) {dir.create(graph_fd, recursive=TRUE)}

for (i in 1:3) {
    y <- get(paste0("df", i, ".q"))
    p1 <- plotQTLStats(SNPset = y, var = "negLog10Pval")
    p2 <- plotQTLStats(SNPset = y, var = "Gprime", plotThreshold = TRUE, q = 0.01)
    p3 <- plotQTLStats(SNPset = y, var = "nSNPs")
    p4 <- plotGprimeDist(SNPset = y, outlierFilter = "Hampel", filterThreshold = 0.1)
    pdf(paste0(graph_fd, "QTLseqr ", i, ".pdf"), width = 15)
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    dev.off()
}
