#!/usr/bin/env Rscript
# Title: X-QTL_module.R
# Version: 1.3
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-04-13
# Modified in: 2023-04-27



#==========#
# Comments #
#==========#

# Module performing X-QTL analysis. The module must be called from the main script.



#==========#
# Versions #
#==========#

# v1.3 - 2023-04-27: update path to plotting function due to new genome version
# v1.2 - 2023-03-22: update plotting section due to new genome version
# v1.1 - 2021-08-26: save filtered table / remove unnecessary graph section
# v1.0 - 2021-08-24: rename script / turn it into a module / clean code
# v0.3 - 2019-05-15: functions split in another file
# v0.2 - 2018-05-01: update source data to use VCF file directly
# v0.1 - 2018-03-15: update email address / remove unnecessary parts
# v0.0 - 2015-04-13: creation



#======================#
# Packages and options #
#======================#

cat("Loading packages\n")

suppressMessages({
    library("rtracklayer")
    library("magrittr")
    library("vcfR")
    library("parallel")
    library("doParallel")
})

options(stringsAsFactors = FALSE)



#===========#
# Functions #
#===========#

source("X-QTL_func.R")

source("functions/coord_intersect.R")

#-------------------#
# Sm exome plotting #
#-------------------#

# Chromosome renaming
source("functions/rename_chr.R")

# Plotting function
## For v7 genome
source("functions/Sm.matplot.data.R")
data.order_v7 <- data.order
matplot.data_v7 <- matplot.data

## For v9 genome
source("functions/Sm.matplot.data_v5.0.R")



#===========#
# Variables #
#===========#

# Check files
if (missing(myvcf_file) || ! file.exists(myvcf_file)) {stop("File path in myvcf_file is missing or does not exist.", call.=FALSE)}
if (missing(mygff_file) || ! file.exists(mygff_file)) {stop("File path in mygff_file is missing or does not exist.", call.=FALSE)}
if (missing(myann_file) || ! file.exists(myann_file)) {stop("File path in myann_file is missing or does not exist.", call.=FALSE)}

# Check folders
if (missing(graph_fd) ) {stop("File path in graph_fd is missing or does not exist.",  call.=FALSE)}
if (missing(result_fd)) {stop("File path in result_fd is missing or does not exist.", call.=FALSE)}

# Check variables
if (missing(myexpr))     {stop("Variable myexpr is missing or does not exist.",     call.=FALSE)}
if (missing(myexpr.cln)) {stop("Variable myexpr.cln is missing or does not exist.", call.=FALSE)}
if (missing(myexpr.nm))  {stop("Variable myexpr.nm is missing or does not exist.",  call.=FALSE)}

if (missing(mylib))     {stop("Variable mylib is missing or does not exist.",     call.=FALSE)}
if (missing(mylib.ind)) {stop("Variable mylib.ind is missing or does not exist.", call.=FALSE)}
if (missing(mylib.rep)) {stop("Variable mylib.rep is missing or does not exist.", call.=FALSE)}

if (missing(mypol.all))  {stop("Variable mypol.all is missing or does not exist.",  call.=FALSE)}
if (missing(mycomp))     {stop("Variable mycomp is missing or does not exist.",     call.=FALSE)}
if (missing(mycomp.cb))  {stop("Variable mycomp.cb is missing or does not exist.",  call.=FALSE)}
if (missing(mypv.pat))   {stop("Variable mypv.pat is missing or does not exist.",   call.=FALSE)}
if (missing(myfrq.pat)) {stop("Variable myfreq.pat is missing or does not exist.", call.=FALSE)}

if (! exists("mysfx")) { mysfx <- NULL }
if (! is.null(mysfx)) {
    graph_fd  <- paste0(graph_fd, "/", mysfx, "/")
    result_fd <- paste0(result_fd, "/", mysfx, "/")
}


#=================#
# Data processing #
#=================#

# Set core cluster
registerDoParallel(cores = detectCores())

# Load data
myvcf <- read.vcfR(myvcf_file)
mygff <- readGFF(mygff_file)
myann <- read.csv(myann_file, header=TRUE, sep="\t")

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

# Read depth table
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


#------------------------#
# Allele frequency stats #
#------------------------#

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

## Run loop
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

# Genotype deduction in F2s or groups
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

# Allele frequency polarization
cat("Generation of polarized allele frequencies\n")
AF.trsh <- 0.5
mymat <- foreach(i=1:nrow(mypol.all), .combine='cbind') %dopar% {
    x <- myfreq.data[ ,paste0(mypol.all[i,1],".freq")]
    mypos <-  which(myfreq.data[,paste0(mypol.all[i,2],".freq")] > AF.trsh)
    x[ mypos ] <- 1 - myfreq.data[ mypos , paste0(mypol.all[i,1],".freq") ]
    return(x)
}
colnames(mymat) <- paste0("pol.freq.", mypol.all[,1], "-", mypol.all[,2])
myfreq.data <- cbind(myfreq.data, mymat)

# Allele frequency subtraction
mymat <- foreach(i=1:nrow(mycomp), .combine='cbind') %dopar% {
    y <- rep(1, nrow(myfreq.data))
    y[myfreq.data[,paste0(mycomp[i,1],".freq")] > AF.trsh] <- -1
    x <- (myfreq.data[ , paste0(mycomp[i,1], ".freq")] - myfreq.data [ , paste0(mycomp[i,2], ".freq")]) * y
}

colnames(mymat) <- paste0("sub.freq.", mycomp[,1], "-", mycomp[,2])
myfreq.data <- cbind(myfreq.data, mymat)


#---------------#
# Z-score stats #
#---------------#

# Z-score calculation
cat("Z-score calulation\n")

mymat <- foreach(i=1:nrow(mycomp), .combine='cbind') %dopar% {

    myctrl.t.r <- myfreq.data[ ,grep(paste(mycomp[i,2],"tot",sep="."), colnames(myfreq.data))]
    myctrl.af  <- myfreq.data[ ,grep(paste(mycomp[i,2],"freq",sep="."), colnames(myfreq.data))]

    mytreat.t.r <- myfreq.data[ ,grep(paste(mycomp[i,1],"tot",sep="."), colnames(myfreq.data))]
    mytreat.af  <- myfreq.data[ ,grep(paste(mycomp[i,1],"freq",sep="."), colnames(myfreq.data))]

    myctrl.ind  <- mylib.ind[grep(mycomp[i,2],mylib)]
    mytreat.ind <- mylib.ind[grep(mycomp[i,1],mylib)]

    # Formula:
    x1 <- (mytreat.af - myctrl.af) / sqrt( rowMeans(cbind(mytreat.af,myctrl.af)) * (1-rowMeans(cbind(mytreat.af,myctrl.af))) * ( 1/(2*myctrl.ind) + 1/(2*mytreat.ind) + 1/myctrl.t.r + 1/mytreat.t.r ) )

    # p-value calculation
    x2 <- pnorm(-abs(x1))

    x <- cbind(x1, x2)
    colnames(x) <- c(paste0("z.",mycomp[i,1],"-",mycomp[i,2]), paste0("pv.",mycomp[i,1],"-",mycomp[i,2]))

    return(x)

}
myfreq.data <- cbind(myfreq.data, mymat)

# Z-scores combination
cat("Z-score combination\n")
for (i in 1:nrow(mycomp.cb)) {
    x <- myfreq.data[ ,grep(paste("z.",mycomp.cb[i,1],sep=""),colnames(myfreq.data))]
    y <- myfreq.data[ ,grep(paste("z.",mycomp.cb[i,2],sep=""),colnames(myfreq.data))]
    myfreq.data[,ncol(myfreq.data)+1] <- combine.z.pv(x,y)
    colnames(myfreq.data)[ncol(myfreq.data)] <- paste("pv.",mycomp.cb[i,1],"-",mycomp.cb[i,2],sep="")
}


#----------------#
# Data filtering #
#----------------#

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


#--------#
# Report #
#--------#

if (! dir.exists(result_fd)) {dir.create(result_fd, recursive=TRUE)}
myext <- paste0("pv.sd-", sd.trsh, ".gq-", gq.trsh, ".rd-", rd.trsh)

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
write(report, paste0(result_fd, "snp_report"))

# Table
pv.vec <- grep(paste0("pv", mypv.pat), colnames(myfreq.data.fltr))
pv.print  <- grep("pv.*", colnames(myfreq.data.fltr))
gt.vec <- grep("GT", colnames(myfreq.data.fltr))
frq.vec <- grep(paste0("pol.freq", myfrq.pat), colnames(myfreq.data.fltr))

bf.cor <- 0.05/nrow(myfreq.data.fltr[is.finite(rowSums(as.data.frame(myfreq.data.fltr[,pv.vec]))), ])

pk.rpt(myfreq.data.fltr[is.finite(rowSums(as.data.frame(myfreq.data.fltr[,pv.vec]))),pv.vec], "<=", bf.cor, cln.print=c(1:8,gt.vec,pv.print), cln.nm.ext=myext, res.folder=result_fd)

## Without Bonferroni correction
pk.rpt(myfreq.data.fltr[is.finite(rowSums(as.data.frame(myfreq.data.fltr[,pv.vec]))),pv.vec], ">=", 0, cln.print=c(1:8,gt.vec,pv.print), cln.nm.ext=paste0(myext, "-no_bf"), res.folder=result_fd)

# Saving object
save(myfreq.data.fltr, file = paste0(result_fd, "myfreq.data.fltr.RData"))


#----------------#
# QTL boundaries #
#----------------#

cat("Delimiting QTL boundaries\n")

myfreq.data.fltr.clean <- myfreq.data.fltr[ ! is.na(myfreq.data.fltr[,pv.vec]), ]
myfreq.data.fltr.bf    <- myfreq.data.fltr.clean[ myfreq.data.fltr.clean[, pv.vec] <= bf.cor, ]
mychr <- (sapply(myfreq.data.fltr.bf[,1] %>% unique(), function(x) (myfreq.data.fltr.bf[,1] == x) %>% sum()) > nrow(myfreq.data.fltr.bf) * 0.1) %>% which() %>% names() %>% sort

bf.lim <- as.data.frame(matrix(rep(NA, 3 * length(mychr)), nrow = length(mychr)))

for (c in mychr) {
    my.qtl <- na.omit(myfreq.data.fltr[ myfreq.data.fltr[,pv.vec] <= bf.cor & myfreq.data.fltr[,1] == c, 2])

    # Determining boundaries of the QTL
    ## significant p-values must be withing 500,000 bp to be in the same QTL
    mystretch <- ( (diff(my.qtl) < 500000)) %>% split(., cumsum(c(TRUE, diff(.) < 0)))
    lg.stretch <- mystretch %>% lapply(., sum) %>% which.max()

    if (lg.stretch == 1) {
        mystretch.start <- 1
        mystretch.end   <- length(mystretch[[1]])
    } else {
        mystretch.start <- lapply(mystretch[1:(lg.stretch-1)], length) %>% unlist %>% sum() +  sum( (! mystretch[[lg.stretch]]) %>% as.integer()) + 1
        mystretch.end   <- lapply(mystretch[1:lg.stretch], length) %>% unlist %>% sum()
    }

    # Store QTL (Bonferroni) boundaries
    bf.lim[ match(c, mychr), 1]  <- c
    bf.lim[ match(c, mychr), -1] <- c(my.qtl[mystretch.start], my.qtl[mystretch.end])
}

# Write bed file format of the boundaries
write.table(bf.lim, file = paste0(result_fd, "QTL.bed"), col.names = FALSE, row.names = FALSE, quote = FALSE)


#----------------#
# Mutation table #
#----------------#

cat("Generating mutation and gene tables related to QTL\n")

mygff.genes <- mygff[ mygff[,3] == "gene", ]

# Reorder chromosomes
chr.order <- sort(levels(mygff.genes[,1]))
mygff.genes[,1] <- ordered(mygff.genes[,1], chr.order)
mygff.genes <- mygff.genes[order(mygff.genes[,1]), ]

levels(mygff.genes[,1]) <- sort(levels(mygff.genes[,1]))
levels(mygff.genes[,1]) <- sort(levels(mygff.genes[,1]))
mygff.genes <- mygff.genes[order(mygff.genes[,1], mygff.genes[,4]),]

# Expression column
if(! exists("myexpr.nm")) { myexpr.nm <- colnames(myexpr)[myexpr.cln] }

myfreq.data.bed <- coord_intersect(myfreq.data, bf.lim)
myfreq.data.fltr.bed <- coord_intersect(myfreq.data.fltr, bf.lim)

# Select combined Z scores
mypv <- pv.vec[ colnames(myfreq.data)[pv.vec] %>% nchar() %>% which.max() ]

if (length(mypv) == 1) { myrows <- myfreq.data.fltr[,mypv] <= bf.cor } else { myrows <- rowSums(myfreq.data.fltr[,mypv] <= bf.cor) > 0 }

# Generating tables
for (i in 1:length(mychr)) {
    mytb <- qtl.tb(x=myfreq.data.bed[[i]], x.fltr=myfreq.data.fltr.bed[[i]], chr=mychr[[i]], pv.cln=mypv, bf.cor, ann.tb=myann, expr.tb=myexpr, expr.cln=myexpr.cln, expr.nm=myexpr.nm, cln.print=c(gt.vec,frq.vec,pv.print), res.folder=result_fd)

    if (! is.null(mytb)) {
        qtl.tb.sum(mytb, ann.tb=myann, expr.tb=myexpr, expr.cln=myexpr.cln, expr.nm=myexpr.nm, gff=mygff.genes, baits.bed=mybaits, res.folder=result_fd)
    }
}



#=========#
# Figures #
#=========#

cat("Generating figures\n")

if (! dir.exists(graph_fd)) {dir.create(graph_fd, recursive=TRUE)}

# Replicate in allele frequency plot
myaf.plot(mydata, mylib.rep, myfn, graph.folder=graph_fd)

# Allele frequency plot from filtered data
vec <- grep("Sm", mylib.rep, value=TRUE)
myaf.plot(mydata[sd.rows,], vec, myfn, ext=".freq-fltr-sd.png", graph.folder=graph_fd)

# Z-score distribution plot
z.vec <- grep("z", colnames(myfreq.data.fltr), value=T)
myz.qq(myfreq.data.fltr, z.vec, myfn, graph.folder=graph_fd)

if (grepl("SM_V7_", myfreq.data.fltr[1,1])) {
    myfreq.data.fltr.r <- rename_chr_SmV7(myfreq.data.fltr, 1)
    myfreq.data.fltr.r <- data.order_v7(myfreq.data.fltr.r)
    g_vers <- 7
} else {
    myfreq.data.fltr.r <- myfreq.data.fltr
    g_vers <- 10
}

# P-value plot
runmed.vec <- c(NA, 5,11,21)
mydef <- "HD"

pv.vec <- grep("pv", colnames(myfreq.data.fltr.r))

f.vec <- grep("freq", colnames(myfreq.data.fltr.r))

sf.vec <- grep("sub.freq", colnames(myfreq.data.fltr.r))

pf.vec <- grep("pol.freq", colnames(myfreq.data.fltr.r))

myfreq.data.fltr.r <- myfreq.data.fltr.r[ is.finite(rowSums(myfreq.data.fltr.r[,pv.vec])), ]

for (j in runmed.vec) {
    mygraph(myfreq.data.fltr.r, pv.vec, j, "pvalue", myfn, def=mydef, graph.folder=graph_fd, version=g_vers)

#   mygraph(myfreq.data.fltr.r, f.vec, j, "freq", myfn, graph.folder=graph_fd, version=g_vers)

    mygraph(myfreq.data.fltr.r, sf.vec, j, "freq", myfn, -0.75, 0.75, ylab="Difference in allele frequency", def=mydef, graph.folder=graph_fd, version=g_vers)

    mygraph(myfreq.data.fltr.r, pf.vec, j, "freq", myfn, def=mydef, graph.folder=graph_fd, version=g_vers)
}
