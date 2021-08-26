#!/usr/bin/env Rscript
# Title: PZQ_analysis
# Version: 0.3
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-05-15
# Modified in: 2018-05-01



#==========#
# Comments #
#==========#

# v0.3 - 2019-05-15: functions split in another file
# v0.2 - 2018-05-01: update source data to use VCF file directly
# v0.1 - 2018-03-15: update email address / remove unnecessary parts



#======================#
# Packages and options #
#======================#

library("rtracklayer")
library("magrittr")
library("vcfR")
library("parallel")
library("doParallel")

options(stringsAsFactors = FALSE)



#===========#
# Functions #
#===========#

source("PZQ_analysis_func.R")

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

source("PZQ_analysis_config.R")



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


# Build read depth table
mydata.ad     <- extract.gt(myvcf, "AD")
mydata.gq     <- extract.gt(myvcf, "GQ")


mydata <- as.data.frame(cbind( getFIX(myvcf), getINFO(myvcf), matrix(NA, nrow=nrow(mydata.ad), ncol=ncol(mydata.ad)*4) ))

#mydata[,9:ncol(mydata)]  <- foreach(i=1:ncol(mydata.ad), .combine='cbind') %dopar% {
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


# Allele frequency accross all replicates

## Build primary table
myfreq.data <- mydata[,1:8]

### Run loop
mymat <- foreach(i=mylib.rep, .combine='cbind') %dopar% {
    
    # Method recalculating AF accross all the samples
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

# Genotype deduction in F2s
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
    mypos <-  which(myfreq.data[,paste0(mypol.all[i,2],".freq")] > AF.trsh)
    x[ mypos ] <- 1 - myfreq.data[ mypos , paste0(mypol.all[i,1],".freq") ]
    return(x)
}
colnames(mymat) <- paste0("pol.freq.", mypol.all[,1], "-", mypol.all[,2])
myfreq.data <- cbind(myfreq.data, mymat)

# Allele frequency subtraction
# mymat <- foreach(i=1:nrow(mycomp), .combine='cbind') %dopar% {
#     x <- myfreq.data[ ,grep(paste0("pol.freq.",mycomp[i,1]),colnames(myfreq.data))] - myfreq.data [ ,grep(paste0("pol.freq.",mycomp[i,2]),colnames(myfreq.data))]
# }
mymat <- foreach(i=1:nrow(mycomp), .combine='cbind') %dopar% {
    y <- rep(1, nrow(myfreq.data))
    y[myfreq.data[,paste0(mycomp[i,1],".freq")] > AF.trsh] <- -1
    x <- (myfreq.data[ , paste0(mycomp[i,1], ".freq")] - myfreq.data [ , paste0(mycomp[i,2], ".freq")]) * y
}

colnames(mymat) <- paste0("sub.freq.", mycomp[,1], "-", mycomp[,2])
myfreq.data <- cbind(myfreq.data, mymat)


#~~~~~~~~~~~~~~~#
# Z-score stats #
#~~~~~~~~~~~~~~~#

# Z-score calculation
cat("Z-score calulation\n")
## Original formula:
## mydata[i,mycol.6] <- (mydata[i,LGMT.AF]-mydata[i,LGMNT.AF])/sqrt(mean(c(mydata[i,LGMT.AF],mydata[i,LGMNT.AF]))*(1-mean(c(mydata[i,LGMT.AF],mydata[i,LGMNT.AF])))*(1/(2*147) + 1/(2*16) + 1/(mydata[i,LGMNT.R.alt]+mydata[i,LGMNT.R.ref]) + 1/(mydata[i,LGMT.R.alt]+mydata[i,LGMT.R.ref])))

mymat <- foreach(i=1:nrow(mycomp), .combine='cbind') %dopar% {
    
#    myctrl.nr.r <- myfreq.data[ ,grep(paste(mycomp[i,2],"alt",sep="."), colnames(myfreq.data))]
    myctrl.t.r <- myfreq.data[ ,grep(paste(mycomp[i,2],"tot",sep="."), colnames(myfreq.data))] 
    myctrl.af  <- myfreq.data[ ,grep(paste(mycomp[i,2],"freq",sep="."), colnames(myfreq.data))]
    
#   mytreat.nr.r <- myfreq.data[ ,grep(paste(mycomp[i,1],"alt",sep="."), colnames(myfreq.data))]
    mytreat.t.r <- myfreq.data[ ,grep(paste(mycomp[i,1],"tot",sep="."), colnames(myfreq.data))] 
    mytreat.af  <- myfreq.data[ ,grep(paste(mycomp[i,1],"freq",sep="."), colnames(myfreq.data))]

    myctrl.ind  <- mylib.ind[grep(mycomp[i,2],mylib)]
    mytreat.ind <- mylib.ind[grep(mycomp[i,1],mylib)]

    # Forumla:
#    myfreq.data.col <- ncol(myfreq.data)+1
#    for (j in 1:length(my.ctrl.af)) {
#        myfreq.data[j,myfreq.data.col] <- (mytreat.af[j] - myctrl.af[j]) / sqrt( mean(c(mytreat.af[j], myctrl.af[j])) * (1-mean(c(mytreat.af[j], myctrl.af[j]))) * ( 1/(2*myctrl.ind) + 1/(2*mytreat.ind) + 1/(myctrl.nr.r[j]+myctrl.r.r[j]) + 1/(mytreat.nr.r[j]+mytreat.r.r[j]) ) )
#    }
    
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

#myfreq.data.fltr <- myfreq.data[ apply( cbind(snp.rows, bial.rows, sd.rows, gq.rows, rd.rows), 1, all), ]
myfreq.data.fltr <- myfreq.data[ apply( cbind(sd.rows, gq.rows, rd.rows), 1, all), ]


#~~~~~~~~#
# Report #
#~~~~~~~~#

res.folder <- paste0(result_fd, "/sd-", sd.trsh, ".gq-", gq.trsh, ".rd-", rd.trsh, "/")
if (file.exists(res.folder) == FALSE) {dir.create(res.folder, recursive=TRUE)}

# Summary
#report <- c(
#    paste("Nb of variable sites before filtering:",                         nb.var                                                                     ),
#    paste("Nb of variable sites that are SNPs:",                            sum(snp.rows)                                                              ),
#    paste("Nb of variable sites that are biallelic:",                       sum(apply( cbind(snp.rows, bial.rows), 1, all))                            ),
#    paste("Nb of variable sites after filtering on sd between replicates:", sum(apply( cbind(snp.rows, bial.rows, sd.rows), 1, all))                   ),
#    paste("Nb of variable sites after filtering on GQ:",                    sum(apply( cbind(snp.rows, bial.rows, sd.rows, gq.rows), 1, all))          ),
#    paste("Nb of variable sites after filtering on read depth:",            sum(apply( cbind(snp.rows, bial.rows, sd.rows, gq.rows, rd.rows), 1, all)) )  #,
#)
report <- c(
    paste("Nb of variable sites before filtering:",                         nb.var      ),
    paste("Nb of variable sites that are SNPs:",                            snp.rows    ),
    paste("Nb of variable sites that are biallelic:",                       bial.rows   ),
    paste("Nb of variable sites after filtering on sd between replicates:", sum(sd.rows)),
    paste("Nb of variable sites after filtering on read depth:",            sum(rd.rows)),
    paste("Nb of variable sites retained:",                       nrow(myfreq.data.fltr))  #,
)

cat(" ", report, " ", sep="\n")
write(report, paste(res.folder,"report",sep=""))

# Table 
pv.vec <- grep(paste0("pv", mypv.pat), colnames(myfreq.data.fltr))
pv.print  <- grep("pv.*", colnames(myfreq.data.fltr))
gt.vec <- grep("GT", colnames(myfreq.data.fltr))
frq.vec <- grep(paste0("pol.freq", myfrq.pat), colnames(myfreq.data.fltr))


bf.cor <- 0.05/nrow(myfreq.data.fltr[is.finite(rowSums(as.data.frame(myfreq.data.fltr[,pv.vec]))), ])

pk.rpt(myfreq.data.fltr[is.finite(rowSums(as.data.frame(myfreq.data.fltr[,pv.vec]))),pv.vec], "<=", bf.cor, cln.print=c(1:8,gt.vec,pv.print), cln.nm.ext=paste0("AB4-freq-pv.sd-",sd.trsh), res.folder=res.folder)

## Without bonferroni correction
pk.rpt(myfreq.data.fltr[is.finite(rowSums(as.data.frame(myfreq.data.fltr[,pv.vec]))),pv.vec], ">=", 0, cln.print=c(1:8,gt.vec,pv.print), cln.nm.ext=paste0("AB4-freq-pv.sd-",sd.trsh,"-bf=0"), res.folder=res.folder)
#pk.rpt(myfreq.data.fltr[is.finite(rowSums(myfreq.data.fltr[,pv.vec])),pv.vec], ">=", 0, cln.print=c(1:8,gt.vec,pv.vec,frq.vec), cln.nm.ext=paste("AB4-freq-pv.sd-",sd.trsh,"-bf=0",sep=""), res.folder=res.folder)




a <- myfreq.data.fltr[ ! is.na(myfreq.data.fltr[,pv.vec]), ]
b <- a[ a[, pv.vec] <= bf.cor, ]
mychr <- (sapply(b[,1] %>% unique(), function(x) (b[,1] == x) %>% sum()) > nrow(b) * 0.1) %>% which() %>% names() %>% sort

bf.lim <- as.data.frame(matrix(rep(NA, 3 * length(mychr)), nrow = length(mychr)))

for (c in mychr) {
    my.qtl <- na.omit(myfreq.data.fltr[ myfreq.data.fltr[,pv.vec] <= bf.cor & myfreq.data.fltr[,1] == c, 2])

    # Determining boundaries of the QTL
    mystretch <- ( (diff(my.qtl) < 500000)) %>% split(., cumsum(c(TRUE, diff(.) < 0)))
    lg.stretch <- mystretch %>% lapply(., sum) %>% which.max()

    if (lg.stretch == 1) {
        mystretch.start <- 1
        mystretch.end   <- length(mystretch[[1]])
    } else {
        mystretch.start <- lapply(mystretch[1:(lg.stretch-1)], length) %>% unlist %>% sum() +  sum( (! mystretch[[lg.stretch]]) %>% as.integer()) + 1
        mystretch.end   <- lapply(mystretch[1:lg.stretch], length) %>% unlist %>% sum()
    }


    bf.lim[ match(c, mychr), ] <- c(c, my.qtl[mystretch.start], my.qtl[mystretch.end])
}


write.table(bf.lim, file = paste0(result_fd, "QTL.bed"), col.names = FALSE, row.names = FALSE, quote = FALSE)

#
# GFF
#

mygff.genes <- mygff[ mygff[,3] == "gene", ]

# Reorder chromosomes
chr.order <- sort(levels(mygff.genes[,1]))
mygff.genes[,1] <- ordered(mygff.genes[,1], chr.order)
mygff.genes <- mygff.genes[order(mygff.genes[,1]), ]

levels(mygff.genes[,1]) <- sort(levels(mygff.genes[,1]))
levels(mygff.genes[,1]) <- sort(levels(mygff.genes[,1]))
mygff.genes <- mygff.genes[order(mygff.genes[,1], mygff.genes[,4]),]





#myann.tb <- myfreq.data.fltr
#myann.tb[,ncol(myann.tb)+1] <- NA

#myann.tb <- rep(NA, nrow(myfreq.data.fltr))
#
#for (i in 1:nrow(myfreq.data.fltr)) {
#    mychr <- myfreq.data.fltr[i,1]
#    mypos <- myfreq.data.fltr[i,2]
#    mygff.genes.sub <- mygff.genes[ mygff.genes[,1] == mychr, ]
##    myann.tb[i, ncol(myann.tb)] <- apply(mygff.genes[ mygff.genes[,1] == mychr,4:5], 1, function(x) findInterval(mypos, as.vector(x)) == 1 )
#    myrows <- apply(mygff.genes.sub[,4:5], 1, function(x) findInterval(mypos, as.vector(x)) == 1 )
#    if (any(myrows)) {
#        if (sum(myrows) == 1) {
#            myann.tb[i] <- mygff.genes.sub[myrows, 9]
#        } else {
#            myann.tb[i] <- mygff.genes.sub[myrows, 9][1]
#            warning("The variant ", i, " belongs to ", sum(myrows), "genes. Only the first one is kept. The GFF file is probably malformated.")
#        }
#    }
#}


#--------------#
# Detail table #
#--------------#

if(! exists("myexpr.nm")) { myexpr.nm <- colnames(myexpr)[myexpr.cln] }

k <- length(myexpr.cln)

#myinfo <- strsplit(myfreq.data.fltr[,8], ";", fixed=T)
myinfo <- strsplit(myfreq.data[,8], ";", fixed=T)

myinfo.lg <- 1:length(myinfo)
n <- detectCores()-1
myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))

#myinfo.ann <- lapply(myinfo, function(x) unlist(grep("ANN",x, value=T) %>% strsplit(., "|", fixed=T)))
myinfo.ann <- foreach(i=1:length(myinfo.lg), .combine='c') %dopar% {
    lapply(myinfo[myinfo.lg[[i]]], function(x) unlist(grep("ANN",x, value=T) %>% strsplit(., "|", fixed=T)))
}

mygenes.ls <- vector("list", length(myinfo.ann))
cln.nm <- c("Gene", "GFF_annotation", "HHPred_annotation", myexpr.nm)
mymain.genes.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(myinfo.ann)))
colnames(mymain.genes.tb) <- cln.nm
#for (i in 1:length(myinfo.ann)) {
#    x <- myinfo.ann[[i]]
#    todo <- TRUE
#
#    mycln.tmp <- seq(0, by=15, length.out=round(length(x)/15))
#    
#    if (length(mycln.tmp) > 0) {
#        mymtx <- sapply(mycln.tmp, function(y) x[y+c(2,3,4,7,8,10,11)])
#
#        mymain.gene <- names(sort(table(mymtx[3,]), decreasing=T)[1]) 
#   
#        myRNA <- mymtx[4,]  
#        myRNA.tb <- matrix(NA, nrow=2+length(myexpr.cln), ncol=ncol(mymtx))
#
#        for (r in 1:length(myRNA)) {
#            
#            # Corresponding gene
#            mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
#            if (is.na(mygene)) next 
#            
#            # Transcript data
#            myann.tmp <- myann[ myann[,1] == myRNA[r] , ]
#            # If no transcript detected
#            if (length(myann.tmp) == 0) { myann.tmp <- myann[ myann[,1] == mygene , ] }
#            # If really nothing detected
#            if (dim(myann.tmp)[1] == 0) { myann.tmp <- matrix(NA, ncol=3, nrow=1) }
#            # Add data to the matrix
#            myRNA.tb[1:2, r] <- t(as.matrix(myann.tmp[,2:3]))
#
#            # Expression data
#            myexpr.tmp <- myexpr[ myexpr[,1] == myRNA[r] , myexpr.cln ]
#            # If no expression found
#            if (length(myexpr.tmp) == 0 | dim(myexpr.tmp)[1] == 0) { myexpr.tmp <- rep(NA, length(myexpr.cln)) }
#            # Add data to matrix
#            myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(myexpr.tmp)
#
#            # Add info for main gene table
#            if (mygene == mymain.gene & isTRUE(todo)) {
#                mymain.genes.tb[i,] <- cbind(mygene, t(myRNA.tb[,r]))
#
#                # If expression is NA, is there any other expression data for this given gene?
#                if (all(is.na(myexpr.tmp)) | sum(myexpr.tmp) == 0) {
#                    #myexpr.tmp <- myexpr[ grep(myRNA[r], myexpr[,1]) , myexpr.cln ]
#                    myexpr.tmp <- myexpr[ grep(paste0(mygene,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ]
#                    #if (! all(is.na(myexpr.tmp))) { mymain.genes.tb[i,4] <- max(myexpr[ grep(paste0(mygene,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ], na.rm=TRUE) }
#                    if (! all(is.na(myexpr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- apply(myexpr[ grep(paste0(mygene,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
#                }
#
#                # Don't do anything else for this gene in the future
#                todo <- FALSE
#            }
#        }
#
#        # Combine matrix
#        mymtx <- rbind(mymtx, myRNA.tb)
#        
#        # Reorder rows
#        mymtx <- mymtx[c(3, 8:(9+k), 4, 1:2, 5:7), ]
#        
#        mygenes.ls[[i]] <- mymtx
#    }
#    
#}


myinfo.lg <- 1:length(myinfo)
n <- round(20*log(length(myinfo))) 
myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))

#mygenes.ls <- foreach(j=1:length(myinfo.lg), .combine='c') %dopar% {
mygenes.ls2 <- foreach(j=1:length(myinfo.lg), .combine='c') %dopar% {

    myidx <- myinfo.lg[[j]]

    mygenes.ls.tmp <- mygenes.ls[myidx]

    for (i in 1:length(myidx)) {
        x <- myinfo.ann[[myidx[i]]]
        todo <- TRUE

        mycln.tmp <- seq(0, by=15, length.out=round(length(x)/15))
        
        if (length(mycln.tmp) > 0) {
            mymtx <- sapply(mycln.tmp, function(y) x[y+c(2,3,4,7,8,10,11)])

            mymain.gene <- names(sort(table(mymtx[3,]), decreasing=T)[1]) 
       
            myRNA <- mymtx[4,]  
            myRNA.tb <- matrix(NA, nrow=2+length(myexpr.cln), ncol=ncol(mymtx))

            for (r in 1:length(myRNA)) {
                
                # Corresponding gene
                mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
                if (is.na(mygene)) next 
                
                # Transcript data
                myann.tmp <- myann[ myann[,1] == myRNA[r] , ]
                # If no transcript detected
                if (length(myann.tmp) == 0) { myann.tmp <- myann[ myann[,1] == mygene , ] }
                # If really nothing detected
                if (dim(myann.tmp)[1] == 0) { myann.tmp <- matrix(NA, ncol=3, nrow=1) }
                # Add data to the matrix
                myRNA.tb[1:2, r] <- t(as.matrix(myann.tmp[,2:3]))

                # Expression data
                myexpr.tmp <- myexpr[ myexpr[,1] == myRNA[r] , myexpr.cln ]
                # If no expression found
                if (length(myexpr.tmp) == 0 | dim(myexpr.tmp)[1] == 0) { myexpr.tmp <- rep(NA, length(myexpr.cln)) }
                # Add data to matrix
                myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(myexpr.tmp)

                # Add info for main gene table
                if (mygene == mymain.gene & isTRUE(todo)) {
                    mymain.genes.tb[i,] <- cbind(mygene, t(myRNA.tb[,r]))

                    # If expression is NA, is there any other expression data for this given gene?
                    if (all(is.na(myexpr.tmp)) | sum(myexpr.tmp) == 0) {
                        #myexpr.tmp <- myexpr[ grep(myRNA[r], myexpr[,1]) , myexpr.cln ]
                        myexpr.tmp <- myexpr[ grep(paste0(mygene,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ]
                        #if (! all(is.na(myexpr.tmp))) { mymain.genes.tb[i,4] <- max(myexpr[ grep(paste0(mygene,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ], na.rm=TRUE) }
                        if (! all(is.na(myexpr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- apply(myexpr[ grep(paste0(mygene,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
                    }

                    # Don't do anything else for this gene in the future
                    todo <- FALSE
                }
            }

            # Combine matrix
            mymtx <- rbind(mymtx, myRNA.tb)
            
            # Reorder rows
            mymtx <- mymtx[c(3, 8:(9+k), 4, 1:2, 5:7), ]
            
            mygenes.ls.tmp[[i]] <- mymtx
        }
        
    }

    return(mygenes.ls.tmp)
}

mygenes.ls[ sapply(mygenes.ls, is.null) ] <- NA
mymax      <- max(sapply(mygenes.ls, length))
mygenes.tb <- lapply(mygenes.ls, function(x) x[1:mymax])
mygenes.tb <- t(do.call(cbind, mygenes.tb))
mygenes.tb[ mygenes.tb == "" ] <- NA
mygenes.tb <- as.data.frame(mygenes.tb)

cln.nm.b <- c("Gene (hit ", "GFF annotation (hit ", "HHPred annotation (hit ", paste(myexpr.nm, "(hit "), "Transcript (hit ", "Region (hit ", "Impact (hit ", "Type (hit ", "DNA mutation (hit ", "Protein mutation (hit ")
cln.nm <- NULL
for (i in 1:(mymax/length(cln.nm.b))) { cln.nm <- c(cln.nm, paste0(cln.nm.b, i, ")")) }
colnames(mygenes.tb) <- cln.nm

# Final table
#myfinal.tb <- cbind(myfreq.data.fltr[,c(1,2,4,5)], mymain.genes.tb, myfreq.data.fltr[,c(gt.vec,frq.vec,pv.print)], mygenes.tb)
myfinal.tb <- cbind(myfreq.data[,c(1,2,4,5)], mymain.genes.tb, myfreq.data[,c(gt.vec,frq.vec,pv.print)], mygenes.tb)

## Select rows in the QTL region
myrows <- na.omit(rownames(myfinal.tb)[myfinal.tb[,1] == "Chr_3" & myfinal.tb[,2] >= bf.lim[1] & myfinal.tb[,2] <= bf.lim[2]])
myfinal.tb.qtl <- myfinal.tb[ myrows, ]

## Select only informative columns then
myclns <- ! apply(myfinal.tb.qtl, 2, function(x) all(is.na(x)))
myfinal.tb.qtl <- myfinal.tb.qtl[ , myclns ]

write.table(myfinal.tb.qtl, paste0(res.folder,"/QTL_table_details.tsv"), row.names=FALSE, quote=FALSE, sep="\t")


##---------------#
## Summary table #
##---------------#
#
#mygff.genes.qtl <- mygff.genes[mygff.genes[,1] == "Chr_3" & mygff.genes[,5] >= bf.lim[1] & mygff.genes[,4] <= bf.lim[2],]
#
#mygenes.occ <- unique(mygff.genes.qtl[,9])
#mygenes.occ <- na.omit(unique(c(mygenes.occ, myfinal.tb.qtl[,5])))
#
#myimpct.cln <- grep("Impact", colnames(myfinal.tb))
#myimpct.tb  <- as.data.frame(matrix(c("LOW", "MODERATE", "HIGH", "MODIFIER", 0.3, 0.6, 0.9, 0.3), ncol=2))
#
#myfreq.cln  <- grep("pol.freq",colnames(myfinal.tb))
#
#cln.nm   <- c("Chr.", "Gene ID", "v5 gene ID", "Start", "End", "GFF annotation", "HHPred annotation", myexpr.nm, "Captured", "Nb of variable sites", "Impact score", "Weighted impact score", "Mean polarized all. freq.", "Median polarized all. freq.", "Highest polarized all. freq.", "Mean p-value", "Median p-value", "Lowest p-value", "Mean -log10(p-value)", "Median -log10(p-value)", "Highest -log10(p-value)")
#mysum.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(mygenes.occ)))
#for (g in mygenes.occ) {
#    myidx <- match(g, mygenes.occ)
#    mypv  <- rev(grep("pv.", colnames(myfinal.tb), fixed=TRUE))[1]
#
#    if (! is.null(mybaits)) {
#        # Is the gene in the capture array?
#        g.crdt <- as.data.frame(mygff.genes[ mygff.genes[,9] == g, c(1,4:5) ])
#        if (dim(g.crdt)[1] > 0) {
#            capture <- any( findInterval(mybaits[ mybaits[,1] == g.crdt[,1], 3], g.crdt[,2:3]) == 1 )
#        } else {
#            capture <- FALSE
#        }
#        if (capture) { capture <- "yes" } else { capture <- "no" }
#    } else {
#        capture <- NA
#    }
#    
#    # Is there any GFF v5 ID(s) associated to the current gene?
#    old.id <- as.data.frame(mygff.genes[ mygff.genes[,9] == g, 15 ])[,3]
#    if (length(old.id) == 0) { 
#        old.id <- NA 
#    } else {
#        old.id <- paste(old.id, collapse=", ")
#    }
#    
#    # Pull out information if the gene is in the variant data table
#    if (any(myfinal.tb[,5] == g, na.rm=TRUE)) {
#        mytb.tmp <- myfinal.tb[ ! is.na(myfinal.tb[,5]) & myfinal.tb[,5] == g, ]
#    
#        mypos <- as.data.frame(mygff.genes[ mygff.genes[,9] == g, c(4,5) ])
#        if (dim(mypos)[1] == 0) { mypos <- c(NA,NA) }
#
#        nb.var <- nrow(mytb.tmp)
#
#        mytb.impct <- mytb.tmp[,myimpct.cln]
#        for (i in 1:nrow(myimpct.tb)) { mytb.impct[ mytb.impct == myimpct.tb[i,1] ] <- myimpct.tb[i,2] }
#        myimpct.vec <- na.omit(as.numeric(unlist((mytb.impct))))
#        myimpct.sc  <- sum(myimpct.vec)
#        myimpct.sc2 <- round(sum(myimpct.vec)/(length(myimpct.vec)/nb.var), digits=1)
#
#        
#        mysum.tb[myidx,1:2]  <- mytb.tmp[1,c(1,5)] 
#        mysum.tb[myidx,3]    <- old.id
#        mysum.tb[myidx,4:5]  <- mypos 
#        mysum.tb[myidx,6:7]  <- mytb.tmp[1,6:7]
#        mysum.tb[myidx,8:(7+k)]  <- as.matrix(mytb.tmp[1,8:(7+k)])
#        mysum.tb[myidx,8+k]  <- capture
#        mysum.tb[myidx,9+k]  <- nb.var
#        mysum.tb[myidx,10+k] <- myimpct.sc
#        mysum.tb[myidx,11+k] <- myimpct.sc2
#        mysum.tb[myidx,12+k] <- mean  (unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
#        mysum.tb[myidx,13+k] <- median(unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
#        mysum.tb[myidx,14+k] <- max   (unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
#        if (! all(is.na(mytb.tmp[,mypv]))) {
#            mysum.tb[myidx,15+k] <- mean  (mytb.tmp[,mypv], na.rm=TRUE)
#            mysum.tb[myidx,16+k] <- median(mytb.tmp[,mypv], na.rm=TRUE)
#            mysum.tb[myidx,17+k] <- min   (mytb.tmp[,mypv], na.rm=TRUE)
#            mysum.tb[myidx,18+k] <- mean  (-log10(mytb.tmp[,mypv]), na.rm=TRUE)
#            mysum.tb[myidx,19+k] <- median(-log10(mytb.tmp[,mypv]), na.rm=TRUE)
#            mysum.tb[myidx,20+k] <- max   (-log10(mytb.tmp[,mypv]), na.rm=TRUE)
#        }
#
#    } else {
#        myexpr.tmp <- myexpr[ grep(g, myexpr[,1]) , myexpr.cln ]
#        if (! all(is.na(myexpr.tmp))) { myexpr.tmp <- apply(myexpr[ grep(paste0(g,"."), myexpr[,1], fixed=TRUE) , myexpr.cln ], 2, function(x) max(x, na.rm=TRUE)) } else { myexpr.tmp <- rep(NA, k) }
#
#        mysum.tb[myidx,1]        <- as.character(mygff.genes.qtl[ grep(g, mygff.genes.qtl[,9]), 1 ])
#        mysum.tb[myidx,c(2,4,5)] <- as.data.frame(mygff.genes.qtl[ grep(g, mygff.genes.qtl[,9]), c(9,4:5) ])
#        mysum.tb[myidx,3]        <- old.id
#        mysum.tb[myidx,6:7]      <- myann[ grep(g,myann[,1]), 2:3 ][1,]
#        mysum.tb[myidx,8:(7+k)]  <- myexpr.tmp
#        mysum.tb[myidx,8+k]      <- capture
#
#    }
#
#}
#
## Remove empty rows
#mysum.tb <- as.data.frame(mysum.tb[ ! rowSums(is.na(mysum.tb)) == ncol(mysum.tb), ])
#
## Name column
#colnames(mysum.tb) <- cln.nm 
#
#myrows <- na.omit(rownames(mysum.tb)[mysum.tb[,1] == "Chr_3" & mysum.tb[,5] >= bf.lim[1] & mysum.tb[,4] <= bf.lim[2]])
#mysum.tb.qtl <- mysum.tb[ myrows, ]
#
#write.table(mysum.tb.qtl, paste0(res.folder,"/QTL_table_summary.tsv"), row.names=FALSE, quote=FALSE, sep="\t")



mypv <- pv.vec[ colnames(myfreq.data)[pv.vec] %>% nchar() %>% which.max() ]

if (length(mypv) == 1) { myrows <- myfreq.data.fltr[,mypv] <= bf.cor } else { myrows <- rowSums(myfreq.data.fltr[,mypv] <= bf.cor) > 0 }

mychr <- myfreq.data.fltr[ myrows , 1 ] %>% sort() %>% unique()

for (i in mychr) {
    a <- qtl.tb(x=myfreq.data, x.fltr=myfreq.data.fltr, chr=i, pv.cln=mypv, bf.cor, ann.tb=myann, expr.tb=myexpr, expr.cln=myexpr.cln, expr.nm=myexpr.nm, cln.print=c(gt.vec,frq.vec,pv.print))
    qtl.tb.sum(a, ann.tb=myann, expr.tb=myexpr, expr.cln=myexpr.cln, expr.nm=myexpr.nm, gff=mygff.genes, baits.bed=mybaits)
}


#=========#
# Figures #
#=========#

# Replicate in allele frequency plot
myaf.plot(mydata, mylib.rep, myfn)

# Allele frequency plot from filtered data
vec <- grep("Sm", mylib.rep, value=TRUE)
myaf.plot(mydata[sd.rows,], vec, myfn, ext=".freq-fltr-sd.png")

# Z-score distribution plot
z.vec <- grep("z", colnames(myfreq.data.fltr), value=T)
myz.qq(myfreq.data.fltr, z.vec, myfn)

#--------
# Not sure that is a good idea (see above)

## p-value distribution plot
#p.vec <- grep("pv", colnames(myfreq.data.fltr), value=T)
#myp.qq(myfreq.data.fltr, p.vec, paste(myfn,"pv",sep="."))

#--------

myfreq.data.fltr.r <- rename_chr_SmV7(myfreq.data.fltr, 1)


# P-value plot
runmed.vec <- c(NA, 5,11,21)
mydef <- "HD"

pv.vec <- grep("pv", colnames(myfreq.data.fltr.r))

f.vec <- grep("freq", colnames(myfreq.data.fltr.r))

sf.vec <- grep("sub.freq", colnames(myfreq.data.fltr.r))

pf.vec <- grep("pol.freq", colnames(myfreq.data.fltr.r))

myfreq.data.fltr.r <- myfreq.data.fltr.r[ is.finite(rowSums(myfreq.data.fltr.r[,pv.vec])), ]

for (j in runmed.vec) {
    mygraph(myfreq.data.fltr.r, pv.vec, j, "pvalue", myfn, def=mydef)
    
#   mygraph(myfreq.data.fltr.r, f.vec, j, "freq", myfn)

    mygraph(myfreq.data.fltr.r, sf.vec, j, "freq", myfn, -0.75, 0.75, ylab="Difference in allele frequency", def=mydef)
    
    mygraph(myfreq.data.fltr.r, pf.vec, j, "freq", myfn, def=mydef)
}


j <- 1
pv.vec2 <- 47
sf.vec2 <- 29
pf.vec2 <- 25
mygraph(myfreq.data.fltr.r, pv.vec, j, "pvalue", myfn, def=mydef)
mygraph(myfreq.data.fltr.r, sf.vec, j, "freq", myfn, -0.75, 0.75, ylab="Difference in allele frequency", def=mydef)
    
mygraph(myfreq.data.fltr.r, pv.vec, j, "pvalue", myfn, def=mydef)
mygraph(myfreq.data.fltr.r, pf.vec, j, "freq", myfn, def=mydef)


mybc <- -log10(0.05/nrow(myfreq.data.fltr.r))
png(paste0("graphs_HD/", myfn, ".pv-sf.png"), width=72*12, height=72*10)
layout(matrix(1:2, ncol = 1))
par(mar=c(5,4,1,1)+0.1) # For small size
matplot.data(myfreq.data.fltr.r, pv.vec2, "pvalue", abline.h=NULL, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), by.pos=TRUE, type="p")
matplot.data(myfreq.data.fltr.r, sf.vec2, "freq", ylim.min=-1, ylim.max=1, abline.h=mybc, abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), ylab="Difference in allele frequency", by.pos=TRUE, type="p")
dev.off()


mybc <- -log10(0.05/nrow(myfreq.data.fltr.r))
png(paste0("graphs_HD/", myfn, ".pv-pf.png"), width=72*12, height=72*10)
layout(matrix(1:2, ncol = 1))
par(mar=c(5,4,1,1)+0.1) # For small size
matplot.data(myfreq.data.fltr.r, pv.vec2, "pvalue", abline.h=mybc, abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), by.pos=TRUE, type="p")
matplot.data(myfreq.data.fltr.r, pf.vec2, "freq", ylim.min=0, ylim.max=1, abline.h=mybc, abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), ylab="Polarized allele frequency", by.pos=TRUE, type="p")
dev.off()
