#!/usr/bin/env Rscript
# Title: PZQ_analysis_func.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-05-15
# Modified in: 



#==========#
# Comments #
#==========#

# v0.0 - 2019-05-15: store functions in a separate file from the inital script



#======================#
# Packages and options #
#======================#

library("rtracklayer")
library("magrittr")
library("vcfR")

options(stringsAsFactors = FALSE)



#===========#
# Functions #
#===========#

#------------#
# Save plots #
#------------#

save.myplot <- function (filename, width=8, height=8) {
	# Folder creation for graphic output (if does not exist)
	if (file.exists("Graphs/") == FALSE) {dir.create("Graphs")}

	# Graphic saving (png and eps formats)
	dev.copy(png,paste("Graphs/",filename,".png",sep=""), width=width*72, height=height*72, bg="white")
	dev.off()
	dev.copy2pdf(file=paste("Graphs/",filename,".pdf",sep=""), width=width, height=height)
}


#-------------#
# Save tables #
#-------------#

save.mytable <- function (x, filename, ext=".tsv", sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE){
	# Folder creation for table output (if does not exist)
	if (file.exists("Results/") == FALSE) {dir.create("Results")}
	
	# Data export
	write.table(x, file=paste("Results/",filename, ext, sep=""), col.names=col.names, row.names=row.names, append=append, sep=sep)
}


#------------------------------------------#
# Function to combine p-value from Z score #
#------------------------------------------#

# source: http://en.wikipedia.org/wiki/Fisher%27s_method
combine.z.pv <- function (x, y) {
    mysum <- rowSums(cbind(x,y))
    pnorm(-abs(mysum/sqrt(2)))
}


#----------#
# Filename #
#----------#

# Improvement of file_path_sans_ext from the tools package
filename <- function (x, compression=FALSE, length=TRUE) {
    # Test if compression extension
    if (compression) {x <- sub("[.](gz|bz2|xz)$", "", x)}
    # Split name
    myname <- strsplit(x,"\\.")[[1]]

    # Set length of the filename to keep if several point
    if (isTRUE(length)) {
        length <- length(myname)-1
    } else if (! is.numeric(length)) {
        stop("Length must be numeric.")
    } else if (length > length(myname)) {
        stop("Specified lenght is greater than the actual length.")
    }

    basename(paste(myname[0:length], collapse='.'))
}



#===========================#
# Data processing functions #
#===========================#

#--------------#
# Peaks report #
#--------------#

# Peaks report
pk.rpt <- function ( x, sign, trsh, cln.print=NULL, cln.nm.ext=NULL, res.folder, filename=paste("peaks_report.",x.nm,".",cln.nm.ext,sep="") ) {
    # Usage
    ## x: dataframe to analyse
    ## sign: >, <, >=, <= or ==
    ## trsh: trehsold value
    ## cln.print: vector containing the column number or name to print in the report (default all)
    ## cln.nm.ext: name to be printed in the filename of the report (default, name from cln.print)
    ## res.folder: folder to store tables
    ## filename: filename of the table
    
    # Check steps
    ## Check mandatory objects
    if (missing(x) || missing(sign) || missing(trsh) || missing(res.folder)) {stop("x, sign, trsh and res.folder are mandatory.", call.=FALSE)}
    ## Check class of objects
#    if (! is.data.frame(x)) {stop("x must be a data.frame.", call.=FALSE)}
    if (! is.numeric(trsh) || length(trsh) > 1) {stop("trsh must be a single numerical value.", call.=FALSE)}
    if (! is.null(cln.print) && ! is.vector(cln.print)) {stop("cln.print must be a vector.", call.=FALSE)}
#    if (! is.null(cln.nm.ext) && ! (is.vector(cln.nm.ext) || is.character(cln.nm.ext))) {stop("cln.nm.ext must be a vector.", call.=FALSE)}
    if (! is.null(cln.nm.ext) && ! is.character(cln.nm.ext) && ! length(cln.nm.ext) > 1) {stop("cln.nm.ext must be a single character string.", call.=FALSE)}
    if (! is.character(res.folder) || length(res.folder) > 1) {stop("res.foler must be a character string of length 1.")}

    ## Check signs
    if (! is.character(sign) || ! (sign == ">" || sign == "<" || sign == "<=" || sign == ">=" || sign == "==")) {stop("sign must be >, <, >=, <= or ==.", call.=FALSE)}
    
    ## Check for directory
    if (file.exists(res.folder) == FALSE) {dir.create(res.folder, recursive=TRUE)}

    # x name 
    x.splt <- unlist(strsplit(deparse(substitute(x)), "\\[|,|\\]$"))[1]
    x.nm <- x.splt[1]

    # Original data.frame
    x.data <- get(x.nm)

    # Column range and names
    if (is.null(cln.print)) { cln.print <- colnames(get(x.nm)) }

    # Table generation
    if (! is.data.frame(x)) { x <- as.data.frame(x) }
    myrows <- apply(x, 1, function(z) any(do.call (sign, list(z, trsh))))
    report <- x.data[ myrows , cln.print ]

    # Table export
    write.table(report, paste0(res.folder,"/",filename), row.names=FALSE, quote=FALSE, sep="\t")
}


##-------------------------------------------------------------#
## Detail and summary tables listing genes and variants in QTL #
##-------------------------------------------------------------#
#
#qtl.tb <- function(x, x.fltr, chr, pv.cln, bf.cor, ann.tb, expr.tb, expr.cln, expr.nm=NULL, cln.print, gff.genes, baits.bed=NULL) {
#
#    # Usage
#    ## x            a table with allele frequencies, genotypes, and p-values
#    ## x.flt        filtered x table
#    ## chr          chromosome to filter on
#    ## pv.cln       column(s) that contains p-values to filter on       
#    ## bf.cor       Bonferroni correction
#    ## ann.tb       annotation table
#    ## expr.tb      gene expression table
#    ## expr.cln     column of the gene expression table to include
#    ## expr.nm      name of expression stage
#    ## cln.print    column to print in the final table
#    ## gff.genes    GFF file with only gene information
#    ## baits.bed    BED coordinated of the baits
#
#    cat("QTL tables: processing ", chr, "...\n", sep="")
#
#    if (length(pv.cln) == 1) {
#        my.qtl <- na.omit(x.fltr[ x.fltr[,pv.cln] <= bf.cor & x.fltr[,1] == chr, 2])
#    } else {
#        my.qtl <- x.fltr[ rowSums(x.fltr[,pv.cln] <= bf.cor, na.rm=TRUE) > 0 & x.fltr[,1] == chr, 2]
#    }
#
#    
#    bf.lim <- c(my.qtl[1], my.qtl[length(my.qtl)])
#    
#    # Select rows in the QTL region
#    x <- x[ x[,1] == chr & x[,2] >= bf.lim[1] & x[,2] <= bf.lim[2], ]
#    x.fltr <- x.fltr[ x.fltr[,1] == chr & x.fltr[,2] >= bf.lim[1] & x.fltr[,2] <= bf.lim[2], ]
#
#    if(is.null(expr.nm)) { expr.nm <- colnames(expr)[expr.cln] }
#
#    k <- length(expr.cln)
#
#    #--------------#
#    # Detail table #
#    #--------------#
#
#    myinfo <- strsplit(x[,8], ";", fixed=T)
#    
#    # Prepare parallelization
#    myinfo.lg <- 1:length(myinfo)
#    n <- detectCores()-1
#    myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))
#
##    myinfo.ann <- lapply(myinfo, function(x) unlist(grep("ANN",x,value=T) %>% strsplit(., "|", fixed=T)))
#    myinfo.ann <- foreach(i=1:length(myinfo.lg), .combine='c') %dopar% {
#        lapply(myinfo[myinfo.lg[[i]]], function(x) unlist(grep("ANN",x, value=T) %>% strsplit(., "|", fixed=T)))
#    }
#
#    mygenes.ls <- vector("list", length(myinfo.ann))
#    cln.nm <- c("Gene", "GFF_annotation", "HHPred_annotation", expr.nm)
#    mymain.genes.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(myinfo.ann)))
#    colnames(mymain.genes.tb) <- cln.nm
#
#
##    myinfo.lg <- 1:length(myinfo)
##    n <- round(20*log(length(myinfo))) 
##    myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))
#
#    old.workers <- getDoParWorkers()
#    registerDoParallel(4)
#    
#    mygenes.ls <- foreach(j=1:length(myinfo.lg), .combine='c') %dopar% {
#        myidx <- myinfo.lg[[j]]
#
#        mygenes.ls.tmp <- mygenes.ls[myidx]
#        
#        for (i in 1:length(myinfo.ann)) {
##            ann <- myinfo.ann[[i]]
#            ann <- myinfo.ann[[myidx[i]]]
#            todo <- TRUE
#
#            mycln.tmp <- seq(0, by=15, length.out=round(length(ann)/15))
#            
#            if (length(mycln.tmp) > 0) {
#                mymtx <- sapply(mycln.tmp, function(y) ann[y+c(2,3,4,7,8,10,11)])
#
#                mymain.gene <- names(sort(table(mymtx[3,]), decreasing=T)[1]) 
#           
#                myRNA <- mymtx[4,]  
#                myRNA.tb <- matrix(NA, nrow=2+length(expr.cln), ncol=ncol(mymtx))
#
#                for (r in 1:length(myRNA)) {
#                    
#                    # Corresponding gene
#                    mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
#                    if (is.na(mygene)) next 
#                    
#                    # Transcript data
#                    ann.tmp <- ann.tb[ ann.tb[,1] == myRNA[r] , ]
#                    # If no transcript detected
#                    if (length(ann.tmp) == 0) { ann.tmp <- ann[ ann[,1] == mygene , ] }
#                    # If really nothing detected
#                    if (dim(ann.tmp)[1] == 0) { ann.tmp <- matrix(NA, ncol=3, nrow=1) }
#                    # Add data to the matrix
#                    myRNA.tb[1:2, r] <- t(as.matrix(ann.tmp[,2:3]))
#
#                    # Expression data
#                    expr.tmp <- expr.tb[ expr.tb[,1] == myRNA[r] , expr.cln ]
#                    # If no expression found
#                    if (k == 1 && length(expr.tmp) == 0 ) { expr.tmp <- rep(NA, length(expr.cln)) }
#                    if (k > 1  && (length(expr.tmp) == 0 | dim(expr.tmp)[1] == 0)) { expr.tmp <- rep(NA, length(expr.cln)) }
#                    # Add data to matrix
#                    myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(expr.tmp)
#
#                    # Add info for main gene table
#                    if (mygene == mymain.gene & isTRUE(todo)) {
#                        mymain.genes.tb[i,] <- cbind(mygene, t(myRNA.tb[,r]))
#
#                        # If expression is NA, is there any other expression data for this given gene?
#                        if (all(is.na(expr.tmp)) | sum(expr.tmp) == 0) {
#                            expr.tmp <- expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ]
#                            if (k == 1 && ! all(is.na(expr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- max(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE) }
#                            if (k > 1  && ! all(is.na(expr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- apply(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
#                        }
#
#                        # Don't do anything else for this gene in the future
#                        todo <- FALSE
#                    }
#                }
#
#                # Combine matrix
#                mymtx <- rbind(mymtx, myRNA.tb)
#                
#                # Reorder rows
#                mymtx <- mymtx[c(3, 8:(9+k), 4, 1:2, 5:7), ]
#                
##                mygenes.ls[[i]] <- mymtx
#                mygenes.ls.tmp[[i]] <- mymtx
#            }
#        }
#
#        return(mygenes.ls.tmp)
#
##        rm(mygenes.ls.tmp)
##        gc(reset=TRUE)
#
#    }
#
#    # Reset initial workers
#    registerDoParallel(old.workers)
#
#    mygenes.ls[ sapply(mygenes.ls, is.null) ] <- NA
#    mymax      <- max(sapply(mygenes.ls, length))
#    mygenes.tb <- lapply(mygenes.ls, function(x) x[1:mymax])
#    mygenes.tb <- t(do.call(cbind, mygenes.tb))
#    mygenes.tb[ mygenes.tb == "" ] <- NA
#    mygenes.tb <- as.data.frame(mygenes.tb)
#
#    cln.nm.b <- c("Gene (hit ", "GFF annotation (hit ", "HHPred annotation (hit ", paste(expr.nm, "(hit "), "Transcript (hit ", "Region (hit ", "Impact (hit ", "Type (hit ", "DNA mutation (hit ", "Protein mutation (hit ")
#    cln.nm <- NULL
#    for (i in 1:(mymax/length(cln.nm.b))) { cln.nm <- c(cln.nm, paste0(cln.nm.b, i, ")")) }
#    colnames(mygenes.tb) <- cln.nm
#
#    # Final table
#    myfinal.tb <- cbind(x[,c(1,2,4,5)], mymain.genes.tb, x[,cln.print], mygenes.tb)
#    myfinal.tb.qtl <- myfinal.tb 
#
##    ## Select rows in the QTL region
##    myrows <- na.omit(rownames(myfinal.tb)[myfinal.tb[,1] == chr & myfinal.tb[,2] >= bf.lim[1] & myfinal.tb[,2] <= bf.lim[2]])
##    myfinal.tb.qtl <- myfinal.tb[ myrows, ]
#
#    ## Select only informative columns then
#    myclns <- ! apply(myfinal.tb.qtl, 2, function(x) all(is.na(x)))
#    myfinal.tb.qtl <- myfinal.tb.qtl[ , myclns ]
#
#    if(! exists("res.folder")) {res.folder <- "results/"}
#
#    write.table(myfinal.tb.qtl, paste0(res.folder,"/",chr,"_QTL_table_details.tsv"), row.names=FALSE, quote=FALSE, sep="\t")
#
#
#    #---------------#
#    # Summary table #
#    #---------------#
#
#    gff.genes.qtl <- gff.genes[gff.genes[,1] == chr & gff.genes[,5] >= bf.lim[1] & gff.genes[,4] <= bf.lim[2],]
#
#    mygenes.occ <- unique(gff.genes.qtl[,9])
#    mygenes.occ <- na.omit(unique(c(mygenes.occ, myfinal.tb.qtl[,5])))
#
#    myimpct.cln <- grep("Impact", colnames(myfinal.tb))
#    myimpct.tb  <- as.data.frame(matrix(c("LOW", "MODERATE", "HIGH", "MODIFIER", 0.3, 0.6, 0.9, 0.3), ncol=2))
#
#    myfreq.cln  <- grep("pol.freq",colnames(myfinal.tb))
#
#    cln.nm   <- c("Chr.", "Gene ID", "v5 gene ID", "Start", "End", "GFF annotation", "HHPred annotation", expr.nm, "Captured", "Nb of variable sites", "Impact score", "Weighted impact score", "Mean polarized all. freq.", "Median polarized all. freq.", "Highest polarized all. freq.", "Mean p-value", "Median p-value", "Lowest p-value", "Mean -log10(p-value)", "Median -log10(p-value)", "Highest -log10(p-value)")
#    mysum.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(mygenes.occ)))
#    for (g in mygenes.occ) {
#        myidx <- match(g, mygenes.occ)
#        mypv  <- rev(grep("pv.", colnames(myfinal.tb), fixed=TRUE))[1]
#
#        if (! is.null(mybaits)) {
#            # Is the gene in the capture array?
#            g.crdt <- as.data.frame(mygff.genes[ mygff.genes[,9] == g, c(1,4:5) ])
#            if (dim(g.crdt)[1] > 0) {
#                capture <- any( findInterval(mybaits[ mybaits[,1] == g.crdt[,1], 3], g.crdt[,2:3]) == 1 )
#            } else {
#                capture <- FALSE
#            }
#            if (capture) { capture <- "yes" } else { capture <- "no" }
#        } else {
#            capture <- NA
#        }
#        
#        # Is there any GFF v5 ID(s) associated to the current gene?
#        old.id <- as.data.frame(gff.genes[ gff.genes[,9] == g, 15 ])[,3]
#        if (length(old.id) == 0) { 
#            old.id <- NA 
#        } else {
#            old.id <- paste(old.id, collapse=", ")
#        }
#        
#        # Pull out information if the gene is in the variant data table
#        if (any(myfinal.tb[,5] == g, na.rm=TRUE)) {
#            mytb.tmp <- myfinal.tb[ ! is.na(myfinal.tb[,5]) & myfinal.tb[,5] == g, ]
#        
#            mypos <- as.data.frame(gff.genes[ gff.genes[,9] == g, c(4,5) ])
#            if (dim(mypos)[1] == 0) { mypos <- c(NA,NA) }
#
#            nb.var <- nrow(mytb.tmp)
#
#            mytb.impct <- mytb.tmp[,myimpct.cln]
#            for (i in 1:nrow(myimpct.tb)) { mytb.impct[ mytb.impct == myimpct.tb[i,1] ] <- myimpct.tb[i,2] }
#            myimpct.vec <- na.omit(as.numeric(unlist((mytb.impct))))
#            myimpct.sc  <- sum(myimpct.vec)
#            myimpct.sc2 <- round(sum(myimpct.vec)/(length(myimpct.vec)/nb.var), digits=1)
#
#            
#            mysum.tb[myidx,1:2]  <- mytb.tmp[1,c(1,5)] 
#            mysum.tb[myidx,3]    <- old.id
#            mysum.tb[myidx,4:5]  <- mypos 
#            mysum.tb[myidx,6:7]  <- mytb.tmp[1,6:7]
#            mysum.tb[myidx,8:(7+k)]  <- as.matrix(mytb.tmp[1,8:(7+k)])
#            mysum.tb[myidx,8+k]  <- capture
#            mysum.tb[myidx,9+k]  <- nb.var
#            mysum.tb[myidx,10+k] <- myimpct.sc
#            mysum.tb[myidx,11+k] <- myimpct.sc2
#            if (! all(is.na(mytb.tmp[,myfreq.cln]))) {
#                mysum.tb[myidx,12+k] <- mean  (unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
#                mysum.tb[myidx,13+k] <- median(unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
#                mysum.tb[myidx,14+k] <- max   (unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
#            }
#            if (! all(is.na(mytb.tmp[,mypv]))) {
#                mysum.tb[myidx,15+k] <- mean  (mytb.tmp[,mypv], na.rm=TRUE)
#                mysum.tb[myidx,16+k] <- median(mytb.tmp[,mypv], na.rm=TRUE)
#                mysum.tb[myidx,17+k] <- min   (mytb.tmp[,mypv], na.rm=TRUE)
#                mysum.tb[myidx,18+k] <- mean  (-log10(mytb.tmp[,mypv]), na.rm=TRUE)
#                mysum.tb[myidx,19+k] <- median(-log10(mytb.tmp[,mypv]), na.rm=TRUE)
#                mysum.tb[myidx,20+k] <- max   (-log10(mytb.tmp[,mypv]), na.rm=TRUE)
#            }
#
#        } else {
#            myexpr.tmp <- expr.tb[ grep(g, expr.tb[,1]) , expr.cln ]
#            if (k == 1 && ! all(is.na(myexpr.tmp))) {
#                myexpr.tmp <- max(expr.tb[ grep(paste0(g,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE)
#            } else if (k > 1  && ! all(is.na(myexpr.tmp))) {
#                myexpr.tmp <- apply(expr.tb[ grep(paste0(g,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE))
#            } else {
#                myexpr.tmp <- rep(NA, k)
#            }
#
#            mysum.tb[myidx,1]        <- as.character(gff.genes.qtl[ grep(g, gff.genes.qtl[,9]), 1 ])
#            mysum.tb[myidx,c(2,4,5)] <- as.data.frame(gff.genes.qtl[ grep(g, gff.genes.qtl[,9]), c(9,4:5) ])
#            mysum.tb[myidx,3]        <- old.id
#            mysum.tb[myidx,6:7]      <- ann.tb[ grep(g,ann.tb[,1]), 2:3 ][1,]
#            mysum.tb[myidx,8:(7+k)]  <- myexpr.tmp
#            mysum.tb[myidx,8+k]      <- capture
#
#        }
#
#    }
#
#    # Remove empty rows
#    mysum.tb <- as.data.frame(mysum.tb[ ! rowSums(is.na(mysum.tb)) == ncol(mysum.tb), ])
#
#    # Name column
#    colnames(mysum.tb) <- cln.nm 
#
#    myrows <- na.omit(rownames(mysum.tb)[mysum.tb[,1] == chr & mysum.tb[,5] >= bf.lim[1] & mysum.tb[,4] <= bf.lim[2]])
#    mysum.tb.qtl <- mysum.tb[ myrows, ]
#
#    write.table(mysum.tb.qtl, paste0(res.folder,"/",chr,"_QTL_table_summary.tsv"), row.names=FALSE, quote=FALSE, sep="\t")
#}


#-------------------------------------------------------------#
# Detail and summary tables listing genes and variants in QTL #
#-------------------------------------------------------------#

qtl.tb <- function(x, x.fltr, chr, pv.cln, bf.cor, ann.tb, expr.tb, expr.cln, expr.nm=NULL, cln.print) {

    # Usage
    ## x            a table with allele frequencies, genotypes, and p-values
    ## x.flt        filtered x table
    ## chr          chromosome to filter on
    ## pv.cln       column(s) that contains p-values to filter on       
    ## bf.cor       Bonferroni correction
    ## ann.tb       annotation table
    ## expr.tb      gene expression table
    ## expr.cln     column of the gene expression table to include
    ## expr.nm      name of expression stage
    ## cln.print    column to print in the final table

    cat("QTL tables: processing ", chr, "...\n", sep="")

    if (length(pv.cln) == 1) {
        my.qtl <- na.omit(x.fltr[ x.fltr[,pv.cln] <= bf.cor & x.fltr[,1] == chr, 2])
    } else {
        my.qtl <- x.fltr[ rowSums(x.fltr[,pv.cln] <= bf.cor, na.rm=TRUE) > 0 & x.fltr[,1] == chr, 2]
    }

    
    bf.lim <- c(my.qtl[1], my.qtl[length(my.qtl)])
    
    # Select rows in the QTL region
    x <- x[ x[,1] == chr & x[,2] >= bf.lim[1] & x[,2] <= bf.lim[2], ]
    x.fltr <- x.fltr[ x.fltr[,1] == chr & x.fltr[,2] >= bf.lim[1] & x.fltr[,2] <= bf.lim[2], ]

    if(is.null(expr.nm)) { expr.nm <- colnames(expr)[expr.cln] }

    k <- length(expr.cln)

    #--------------#
    # Detail table #
    #--------------#

    myinfo <- strsplit(x[,8], ";", fixed=T)
    
    # Prepare parallelization
    myinfo.lg <- 1:length(myinfo)
    n <- detectCores()-1
    myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))

#    myinfo.ann <- lapply(myinfo, function(x) unlist(grep("ANN",x,value=T) %>% strsplit(., "|", fixed=T)))
    myinfo.ann <- foreach(i=1:length(myinfo.lg), .combine='c') %dopar% {
        lapply(myinfo[myinfo.lg[[i]]], function(x) unlist(grep("ANN",x, value=T) %>% strsplit(., "|", fixed=T)))
    }

    mygenes.ls <- vector("list", length(myinfo.ann))


#    myinfo.lg <- 1:length(myinfo)
#    n <- round(20*log(length(myinfo))) 
#    myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))

    old.workers <- getDoParWorkers()
    registerDoParallel(4)
    
    mygenes.ls <- foreach(j=1:length(myinfo.lg), .combine='c') %dopar% {
        myidx <- myinfo.lg[[j]]

        mygenes.ls.tmp <- mygenes.ls[myidx]
        
        for (i in 1:length(myidx)) {
#            ann <- myinfo.ann[[i]]
            ann <- myinfo.ann[[myidx[i]]]
#            todo <- TRUE

            mycln.tmp <- seq(0, by=15, length.out=round(length(ann)/15))
            
            if (length(mycln.tmp) > 0) {
                mymtx <- sapply(mycln.tmp, function(y) ann[y+c(2,3,4,7,8,10,11)])

                mymain.gene <- names(sort(table(mymtx[3,]), decreasing=T)[1]) 
           
                myRNA <- mymtx[4,]  
                myRNA.tb <- matrix(NA, nrow=2+length(expr.cln), ncol=ncol(mymtx))

                for (r in 1:length(myRNA)) {
                    
                    # Corresponding gene
                    mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
                    if (is.na(mygene)) next 
                    
                    # Transcript data
                    ann.tmp <- ann.tb[ ann.tb[,1] == myRNA[r] , ]
                    # If no transcript detected
                    if (length(ann.tmp) == 0) { ann.tmp <- ann[ ann[,1] == mygene , ] }
                    # If really nothing detected
                    if (dim(ann.tmp)[1] == 0) { ann.tmp <- matrix(NA, ncol=3, nrow=1) }
                    # Add data to the matrix
                    myRNA.tb[1:2, r] <- t(as.matrix(ann.tmp[,2:3]))

                    # Expression data
                    expr.tmp <- expr.tb[ expr.tb[,1] == myRNA[r] , expr.cln ]
                    # If no expression found
                    if (k == 1 && length(expr.tmp) == 0 ) { expr.tmp <- rep(NA, length(expr.cln)) }
                    if (k > 1  && (length(expr.tmp) == 0 | dim(expr.tmp)[1] == 0)) { expr.tmp <- rep(NA, length(expr.cln)) }
                    # Add data to matrix
                    myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(expr.tmp)

#                    # Add info for main gene table
#                    if (mygene == mymain.gene & isTRUE(todo)) {
#                        mymain.genes.tb[i,] <- cbind(mygene, t(myRNA.tb[,r]))
#
#                        # If expression is NA, is there any other expression data for this given gene?
#                        if (all(is.na(expr.tmp)) | sum(expr.tmp) == 0) {
#                            expr.tmp <- expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ]
#                            if (k == 1 && ! all(is.na(expr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- max(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE) }
#                            if (k > 1  && ! all(is.na(expr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- apply(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
#                        }
#
#                        # Don't do anything else for this gene in the future
#                        todo <- FALSE
#                    }
                }

                # Combine matrix
                mymtx <- rbind(mymtx, myRNA.tb)
                
                # Reorder rows
                mymtx <- mymtx[c(3, 8:(9+k), 4, 1:2, 5:7), ]
                
#                mygenes.ls[[i]] <- mymtx
                mygenes.ls.tmp[[i]] <- mymtx
            }
        }

        return(mygenes.ls.tmp)

#        rm(mygenes.ls.tmp)
#        gc(reset=TRUE)

    }
    
    cln.nm <- c("Gene", "GFF_annotation", "HHPred_annotation", expr.nm)
    mymain.genes.tb <- foreach(j=1:length(myinfo.lg), .combine='rbind') %dopar% {
        myidx <- myinfo.lg[[j]]

        mymain.genes.tb.tmp <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(myinfo.ann)))
        
        for (i in 1:length(myidx)) {
#            ann <- myinfo.ann[[i]]
            ann <- myinfo.ann[[myidx[i]]]
#            todo <- TRUE

            mycln.tmp <- seq(0, by=15, length.out=round(length(ann)/15))
            
            if (length(mycln.tmp) > 0) {
                mymtx <- sapply(mycln.tmp, function(y) ann[y+c(2,3,4,7,8,10,11)])

                mymain.gene <- names(sort(table(mymtx[3,]), decreasing=T)[1]) 
           
                myRNA <- mymtx[4,]  
                myRNA.tb <- matrix(NA, nrow=2+length(expr.cln), ncol=ncol(mymtx))

                for (r in 1:length(myRNA)) {
                    
                    # Corresponding gene
                    mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
                    if (is.na(mygene)) next 
                    
                    # Transcript data
                    ann.tmp <- ann.tb[ ann.tb[,1] == myRNA[r] , ]
                    # If no transcript detected
                    if (length(ann.tmp) == 0) { ann.tmp <- ann[ ann[,1] == mygene , ] }
                    # If really nothing detected
                    if (dim(ann.tmp)[1] == 0) { ann.tmp <- matrix(NA, ncol=3, nrow=1) }
                    # Add data to the matrix
                    myRNA.tb[1:2, r] <- t(as.matrix(ann.tmp[,2:3]))

                    # Expression data
                    expr.tmp <- expr.tb[ expr.tb[,1] == myRNA[r] , expr.cln ]
                    # If no expression found
                    if (k == 1 && length(expr.tmp) == 0 ) { expr.tmp <- rep(NA, length(expr.cln)) }
                    if (k > 1  && (length(expr.tmp) == 0 | dim(expr.tmp)[1] == 0)) { expr.tmp <- rep(NA, length(expr.cln)) }
                    # Add data to matrix
                    myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(expr.tmp)

                    # Add info for main gene table
                    if (mygene == mymain.gene) {
                        mymain.genes.tb.tmp[i,] <- cbind(mygene, t(myRNA.tb[,r]))

                        # If expression is NA, is there any other expression data for this given gene?
                        if (all(is.na(expr.tmp)) | sum(expr.tmp) == 0) {
                            expr.tmp <- expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ]
                            if (k == 1 && ! all(is.na(expr.tmp))) { mymain.genes.tb.tmp[i,4:ncol(mymain.genes.tb.tmp)] <- max(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE) }
                            if (k > 1  && ! all(is.na(expr.tmp))) { mymain.genes.tb.tmp[i,4:ncol(mymain.genes.tb.tmp)] <- apply(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
                        }

                        # Don't do anything else for this gene in the future
#                       todo <- FALSE
                        break
                    }
                }

            }
        }

        return(mymain.genes.tb.tmp)

#        rm(mygenes.ls.tmp)
#        gc(reset=TRUE)

    }
    
#    mymain.genes.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(myinfo.ann)))
    colnames(mymain.genes.tb) <- cln.nm

    # Reset initial workers
    registerDoParallel(old.workers)

    mygenes.ls[ sapply(mygenes.ls, is.null) ] <- NA
    mymax      <- max(sapply(mygenes.ls, length))
    mygenes.tb <- lapply(mygenes.ls, function(x) x[1:mymax])
    mygenes.tb <- t(do.call(cbind, mygenes.tb))
    mygenes.tb[ mygenes.tb == "" ] <- NA
    mygenes.tb <- as.data.frame(mygenes.tb)

    cln.nm.b <- c("Gene (hit ", "GFF annotation (hit ", "HHPred annotation (hit ", paste(expr.nm, "(hit "), "Transcript (hit ", "Region (hit ", "Impact (hit ", "Type (hit ", "DNA mutation (hit ", "Protein mutation (hit ")
    cln.nm <- NULL
    for (i in 1:(mymax/length(cln.nm.b))) { cln.nm <- c(cln.nm, paste0(cln.nm.b, i, ")")) }
    colnames(mygenes.tb) <- cln.nm

    # Final table
    myfinal.tb <- cbind(x[,c(1,2,4,5)], mymain.genes.tb, x[,cln.print], mygenes.tb)
    myfinal.tb.qtl <- myfinal.tb 

#    ## Select rows in the QTL region
#    myrows <- na.omit(rownames(myfinal.tb)[myfinal.tb[,1] == chr & myfinal.tb[,2] >= bf.lim[1] & myfinal.tb[,2] <= bf.lim[2]])
#    myfinal.tb.qtl <- myfinal.tb[ myrows, ]

    ## Select only informative columns then
    myclns <- ! apply(myfinal.tb.qtl, 2, function(x) all(is.na(x)))
    myfinal.tb.qtl <- myfinal.tb.qtl[ , myclns ]

    if(! exists("res.folder")) {res.folder <- "results/"}

    write.table(myfinal.tb.qtl, paste0(res.folder,"/",chr,"_QTL_table_details.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

    return(invisible(myfinal.tb.qtl))

}

#---------------#
# Summary table #
#---------------#

qtl.tb.sum <- function(x, ann.tb, expr.tb, expr.cln, expr.nm=NULL, gff.genes, baits.bed=NULL) {

    # Usage
    ## x            QTL table from qtl.tb
    ## gff.genes    GFF file with only gene information
    ## baits.bed    BED coordinated of the baits
    
    chr <- unique(x[,1])
    bf.lim <- c(x[1,2], x[nrow(x),2])

    cat("Summary QTL table: processing ", chr, "...\n", sep="")

    gff.genes.qtl <- gff.genes[gff.genes[,1] == chr & gff.genes[,5] >= bf.lim[1] & gff.genes[,4] <= bf.lim[2],]

    mygenes.occ <- unique(gff.genes.qtl[,9])
    mygenes.occ <- na.omit(unique(c(mygenes.occ, x[,5])))

    myimpct.cln <- grep("Impact", colnames(x))
    myimpct.tb  <- as.data.frame(matrix(c("LOW", "MODERATE", "HIGH", "MODIFIER", 0.3, 0.6, 0.9, 0.3), ncol=2))

    myfreq.cln  <- grep("pol.freq",colnames(x))

    cln.nm   <- c("Chr.", "Gene ID", "v5 gene ID", "Start", "End", "GFF annotation", "HHPred annotation", expr.nm, "Captured", "Nb of variable sites", "Impact score", "Weighted impact score", "Mean polarized all. freq.", "Median polarized all. freq.", "Highest polarized all. freq.", "Mean p-value", "Median p-value", "Lowest p-value", "Mean -log10(p-value)", "Median -log10(p-value)", "Highest -log10(p-value)")
    mysum.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(mygenes.occ)))
    for (g in mygenes.occ) {
        myidx <- match(g, mygenes.occ)
        mypv  <- rev(grep("pv.", colnames(x), fixed=TRUE))[1]

        if (! is.null(baits.bed)) {
            # Is the gene in the capture array?
            g.crdt <- as.data.frame(mygff.genes[ mygff.genes[,9] == g, c(1,4:5) ])
            if (dim(g.crdt)[1] > 0) {
                capture <- any( findInterval(baits.bed[ baits.bed[,1] == g.crdt[,1], 3], g.crdt[,2:3]) == 1 )
            } else {
                capture <- FALSE
            }
            if (capture) { capture <- "yes" } else { capture <- "no" }
        } else {
            capture <- NA
        }
        
        # Is there any GFF v5 ID(s) associated to the current gene?
        old.id <- as.data.frame(gff.genes[ gff.genes[,9] == g, 15 ])[,3]
        if (length(old.id) == 0) { 
            old.id <- NA 
        } else {
            old.id <- paste(old.id, collapse=", ")
        }
        
        # Pull out information if the gene is in the variant data table
        if (any(x[,5] == g, na.rm=TRUE)) {
            mytb.tmp <- x[ ! is.na(x[,5]) & x[,5] == g, ]
        
            mypos <- as.data.frame(gff.genes[ gff.genes[,9] == g, c(4,5) ])
            if (dim(mypos)[1] == 0) { mypos <- c(NA,NA) }

            nb.var <- nrow(mytb.tmp)

            mytb.impct <- mytb.tmp[,myimpct.cln]
            for (i in 1:nrow(myimpct.tb)) { mytb.impct[ mytb.impct == myimpct.tb[i,1] ] <- myimpct.tb[i,2] }
            myimpct.vec <- na.omit(as.numeric(unlist((mytb.impct))))
            myimpct.sc  <- sum(myimpct.vec)
            myimpct.sc2 <- round(sum(myimpct.vec)/(length(myimpct.vec)/nb.var), digits=1)

            
            mysum.tb[myidx,1:2]  <- mytb.tmp[1,c(1,5)] 
            mysum.tb[myidx,3]    <- old.id
            mysum.tb[myidx,4:5]  <- mypos 
            mysum.tb[myidx,6:7]  <- mytb.tmp[1,6:7]
            mysum.tb[myidx,8:(7+k)]  <- as.matrix(mytb.tmp[1,8:(7+k)])
            mysum.tb[myidx,8+k]  <- capture
            mysum.tb[myidx,9+k]  <- nb.var
            mysum.tb[myidx,10+k] <- myimpct.sc
            mysum.tb[myidx,11+k] <- myimpct.sc2
            if (! all(is.na(mytb.tmp[,myfreq.cln]))) {
                mysum.tb[myidx,12+k] <- mean  (unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
                mysum.tb[myidx,13+k] <- median(unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
                mysum.tb[myidx,14+k] <- max   (unlist(mytb.tmp[,myfreq.cln]), na.rm=TRUE)
            }
            if (! all(is.na(mytb.tmp[,mypv]))) {
                mysum.tb[myidx,15+k] <- mean  (mytb.tmp[,mypv], na.rm=TRUE)
                mysum.tb[myidx,16+k] <- median(mytb.tmp[,mypv], na.rm=TRUE)
                mysum.tb[myidx,17+k] <- min   (mytb.tmp[,mypv], na.rm=TRUE)
                mysum.tb[myidx,18+k] <- mean  (-log10(mytb.tmp[,mypv]), na.rm=TRUE)
                mysum.tb[myidx,19+k] <- median(-log10(mytb.tmp[,mypv]), na.rm=TRUE)
                mysum.tb[myidx,20+k] <- max   (-log10(mytb.tmp[,mypv]), na.rm=TRUE)
            }

        } else {
            myexpr.tmp <- expr.tb[ grep(g, expr.tb[,1]) , expr.cln ]
            if (k == 1 && ! all(is.na(myexpr.tmp))) {
                myexpr.tmp <- max(expr.tb[ grep(paste0(g,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE)
            } else if (k > 1  && ! all(is.na(myexpr.tmp))) {
                myexpr.tmp <- apply(expr.tb[ grep(paste0(g,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE))
            } else {
                myexpr.tmp <- rep(NA, k)
            }

            mysum.tb[myidx,1]        <- as.character(gff.genes.qtl[ grep(g, gff.genes.qtl[,9]), 1 ])
            mysum.tb[myidx,c(2,4,5)] <- as.data.frame(gff.genes.qtl[ grep(g, gff.genes.qtl[,9]), c(9,4:5) ])
            mysum.tb[myidx,3]        <- old.id
            mysum.tb[myidx,6:7]      <- ann.tb[ grep(g,ann.tb[,1]), 2:3 ][1,]
            mysum.tb[myidx,8:(7+k)]  <- myexpr.tmp
            mysum.tb[myidx,8+k]      <- capture

        }

    }

    # Remove empty rows
    mysum.tb <- as.data.frame(mysum.tb[ ! rowSums(is.na(mysum.tb)) == ncol(mysum.tb), ])

    # Name column
    colnames(mysum.tb) <- cln.nm 

    myrows <- na.omit(rownames(mysum.tb)[mysum.tb[,1] == chr & mysum.tb[,5] >= bf.lim[1] & mysum.tb[,4] <= bf.lim[2]])
    mysum.tb.qtl <- mysum.tb[ myrows, ]

    write.table(mysum.tb.qtl, paste0(res.folder,"/",chr,"_QTL_table_summary.tsv"), row.names=FALSE, quote=FALSE, sep="\t")

    return(invisible(mysum.tb.qtl))
}



#=========#
# Figures #
#=========#

#~~~~~~~~~~~~~~~~~~~~~~~#
# Allele frequency plot #
#~~~~~~~~~~~~~~~~~~~~~~~#

myaf.plot <- function(x, lib.vec, plot.filename, ext=".png", def="HD") {

    # Picture size
    if (def == "HD") {
        mywidth <- 12
        myheight <- 8
    } else if (def == "LD") {
        mywidth <- 6
        myheight <- 4
    } else if (def != "HD" && def != "LD") {
        stop("def can be HD or LD only")
    }
    
    graph.folder <- paste("graphs_",def,"/","sd-", sd.trsh, ".gq-", gq.trsh, ".rd-", rd.trsh, "/", sep="")
    if (file.exists(graph.folder) == FALSE) {dir.create(graph.folder, recursive=TRUE)}

    Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
    
    for (i in lib.vec) {
        # Remove NA
        mycln <- grep(paste(i,".*.freq",sep=""),colnames(x))
        if (length(mycln) > 1) {
            x <- x[ is.finite(rowSums(x[,mycln])), ]
        } else {
             x <- x[ ! is.na(x[,mycln]), ]
        }

        # Start recording graphs
        ## Layout
        nb.plot <- length(mycln)*2
        nb.plot.round <- ceiling(nb.plot/2)*2
        
        png(paste(graph.folder,plot.filename,".afplot.",i,ext,sep=""), width=72*myheight*nb.plot.round/2)

        layout(matrix(1:nb.plot.round, ncol=nb.plot.round/2))
        
        ## Pairwise comparison
        mycbm <- combn(mycln,2)

        ## Plot
        for (j in 1:ncol(mycbm)) {
            ## Linear regression
            mylm <- lm(x[,mycbm[1,j]] ~ x[,mycbm[2,j]])

            plot(x[,mycbm[1,j]], x[,mycbm[2,j]], col=rgb(0,100,0,50,maxColorValue=255), xlab=colnames(x)[mycbm[1,j]], ylab=colnames(x)[mycbm[2,j]])
            abline(mylm, col="red")
            text(0, (par("usr")[4]+par("usr")[4]*0.02), paste("R²=",round(mylm$coefficients[2],2),sep=""), adj = c(0,0), xpd=TRUE)
            
            smoothScatter(x[,mycbm[1,j]], x[,mycbm[2,j]], xlab=colnames(x)[mycbm[1,j]], ylab=colnames(x)[mycbm[2,j]], colramp = Lab.palette)
            abline(mylm, col="red")
        }

        dev.off()
        
    }
}


#~~~~~~~~~~~~~~~~~~~~~#
# QQplot for Z-scores #
#~~~~~~~~~~~~~~~~~~~~~#

myz.qq <- function(x, z.vec, plot.filename, def="HD") {
    
    # Picture size
    if (def == "HD") {
        mywidth <- 12
        myheight <- 8
    } else if (def == "LD") {
        mywidth <- 6
        myheight <- 4
    } else if (def != "HD" && def != "LD") {
        stop("def can be HD or LD only")
    }
    
    graph.folder <- paste("graphs_",def,"/","sd-", sd.trsh, ".gq-", gq.trsh, ".rd-", rd.trsh, "/", sep="")
    if (file.exists(graph.folder) == FALSE) {dir.create(graph.folder, recursive=TRUE)}
    
    for (i in z.vec) {
#        mycln <- grep(i,colnames(x))
        
        png(paste(graph.folder,plot.filename,".qq-",i,".png",sep=""), width=72*mywidth, height=72*myheight)
        
        myp <- qqnorm(x[,i], main=paste0("Normal Q-Q Plot of ",i), type="n")
        smoothScatter(myp, col="darkblue", add=TRUE)
        qqline(x[,i], col="red")

        dev.off()
    }
}


#-------------

# NOT SURE that is a good idea because p-values are not independant and inflation will be likely to observed

myp.qq <- function(x, p.vec, plot.filename, def="HD") {
    
    # Picture size
    if (def == "HD") {
        mywidth <- 12
        myheight <- 8
    } else if (def == "LD") {
        mywidth <- 6
        myheight <- 4
    } else if (def != "HD" && def != "LD") {
        stop("def can be HD or LD only")
    }
    
    graph.folder <- paste("graphs_",def,"/","sd-", sd.trsh, ".gq-", gq.trsh, ".rd-", rd.trsh, "/", sep="")
    if (file.exists(graph.folder) == FALSE) {dir.create(graph.folder, recursive=TRUE)}
    
    for (i in p.vec) {
        png(paste(graph.folder,plot.filename,".qq-",i,".png",sep=""), width=72*mywidth, height=72*myheight)
        
#        # http://GettingGeneticsDone.blogspot.com/
#        # See http://gettinggeneticsdone.blogspot.com/p/copyright.html
#        o = -log10(sort(x[,i],decreasing=F))
#        e = -log10( 1:length(o)/length(o) )
#        smoothScatter(e,o, col="darkblue", xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)))
#        lines(e,e,col="red")

        qqunif.plot(x[,i])

        dev.off()
    }
}



library(lattice)
qqunif.plot<-function(pvalues, 
    should.thin=T, thin.obs.places=2, thin.exp.places=2, 
    xlab=expression(paste("Expected (",-log[10], " p-value)")),
    ylab=expression(paste("Observed (",-log[10], " p-value)")), 
    draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
    already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
    par.settings=list(superpose.symbol=list(pch=pch)), ...) {
    
    
    #error checking
    if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
    if(!(class(pvalues)=="numeric" || 
        (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
        stop("pvalue vector is not numeric, can't draw plot")
    if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
    if (already.transformed==FALSE) {
        if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
    } else {
        if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
    }
    
    
    grp<-NULL
    n<-1
    exp.x<-c()
    if(is.list(pvalues)) {
        nn<-sapply(pvalues, length)
        rs<-cumsum(nn)
        re<-rs-nn+1
        n<-min(nn)
        if (!is.null(names(pvalues))) {
            grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
            names(pvalues)<-NULL
        } else {
            grp=factor(rep(1:length(pvalues), nn))
        }
        pvo<-pvalues
        pvalues<-numeric(sum(nn))
        exp.x<-numeric(sum(nn))
        for(i in 1:length(pvo)) {
            if (!already.transformed) {
                pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
            } else {
                pvalues[rs[i]:re[i]] <- pvo[[i]]
                exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
            }
        }
    } else {
        n <- length(pvalues)+1
        if (!already.transformed) {
            exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
            pvalues <- -log10(pvalues)
        } else {
            exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
        }
    }


    #this is a helper function to draw the confidence interval
    panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
        require(grid)
        conf.points = min(conf.points, n-1);
        mpts<-matrix(nrow=conf.points*2, ncol=2)
            for(i in seq(from=1, to=conf.points)) {
                    mpts[i,1]<- -log10((i-.5)/n)
                    mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
                    mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
                    mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
            }
            grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
        }

    #reduce number of points to plot
    if (should.thin==T) {
        if (!is.null(grp)) {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                exp.x = round(exp.x, thin.exp.places),
                grp=grp))
            grp = thin$grp
        } else {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                exp.x = round(exp.x, thin.exp.places)))
        }
        pvalues <- thin$pvalues
        exp.x <- thin$exp.x
    }
    gc()
    
    prepanel.qqunif= function(x,y,...) {
        A = list()
        A$xlim = range(x, y)*1.02
        A$xlim[1]=0
        A$ylim = A$xlim
        return(A)
    }

    #draw the plot
    xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
        prepanel=prepanel, scales=list(axs="i"), pch=pch,
        panel = function(x, y, ...) {
            if (draw.conf) {
                panel.qqconf(n, conf.points=conf.points, 
                    conf.col=conf.col, conf.alpha=conf.alpha)
            };
            panel.xyplot(x,y, ...);
            panel.abline(0,1);
        }, par.settings=par.settings, ...
    )
}


#-------------------


#~~~~~~~~~~~~~~#
# P-value plot #
#~~~~~~~~~~~~~~#

mygraph <- function(data, col.vec, runmed=NULL, type="pvalue", plot.filename, ylim.min=NULL, ylim.max=NULL, ylab=NULL, def="HD") {
    
    # Picture size
    if (def == "HD") {
        mywidth <- 12
        myheight <- 8
    } else if (def == "LD") {
        mywidth <- 6
        myheight <- 4
    } else if (def != "HD" && def != "LD") {
        stop("def can be HD or LD only")
    }
    
    graph.folder <- paste("graphs_",def,"/","sd-", sd.trsh, ".gq-", gq.trsh, ".rd-", rd.trsh, "/", sep="")
    if (file.exists(graph.folder) == FALSE) {dir.create(graph.folder, recursive=TRUE)}

    # Bonferroni correction 
    mybc <- -log10(0.05/nrow(data))

    if (is.null(runmed) || is.na(runmed)) {runmed <- NULL ; myrmd <- 0} else {myrmd <- runmed}

    for (i in 1:length(col.vec)) {

        png(paste(graph.folder,plot.filename,".",colnames(data)[col.vec[i]],".rmd-",myrmd,".png",sep=""), width=72*mywidth, height=72*myheight)
        
        par(mar=c(5,4,1,1)+0.1) # For small size
        matplot.data(data, col.vec[i], type, myrunmed=runmed, ylim.min=ylim.min, ylim.max=ylim.max, abline.h=mybc, abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), ylab=ylab, by.pos=TRUE, type="p")
    
        dev.off()
    }

}

