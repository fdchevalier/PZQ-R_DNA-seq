#==============#
# Dependencies #
#==============#

# Mandatory object
#if (! exists("call.tb.ls")) { stop("call.tb.ls not found.") }

# Package
library("magrittr")
library("GenomicFeatures")
library("rtracklayer")
library("Sushi")

# v3.2
source("functions/Sm.matplot.data.R")
source("functions/rename_chr.R")


#====================#
# Data and variables #
#====================#

# Load data
mydata <- read.csv("../results/1-QTL/sd-0.1.gq-0.rd-10/peaks_report.myfreq.data.fltr.AB4-freq-pv.sd-0.1-bf=0", header=TRUE, sep="\t")

# mychr <- "Chr_3"
mychr <- "SM_V7_3"
exp.tag <- c("Exp1","Exp2")

# Load GFF
# mygff.fl <- "~/data/sm_Gene_table/Sm_v7.1_renamed.gff"
mygff.fl <- "~/data/sm_Gene_table/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3"
mygff <- readGFF(mygff.fl)
mygff.genes <- mygff[ mygff[,3] == "gene", ]

# Gene expression
# myexpr <- read.csv("~/data/sm_Gene_table/TPM_isoforms_Sm_v7.1.tsv", header=TRUE, sep="\t")
myexpr <- read.csv("~/analyses/00 - Misc/2018-04-19 v7 genome gene expression/TPM_isoforms_Sm_2020-08-07.tsv", header=TRUE, sep="\t")


myseed <- 553309


#=================#
# Data processing #
#=================#

#------------------#
# Exons extraction #
#------------------#

c <- makeTxDbFromGFF(mygff.fl)
d <- exonsBy(c, use.names=T)

e <- as.data.frame(d)

# Add fake score
e[,ncol(e)+1] <- "."

# Rename transcripts into genes
e[,2] <- lapply(strsplit(e[,2], ".", fixed=T), function(x) x[1]) %>% unlist()
e <- unique(e)

# Remove first column
e <- e[,-1]

# Reorder data
e <- e[,c(2:4,1,ncol(e),6)]

#levels(e[,6])[grep("+", levels(e[,6]), fixed=T)] <- "+1"
#levels(e[,6])[grep("-", levels(e[,6]), fixed=T)] <- "-1"

e[,6] <- as.numeric(e[,6])
e[e[,6] == 2, 6] <- -1
e[e[,6] == 3, 6] <- 0

e[,1] <- as.character(e[,1])

colnames(e) <- c("chrom", "start", "stop", "gene", "score", "strand")


#----------------#
# Table redesign #
#----------------#

# Get combined p.values and build table
myp.val.cln <- rev(grep(paste0(".*",paste0(exp.tag,sep=".*", collapse="")), colnames(mydata), value=T))[1]
myp.val.tb  <- cbind(mydata[,1:2], mydata[,myp.val.cln])
myp.val.tb  <- myp.val.tb[! is.na(myp.val.tb[, 3]), ]

# Bonferroni correction
# bf.cor <- -log10(0.05/nrow(mydata))
bf.cor <- 0.05/nrow(myp.val.tb)
my.qtl <- myp.val.tb[myp.val.tb[,3] <= bf.cor & myp.val.tb[,1] == mychr, 2]
# bf.lim <- c(my.qtl[1], my.qtl[length(my.qtl)])
bf.lim <- c(my.qtl[1], my.qtl[which(diff(my.qtl) > 1e6)[1]])

# Sub-sampling
# set.seed(myseed)
# myrow <- rownames(myp.val.tb[myp.val.tb[,3] > bf.cor, ]) %>% sample(., dim(myp.val.tb)[1] * 0.1)
# myrow <- rownames(myp.val.tb[myp.val.tb[,3] <= bf.cor, ]) %>% c(., myrow) %>% as.numeric() %>% sort() %>% as.character()

# myp.val.tb <- myp.val.tb[myrow, ]
myp.val.tb  <- myp.val.tb[! is.na(myp.val.tb[, 3]), ]


#-------#
# Graph #
#-------#


my.x1 <- 1
my.x2 <- 5e6

my.x1b <- 6.5e5
my.x2b <- 1.25e6
#my.x2b <- 7.5e5

# mygff.genes <- mygff.genes[ mygff.genes[,1] == mychr & mygff.genes[,4] >= bf.lim[1] & mygff.genes[,5] <= bf.lim[2] , ]
mygff.genes <- mygff.genes[ mygff.genes[,1] == mychr & mygff.genes[,4] >= my.x1 & mygff.genes[,5] <= my.x2 , ]

# Gene expression
mygenes <- mygff.genes[, "Name"]
mycol   <- NULL
for (i in mygenes) {
    # exp.gene <- (myexpr[ myexpr[,1] == i, 11:12] > 0) %>% sum(.) > 0
    exp.gene <- (myexpr[ grepl(paste0(i, "\\."), myexpr[,1]), 14:15] > 0) %>% sum(.) > 0
    if (exp.gene) {
        mycol[i] <- "black"
    } else {
        mycol[i] <- "gray"
    }
}

# Data restricted to QTL
myp.val.tb.chr     <- myp.val.tb[ myp.val.tb[,1] == mychr & myp.val.tb[,2] <= my.x2, ]
myp.val.tb.chr[,3] <- -log10(myp.val.tb.chr[,3])

myp.val.tb <- rename_chr_SmV7(myp.val.tb, 1)

#png(paste0("graphs/zoom_",mychr,"_QTL.png"), width=13*72, height=16*72)
#pdf(paste0("graphs/zoom_",mychr,"_QTL.pdf"), width=8, height=10, useDingbats=FALSE)
# pdf(paste0("../graphs/zoom_",mychr,"_QTL_v3.pdf"), width=8, height=10, useDingbats=FALSE)
png(paste0("../graphs/zoom_",mychr,"_QTL_v3_2.png"), width=8*72, height=10*72)
    
#    layout(matrix(1:4, ncol=1), height=c(1,0.7,0.3,1))
    layout(matrix(1:3, ncol=1), height=c(1,0.7,0.3))
    
    par(mar=c(5,5,2,2)+0.1)
    mycex.axis <- 1.2
    mycex <- 1.5
    
    # Plot 1: the complete data
    matplot.data(myp.val.tb, ncol(myp.val.tb), "pvalue", ylim.min=0, ylim.max=22, abline.h=-log10(bf.cor), abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), by.pos=TRUE, type="p", cex.axis=mycex.axis) #, myrunmed=5) #, cex=mycex)

    ## Compute shift to put box (length of the previous chr)
    myshift <- 0
    for (i in c("Chr_1", "Chr_2", "Chr_2.un.SC_H001")) { myshift <- sum(myshift, myp.val.tb[ tail(which(myp.val.tb[,1] == i),1), 2 ]) }

    zoomsregion(c(my.x1, my.x2) + myshift)


    # Plot 2: p-values
    par(mar=c(0,5,4,2)+0.1)


    plot(myp.val.tb.chr[,2:3], xlim=c(my.x1,my.x2), type="n", bty="n", axes=FALSE, ann=FALSE)
#    rect(bf.lim[1],par("usr")[3],bf.lim[2],par("usr")[4], col="grey", border=NA)
    points(myp.val.tb.chr[,2:3], col="red", pch=19)

    ablineclip(h=-log10(bf.cor), x1=0, x2=my.x2,, col="blue", lty=3, lwd=2)

#    rect(my.x1b, par("usr")[3]-10, my.x2b, par("usr")[4], lty=2)
    
    ## Annotations
    axis (2, at=format(seq(0,21,length=5),digits=2), cex.axis=mycex.axis)
    mtext(2, text=expression(-Log[10]*" p-value"), line=3) #, cex=mycex)

    # Plot 3: genes
    par(mar=c(5,5,0,2)+0.1)
    plot(0, xlim=c(my.x1,my.x2), ylim=c(-.25, .25), type="n", bty="n", axes=FALSE, ann=FALSE)
    my.y <- par("usr")[4]*0.2

    for(g in 1:nrow(mygff.genes)) {
        rect(mygff.genes[g,4], -my.y, mygff.genes[g,5], my.y, col=mycol[g], border=NA)
    }

    axis (1, at=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1), labels=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1)/1e6, cex.axis=mycex.axis)
    mtext(1, text="Position (Mb)", line=3) #, cex=mycex)


##    zoomsregion(c(my.x1b, my.x2b), extend=c(10,0))
#    zoomsregion(c(my.x1b, my.x2b))
#
#    # Plot 4: genes (and not isoforms!!)
#    par(mar=c(5,5,2,2)+0.1)
#    e.tmp <- e[ e[,1] == mychr & (e[,2] >= my.x1b & e[,3] <= my.x2b) , ]
#    
#    ## Add expression data from adult males
#    for (i in unique(e.tmp[,4])) { 
#        myexpr.tmp <- max(myexpr[ grepl(i, myexpr[,1]) , 11], na.rm=T)
#        if (length(myexpr.tmp) > 0) { e.tmp[ e.tmp[,4] == i , 5 ] <- log10(myexpr.tmp) }
#    }
#
#    e.tmp[,5] <- as.numeric(e.tmp[,5])
#
#    e.tmp <- e.tmp[order(e.tmp[,2],e.tmp[,4]),]
#    
#    pg <- plotGenes(e.tmp, mychr, my.x1b, my.x2b, labeltext=TRUE, colorby=e.tmp[,5], colorbycol=SushiColors(5), colorbyrange=c(0,max(e.tmp[,5])), fontsize=0.5, plotgenetype="arrow", bentline=FALSE, wigglefactor=0.1)
#    addlegend(pg[[1]], palette=pg[[2]], title=expression(Log[10] * "TPM"), side="left", xoffset=-0.1, bottominset=0.5)
#    
#    axis (1, at=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1), labels=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1)/1e6)
#    mtext(1, text="Position (Mb)", line=3)
 
dev.off()

