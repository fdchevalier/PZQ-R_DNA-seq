
##-------------------
## Packages
##-------------------

library("gplots")
library("plotrix")

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
source("functions/line2user.R")



#====================#
# Data and variables #
#====================#

# Load data
mydata1 <- read.csv("../data/phenotype/Lactate assay_worms_media_LE_PZQ-R_association_Exp1.csv")
mydata2 <- read.csv("../data/phenotype/Lactate assay_worms_media_LE_PZQ-R_association_Exp2.csv")

mydata <- read.csv("../results/1-QTL/GWASsd-0.1.gq-0.rd-10/peaks_report.myfreq.data.fltr.AB4-freq-pv.sd-0.1-bf=0", header=TRUE, sep="\t")

# mychr <- "Chr_3"
mychr <- "SM_V7_3"
exp.tag <- c("Exp1","Exp2")

# Load GFF
mygff.fl <- "../data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3"
mygff <- readGFF(mygff.fl)
mygff.genes <- mygff[ mygff[,3] == "gene", ]

# Gene expression
myexpr <- read.csv("../data/genome/Sm_transcript_table_gff-hhpred_2020-09-21.tsv", header=TRUE, sep="\t")

# Gene of interest
goi <- "Smp_246790"

myseed <- 553309

mycoltop <- "#d8c99c"
mycolbot <- "#8cd6cf"


#=================#
# Data processing #
#=================#

# Add unique ID
mydata1[, ncol(mydata1)+1] <- paste0("P", mydata1[,4], " - ", mydata1[,5])
mydata2[, ncol(mydata2)+1] <- paste0("P", mydata2[,4], " - ", mydata2[,5])

# Remove individuals with too low quantity of lactate before treatment 
mydata1 <- mydata1[ mydata1[,6] >= 40, ]
mydata2 <- mydata2[ mydata2[,6] >= 40, ]

## Relative difference (normalized by starting [lactate])
#pzq.lactacte <- (mydata[,7]- mydata[,6])/mydata[,6]

# Ordering data on the lactate quantity after PZQ treatment
mydata1 <- mydata1[order(mydata1[,7]),]
mydata2 <- mydata2[order(mydata2[,7]),]

# Number of individuals corresponding to 20% cut-off
mynb1 <- round(nrow(mydata1)*0.2)
mynb2 <- round(nrow(mydata2)*0.2)

# Get individuals
mybottom1 <- mydata1[1:mynb1,]
mytop1	  <- mydata1[(nrow(mydata1)-mynb1+1):nrow(mydata1),]

mybottom2 <- mydata2[1:mynb2,]
mytop2	  <- mydata2[(nrow(mydata2)-mynb2+1):nrow(mydata2),]


##---------------------------------------------------------------------------------
## Export list of chosen worms for preparing gDNA pools (contracted vs. recovered) 
##---------------------------------------------------------------------------------

write(mybottom1[,8], "Exp1 - bottom.txt")
write(mytop1[,8],    "Exp1 - top.txt")

write(mybottom2[,8], "Exp2 - bottom.txt")
write(mytop2[,8],    "Exp2 - top.txt")




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
bf.cor <- 0.05/nrow(myp.val.tb)
my.qtl <- myp.val.tb[myp.val.tb[,3] <= bf.cor & myp.val.tb[,1] == mychr, 2]
bf.lim <- c(my.qtl[1], my.qtl[which(diff(my.qtl) > 1e6)[1]])

myp.val.tb  <- myp.val.tb[! is.na(myp.val.tb[, 3]), ]


#-------#
# Graph #
#-------#


my.x1 <- 1
my.x2 <- 5e6

mygff.genes <- mygff.genes[ mygff.genes[,1] == mychr & mygff.genes[,4] >= my.x1 & mygff.genes[,5] <= my.x2 , ]

# Gene expression
mygenes <- mygff.genes[, "Name"]
mycol   <- NULL
for (i in mygenes) {
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

# pdf(paste0("../graphs/zoom_",mychr,"_QTL_v3.pdf"), width=8, height=10, useDingbats=FALSE)
png(paste0("../graphs/zoom_",mychr,"_QTL_v3_2.png"), width=8*72, height=10*72)
# pdf("PZQ association study_2.pdf", width=10, height=6, useDingbats=FALSE)

    layout(matrix(c(1,2,3,2),2,2),heights= c(0.99,0.1), widths = c(0.99,0.99))
    layout(matrix(c(1:3, 3, 4, 4, 5, 5, 6, 6), ncol = 2), height=c(1, 0.1, 1, 0.7, 0.3))

    #---Experiment 1 Histogram---#
    par(mar=c(2,5,2,1))

	hist(mydata1[,7], breaks=20, xlim=c(0,120), ylab="Number of male worms", cex.lab=1.5, main="Experiment 1")
	hist(mytop1[,7],  breaks=20, col=mycoltop, add=T)
	hist(mybottom1[,7], col=mycolbot, add=T)

    #---Add letters on the figure panel---#
    mtext("A", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1)


    #---Experiment 2 Histogram---#

    par(mar=c(2,2,2,2))

	hist(mydata2[,7],   breaks=20, xlim=c(0,120), ylab=" ", main="Experiment 2")
	hist(mytop2[,7],    breaks=20, col=mycoltop, add=T)
	hist(mybottom2[,7], breaks=5,  col=mycolbot, add=T)

    # Empty plot to add the x axis label
    par(mar = c(0,0,0,0))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    text(x = 0.5, y = 0.5, "Quantity of L-lactate produced (nmol/h)", cex=1.5)


    
    
    par(mar=c(5,5,2,2)+0.1)
    mycex.axis <- 1.2
    mycex <- 1.5
    
    # Plot 1: the complete data
    matplot.data(myp.val.tb, ncol(myp.val.tb), "pvalue", ylim.min=0, ylim.max=22, abline.h=-log10(bf.cor), abline.lwd=2, xlab.axis=c("1","2","3","4","5","6","7","Z","Unass. sc."), by.pos=TRUE, type="p", cex.axis=mycex.axis) #, myrunmed=5) #, cex=mycex)

    ## Compute shift to put box (length of the previous chr)
    myshift <- 0
    for (i in c("Chr_1", "Chr_2", "Chr_2.un.SC_H001")) { myshift <- sum(myshift, myp.val.tb[ tail(which(myp.val.tb[,1] == i),1), 2 ]) }

    zoomsregion(c(my.x1, my.x2) + myshift)
    
    mtext("B", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

    # Plot 2: p-values
    par(mar=c(0,5,4,2)+0.1)


    plot(myp.val.tb.chr[,2:3], xlim=c(my.x1,my.x2), type="n", bty="n", axes=FALSE, ann=FALSE)
    points(myp.val.tb.chr[,2:3], col="red", pch=19)

    ablineclip(h=-log10(bf.cor), x1=0, x2=my.x2,, col="blue", lty=3, lwd=2)

    
    ## Annotations
    axis (2, at=format(seq(0,21,length=5),digits=2), cex.axis=mycex.axis)
    mtext(2, text=expression(-Log[10] ~ p-value), line=3) #, cex=mycex)
    mtext("C", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

    # Plot 3: genes
    par(mar=c(5,5,0,2)+0.1)
    plot(0, xlim=c(my.x1,my.x2), ylim=c(-.25, .25), type="n", bty="n", axes=FALSE, ann=FALSE)
    my.y <- par("usr")[4]*0.2

    for(g in 1:nrow(mygff.genes)) {
        rect(mygff.genes[g,4], -my.y, mygff.genes[g,5], my.y, col=mycol[g], border=NA)

        if (mygff.genes[g, "Name"] == goi) {
            old_xpd <- par("xpd")
            par(xpd = NA)

            bx_y <- ((par("usr")[4] - par("usr")[3]) * 0.05) %>% abs()
            bx_x <- bx_y * 10^floor(log10(diff(par("usr")[1:2])))
            
            rect(mygff.genes[g,4] - bx_x, -my.y - bx_y, mygff.genes[g,5] + bx_x, 2.25, lty = 2)

            par(xpd = old_xpd)
        }
    }



    axis (1, at=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1), labels=seq(par("xaxp")[1], par("xaxp")[2], length=par("xaxp")[3]+1)/1e6, cex.axis=mycex.axis)
    mtext(1, text="Position (Mb)", line=3) #, cex=mycex)


dev.off()

