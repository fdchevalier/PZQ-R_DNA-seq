#!/usr/bin/env Rscript
# Title: FigS1_Single_worm_lactate_data_sorted.R
# Version: 0.1
# Author: Winka Le Clec'h <winkal@txbiomed.org>
# Modified in: 2017
# Created in: 2021-08-09


#-------------------
# Packages
#-------------------

suppressMessages({
    library("gplots")
    library("plotrix")
})

#----------------
# Datas
#----------------

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
pheno_fd <- "../data/phenotypes/1-Phenotyping_data/"
graph_fd <- "../graphs/"

mydata <- read.csv(paste0(pheno_fd, "2-Single_worm_lactate_data_sorted.csv"), dec=".", sep=",", header=T)

#---------
# Figures
#---------

pdf(file=paste0(graph_fd, "Supp. Fig. 1 - Lactate_sorting_data.pdf"), width=6, height=6,useDingbats=FALSE)

#layout margins
par(mar=c(5,5,1,2)) 

#sorting then plotting the data according to lactate production
mydata_sorted <- mydata[order(mydata[,11]),]

#color of the points according to the group of worms
mycol <- rep("orange1", nrow(mydata_sorted))
mycollist <- mydata_sorted[,6] == "PZQ- died"
mycollist2 <- mydata_sorted[,6] == "PZQ-recovered"
mycollist3 <- mydata_sorted[,6] == "C+ killed worm"
mycollist4 <- mydata_sorted[,6] == "DMSO C-"
mycol[mycollist] <- "red"
mycol[mycollist2] <- "springgreen"
mycol[mycollist3] <- "chocolate4"
mycol[mycollist4] <- "turquoise3"

# Sorted plot
plot(sort(mydata[,11]), xlab="Samples ID", ylab="Lactate production (nmol/h)", ylim=c(0,80), col=mycol, pch=16, bty="n", axes=FALSE, cex=1.5, cex.lab=1.5)

#axes
axis(1, at=c(1:30), cex.axis=0.8, las=2)
axis(2, at=c(0,10,20,30,40,50,60,70,80), cex.axis=1)

# Legend of the figure
legend("bottomright", legend= c("Heat-killed worms - no movement","Contracted worms - no movement", "Contracted worms - some movements","Recovered worms - moving", "DMSO 0.3% treated worms - moving"), pch=16, col=c("chocolate4","red","orange1","springgreen","turquoise3"), bty="n", hor=F, cex=0.8)


#-----------------------------------
# Saving plot
#-----------------------------------

dev.off()


