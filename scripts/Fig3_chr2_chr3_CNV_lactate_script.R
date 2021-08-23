#!/usr/bin/env Rscript
# Title: Fig3_chr2_chr3_CNV_lactate_script.R 
# Version: 0.1
# Author: Winka Le Clec'h <winkal@txbiomed.org>
# Created in: 2020
# Modified in: 2021-08-23


#-------------------
# Packages
#-------------------

suppressMessages({
    library("gplots")
    library("plotrix")
    library("doBy")
    library("dplyr")
})

#-----------------------------
# Loading dataset
#-----------------------------

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
geno_fd  <- "../data/phenotypes/2-Genotyping_data/"
graph_fd <- "../graphs/"

mydata <- read.csv(paste0(geno_fd, "1-TRP_chr2_chr3_CNV_indiv_worm_genotyping.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")

#===========#
# Functions #
#===========#

source("functions/line2user.R")

#-----------------------------
# Datas processing
#-----------------------------

mycolor <- c("#d8c99c", "#cdb5de", "#8cd6cf")

#----------------------------------#
#-------chromosome 2 locus---------#
#----------------------------------#

##Generate the table with means and standard deviation of lactate production for each treatment
sum_table_chr2 <- summaryBy(Lactate_variation ~ Chr2_genotype, data=mydata, FUN=c(length,mean,sd))
sum_table_chr2 <- na.omit(sum_table_chr2)

    
##Rename column lactate_prod.length.length to just N
names(sum_table_chr2)[names(sum_table_chr2)=="Lactate_variation.length"] <- "N"

## Calculate standard err_tmp of the mean
sum_table_chr2$Lactate_prod.se <- sum_table_chr2$Lactate_variation.sd / sqrt(sum_table_chr2$N)


#----------------------------------#
#-------chromosome 3 locus---------#
#----------------------------------#

##Generate the table with means and standard deviation of lactate production for each treatment
sum_table_chr3 <- summaryBy(Lactate_variation ~ Chr3_TRP_genotype, data=mydata, FUN=c(length,mean,sd))
sum_table_chr3 <- na.omit(sum_table_chr3)

##Rename column lactate_prod.length.length to just N
names(sum_table_chr3)[names(sum_table_chr3)=="Lactate_variation.length"] <- "N"

## Calculate standard err_tmp of the mean
sum_table_chr3$Lactate_prod.se <- sum_table_chr3$Lactate_variation.sd / sqrt(sum_table_chr3$N)


#----------------------------------#
#-------CNV genotyping-------------#
#----------------------------------#

##Generate the table with means and standard deviation of lactate production for each treatment
sum_table_CNV <- summaryBy(Lactate_variation ~ CNV_genotype, data=mydata, FUN=c(length,mean,sd))
sum_table_CNV <- na.omit(sum_table_CNV)

##Rename column lactate_prod.length.length to just N
names(sum_table_CNV)[names(sum_table_CNV)=="Lactate_variation.length"] <- "N"

## Calculate standard err_tmp of the mean
sum_table_CNV$Lactate_prod.se <- sum_table_CNV$Lactate_variation.sd / sqrt(sum_table_CNV$N)

## Rearrange the rows of the table
sum_table_CNV <- sum_table_CNV %>% arrange(factor(CNV_genotype, levels = c("Deletion", "1 copy", "2 copies")))


#---------
# Figures
#---------

##----Output files names and extension--------#

pdf(file=paste0(graph_fd, "Fig. 3 - chr2_chr3_CNV_indiv_worm_genotyping.pdf"), width=10, height=5, useDingbats=FALSE)

#layout(matrix(c(1,2,3),1,3))
layout(matrix(c(1,2,3,2,4,2),2,3), heights= c(0.99,0.1), widths = c(0.99,0.99,0.99))

#-------------------------------------------#
#-------Graph chromosome 2 locus------------#
#-------------------------------------------#

par(mar=c(5,6,3,3))

barplot2(sum_table_chr2[,3], col=mycolor, ylab="Change in lactate production (%)", names.arg= c("R/R", "R/S", "S/S"), space=0.5, cex.axis=1.5, cex.lab=2, cex=1.5, plot.ci=TRUE, ci.u=sum_table_chr2[,3]+sum_table_chr2[,5], ci.l=sum_table_chr2[,3]-sum_table_chr2[,5], ci.width=0.1, ci.lwd=1.4, lwd=1.4, xpd=TRUE)


##---Add statistical analysis on the plot---#
mypos <- barplot2(sum_table_chr2[,3], plot=FALSE)
x0.seg.chr2.fig <- 0.5

x1.seg.chr2.fig <- 4.5

x.chr2.text <- mean(c(x0.seg.chr2.fig, x1.seg.chr2.fig))

y.maxs.chr2.fig <- par("usr")[3] - diff(par("usr")[3:4]) * 0.15


text(x.chr2.text, y.maxs.chr2.fig, "NS", xpd=TRUE, cex=2)
segments(x0.seg.chr2.fig, y.maxs.chr2.fig * 0.95, x1.seg.chr2.fig, lwd=2, xpd=TRUE)

#---Add letters on the figure panel---#
mtext("A", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1, xpd=TRUE)

#---------------------------#
#------- X axis label-------#
#---------------------------#
# Empty plot to add the x axis label
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, "Worm genotype", cex=2)


#-------------------------------------------#
#-------Graph chromosome 3 locus------------#
#-------------------------------------------#

par(mar=c(5,2,3,3))

barplot2(sum_table_chr3[,3], col=mycolor, ylab=" ", names.arg= c("R/R", "R/S", "S/S"), space=0.5, cex.axis=1.5, cex.lab=2, cex=1.5, plot.ci=TRUE, ci.u=sum_table_chr3[,3]+sum_table_chr3[,5], ci.l=sum_table_chr3[,3]-sum_table_chr3[,5], ci.width=0.1, ci.lwd=1.4, lwd=1.4, xpd=TRUE)


##---Add statistical analysis on the plot---#
mypos <- barplot2(sum_table_chr3[,3], plot=FALSE)
x0.seg.chr3.fig <- 2

x1.seg.chr3.fig <- 4.5

x.chr3.text <- mean(c(x0.seg.chr3.fig, x1.seg.chr3.fig))

y.maxs.chr3 <- c()  
for (i in 1:nrow(sum_table_chr3)) {   
   y.maxs.chr3[i] <- (sum_table_chr3[i,3] - sum_table_chr3[i,5]) + (sum_table_chr3[i,3] - sum_table_chr3[i,5]) * 0.2
}

y.maxs.chr3.fig <- par("usr")[3] - diff(par("usr")[3:4]) * 0.15

text(1, y.maxs.chr3[1], "***", xpd=TRUE, cex=2)
text(x.chr3.text, y.maxs.chr3.fig, "NS", xpd=TRUE, cex=2)
segments(x0.seg.chr3.fig, y.maxs.chr3.fig * 0.95, x1.seg.chr3.fig, lwd=2, xpd=TRUE)

#---Add letters on the figure panel---#
mtext("B", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1, xpd=TRUE)

#-------------------------------------------#
#-------Graph chromosome CNV------------#
#-------------------------------------------#

par(mar=c(5,2,3,3))

barplot2(sum_table_CNV[,3], col=mycolor, ylab=" ", names.arg= c("Deletion", "1 copy", "2 copies"), space=0.5, cex.axis=1.5, cex.lab=2, cex=1.5, plot.ci=TRUE, ci.u=sum_table_CNV[,3]+sum_table_CNV[,5], ci.l=sum_table_CNV[,3]-sum_table_CNV[,5], ci.width=0.1, ci.lwd=1.4, lwd=1.4, xpd=TRUE)


##---Add statistical analysis on the plot---#
mypos <- barplot2(sum_table_CNV[,3], plot=FALSE)
x0.seg.CNV.fig <- 2

x1.seg.CNV.fig <- 4.5

x.CNV.text <- mean(c(x0.seg.CNV.fig, x1.seg.CNV.fig))

y.maxs.CNV <- c()  
for (i in 1:nrow(sum_table_CNV)) {   
   y.maxs.CNV[i] <- (sum_table_CNV[i,3] - sum_table_CNV[i,5]) + (sum_table_CNV[i,3] - sum_table_CNV[i,5]) * 0.2
}

y.maxs.CNV.fig <- par("usr")[3] - diff(par("usr")[3:4]) * 0.15

text(1, y.maxs.CNV[1], "**", xpd=TRUE, cex=2)
text(x.CNV.text, y.maxs.CNV.fig, "NS", xpd=TRUE, cex=2)
segments(x0.seg.CNV.fig, y.maxs.CNV.fig * 0.95, x1.seg.CNV.fig, lwd=2, xpd=TRUE)


#---Add letters on the figure panel---#
mtext("C", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1, xpd=TRUE)

##-----------------------------------
## Saving plots
##-----------------------------------


dev.off()



