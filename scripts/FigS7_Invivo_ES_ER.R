#!/usr/bin/env Rscript
# Title: FigS7_Invivo_ES_ER.R
# Version: 0.1
# Author: Winka Le Clec'h <winkal@txbiomed.org>
# Created in: 2020
# Modified in: 2021-08-23


#-------------------
# Loading packages
#-------------------

suppressMessages({
    library("gplots")
    library("plotrix")
})


#-----------------------------
# Loading dataset
#-----------------------------

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
pheno_fd <- "../data/phenotypes/1-Phenotyping_data/"
graph_fd <- "../graphs/"

data_mice <-read.table(paste0(pheno_fd, "6-In_vivo_PZQ_resistance_exp.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")


#===========#
# Functions #
#===========#

source("functions/line2user.R")

#-----------------------------
# Data processing
#-----------------------------

# Create a new column in the table
data_mice[,ncol(data_mice)+1] <- (data_mice[,18]/data_mice[,16])*100
colnames(data_mice)[ncol(data_mice)] <- "Normalized spleen weight"

data_mice$Normalized_spleen_weight <- data_mice[,18]/data_mice[,16]

data.mice.ES.PZQ <- data_mice[data_mice[,1] == "SmLE-PZQ-ES" & data_mice[,10] == 120, ] 
data.mice.ES.C <- data_mice[data_mice[,1] == "SmLE-PZQ-ES" & data_mice[,10] == 0,  ]

data.mice.ER.PZQ <- data_mice[data_mice[,1] == "SmLE-PZQ-ER" & data_mice[,10] == 120, ] 
data.mice.ER.C <- data_mice[data_mice[,1] == "SmLE-PZQ-ER" & data_mice[,10] == 0, ] 

#---------
# Figures
#---------

pdf(file=paste0(graph_fd, "Supp. Fig. 7 - Invivo_ES_ER_4.pdf"), width=8, height=8, useDingbats=FALSE)

layout(matrix(c(1,2,3,4),2,2))

c.type <- c("ES-PZQ","ES-Ctrl", "ER-PZQ", "ER-Ctrl")
mycolor <- c("#3bc188", "white", "#aa7b4d", "white")
myborder <- c("black", "#3bc188", "white", "#aa7b4d")

#-------------------------------#
#-----Total worms comparison----#
#-------------------------------#

par(mar=c(5,5,3,2))
boxplot(data.mice.ES.PZQ[,26], data.mice.ES.C[,26], data.mice.ER.PZQ[,26], data.mice.ER.C[,26], ylim=c(0,40), ylab= "Total number of worms", boxwex=0.6, cex.axis=1.3, col=mycolor, border=myborder, names=c.type, cex.names=1, cex.lab=1.5, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2,3,4),labels=c.type, las=1, tcl=-0.3, cex.axis=1)
axis(2, tcl=-0.3, cex.axis=1.5)

#---Add statistical analysis on the plot---#

y.maxs.ES <- c()  
y.maxs.ES <- boxplot(data.mice.ES.PZQ[,26], data.mice.ES.C[,26], plot=FALSE)[[1]][5,]

y.maxs.ER <- c()  
y.maxs.ER <- boxplot(data.mice.ER.PZQ[,26], data.mice.ER.C[,26], plot=FALSE)[[1]][5,]

y.maxs.grp.ES <- sapply(seq(1, length(y.maxs.ES)), function(x) max(y.maxs.ES[x:(x+length(y.maxs.ES)-1)]))
text(1.5, y.maxs.grp.ES*1.1, "**", xpd=TRUE, cex=1.5)
segments(1.10, y.maxs.grp.ES*1.05, 1.90)

y.maxs.grp.ER <- sapply(seq(1, length(y.maxs.ER)), function(x) max(y.maxs.ER[x:(x+length(y.maxs.ER)-1)]))
text(3.5, y.maxs.grp.ER*1.2, "NS", xpd=TRUE, cex=1.5)
segments(3.10, y.maxs.grp.ER*1.1, 3.90)


#---Add letters on the figure panel---#
mtext("A", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#-------------------------------#
#-----Total males comparison----#
#-------------------------------#

par(mar=c(5,5,3,2))
boxplot(data.mice.ES.PZQ[,24], data.mice.ES.C[,24], data.mice.ER.PZQ[,24], data.mice.ER.C[,24], ylim=c(0,40), ylab= "Number of male worms", boxwex=0.6, cex.axis=1.3, col=mycolor, border=myborder, names=c.type, cex.names=1, cex.lab=1.5, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2,3,4),labels=c.type, las=1, tcl=-0.3, cex.axis=1)
axis(2, tcl=-0.3, cex.axis=1.5)

#---Add statistical analysis on the plot---#

y.maxs.ES <- c()  
y.maxs.ES <- boxplot(data.mice.ES.PZQ[,24], data.mice.ES.C[,24], plot=FALSE)[[1]][5,]

y.maxs.ER <- c()  
y.maxs.ER <- boxplot(data.mice.ER.PZQ[,24], data.mice.ER.C[,24], plot=FALSE)[[1]][5,]

y.maxs.grp.ES <- sapply(seq(1, length(y.maxs.ES)), function(x) max(y.maxs.ES[x:(x+length(y.maxs.ES)-1)]))
text(1.5, y.maxs.grp.ES*1.25, "NS", xpd=TRUE, cex=1.5)
segments(1.10, y.maxs.grp.ES*1.1, 1.90)

y.maxs.grp.ER <- sapply(seq(1, length(y.maxs.ER)), function(x) max(y.maxs.ER[x:(x+length(y.maxs.ER)-1)]))
text(3.5, y.maxs.grp.ER*1.3, "NS", xpd=TRUE, cex=1.5)
segments(3.10, y.maxs.grp.ER*1.1, 3.90)

#---Add letters on the figure panel---#
mtext("C", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)


#---------------------------------------------------------------#
#-----Worms reduction after PZQ treatment compare to control----#
#---------------------------------------------------------------#

par(mar=c(5,5,3,2))
boxplot(data.mice.ES.PZQ[,29], data.mice.ER.PZQ[,29], ylim=c(-100,100), ylab= "Worms reduction (%)", boxwex=0.6, cex.axis=1.3, col=c("#3bc188", "#aa7b4d"), names=c("ES-PZQ","ER-PZQ"), cex.names=1, cex.lab=1.5, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2),labels=c("ES-PZQ","ER-PZQ"), las=1, tcl=-0.3, cex.axis=1)
axis(2, tcl=-0.3, cex.axis=1.5)

#---Add statistical analysis on the plot---#

y.maxs <- c()  
y.maxs <- boxplot(data.mice.ES.PZQ[,29], data.mice.ER.PZQ[,29], plot=FALSE)[[1]][5,]

y.maxs.grp <- sapply(seq(1, length(y.maxs)), function(x) max(y.maxs[x:(x+length(y.maxs)-1)]))
text(1.5, y.maxs.grp*1.3, "**", xpd=TRUE, cex=1.5)
segments(1.10, y.maxs.grp*1.1, 1.90)

#---Add letters on the figure panel---#
mtext("B", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#---------------------------------#
#-----Total females comparison----#
#---------------------------------#

par(mar=c(5,5,3,2))
boxplot(data.mice.ES.PZQ[,25], data.mice.ES.C[,25], data.mice.ER.PZQ[,25], data.mice.ER.C[,25], ylim=c(0,40), ylab= "Number of female worms", boxwex=0.6, cex.axis=1.3, col=mycolor, border=myborder, names=c.type, cex.names=1, cex.lab=1.5, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2,3,4),labels=c.type, las=1, tcl=-0.3, cex.axis=1)
axis(2, tcl=-0.3, cex.axis=1.5)

#---Add statistical analysis on the plot---#

y.maxs.ES <- c()  
y.maxs.ES <- boxplot(data.mice.ES.PZQ[,25], data.mice.ES.C[,25], plot=FALSE)[[1]][5,]

y.maxs.ER <- c()  
y.maxs.ER <- boxplot(data.mice.ER.PZQ[,25], data.mice.ER.C[,25], plot=FALSE)[[1]][5,]

y.maxs.grp.ES <- sapply(seq(1, length(y.maxs.ES)), function(x) max(y.maxs.ES[x:(x+length(y.maxs.ES)-1)]))
text(1.5, y.maxs.grp.ES*1.2, "**", xpd=TRUE, cex=1.5)
segments(1.10, y.maxs.grp.ES*1.1, 1.90)

y.maxs.grp.ER <- sapply(seq(1, length(y.maxs.ER)), function(x) max(y.maxs.ER[x:(x+length(y.maxs.ER)-1)]))
text(3.5, y.maxs.grp.ER*1.25, "NS", xpd=TRUE, cex=1.5)
segments(3.10, y.maxs.grp.ER*1.1, 3.90)

#---Add letters on the figure panel---#
mtext("D", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

##---------------------------------#
##-----Spleen weight comparison----#
##---------------------------------#

#par(mar=c(5,5,3,2))
#boxplot(data.mice.ES.PZQ[,30], data.mice.ES.C[,30], data.mice.ER.PZQ[,30], data.mice.ER.C[,30], ylim=c(0,5), ylab= "Normalized spleen weight", main="Normalized spleen weight", boxwex=0.6, cex.axis=1.3, col=mycolor, border=myborder, names=c.type, cex.names=1, cex.lab=2, outline=FALSE, axes=FALSE)

#axis(1,at=c(1,2,3,4),labels=c("ES-PZQ","ES-Control", "ER-PZQ", "ER-Control"), las=1, tcl=-0.3, cex.axis=1.3)
#axis(2, tcl=-0.3, cex.axis=1.5)

##---------------------------#
##----Sexratio comparison----#
##---------------------------#

#par(mar=c(5,5,3,2))
#boxplot(data.mice.ES.PZQ[,27], data.mice.ES.C[,27], data.mice.ER.PZQ[,27], data.mice.ER.C[,27], ylim=c(0,250), ylab= "Sexratio", main="Sexratio", boxwex=0.6, cex.axis=1.3, col=mycolor, border=myborder, names=c.type, cex.names=1, cex.lab=2, outline=FALSE, axes=FALSE)

#axis(1,at=c(1,2,3,4),labels=c("ES-PZQ","ES-Control", "ER-PZQ", "ER-Control"), las=1, tcl=-0.3, cex.axis=1.3)
#axis(2, tcl=-0.3, cex.axis=1.5)

##------------------------------------------------------------#
##-----Total worms normalized by liver weight/mouse weight----#
##------------------------------------------------------------#

#par(mar=c(5,5,3,2))
#boxplot(data.mice.ES.PZQ[,28], data.mice.ES.C[,28], data.mice.ER.PZQ[,28], data.mice.ER.C[,28], ylim=c(0,450), ylab= "Normalized number of worms", boxwex=0.6, cex.axis=1.3, col=mycolor, border=myborder, names=c.type, cex.names=1, cex.lab=2, outline=FALSE, axes=FALSE)

#axis(1,at=c(1,2,3,4),labels=c("ES-PZQ","ES-Control", "ER-PZQ", "ER-Control"), las=1, tcl=-0.3, cex.axis=1.3)
#axis(2, tcl=-0.3, cex.axis=1.5)

#-----------------------------------
# Plot saving
#-----------------------------------

dev.off()


