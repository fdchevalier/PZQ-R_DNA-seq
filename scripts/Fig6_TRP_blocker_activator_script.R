#!/usr/bin/env Rscript
# Title: Fig6_TRP_blocker_activator_script.R
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
pheno_fd <- "../data/phenotypes/1-Phenotyping_data/"
graph_fd <- "../graphs/"

mydata <- read.csv(paste0(pheno_fd, "7-TRP_MB2_blocker_MV1_activator_ER_ES.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")

#===========#
# Functions #
#===========#

source("functions/line2user.R")

#-----------------------------
# Datas processing
#-----------------------------

data.ES <- mydata[mydata[,1] == "SmLE-PZQ-ES", ] 
data.ER <- mydata[mydata[,1] == "SmLE-PZQ-ER", ]

#mypop <- c("SmLE-PZQ-ES", "SmLE-PZQ-ER")
c.type <- c("PZQ","DMSO")
c.treat.colnames <- c("Control", "TRP blocker", "TRP activator (10uM)", "TRP activator (25uM)", "TRP activator (50uM)")
mycolor.type <- c("grey", "white")

## Create an empty list to stock results
res.ES <- vector("list", length(c.type)) 
names(res.ES) <- c.type

for (i in c.type){
        mydata_tmp <- data.ES[data.ES[,6] == i,]
    
        ##Generate the table with means and standard deviation of lactate production for each treatment
        sum_table_tmp <- summaryBy(Lactate_production ~ TRP_treatment, data=mydata_tmp, FUN=c(length,mean,sd))
    
        ## Rename column Lactate_production.length.length to just N
        names(sum_table_tmp)[names(sum_table_tmp)=="Lactate_production.length"] <- "N"

        ## Calculate standard err_tmp of the mean
        sum_table_tmp$Lactate_production.se <- sum_table_tmp$Lactate_production.sd / sqrt(sum_table_tmp$N)

        # Create the list with 3 slots
        res_tmp <- vector("list", 1)
    
        #Extraction of the SE of the means from the summary table _tmp
        res_tmp[[1]] <- sum_table_tmp
        
        res.ES[[i]] <- res_tmp
    
    
}

ES_PZQ_DMSO <- sapply(res.ES, function(x) x[[1]][,3]) %>% t()
colnames(ES_PZQ_DMSO) <- c.treat.colnames

myse_ES_PZQ_DMSO <- sapply(res.ES, function(x) x[[1]][,5]) %>% t()
colnames(myse_ES_PZQ_DMSO) <- c.treat.colnames
 
## Create an empty list to stock results
res.ER <- vector("list", length(c.type)) 
names(res.ER) <- c.type
for (i in c.type){
        mydata_tmp <- data.ER[data.ER[,6] == i,]
    
        ##Generate the table with means and standard deviation of lactate production for each treatment
        sum_table_tmp <- summaryBy(Lactate_production ~ TRP_treatment, data=mydata_tmp, FUN=c(length,mean,sd))
    
        ## Rename column Lactate_production.length.length to just N
        names(sum_table_tmp)[names(sum_table_tmp)=="Lactate_production.length"] <- "N"

        ## Calculate standard err_tmp of the mean
        sum_table_tmp$Lactate_production.se <- sum_table_tmp$Lactate_production.sd / sqrt(sum_table_tmp$N)

        # Create the list with 3 slots
        res_tmp <- vector("list", 1)
    
        #Extraction of the SE of the means from the summary table _tmp
        res_tmp[[1]] <- sum_table_tmp
        
        res.ER[[i]] <- res_tmp
    
    
}

ER_PZQ_DMSO <- sapply(res.ER, function(x) x[[1]][,3]) %>% t()
colnames(ER_PZQ_DMSO) <- c.treat.colnames

myse_ER_PZQ_DMSO <- sapply(res.ER, function(x) x[[1]][,5]) %>% t()
colnames(myse_ER_PZQ_DMSO) <- c.treat.colnames

#---------
# Figures
#---------

##----Output files names and extension--------#

pdf(file=paste0(graph_fd, "Fig. 6 - Barplot_lactate_MB2_MV1_ES_ER_3.pdf"), width=8, height=10, useDingbats=FALSE)

layout(matrix(c(1,2),2,1), heights= c(0.5,0.5))

#-------------------------------------------#
#-------Graph ES grouped by treatment-------#
#-------------------------------------------#

par(mar=c(5,5,2,2))

barplot2(ES_PZQ_DMSO, col=mycolor.type, main= "PZQ-ES", beside=TRUE, ylab="Lactate production (nmol/h)", ylim=c(0,60), cex.axis=1, cex.main=1.5, cex.lab=1.2, cex=1, las=2, names.arg=c("", "", "", "", ""), plot.ci=TRUE, ci.u=ES_PZQ_DMSO+myse_ES_PZQ_DMSO, ci.l=ES_PZQ_DMSO-myse_ES_PZQ_DMSO, ci.width=0.1, ci.lwd=1.4, lwd=1.4, legend.text = c("+PZQ", "+DMSO")) 


#---Add statistical analysis on the plot---#
x0.seg.ES <- barplot2(ES_PZQ_DMSO, beside=TRUE, plot=FALSE)
x0.seg.ES <- x0.seg.ES[1, ] - 0.2

x1.seg.ES <- barplot2(ES_PZQ_DMSO, beside=TRUE, plot=FALSE)
x1.seg.ES <- x1.seg.ES[2, ] + 0.2

x.ES.text <- c(2,5,8,11,14)

y.maxs.ES <- c()  
for (i in 1:length(c.treat.colnames)) {   
   y.maxs.ES[i] <- max(ES_PZQ_DMSO[,i]) + max(myse_ES_PZQ_DMSO[,i]) 
}

text(x.ES.text, y.maxs.ES+3.5, c("***", "NS", "***", "***", "***"), xpd=TRUE, cex=1.5)
segments(x0.seg.ES, y.maxs.ES+1.5, x1.seg.ES, xpd=TRUE)

#---Add letters on the figure panel---#
mtext("A", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

#-------------------------------------------#
#-------Graph ER grouped by treatment-------#
#-------------------------------------------#

par(mar=c(9,5,1.5,2))

barplot2(ER_PZQ_DMSO, col=mycolor.type, main= "PZQ-ER", beside=TRUE, ylab="Lactate production (nmol/h)", ylim=c(0,60), cex.main=1.5, cex.axis=1, cex.lab=1.2, cex=1, las=2, plot.ci=TRUE, ci.u=ER_PZQ_DMSO+myse_ER_PZQ_DMSO, ci.l=ER_PZQ_DMSO-myse_ER_PZQ_DMSO, ci.width=0.1, ci.lwd=1.4, lwd=1.4, xpd=TRUE)


#---Add statistical analysis on the plot---#
x0.seg.ER <- barplot2(ER_PZQ_DMSO, beside=TRUE, plot=FALSE)
x0.seg.ER <- x0.seg.ER[1, ] - 0.2

x1.seg.ER <- barplot2(ER_PZQ_DMSO, beside=TRUE, plot=FALSE)
x1.seg.ER <- x1.seg.ER[2, ] + 0.2

x.ER.text <- c(2,5,8,11,14)

y.maxs.ER <- c()  
for (i in 1:length(c.treat.colnames)) {   
   y.maxs.ER[i] <- max(ER_PZQ_DMSO[,i]) + max(myse_ER_PZQ_DMSO[,i]) 
}

text(x.ER.text, y.maxs.ER+3.9, c("NS", "NS", "***", "***", "***"), xpd=TRUE, cex=1.5)
segments(x0.seg.ER, y.maxs.ER+1.5, x1.seg.ER, xpd=TRUE)


#---Add letters on the figure panel---#
mtext("B", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

##-----------------------------------
## Saving plots
##-----------------------------------


dev.off()



