#-------------------
# Loading packages
#-------------------

library("gplots")
library("plotrix")
#library("Hmisc")


#-----------------------------
# Loading dataset
#-----------------------------

data_snails <-read.table("Data SmLE-PZQ-ER-ES_snails_infections.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
data_hamsters <-read.table("Data SmLE-PZQ-ER-ES_hamsters_infections.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")

#===========#
# Functions #
#===========#

# Line in units
# source: https://stackoverflow.com/a/30835971
line2user <- function(line, side) {
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
        `1` = grconvertY(-line * y_off, 'npc', 'user'),
        `2` = grconvertX(-line * x_off, 'npc', 'user'),
        `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
        `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
        stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

#-----------------------------
# Data processing
#-----------------------------

data.snail.ES <- data_snails[data_snails[,1] == "SmLE-PZQ-ES", ] 
#data.snail.ES <- data.snail.ES[data.snail.ES[,3] == "Bg36", ]

data.snail.ER <- data_snails[data_snails[,1] == "SmLE-PZQ-ER", ] 
#data.snail.ER <- data.snail.ER[data.snail.ER[,3] == "Bg36", ]

data.hamsters.ES <- data_hamsters[data_hamsters[,1] == "SmLE-PZQ-ES", ] 
data.hamsters.ER <- data_hamsters[data_hamsters[,1] == "SmLE-PZQ-ER", ] 

#---------
# Figures
#---------

pdf(file="fitness_ES_ER_4.pdf", width=10, height=4, useDingbats=FALSE)

layout(matrix(c(1,2,3),1,3))

c.type <- c("SmLE-PZQ-ES","SmLE-PZQ-ER")
mycolor <- c("#3bc188","#aa7b4d")

#-----------------------#
#-----Snail survival----#
#-----------------------#

par(mar=c(5,5,3,2))
boxplot(data.snail.ES[,8], data.snail.ER[,8], ylim=c(0,100), ylab= "Survival (%)", main="Snail survival", boxwex=0.6, cex.axis=1.3, col=mycolor, names=c.type, cex.names=1, cex.lab=2, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2),labels=c("PZQ-ES", "PZQ-ER"), las=1, tcl=-0.3, cex.axis=1.3)
axis(2, tcl=-0.3, cex.axis=1.5)


#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.type)) {   
       inds <- which(data_snails[,1] == c.type[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(data.snail.ES[,8][inds], data.snail.ER[,8][inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.1, "NS", xpd=TRUE, cex=1.5)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*1.03, 1.90)
}

#---Add letters on the figure panel---#
mtext("A", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2.5, adj=1)


#--------------------------#
#----Snails infectivity----#
#--------------------------#

par(mar=c(5,5,3,2))
boxplot(data.snail.ES[,7], data.snail.ER[,7], ylim=c(0,100), ylab= "Infectivity (%)", main="Snail infectivity", boxwex=0.6, cex.axis=1.3, col=mycolor, names=c.type, cex.names=1, cex.lab=2, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2),labels=c("PZQ-ES", "PZQ-ER"), las=1, tcl=-0.3, cex.axis=1.3)
axis(2, tcl=-0.3, cex.axis=1.5)

#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.type)) {   
       inds <- which(data_snails[,1] == c.type[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(data.snail.ES[,7][inds], data.snail.ER[,7][inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.1, "NS", xpd=TRUE, cex=1.5)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*1.03, 1.90)
}

#---Add letters on the figure panel---#
mtext("B", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

#----------------------------#
#----Hamsters infectivity----#
#----------------------------#

par(mar=c(5,5,3,2))
boxplot(data.hamsters.ES[,7], data.hamsters.ER[,7], ylim=c(0,100), ylab= "Infectivity (%)", main="Hamster infectivity", boxwex=0.6, cex.axis=1.3, col=mycolor, names=c.type, cex.names=1, cex.lab=2, outline=FALSE, axes=FALSE)

axis(1,at=c(1,2),labels=c("PZQ-ES", "PZQ-ER"), las=1, tcl=-0.3, cex.axis=1.3)
axis(2, tcl=-0.3, cex.axis=1.5)

#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.type)) {   
       inds <- which(data_hamsters[,1] == c.type[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(data.hamsters.ES[,7][inds], data.hamsters.ER[,7][inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.01, "NS", xpd=TRUE, cex=1.5)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*0.9, 1.90)
}

#---Add letters on the figure panel---#
mtext("C", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

#-----------------------------------
# Plot saving
#-----------------------------------

dev.off()
