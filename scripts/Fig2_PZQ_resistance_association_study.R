
##-------------------
## Packages
##-------------------

library("gplots")
library("plotrix")

#----------------
# Loading dataset
#----------------

# Load data
mydata1 <- read.csv("Lactate assay_worms_media_LE_PZQ-R_association_Exp1.csv")
mydata2 <- read.csv("Lactate assay_worms_media_LE_PZQ-R_association_Exp2.csv")

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

#------------------
# Datas processing
#------------------

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

mycoltop <- "#d8c99c"
mycolbot <- "#8cd6cf"

##---------
## Figures
##---------

pdf("PZQ association study_2.pdf", width=10, height=6, useDingbats=FALSE)

layout(matrix(c(1,2,3,2),2,2),heights= c(0.99,0.1), widths = c(0.99,0.99))

#---Experiment 1 Histogram---#
par(mar=c(2,5,2,1))

	hist(mydata1[,7], breaks=20, xlim=c(0,120), ylab="Number of male worms", cex.lab=1.5, main="Experiment 1")
	hist(mytop1[,7],  breaks=20, col=mycoltop, add=T)
	hist(mybottom1[,7], col=mycolbot, add=T)

#---Add letters on the figure panel---#
mtext("A", side=3, line=0.2, at=line2user(3,2), cex=par("cex")*2.5, adj=1)

# Empty plot to add the x axis label
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, "Quantity of L-lactate produced (nmol/h)", cex=1.5)

#---Experiment 2 Histogram---#

par(mar=c(2,2,2,2))

	hist(mydata2[,7],   breaks=20, xlim=c(0,120), ylab=" ", main="Experiment 2")
	hist(mytop2[,7],    breaks=20, col=mycoltop, add=T)
	hist(mybottom2[,7], breaks=5,  col=mycolbot, add=T)



dev.off()

##---------------------------------------------------------------------------------
## Export list of chosen worms for preparing gDNA pools (contracted vs. recovered) 
##---------------------------------------------------------------------------------

write(mybottom1[,8], "Exp1 - bottom.txt")
write(mytop1[,8],    "Exp1 - top.txt")

write(mybottom2[,8], "Exp2 - bottom.txt")
write(mytop2[,8],    "Exp2 - top.txt")
