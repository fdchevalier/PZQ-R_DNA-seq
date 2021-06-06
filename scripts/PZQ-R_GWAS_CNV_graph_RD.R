#!/usr/bin/env Rscript
# Title: PZQ-R_GWAS_CNV_graph_RD.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-06-05
# Modified in:



#==========#
# Comments #
#==========#

# Plot read depth of the QTL of chr. 3



#==========#
# Versions #
#==========#

# v0.0 - 2021-06-05: creation



#===========#
# Variables #
#===========#

# Folders
data_fd   <- "../data/libraries/1-GWAS/"
graph_fd  <- "../graphs/"

# Samples to plot
myspl <- dir(data_fd)

# Regions to plot
regions <- list(c(7e5, 1.05e6), c(1.15e6, 1.4e6))

myseed <- 1542268



#=================#
# Data processing #
#=================#

cat("Loading and processing data. This may take a while...\n")

# Load data
mydata <- vector("list", length(myspl))
names(mydata) <- myspl
for (i in myspl) {
    mydata[[i]] <- read.delim(paste0(data_fd,  i, "/", i, "_3_QTL.cov"), header = FALSE)
}

mydata_ls <- vector("list", length(regions))
graph_ls  <- vector("list", length(regions))

for (r in 1:length(regions)) {
    mystart <- regions[[r]][1]
    myend   <- regions[[r]][2]

    # Focusing on the regions of interest
    mydata_tmp <- lapply(mydata, function(x) { x <- x[ x[, 2] > mystart & x[ ,2] < myend, ] ; return(x) })

    # Smoothing data
    mydata_tmp <- lapply(mydata_tmp, function(x) { x[, 3] <- runmed(x[, 3], 2001) ; return(x) })

    # Reducing data by removing low covered regions
    mydata_tmp <- lapply(mydata_tmp, function(x) { x <- x[x[, 3] > 5, ] ; return(x) })

    # Sub-sampling
    set.seed(myseed)
    mydata_tmp <- lapply(mydata_tmp, function(x) x[sort(sample(1:nrow(x), nrow(x)/200)), ])

    # x max position
    x_max <- max(unlist(lapply(mydata_tmp, function(x) max(x[,2]))))

    # x min position
    x_min <- max(unlist(lapply(mydata_tmp, function(x) min(x[,2]))))

    # y max position
    y_max <- max(unlist(lapply(mydata_tmp, function(x) max(x[,3]))))

    # Stroring data
    mydata_ls[[r]] <- mydata_tmp
    graph_ls[[r]]  <- c(x_max, x_min, y_max)
}


#=========#
# Figures #
#=========#

cat("\nDrawing graphs. This may take a while...\n")

png(paste0(graph_fd, "Supp. Fig. 1.png"), width = 5 * length(myspl) / 2, height = 3 * length(regions) * 2, unit = "in", res = 300)
layout(matrix(1:(length(myspl) * length(regions)), ncol = length(myspl) / 2, byrow = TRUE))

for (r in 1:length(regions)) {
    mydata_tmp <- mydata_ls[[r]]

    x_max <- graph_ls[[r]][1]
    x_min <- graph_ls[[r]][2]
    y_max <- graph_ls[[r]][3]

    for (i in myspl) {
        # Point graph
        plot(mydata_tmp[[i]][,3] ~ mydata_tmp[[i]][,2], xlab = "Position on chromosome Z (Mb)", ylab = "Read depth", xlim = c(x_min, x_max), ylim = c(5, y_max), main = i, xaxt = "n", pch = 20, col = "grey", log = "y")
        #axis(1, at = seq(0, 80, 10) * 1e6, labels = seq(0, 80, 10))
        axis(1)

        # Smoothed line
        lw1 <- loess(mydata_tmp[[i]][,3] ~ mydata_tmp[[i]][,2], span=0.05)
        lines(mydata_tmp[[i]][,2], predict(lw1, mydata_tmp[[i]][,2]), col = "red", lwd = 3)
    }
}
dev.off()
