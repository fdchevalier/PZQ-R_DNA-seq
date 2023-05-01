#!/usr/bin/env Rscript
# Title: PZQ-R_GWAS_CNV_graph_RD.R
# Version: 0.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-06-05
# Modified in: 2023-05-01



#==========#
# Comments #
#==========#

# Plot read depth of the QTL of chr. 3



#==========#
# Versions #
#==========#

# v0.3 - 2023-05-01: adapt the PZQ-R_GWAS_CNV_graph_RD.R to S. mansoni v10 genome
# v0.1 - 2021-08-27: add missing code / update file path / load data more efficiently
# v0.0 - 2021-06-05: creation



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    library("magrittr")
    library("data.table")
})



#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd   <- "../data/libraries/1-GWAS/"
graph_fd  <- "../graphs/1-GWAS/"

mysfx     <- "v10"
graph_fd  <- paste0(graph_fd, "/", mysfx, "/")

# Samples to plot
myspl <- dir(data_fd)

# Regions to plot
#regions <- list(c(7e5, 1.05e6), c(1.15e6, 1.4e6))
regions <- list(c(2.7e6, 3e6), c(3.15e6, 3.35e6))

myseed <- 1542268



#=================#
# Data processing #
#=================#

cat("Loading data. This may take a while...\n")

# Load data
mydata <- vector("list", length(myspl))
names(mydata) <- myspl
for (i in myspl) {
    mydata[[i]] <- fread(paste0(data_fd,  i, "/", i, "_", mysfx, "_3_QTL.cov"), header = FALSE, sep = "\t")
}

cat("Processing data. This may take a while...\n")

mydata_ls <- vector("list", length(regions))
graph_ls  <- vector("list", length(regions))

for (r in 1:length(regions)) {
    mystart <- regions[[r]][1]
    myend   <- regions[[r]][2]

    # Focusing on the regions of interest
    mydata_tmp <- lapply(mydata, function(x) { x <- x[ V2 > mystart & V2 < myend ] ; return(x) })

    # Smoothing data
    mydata_tmp <- lapply(mydata_tmp, function(x) { x[, 3] <- runmed(x[, V3], 2001) ; return(x) })

    # Reducing data by removing low covered regions
    mydata_tmp <- lapply(mydata_tmp, function(x) { x <- x[ V3 > 5, ] ; return(x) })

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

png(paste0(graph_fd, "CNV_panel.png"), width = 5 * length(myspl) / 2, height = 3 * length(regions) * 2, unit = "in", res = 300)
layout(matrix(1:(length(myspl) * length(regions)), ncol = length(myspl) / 2, byrow = TRUE))

for (r in 1:length(regions)) {
    mydata_tmp <- mydata_ls[[r]]

    x_max <- graph_ls[[r]][1]
    x_min <- graph_ls[[r]][2]
    y_max <- graph_ls[[r]][3]

    for (i in myspl) {
        # Point graph
        plot(mydata_tmp[[i]][,V3] ~ mydata_tmp[[i]][,V2], xlab = "Position on chromosome Z (Mb)", ylab = "Read depth", xlim = c(x_min, x_max), ylim = c(5, y_max), main = i, xaxt = "n", pch = 20, col = "grey", log = "y")
        #axis(1, at = seq(0, 80, 10) * 1e6, labels = seq(0, 80, 10))
        axis(1)

        # Smoothed line
        lw1 <- loess(mydata_tmp[[i]][,V3] ~ mydata_tmp[[i]][,V2], span=0.05)
        lines(mydata_tmp[[i]][,V2], predict(lw1, mydata_tmp[[i]][,V2]), col = "red", lwd = 3)
    }
}
dev.off()
