#!/usr/bin/env Rscript
# Title: PZQ-R_ER-ES_CNV_graph_RD.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-06-05
# Modified in: 2021-08-10



#==========#
# Comments #
#==========#

# Plot read depth of the QTL of chr. 3



#==========#
# Versions #
#==========#

# v0.1 - 2021-08-10: add packages / load data more efficiently / redesign graph layout
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
data_fd   <- "../data/libraries/2-ER-ES_populations/"
graph_fd  <- "../graphs/"

# Samples to plot
myspl <- dir(data_fd)

# Regions to plot
regions <- list(c(7e5, 1e6), c(1.15e6, 1.5e6))

myseed <- 1542268



#=================#
# Data processing #
#=================#

cat("Loading data. This may take a while...\n")

# Load data
mydata <- vector("list", length(myspl))
names(mydata) <- myspl
for (i in myspl) {
    mydata[[i]] <- fread(paste0(data_fd,  i, "/", i, "_3_QTL.cov"), header = FALSE, sep = "\t")
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

cat("\nDrawing graphs...\n")

pdf(paste0(graph_fd, "Supp. Fig. 4.pdf"), width = 15, height = 8)

mycex    <- 1.2
mycex_ax <- 1.2

# Layout
plot_nb <- length(myspl) * length(regions)
mymat <- matrix(1:plot_nb, ncol = length(myspl) / 2, byrow = TRUE)

## Add rows
mymat <- c(rep(max(mymat) + 1, length(myspl) / 4), rep(max(mymat) + 2, length(myspl) / 4)) %>% rbind(., mymat, deparse.level = 0)
mymat <- rep(max(mymat) + 1, length(myspl) / 2) %>% rbind(mymat, ., deparse.level = 0)

## Add columns
mymat <- c(0, rep(max(mymat) + 1, length(regions)), rep(max(mymat) + 2, length(regions)), 0) %>% cbind(., mymat, deparse.level = 0)
mymat <- c(0, max(mymat) + 1:4, 0) %>% cbind(mymat, ., deparse.level = 0)

layout(mymat, widths = c(0.1, rep(1, length(myspl) / 2), 0.1), heights = c(0.2, rep(1, length(regions) * 2), 0.1))

# Adjust margins
par(mar = c(3, 3, 1, 1))

# Plot regions
for (r in 1:length(regions)) {
    mydata_tmp <- mydata_ls[[r]]

    x_max <- graph_ls[[r]][1]
    x_min <- graph_ls[[r]][2]
    y_max <- graph_ls[[r]][3]

    for (i in 1:length(myspl)) {
        # Title
        mytitle <- NULL
        if (i == 1 | i == 4) { mytitle <- "Replicate 1" }
        if (i == 2 | i == 5) { mytitle <- "Replicate 2" }
        if (i == 3 | i == 6) { mytitle <- "Replicate 3" }

        # Point graph
        plot(mydata_tmp[[i]][, V3] ~ mydata_tmp[[i]][, V2], xlab = "", ylab = "", xlim = c(x_min, x_max), ylim = c(5, y_max), axes = FALSE, pch = 20, col = "grey", log = "y", frame.plot = TRUE)

        mtext(mytitle, cex = 0.7)

        # Axes
        if (i %in% c(7:12, 19:24))  { 
            mymarks <- seq(par("xaxp")[1], par("xaxp")[2], length.out=par("xaxp")[3])
            axis(1, at = mymarks, labels = format(mymarks, scientific = TRUE, digits = 2))
        }
        if (i %in% c(1, 7, 13, 19)) { axis(2) }        

        # Smoothed line
        lw1 <- loess(mydata_tmp[[i]][, V3] ~ mydata_tmp[[i]][, V2], span=0.05)
        lines(mydata_tmp[[i]][, V2], predict(lw1, mydata_tmp[[i]][, V2]), col = "red", lwd = 3)
    }
}

# Adjust margins
par(mar = c(0, 3, 0, 1))

# Title
mysex <- c("Male", "Female")
for (i in mysex) {
    plot.new()
    lines(par("usr")[1:2], c(0.12, 0.12))
    text(0.5, 0.5, labels = i, cex = mycex * 1.5, font = 2)
}

# Adjust margins
par(mar = rep(0, 4))

# for (i in 1:length(mysex)) {
#     plot.new()
#     text(0.5, 0.5, labels = "Position on chromosome 3 (bp)", cex = mycex)
# }

# X axis label
plot.new()
text(0.5, 0.5, labels = "Position on chromosome 3 (bp)", cex = mycex)

# Y axis labels
for (i in 1:length(regions)) {
    plot.new()
    text(0.5, 1, labels = LETTERS[i], cex = mycex * 1.5)
    text(0.5, 0.5, labels = "Read depth", cex = mycex, srt = 90)
}

# Population tags
for (i in 1:length(regions)) {
    plot.new()
    text(0.5, 0.5, labels = "ER", cex = mycex, font = 2)

    plot.new()
    text(0.5, 0.5, labels = "ES", cex = mycex, font = 2)
}

# Separator
par(xpd = NA)
abline(h = 2.15, lwd = 2)

dev.off()
