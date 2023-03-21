# Title: Sm.matplot.data.R
# Version: 4.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2012-08-03
# Modified in: 2023-03-20



#==========#
# Comments #
#==========#

# Plot of values along the Schistosoma mansoni genome.



#==========#
# Versions #
#==========#

# v4.0 - 2023-03-20: adapt function to the v9 genome / add and update y-axis related arguments / correct typos
# v3.4 - 2021-08-11: correct minor bug regarding arrow.head
# v3.3 - 2018-01-03: transform data to respect ylim.max if use
# v3.2 - 2017-08-28: add arrow.head option
# v3.1 - 2017-06-23: by.pos improved to plot empty chromosomes
# v3.0 - 2017-03-07: info about usage added / loess option added / by.pos option added / bug with xlim.max when using abline corrected / code cleaned up
# v2.5 - 2016-07-17: script updated for Sm version 6 genome
# v2.4 - 2016-03-05: data reorder step made optional
# v2.3 - 2015-08-25: check if column of data is a vector / y axis improved for values within [-1;1] / bug in cex of x axis corrected / default value of cex.axis revised / color bug correction in the unass. scaffold plotting (1st scaffold red)
# v2.2 - 2015-05-29: xlab.axis added / abline.lwd added
# v2.1 - 2015-05-04: ylim.min added / y axis improved for values within [-1;1] / error message improved / check steps for input added / abline added
# v2.0 - 2014-11-03: data order function added / y axis designed improved (scientific format when big numbers)
# v1.2 - 2014-10-31: identification of the chromosomes present in the dataset (useful if some missing) / y axis sequence when ylim < 0.01 improved / color setup for line plot improved
# v1.1 - 2014-10-25: runmed calculation done by chromosomes/scaffolds / color setup for point plot dramatically improved
# v1.0 - 2014-04-07: y axis improved (choice in max and label) / x axis improve (label) / type added (line or points) / add option added / choice in color of assembled and unassembled scaffold added
# v0.4 - 2014-02-19: improvement of the y axis design
# v0.3 - 2013-09-24: windows option removed / initial empty plot and y axis improved / mitochondria removed if present
# v0.2 - 2012-10-22: windows option added for runmed
# v0.1 - 2012-10-14: cex options added / initial empty plot and y axis improved
# v0.0 - 2012-08-03: creation


#==============#
# Dependencies #
#==============#

if (! suppressWarnings(require("plotrix", quietly=T))) {stop("plotrix package is missing.", call.=FALSE)}    # for ablineclip function
if (! suppressWarnings(require("shape",   quietly=T))) {stop("shape package is missing.",   call.=FALSE)}    # for Arrowhead function


#===========#
# Functions #
#===========#


#-------------------------------------------------#
# Function to test if the number is a wholenumber #
#-------------------------------------------------#

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#-----------------#
# Data reordering #
#-----------------#

data.order <- function(x) {
    chr.names <- c("SM_V9_1", "SM_V9_2", "SM_V9_3", "SM_V9_4", "SM_V9_5", "SM_V9_6", "SM_V9_7", "SM_V9_PAR1", "SM_V9_ZSR", "SM_V9_WSR", "SM_V9_PAR2")

    data.tmp <- NULL

    for (i in chr.names) {
        data.tmp <- rbind(data.tmp, x[grep(i,x[,1]),])
    }

    return(data.tmp)
}


#------------#
# Arrow head #
#------------#

arr.gene <- function (x, y, col="black") {
   Arrowhead(x0 = x, y0 = y, angle = -90, arr.length = 0.2, arr.type = "triangle", arr.col = col, lty=0, arr.adj = 1)
}


#----------------------------------------------#
# Function to plot the variants by chromosomes #
#----------------------------------------------#

matplot.data <- function (data.tab, column, datatype, myrunmed=NULL, loess.span=NULL, ylim=NULL, ylab=NULL, ymir=FALSE, xlab="Chromosomes", xlab.axis=c("1", "2", "3", "4", "5", "6", "7", "PAR1", "ZSR", "WSR", "PAR2"), col=c("aquamarine4", "grey30"), chr.bg="grey90", cex.axis=0.9, cex=1, lwd=1, type="l", pch=20, abline.h=NULL, abline.col="blue", abline.lty=3, abline.lwd=1, arrow.head=NULL, arrow.head.col="black", data.order=TRUE, by.pos=FALSE, gap=0.2, add=FALSE, verbose=FALSE) {


    #~~~~~~~#
    # Usage #
    #~~~~~~~#

    # data.tab      dataframe containing chr names, position and values to plot
    # column        column of data.tab to plot
    # datatype      can be freq or pvalue
    # myrunmed      value of the runmed to smooth data
    # loess.span    span value of the loess function to smooth data
    # ylim
    # ylab          name of the y axis
    # ymir          mirror y axis (useful for differences in allele frequency)
    # xlab          name of x axis
    # xlab.axis     name of the chromosomes
    # col           color of the line (first: chr, second: scaffold)
    # chr.bg        color of the background rectangle for chr
    # cex.axis
    # cex
    # lwd           line width if type "l"
    # type          "l" for line, "p" for point
    # pch           pch of the point
    # abline.h      y value of the abline
    # abline.col    color of the abline
    # abline.lty    lty for abline
    # abline.lwd    lwd for abline
    # arrow.head        matrix containing coordinate to place arrow head(s)
    # arrow.head.col    color of the arrow head(s)
    # data.order    order data prior plotting
    # by.pos        plot by position
    # gap           add white space between chromosomes
    # add           add data to an existing plot
    # verbose


	# Checking variables
    obj.nm <- deparse(substitute(data.tab))
    if (! is.data.frame(data.tab)) {stop("'data.tab' must be a data.frame.")}
    if (length(column) != 1) {stop("'column' must not be greater than 1.")}
    if (! is.vector(data.tab[,column])) {stop(paste(obj.nm,"[,",column,"] must be a vector.",sep=""))}
    if (! is.null (ylim) & ! (is.numeric(ylim) && length(ylim) == 2)) {stop("ylim must be numeric vector of 2 values.")}
    if (! is.null (myrunmed) & ! is.null(loess.span)) {stop("myrunmed and loess.span cannot be used together.")}
	if (! is.null (myrunmed) & ! is.numeric(myrunmed)) {stop("myrunmed must be numeric.")}
	if (! is.null (loess.span) & (! is.numeric(loess.span))) {
        if (loess.span < 0 | loess.span > 1) {stop("loess.span must a number in [0;1].")}
    }
	if (type != "l" & type != "p") {stop("type must be l or p.")}
    if (! is.null(xlab.axis) && length(xlab.axis) != 11) {stop("xlab.axis must be a vector of 11 values.")}
    if (! is.null(by.pos) & ! ( is.logical(by.pos) | is.data.frame(by.pos) | is.matrix(by.pos) )) {stop("by.pos must be TRUE/FALSE or a dataframe/matrix of chromosome names and positions.")}
    if (! is.null(arrow.head) & ! ( isTRUE(by.pos) | is.data.frame(by.pos) | is.matrix(by.pos) )) {stop("arrow.head must be a dataframe or a matrix of chromosome names and positions and requires by.pos.")}
	if (! is.numeric(gap)) {stop("gap must be numeric.")}
	if (! is.logical(ymir)) {stop("ymir must be logical.")}

    # Data ordering
    if (data.order) {data.tab <- data.order(data.tab)}

	# Smoothing data using runmed
	if (is.null(myrunmed) & is.null(loess.span)) {
		rmdt <- data.tab[,column]
	} else if (! is.null(myrunmed)) {
		rmdt <- unlist(tapply(data.tab[,column], data.tab[,1], function(x) runmed(x,myrunmed)), use.names=FALSE) # runmed(data.tab[,column], myrunmed)
	} else if (! is.null(loess.span)) {
        ids <- 1:nrow(data.tab)
		rmdt <- predict(loess(data.tab[,column] ~ ids, span=loess.span))
	}

    # # Apply limits to the data if present
    # if (! is.null(ylim.min)) { rmdt[ rmdt < ylim.min ] <- ylim.min }
    # if (! is.null(ylim.max)) { rmdt[ rmdt > ylim.max ] <- ylim.max }

	# Type of data
	if (datatype == "freq" & is.null(ylab)) {
		ylab <- "Allele frequency"
	} else if (datatype == "pvalue" & is.null(ylab)){
		rmdt <- -log10(rmdt)
		ylab <- expression(paste(-Log[10]," p-value", sep=""))
	} else if (datatype != "freq" & datatype != "pvalue") {
		stop("datatype must be \"freq\" or \"pvalue\"")
	}

    # Define positions: where to look at coordinate data (within data.tab or within given table)
    if (is.matrix(by.pos) | is.data.frame(by.pos)) {
        data.coord <- by.pos
        # Reset by.pos
        by.pos <- TRUE
    } else {
        data.coord <- data.tab[,1:2]
    }

	# Chromosome parameters use for search and plotting
	chr.names      <- c("SM_V9_1", "SM_V9_2", "SM_V9_3", "SM_V9_4", "SM_V9_5", "SM_V9_6", "SM_V9_7", "SM_V9_PAR1", "SM_V9_ZSR", "SM_V9_WSR", "SM_V9_PAR2")
    chr.names.true <- xlab.axis
    chr.names.all  <- unique(as.vector(data.coord[,1]))

    if (all(! chr.names.all %in% chr.names)) {stop("data.tab does not seem well formated. No chromosome name detectable in the first column of the table.")}

    if (by.pos) {
        # Get the last pb position of each chr
        chr.length <- c(0, sapply(chr.names, function (x) tail(data.coord[ data.coord[,1] == x, 2],1)))
        # Transform the vector to have cumulative positions
        chr.length <- c(0, sapply(2:length(chr.length), function(i) sum(chr.length[1:(i-1)])+chr.length[i]))
        # Get the positions of the chr + scaffold associated
        chr.SNP.length <- c(0, sapply(chr.names, function(x) tail(chr.length[grep(x, names(chr.length))],1)))
    } else {
        # Get the cumulative number of SNP per chr
        chr.SNP.length <- c(0, sapply(chr.names, function(x) tail(grep(x,data.coord[,1]),1)) )
    }

    # Gap factor
    gap <- tail(chr.SNP.length, 1) * gap / 10

    # Positions of the arrow head(s)
    if (! is.null(arrow.head)) {
        arrow.pos <- NULL

        for (i in 1:nrow(arrow.head)) {
            myidx <- grep(paste0("^", arrow.head[i, 1], "$"), chr.names.all)
            if (myidx == 0) {
                arrow.pos <- c(arrow.pos, as.numeric(as.character(arrow.head[i,2])))
            } else {
                arrow.pos <- c(arrow.pos, as.numeric(as.character(arrow.head[i,2])) + chr.length[myidx] + gap * myidx)
            }
        }

    }


    #~~~~~~~~~~~~~~~~~~#
	# Plot preparation #
    #~~~~~~~~~~~~~~~~~~#

    # X axis parameter
    xlim.max <- tail(chr.SNP.length, 1) + (gap * (length(chr.names) - 1))

    # Creating initial plot
	if (add == "FALSE") {
		# Y axis parameter
        if (is.null(ylim)) {
            if (ymir) {
                mymax   <- max(abs(rmdt))
                myrange <- c(-mymax, mymax)
            } else {
                myrange <- range(rmdt)
            }

            ylim <- range(pretty(myrange))
        }

        # Empty plot
		plot(1,1, ylim = ylim, xlim = c(0, xlim.max), xlab = xlab, type = "n", bty = "n", axes = FALSE, ann = FALSE)

		# Axes draw
        mypos_lb <- NULL
        mypos_tk <- NULL
		for (i in 1:length(chr.names)) {
			mypos_lb[i]   <- chr.SNP.length[i] + (chr.SNP.length[i+1]-chr.SNP.length[i]) / 2 + gap * (i - 1)
            mypos_tk[[i]] <- c(chr.SNP.length[i], chr.SNP.length[i+1]) + gap * (i - 1)
		}

        ## Axis lines
        if (gap == 0) {
            axis(1, at = chr.SNP.length, labels=FALSE)
        } else {
            axis(1, at = unlist(mypos_tk), labels = FALSE, col = NA, col.ticks = 1)
            for (i in 1:length(mypos_tk)) {
                lines(c(mypos_tk[[i]][1], mypos_tk[[i]][2]), rep(par("usr")[3], 2))
            }
        }

        ## Axis labels
		axis(1, at=mypos_lb, labels=chr.names.true, cex.axis=cex.axis, tick=FALSE)
		axis(2, cex.axis=cex.axis)

        mtext(xlab, side=1, line=3)
		mtext(ylab, side=2, line=3)

		# Plot the grey rectangles to delimitate the chromosomes
        ymin_r <- ylim[1]
        if (datatype == "pvalue") { ymin_r <- min(rmdt) }
		for (i in 1:(length(chr.names))) {
			if (is.wholenumber(i/2) == FALSE){
				rect(mypos_tk[[i]][1], ymin_r, mypos_tk[[i]][2], ylim[2], col = chr.bg, border = NA)
			}
		}
	}

    # Define x positions
    if (by.pos) {
        cum.pos <- sapply(chr.names.all, function(x) data.tab[ data.tab[,1] == x, 2] + chr.length[match(x,chr.names.all)] + gap * (match(x,chr.names.all) - 1))
    } else {
        # cum.pos <- 1:xlim.max
        cum.pos <- sapply(chr.names.all, function(x) which(data.tab[,1] == x) + gap * (match(x,chr.names.all) - 1))
    }

    # Colors
    if (length(col) != length(xlab.axis)) {
        col <- rep(col, length.out = length(xlab.axis))
    }

	# Point plot
	if (type == "p") {
        col <- unlist(sapply(chr.names.all, function(x) rep(col[match(x,chr.names.all)], sum(data.tab[,1] == x))))
        points(unlist(cum.pos), rmdt, pch = pch, col = col)
	}


	# Line plot
    if (type == "l") {
		for (i in 1:length(chr.names)) {			# Plotting the data

            # Get data
            rmdt.tmp <- rmdt[ which(data.tab[,1] == chr.names[i]) ]

            # Draw
            matplot(cum.pos[[i]], rmdt.tmp, col=col[i], type="l", add=TRUE, lwd=lwd)

			if (verbose) {cat(paste(i, "-", chr.names[i]), "\n")}
		}
	}

    # Abline
    if (! is.null(abline.h)) {
        ablineclip(h=abline.h, x1=0, x2=xlim.max, col=abline.col, lty=abline.lty, lwd=abline.lwd)
    }

    if (! verbose) {cat(paste("\tGraph related to", colnames(data.tab)[column], "created.\n"))}

    # Arrow head
    if (! is.null(arrow.head)) {
        arr.gene(arrow.pos, ylim[2], col=arrow.head.col)
    }

} # End of the function

