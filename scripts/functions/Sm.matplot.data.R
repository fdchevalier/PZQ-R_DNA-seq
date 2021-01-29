# Title: Sm.matplot.data.R
# Version: 3.3
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2012-08-03
# Modified in: 2018-01-03



#==========#
# Comments #
#==========#

# Plot of values along the Schistosoma mansoni genome.



#==========#
# Versions #
#==========#

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
# v1.2 - 2014-10-31: indification of the chromosomes present in the dataset (useful if some missing) / y axis sequence when ylim < 0.01 improved / color setup for line plot improved
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


#------------------------------------------------#
# Fuction to test if the number is a wholenumber #
#------------------------------------------------#

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#----------------#
# Data reodering #
#----------------#

data.order <- function(x) {
    chr.names <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_6", "Chr_7", "Chr_W", "i.SC|^SC")

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


#-----------------------------------------#
# Function to plot the SNP by chromosomes #
#-----------------------------------------#

matplot.data <- function (data.tab, column, datatype, myrunmed=NULL, loess.span=NULL, ylim.min=NULL, ylim.max=NULL, ylab=NULL, xlab="Chromosomes", xlab.axis=c("1", "2", "3", "4", "5", "6", "7", "Z", "Unassembled scaffolds"), col=c("red","black"), chr.bg="grey90", cex.axis=0.9, cex=1, lwd=1, type="l", pch=20, abline.h=NULL, abline.col="blue", abline.lty=3, abline.lwd=1, arrow.head=NULL, arrow.head.col="black", data.order=TRUE, by.pos=FALSE, add=FALSE, verbose=FALSE) {


    #~~~~~~~#
    # Usage #
    #~~~~~~~#

    # data.tab      dataframe containing chr names, position and values to plot
    # column        column of data.tab to plot
    # datatype      can be freq or pvalue
    # myrunmed      value of the runmed to smooth data
    # loess.span    span value of the loess function to smooth data
    # ylim.min
    # ylim.max
    # ylab          name of the y axis
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
    # add           add data to an existing plot
    # verbose


	# Checking variables
    obj.nm <- deparse(substitute(data.tab))
    if (! is.data.frame(data.tab)) {stop("'data.tab' must be a data.frame.")}
    if (length(column) != 1) {stop("'column' must not be greater than 1.")}
    if (! is.vector(data.tab[,column])) {stop(paste(obj.nm,"[,",column,"] must be a vector.",sep=""))}
    if (! is.null (ylim.min) & ! is.numeric(ylim.min)) {stop("ylim.max must be numeric.")}
    if (! is.null (ylim.max) & ! is.numeric(ylim.max)) {stop("ylim.max must be numeric.")}
    if (! is.null (myrunmed) & ! is.null(loess.span)) {stop("myrunmed and loess.span cannot be used together.")}
	if (! is.null (myrunmed) & ! is.numeric(myrunmed)) {stop("myrunmed must be numeric.")}
	if (! is.null (loess.span) & (! is.numeric(loess.span))) {
        if (loess.span < 0 | loess.span > 1) {stop("loess.span must a number in [0;1].")}
    }
	if (type != "l" & type != "p") {stop("type must be l or p.")}
    if (! is.null(xlab.axis) && length(xlab.axis) != 9) {stop("xlab.axis must be a vector of 9 values.")}
    if (! is.null(by.pos) & ! ( is.logical(by.pos) | is.data.frame(by.pos) | is.matrix(by.pos) )) {stop("by.pos must be TRUE/FALSE or a dataframe/matrix of chromosome names and positions.")}
    if (! is.null(arrow.head) & ! ( isTRUE(by.pos) | is.data.frame(by.pos) | is.matrix(by.pos) )) {stop("arrow.head must be a dataframe or a matrix of chromosome names and positions and requires by.pos.")}

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

    # Apply limits to the data if present
    if (! is.null(ylim.min)) { rmdt[ rmdt < ylim.min ] <- ylim.min }
    if (! is.null(ylim.max)) { rmdt[ rmdt > ylim.max ] <- ylim.max }

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
	chr.names.cplt      <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_6", "Chr_7", "Chr_W", "i.SC|^SC")
    chr.names.true.cplt <- xlab.axis

    chr.names      <- NULL
    chr.names.true <- NULL
    chr.names.all  <- unique(as.vector(data.coord[,1]))
    
    for (i in 1:length(chr.names.cplt)) {
        if (any(grepl(chr.names.cplt[i], chr.names.all))) {
            chr.names <- c(chr.names, chr.names.cplt[i])
        	chr.names.true <- c(chr.names.true, chr.names.true.cplt[i])
        }
    }

    if (length(chr.names) == 0) {stop("data.tab does not seem well formated. No chromosome name detectable in the first column of the table.")}

    if (by.pos) {
        # Get the last pb position of each chr
        chr.length <- c(0, sapply(chr.names.all, function (x) tail(data.coord[ data.coord[,1] == x, 2],1)))
        # Transform the vector to have cumulative positions
        chr.length <- c(0, sapply(2:length(chr.length), function(i) sum(chr.length[1:(i-1)])+chr.length[i]))
        # Get the positions of the chr + scaffold associated
        chr.SNP.length <- c(0, sapply(chr.names, function(x) tail(chr.length[grep(x, names(chr.length))],1)))
    } else {
        # Get the cumulative number of SNP per chr
        chr.SNP.length <- c(0, sapply(chr.names, function(x) tail(grep(x,data.coord[,1]),1)) )
    }

    # Positions of the arrow head(s)
    if (! is.null(arrow.head)) {
        arrow.pos <- NULL

        for (i in 1:nrow(arrow.head)) {
            myidx <- grep(arrow.head[i,1], chr.names.all) - 1
            if (myidx == 0) {
                arrow.pos <- c(arrow.pos, as.numeric(as.character(arrow.head[i,2])))
            } else {
                arrow.pos <- c(arrow.pos, as.numeric(as.character(arrow.head[i,2])) + chr.length[myidx])
            }
        }

    }


    #~~~~~~~~~~~~~~~~~~#
	# Plot preparation #
    #~~~~~~~~~~~~~~~~~~#

    # X axis parameter
    xlim.max <- tail(chr.SNP.length, 1)

    # Creating initial plot
	if (add == "FALSE") {
		# Y axis parameter
        if (max(ceiling(rmdt),na.rm=TRUE) <= 1 || min(floor(rmdt),na.rm=TRUE) >= -1) {
#            ylim.max <- ceiling(max(rmdt)*10^ceiling(abs(log10(abs(max(rmdt))))))/10^ceiling(abs(log10(abs(max(rmdt)))))
#            ylim.min <- floor(min(rmdt)*10^ceiling(abs(log10(abs(min(rmdt))))))/10^ceiling(abs(log10(abs(min(rmdt)))))

            pwr <- ceiling(abs(min(c(log10(abs(max(rmdt,na.rm=TRUE)))), log10(abs(min(rmdt,na.rm=TRUE))),na.rm=TRUE)))
            ylim.max.tmp <- ceiling(max(rmdt,na.rm=TRUE)*10^pwr)/10^pwr
            ylim.min.tmp <-   floor(min(rmdt,na.rm=TRUE)*10^pwr)/10^pwr
           
            # When one of the border is 0
            if (is.na(ylim.max.tmp)) { ylim.max.tmp <- max(ceiling(rmdt),na.rm=TRUE) }
            if (is.na(ylim.min.tmp)) { ylim.min.tmp <- min(floor(rmdt),na.rm=TRUE)   }

            # TODO: future option: mirror to have the same max and min => ylim.min/max <- max(abs(c(ylim.min, y.lim.max)))
        }
        
        ## Y min
        if (is.null(ylim.max) && ! exists("ylim.max.tmp")) {
            ylim.max <- max(ceiling(rmdt),na.rm=TRUE)
        } else if (is.null(ylim.max) && exists("ylim.max.tmp")) {
            ylim.max <- ylim.max.tmp
        }
        
        ## Y max
        if (is.null(ylim.min) && ! exists("ylim.min.tmp")) {
            ylim.min <- min(floor(rmdt),na.rm=TRUE)
        } else if (is.null(ylim.min) && exists("ylim.min.tmp")) {
            ylim.min <- ylim.min.tmp
        }


        # Empty plot
		plot(1,1, ylim=c(ylim.min, ylim.max), xlim=c(0, xlim.max), xlab=xlab, type="n", bty="n", axes=FALSE, ann=FALSE)

		# Axes draw
		mypos <- NULL
		for (i in 1:length(chr.names)) {
			mypos[i] <- chr.SNP.length[i]+(chr.SNP.length[i+1]-chr.SNP.length[i])/2
		}
		axis(1, at=chr.SNP.length, labels=FALSE)
		axis(1, at=mypos, labels=chr.names.true, cex.axis=cex.axis, tick=FALSE)
		
        y.max <- ceiling(ylim.max/10^round(log10(ylim.max)-1))*10^round(log10(ylim.max)-1)
        #my.ylabels <- seq(ylim.min, y.max, length=5)
        my.ylabels <- seq(ylim.min, ylim.max, length=5)
        if (log10(max(abs(my.ylabels))) < -2 | log10(max(abs(my.ylabels))) > 4) {sc.write <- TRUE} else {sc.write <- FALSE}
		axis(2, at=format(my.ylabels,scientific=sc.write,digits=2), cex.axis=cex.axis)
		
        mtext(xlab, side=1, line=3)
		mtext(ylab, side=2, line=3)

		# Plot the grey rectangles to delimitate the chromosomes
		for (i in 1:(length(chr.names))) {
			if (is.wholenumber(i/2) == FALSE){
#				rect(head(grep(chr.names[i],data.tab[,1]),1), ylim.min, tail(grep(chr.names[i],data.tab[,1]),1), ylim.max, col=chr.bg, border=NA)
				rect(chr.SNP.length[i], ylim.min, chr.SNP.length[i+1], ylim.max, col=chr.bg, border=NA)
			}
		}
	}

	# Point plot
	if (type == "p") {
        if (by.pos) {
            cum.pos <- unlist(sapply(chr.names.all, function (x) data.tab[ data.tab[,1] == x, 2] + chr.length[match(x,chr.names.all)]))
        } else {
            cum.pos <- 1:xlim.max
        }

        mycolor <- rep(col[1], nrow(data.tab))
        mycolor[grep("SC",data.tab[,1])] <- col[2]
		points(cum.pos, rmdt, pch=pch, col=mycolor)
	}


	# Line plot
    if (type == "l") {
		for (i in chr.names) {			# Plotting the data
            
            j <- grep(i, chr.names.all, value=TRUE)
            # If SC found, concatenate all of them
            if (any(grepl("SC",j))) {
                j.c <- c(grep("SC",j,invert=TRUE, value=TRUE), paste(grep("SC",j, value=TRUE),collapse="$|^"))
            } else {
                j.c <- j
            }

            for (k in j.c) {
                myexp <- paste0("^",k,"$")
            
                # Identify coordinates from the name of the chr or SC and create a sequence from this
                mystart <- head(grep(myexp, data.tab[,1]), 1)
                myend   <- tail(grep(myexp, data.tab[,1]), 1)

                # Get data
                rmdt.tmp <- rmdt[ mystart:myend ]

                # Define x positions
                if (by.pos) {
                    mysc <- grep(myexp, j, value=TRUE)
                    plot.pos <- unlist(lapply(mysc, function (x) data.tab[ data.tab[,1] == x, 2] + chr.length[match(x,chr.names.all)]))
                } else {
                    plot.pos <- seq(mystart, myend)
                }
                
                # Set color depending on scaffold type
                if (any(grepl("SC", k))) { mycolor <- col[2] } else { mycolor <- col[1] }

                # Draw
                matplot(plot.pos, rmdt.tmp, col=mycolor, type="l", add=TRUE, lwd=lwd)
            }
            
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
        arr.gene(arrow.pos, ylim.max, col=arrow.head.col)
    }

} # End of the function

