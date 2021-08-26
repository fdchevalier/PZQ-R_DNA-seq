#!/usr/bin/env Rscript
# Title: line2user.R
# Version: 0.0
# Author: jbaums, Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2018-07-18
# Modified in:



#==========#
# Comments #
#==========#

# Convert line coordinates to physical coordinates for graphic device.
# This function was published on Stackoverflow by jbaums (source: https://stackoverflow.com/a/30835971)



#==========#
# Versions #
#==========#

# v0.0 - 2018-07-18: creation



#==========#
# Function #
#==========#

# Line in units
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
