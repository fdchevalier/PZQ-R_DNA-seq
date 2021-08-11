#!/usr/bin/env Rscript
# Title: rename_chr.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-06-30
# Modified in: 2021-08-11



#==========#
# Comments #
#==========#

# Rename chromosome ID of the V7 S. mansoni genome to the form of the V5 for script backward compatibility.



#==========#
# Versions #
#==========#

# v0.1 - 2021-08-11: correct bug related to some scaffold name
# v0.0 - 2020-06-30: creation



#==========#
# Function #
#==========#

rename_chr_SmV7 <- function(x, cln) {
    y <- as.character(x[,cln])

    y <- gsub("SM_V7_", "", y)
    y <- gsub("(?!^U|[MITO])[A-V,X,Y]+", "\\.un.SC_", y, perl = TRUE)
    y <- gsub("^U", "SC_", y)
    y <- gsub("^(W|\\d\\d+)", "SC_\\1", y, perl = TRUE)
    y <- gsub("^ZW", "W", y)
    y <- gsub("(?!^SC)^", "Chr_", y, perl = TRUE)

    # x <- cbind(original=x, renamed=y)
    x[,cln] <- y

    return(x)
}

