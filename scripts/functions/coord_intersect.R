# Title: coord_intersect.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-05-24
# Modified in:



#==========#
# Comments #
#==========#

# Intersect list of genomic positions with position intervals similarly to the bed intersect command.



#==========#
# Versions #
#==========#

# v0.0 - 2021-05-24: creation



#==========#
# Function #
#==========#

coord_intersect <- function(x, bed) {

    mylist <- vector("list", nrow(bed))

    for (i in 1:nrow(bed)) {
        mylist[[i]] <- x[ x[, 1] == bed[i, 1] & x[, 2] > bed[i, 2] & x[, 2] <= bed[i, 3], ]
    }

    return(mylist)
}
