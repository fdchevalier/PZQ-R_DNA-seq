#-------------------
# Packages
#-------------------

library(doBy)
library(drc)

#----------------
# Datas
#----------------

mydata <- read.table("IC50_LE_PZQR.tsv", header=T, sep="\t", dec=".")


#-----------------
# Data processing
#-----------------

# function from the drc package - IC50 fitting curve
myfct <- "LL.3"
 
mypop <- unique(mydata$type)

## Create a color vector
myclr <- NULL

## Create color matrix
mycolmat <- matrix(c("SmLE", "black",
                     "SmLE-PZQ-R", "red"
                     ), ncol = 2, byrow = TRUE)
                     
for (i in mypop) {
    myclr[i] <- mycolmat[mycolmat[,1] == i, 2]
}

## Create an empty list to stock results
res <- vector("list", length(mypop)) 
names(res) <- mypop

for (i in mypop){
    mydata_tmp <- mydata[mydata[,5] == i,]
    
    ##Generate the table with means and standard deviation of the percentage of worm mortality for each PZQ each dose and each population of worms
    sum_table_tmp <- summaryBy(perc_worm_mortality ~ PZQ_dose, data=mydata_tmp, FUN=c(length,mean,sd))
 
    ## Rename column perc_worm_mortality.length to just N
    names(sum_table_tmp)[names(sum_table_tmp)=="perc_worm_mortality.length"] <- "N"

    ## Calculate standard err_tmp of the mean
    sum_table_tmp$perc_worm_mortality.se <- sum_table_tmp$perc_worm_mortality.sd / sqrt(sum_table_tmp$N)

    #Fitting the data with the model for IC50 curves
    myfit_tmp <- drm(perc_worm_mortality.mean ~ PZQ_dose, data = sum_table_tmp, fct = get(myfct)())
    
    #Fitting the data with the model to ED50 values
    myfit_ED50_tmp <- drm(perc_worm_mortality.mean ~ PZQ_dose, data = sum_table_tmp, fct = LL.3(fixed = c(NA, 100, NA)))
    
    # Create the list with 3 slots
    res_tmp <- vector("list", 3)
    
    #Store the for loop output (fitting data to model) in the list
    res_tmp[[1]] <- myfit_tmp
    
    #Extraction of the SE of the means from the summary table _tmp
    res_tmp[[2]] <- sum_table_tmp
    
    #Store the for loop output (fitting data to model to compute ED50) in the list
    res_tmp[[3]] <- myfit_ED50_tmp
        
    res[[i]] <- res_tmp
    
    
}

# Compute and store IC50 values for all the populations tested
## Create an empty list to stock results
res_IC <- vector("list", length(mypop)) 
names(res_IC) <- mypop

for (i in 1:length(res_IC)) {
        
        #Compute IC50 for each population
        IC_tmp <- ED(res[[i]][[3]], c(50), interval ="delta")
      
        # Create the list with 2 slots
        res_tmp <- vector("list", 1)
        
        #Store the IC50 values output in the list
        res_tmp[[1]] <- IC_tmp
        res_IC[[i]] <- res_tmp

}

#---------
# Figures
#---------

#Plots the IC 50 curves
pdf(file="IC_curves_LE_PZQR.pdf", width=8, height=8,useDingbats=FALSE)

par(mar=c(5,5,4,2)) #layout margins

#myclr <- rainbow(length(res))
mymax <- lapply(res, function(x) {max(x[[2]][,3])}) %>% unlist() %>% max() %>% round()

for (i in 1:length(res)) {
    
    if (i == 1) {
        myadd <- F
    } else {
        myadd <- T
    }
    
    plot(res[[i]][[1]], type = "all", xlab="PZQ concentration (µg/mL)", ylab="% worm mortality", col= myclr[i], pch=16, cex.lab=1.8, axes=FALSE, bty="n", xlim=c(0,100), ylim=c(0, mymax), add = myadd, xpd=T)
    arrows(res[[i]][[2]][,1], res[[i]][[2]][,3]-res[[i]][[2]][,5], res[[1]][[2]][,1], res[[i]][[2]][,3]+res[[i]][[2]][,5], code=3, length=0.02, angle = 90, col= myclr[i], xpd=T)
    
}

#axes
axis(1, at=res[[1]][[2]][,1], cex.axis=1.1)
#axis(1, at=c(0.1,0.3,0.9,2.7,8.1,24.3,72.9), cex.axis=1.1)
axis(2, at=c(0,10,20,30,40,50,60,70,80,90,100), cex.axis=1.1, xpd=T)

#text and line added in the graph
text(45,70,"LE-PZQR",cex=1.8, xpd=TRUE, col="red")
text(45,95,"LE",cex=1.8, xpd=TRUE, col="black")

abline(h=65, lty=2, lwd=0.8, col ="grey46")

#-----------------------------------
# Saving
#-----------------------------------

dev.off()


