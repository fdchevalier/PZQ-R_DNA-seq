#-------------------
# Packages
#-------------------

library(doBy)
library(drc)

#----------------
# Datas
#----------------

mydata <- read.csv("Lactate_prod_LE_PZQ_ER_ES_IC50.csv", header=T, sep=",", dec=".")


#-----------------
# Data processing
#-----------------

# function from the drc package - IC50 fitting curve
myfct <- "LL.4"
 
mypop <- unique(mydata$Parasite)

## Create a color vector
myclr <- NULL

## Create color matrix
mycolmat <- matrix(c("SmLE-PZQ-ER", "#aa7b4d",
                     "SmLE-PZQ-ES", "#3bc188"
                     ), ncol = 2, byrow = TRUE)
                     
for (i in mypop) {
    myclr[i] <- mycolmat[mycolmat[,1] == i, 2]
}

## Create an empty list to stock results
res <- vector("list", length(mypop)) 
names(res) <- mypop

for (i in mypop){
    mydata_tmp <- mydata[mydata[,1] == i,]
    
    ##Generate the table with means and standard deviation of the percentage of worm mortality for each PZQ each dose and each population of worms
    sum_table_tmp <- summaryBy(Ratio_perc ~ PZQ_dose, data=mydata_tmp, FUN=c(length,mean,sd))
 
    ## Rename column perc_worm_mortality.length to just N
    names(sum_table_tmp)[names(sum_table_tmp)=="Ratio_perc.length"] <- "N"

    ## Calculate standard err_tmp of the mean
    sum_table_tmp$Ratio_perc.se <- sum_table_tmp$Ratio_perc.sd / sqrt(sum_table_tmp$N)
    
    ## Variation in lactate production relative to the control dose in percentage
    sum_table_tmp$Ratio_perc_100 <- (100 * sum_table_tmp$Ratio_perc.mean) / sum_table_tmp$Ratio_perc.mean[1]
    
    ## Remove the control dose (PZQ=0ug/mL) from the table
    sum_table_tmp <- sum_table_tmp[ !(sum_table_tmp$PZQ_dose %in% c(0.0)), ]

    #Fitting the data with the model for IC50 curves
    myfit_tmp <- drm(Ratio_perc_100 ~ PZQ_dose, data = sum_table_tmp, fct = get(myfct)())
    
    #Fitting the data with the model to ED50 values
    myfit_ED50_tmp <- drm(Ratio_perc_100 ~ PZQ_dose, data = sum_table_tmp, fct = LL.3(fixed = c(NA,100,NA)))
    
    #Store the for loop output (fitting data to model) in the list
    res[[i]] <- myfit_tmp
        
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


#data_ER <- res[[1]][[1]]
#data_ES <- res[[2]][[1]]
#x <- res[[1]][[2]][,1]
#y <- res[[1]][[2]][,6]
#y2 <- res[[2]][[2]][,6]
 
#sd_ER <- res[[1]][[2]][,5]
#sd_ES <- res[[2]][[2]][,5]

#---------
# Figures
#---------

#Plots the IC 50 curves
pdf(file="IC_curves_ER_ES.pdf", width=8, height=8,useDingbats=FALSE)

par(mar=c(5,5,4,2)) #layout margins

#mymax <- lapply(res, function(x) {max(x[[2]][,3])}) %>% unlist() %>% max() %>% round()

for (i in 1:length(res)) {
    
    if (i == 1) {
        myadd <- F
    } else {
        myadd <- T
    }
    
    plot(res[[i]][[1]], type = "all", xlab="PZQ concentration (µg/mL)", ylab="Lactate production (%)", col= myclr[i], pch=16, cex.lab=1.8, axes=FALSE, bty="n", xlim=c(0,100), ylim=c(0,100), add = myadd)
    arrows(res[[i]][[2]][,1], res[[i]][[2]][,6]-res[[i]][[2]][,5], res[[1]][[2]][,1], res[[i]][[2]][,6]+res[[i]][[2]][,5], code=3, length=0.02, angle = 90, col= myclr[i])
    
}

#plot(data_ER, type = "all", xlab="PZQ concentration (µg/mL)", ylab="Lactate production (%)", ylim=c(0,100), col="#aa7b4d", pch=19, cex.lab=1.5, axes=FALSE, bty="n" )
#arrows(x, y-sd_ER, x, y+sd_ER, code=3, length=0.02, angle = 90, col="#aa7b4d")

#plot(data_ES, type ="all", add=T, col="#3bc188", pch=19, axes=FALSE)
#arrows(x, y2-sd_ES, x, y2+sd_ES, code=3, length=0.02, angle = 90, col="#3bc188")
    

##----Axes---------##
axis(1, at=c(0.1,0.3,0.9,2.7,8.1,24.3,72.9), cex.axis=1.1)
axis(2, at=c(0,10,20,30,40,50,60,70,80,90,100), cex.axis=1.1)

#text and line added in the graph
text(45,80,"SmLE-PZQ-ER",cex=1.2, xpd=TRUE, col=mycolmat[1,2])
text(45,25,"SmLE-PZQ-ES",cex=1.2, xpd=TRUE, col=mycolmat[2,2])

#-----------------------------------
# Saving
#-----------------------------------

dev.off()


