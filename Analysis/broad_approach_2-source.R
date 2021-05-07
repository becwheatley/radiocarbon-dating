#--------------------------------------------------------------------------------------------------------------------------------
# RADIOCARBON DATING PROJECT
# AUSTARCH - BROAD APPROACH USING ONLY SITES THAT HAVE AT LEAST 5 SAMPLES AND THEN ADDING IN OTHER DATA
# Source functions
# Code by Rebecca Wheatley
# Last modified 7 May 2021
#--------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------
# I. GENERATE FREQUENCY DISTRIBUTION
#--------------------------------------------------------------------------------------------------------------------------------

## Create a standardised frequency distribution of dates within bins
## @param data               = our list of uncalibrated data
## @param calibrated_data    = our list of calibrated data
## @param timeRange          = the time range we are interested in sampling over
## @param taphCorrect        = logical specifying whether to taphonomically correct open sites using Williams (2013) (TRUE or FALSE) 
## @param correctForSite     = logical specifying whether to count all radiocarbon dates within a bin or only one per site (TRUE OR FALSE) - note, presently taphCorrect does nothing if this is TRUE
## @param binSize            = size of bin to sort dates into for frequency distribution (200 years used in Williams 2013)
generate_frequency_dist <- function(data, calibrated_data, timeRange, taphCorrect, correctForSite, binSize){
  
 # Get MedianBP values for each calibrated date
 summary <- summary(calibrated_data)
    
 # Get age bins
 data$bin <- cut(summary$MedianBP, 
                 breaks = c(timeRange[2], seq(timeRange[2]+binSize, timeRange[1], by = binSize)), 
                 labels = 1:((timeRange[1] - timeRange[2])/binSize)
                 )
    
 # If we only want to count one occurance of a particular site per bin
 if (correctForSite){
     one_site_per_bin <- unique(data[c("SITE", "bin")])
     count.bins <- as.data.frame(table(one_site_per_bin$bin))$Freq
      
 # Else if we want to count all radiocarbon dates
 } else {
      
   ## Apply taphonomic correction to open sites
   if (taphCorrect){
       open.sites                     <- data %>% subset(Open.or.Closed.Site == "Open")
       count.bins.open                <- as.data.frame(table(open.sites$bin))
       count.bins.open$corrected.freq <- count.bins.open$Freq/(2.107 * 10^7*(median.bin.age + 2754)^(-1.526))
       closed.sites                   <- data %>% subset(Open.or.Closed.Site == "Closed")
       count.bins.closed              <- as.data.frame(as.data.frame(table(closed.sites$bin)))
      
       ## Combine taphonomically corrected open site counts to non-corrected closed site counts to get the corrected frequency distribution in bins
       count.bins                     <- count.bins.open$corrected.freq + count.bins.closed$Freq
   } else {
      
   ## Get uncorrected frequency distribution in bins
   count.bins                     <- as.data.frame(table(data$bin))$Freq
   }
 }
 
 # Standardise the frequency distribution
 # (calculate a new standardised variable that expresses the frequencies as a proportion of the total frequency count)
 standardised_freq = count.bins/sum(count.bins)

 # Return results
 return(standardised_freq)
  
}