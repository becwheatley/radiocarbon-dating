#--------------------------------------------------------------------------------------------------------------------------------
# RADIOCARBON DATING PROJECT
# SIMULATION STUDY - INVESTIGATE THE EFFECT OF SAMPLING BIAS ON COMMON ANALYSIS RESULTS
# Source functions
# Code by Rebecca Wheatley
# Last modified 23 August 2021
#--------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# I. GET EVIDENCE OF HUMAN OCCUPATION
#------------------------------------------------

## Generate some potential evidence of human occupation at some sites with no biases other than (potentially) pull of the recent
## across the specified time range
## @param timeRange  = the time range we are interested in sampling over
## @param no_sites   = the number of sites we want in our data set
## @param no_samples = the number of samples we want to take per site
## @param pop_trend  = the underlying trend in the population we want to mimic
get_available_evidence <- function(timeRange, no_sites, no_samples, pop_trend)
{
  ## Generate the real occupation period for each site (year max, year min in cal BP) using:
  ## 1. A uniform distribution (implies a stable population over time)
  ## 2. An exponential distribution (implies exponential population growth over time)
  ## 3. A triangular distribution (with values below the mode discarded - implies steady population growth over time)
  ## 4. A triangular distribution (implies population growth then decline)
  occupation_history <- matrix(data = NA, nrow = no_sites, ncol = 2)
  evidence.taphloss  <- data.frame(matrix(data = NA, nrow = no_samples * no_sites, ncol = 5))
  evidence.noloss    <- evidence.taphloss
  names(evidence.taphloss) <- names(evidence.noloss) <- c("site", "sample", "age", "error", "Open.Closed")
  
  # For each site:
  for (s in 1:no_sites){
  
    # Establish the occupation range:
    if (pop_trend == "no change"){
      #temp <- extraDistr::rdunif(2, min = timeRange[2]-3000, max = timeRange[1]+3000)
      temp <- extraDistr::rdunif(no_samples*10, min = timeRange[2]-3000, max = timeRange[1]+3000)
      
    } else if (pop_trend == "steady growth"){
      samples <- round(EnvStats::rtri(no_samples*100, min = timeRange[2]-3000, max = timeRange[1]+3000, mode = timeRange[2]))
      samples2 <- unlist(lapply(samples, function(x) {if(x >= timeRange[2]-3000){return(x)}}))
      temp <- samples2[1:(no_samples*10)]
    
    } else if (pop_trend == "exponential growth"){
      temp <- round(truncdist::rtrunc(n = no_samples*10, spec = "exp", a = timeRange[2]-3000, b = timeRange[1]+3000, rate = 0.3/500))
    
    } else if (pop_trend == "growth then decline"){
      temp <- round(EnvStats::rtri(no_samples*10, min = timeRange[2]-3000, max = timeRange[1]+3000, mode = (timeRange[1] - timeRange[2])/2))
    
    #} else if (pop_trend == "growth decline growth"){
      #samples <- round(EnvStats::rtri(100, min = timeRange[2] - (timeRange[1] - timeRange[2])/3 - 1000, max = timeRange[2] + (timeRange[1] - timeRange[2])/3 - 1000, mode = timeRange[2]))
      #samples2 <- unlist(lapply(samples, function(x) {if(x >= timeRange[2]-1000){return(x)}}))
      #temp1 <- samples2[1:2] 
      #temp2 <- round(EnvStats::rtri(2, min = timeRange[2] + (timeRange[1] - timeRange[2])/3, max = timeRange[2] + (timeRange[1] - timeRange[2]) + 1000, mode = timeRange[1] - (timeRange[1] - timeRange[2])/3))
      #temp <- sample(c(temp1, temp2), size = 2, replace = FALSE)
      
    #} else if (pop_trend == "growth then plateau"){
       
    }
    
    # Get rid of samples that are outside the bounds of the region of interest
    samples1 <- unlist(lapply(temp, function(x) {if(x >= timeRange[2] && x <= timeRange[1]){return(x)}}))
    
    # Retain only the number of samples desired (including a buffer for taphonomic loss)
    samples <- samples1[1:no_samples*5]
    
    # Get the occupation history for the site
    occupation_history[s,1] = min(samples)
    occupation_history[s,2] = max(samples)
    
    ## if simulating taphonomic loss, sample from the temp data using weights from an exponential distribution
    ### note that our rate is informed by the results from the taphonomic loss dynamic occupancy model on AustArch
    weights.taphloss <- dexp(samples, rate = 0.04/500)
    data.taphloss    <- sample(x = samples, size = no_samples, replace = FALSE, prob = weights.taphloss)
    
    ## if not simulating taphonomic loss, sample from the temp data uniformly
    weights.noloss <- dunif(samples, min = timeRange[2], max = timeRange[1])
    data.noloss    <- sample(x = samples, size = no_samples, replace = FALSE, prob = weights.noloss)
    
    for (i in 1:no_samples)
    {
      evidence.noloss[(s-1)*no_samples+i, 1] <- s
      evidence.noloss[(s-1)*no_samples+i, 2] <- i
      evidence.noloss[(s-1)*no_samples+i, 3] <- data.noloss[i]
      
      evidence.taphloss[(s-1)*no_samples+i, 1] <- s
      evidence.taphloss[(s-1)*no_samples+i, 2] <- i
      evidence.taphloss[(s-1)*no_samples+i, 3] <- data.taphloss[i]
    }
  }
  
  ## construct normally distributed errors (mean = 100, sd = 50)
  evidence.noloss$error   <- round(rtnorm(no_sites * no_samples, mean = 100, sd = 50, a = 0, b = 500))
  evidence.taphloss$error <- round(rtnorm(no_sites * no_samples, mean = 100, sd = 50, a = 0, b = 500))
  
  ## set all sites to open
  evidence.noloss$Open.Closed <- evidence.taphloss$Open.Closed <- "Open"
  
  ## return evidence
  return(list(evidence.noloss, evidence.taphloss))
}

#--------------------------------------------------------------------------------------------------------------------------------
# II. GENERATE A SAMPLE (FROM SITE EVIDENCE)
#--------------------------------------------------------------------------------------------------------------------------------

## Sample the evidence at each site
## @param evidence        = the possible samples that can be pulled for each site
## @param sampling_method = method of sampling the data (uniform, singleton_ancient, singleton_recent, singleton_random, bracketed)
## @param sampling_effort = the number of samples to take per site (uniform only)
## @param percent_sites   = the percent of sites to apply sampling method to - others sites are exhaustively sampled
get_a_sample <- function(evidence, sampling_method, sampling_effort, percent_sites)
{
  # UNIFORM sampling takes a number of random samples from a uniform distribution for each site
  if (sampling_method == "uniform"){
    
    # So that we are not always sampling from the same sites, assign each site a number from 1 to the total number of sites, 
    # in random order, then arrange the data by the site's random number
    site.RN <- sample(1:length(unique(evidence$site)), length(unique(evidence$site)), replace = FALSE)
    evidence$site_RN <- site.RN[rleid(evidence$site)]
    evidence <- evidence %>% 
      arrange(site_RN)
    
    # Create a new data frame to store our samples in
    samples <- evidence[FALSE,]
    
    # For each site in our evidence set...
    for (i in 1:length(unique(evidence$site))) {
      
      ## get the data for this site only
      just_this_site <- subset(evidence, site_RN == i)
      
      ## if we are in the first percent% of sites, take no_samples from a uniform distribution
      if (i <= (percent_sites/100) * length(unique(evidence$site))){
      
        ## choose which samples to take from this site
        site_samples <- sample(1:nrow(just_this_site), sampling_effort, replace = FALSE) ## this will be no_samples long

        ## this might work      
        for (p in 1:length(site_samples))
        {
          samples <- rbind(samples, just_this_site[site_samples[p],])
        }
        
      ## if we are not in the first percent% of sites, sample exhaustively
      } else {
        samples <- rbind(samples, just_this_site)
      }
    }
    
    
  # SINGLETON ANCIENT sampling takes the most ancient date only from a certain % of sites (the rest are exhaustively sampled)  
  } else if (sampling_method == "singleton_ancient"){
    
    # So that we are not always sampling from the same sites, assign each site a number from 1 to the total number of sites, 
    # in random order, then arrange the data by the site's random number
    site.RN <- sample(1:length(unique(evidence$site)), length(unique(evidence$site)), replace = FALSE)
    evidence$site_RN <- site.RN[rleid(evidence$site)]
    evidence <- evidence %>% 
      arrange(site_RN)
    
    # Create a new data frame to store our samples in
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site_RN[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the last row of the data set (i.e. the last sample of the last site)
        if (i == nrow(evidence)) {
          samples[p,] <- evidence[i,]
          p <- p + 1
        
        ## if the site in this row is NOT the same as the site in the next row (i.e. it is the last sample for the site)
        } else if (evidence$site_RN[i] != evidence$site_RN[i+1]) { 
          
          ## write the line to the new data file (taking the last sample of each site)
          samples[p,] <- evidence[i,]
          p <- p + 1
        }
        
      ## if we are not within the first percent% of sites
      } else {
        
        ## write the row to the sample data frame
        samples[p,] <- evidence[i,]
        p <- p + 1
      }
    }
  
  
  # SINGLETON RANDOM sampling takes one random date from a certain % of sites (the rest are exhaustively sampled)  
  } else if (sampling_method == "singleton_random"){
    
    # So that we are not always sampling from the same sites, assign each site a number from 1 to the total number of sites, 
    # in random order
    site.RN <- sample(1:length(unique(evidence$site)), length(unique(evidence$site)), replace = FALSE)
    evidence$site_RN <- site.RN[rleid(evidence$site)]
    
    # Organise data randomly (to randomise the order of the samples within each site)
    rows <- sample(nrow(evidence))
    evidence <- evidence[rows, ]
    
    # Then arrange the data frame by the site's random number
    evidence <- evidence %>% 
      arrange(site_RN)
    
    # Create a new data frame to store our samples in
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site_RN[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the last row of the data set (i.e. the last sample of the last site)
        if (i == nrow(evidence)) {
          samples[p,] <- evidence[i,]
          p <- p + 1
          
          ## if the site in this row is NOT the same as the site in the next row (i.e. it is the last sample for the site)
        } else if (evidence$site_RN[i] != evidence$site_RN[i+1]) { 
          
          ## write the line to the new data file (taking the last sample of each site)
          samples[p,] <- evidence[i,]
          p <- p + 1
        }
        
        ## if we are not within the first percent% of sites
      } else {
        
        ## write the row to the sample data frame
        samples[p,] <- evidence[i,]
        p <- p + 1
      }
    }

    
  # SINGLETON RECENT sampling takes the most recent date only from a certain % of sites (the rest are exhaustively sampled)      
  } else if (sampling_method == "singleton_recent"){
    
    # So that we are not always sampling from the same sites, assign each site a number from 1 to the total number of sites, 
    # in random order, then arrange the data by the site's random number
    site.RN <- sample(1:length(unique(evidence$site)), length(unique(evidence$site)), replace = FALSE)
    evidence$site_RN <- site.RN[rleid(evidence$site)]
    evidence <- evidence %>% 
      arrange(site_RN)
    
    # Create a new data frame to store our samples in
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site_RN[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the first row of the data set (i.e. the first sample of the first site)
        if (i == 1) {
          samples[p,] <- evidence[i,]
          p <- p + 1
          
        ## if the site in this row is NOT the same as the site in the previous row (i.e. the first sample for a new site)
        } else if (evidence$site_RN[i] != evidence$site_RN[i-1]) {
          
          ## write the line to the new data file (taking the last sample of each site)
          samples[p,] <- evidence[i,]
          p <- p + 1
        }
        
      ## if we are not within the first percent% of sites
      } else {
        
        ## write the row to the sample data frame
        samples[p,] <- evidence[i,]
        p <- p + 1
      }
    }
  
    
  # BRACKETED sampling takes the most recent and most ancient dates only from a certain % of sites (the rest are exhaustively sampled)  
  } else if (sampling_method == "bracketed") {
    
    # So that we are not always sampling from the same sites, assign each site a number from 1 to the total number of sites, 
    # in random order, then arrange the data by the site's random number
    site.RN <- sample(1:length(unique(evidence$site)), length(unique(evidence$site)), replace = FALSE)
    evidence$site_RN <- site.RN[rleid(evidence$site)]
    evidence <- evidence %>% 
      arrange(site_RN)
    
    # Create a new data frame to store our samples in
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site_RN[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the very first row (the first sample of our first site) or the last row (the last sample of our last site),
        ## add it to the data frame
        if (i == 1 || i == nrow(evidence)) { 
          samples[p,] <- evidence[i,] 
          p <- p + 1
        
        ## if we are not in the first row of the data frame, and the site in the current row is not the same as the previous row 
        ## (i.e. it is the first sample for a new site) OR the site in the current row is not the same as the next row (i.e. it is the last
        ## sample for the site), add it to the data frame
        } else if (evidence$site_RN[i] != evidence$site_RN[i-1] || evidence$site_RN[i] != evidence$site_RN[i+1]) {
        
          samples[p,] <- evidence[i,]
          p <- p + 1
        }
        
      ## if we are not within the first percent% of sites
      } else {
        ## write the row to the sample data frame
        samples[p,] <- evidence[i,]
        p <- p + 1
      }
    }
  
  
  }
  
  return(samples)
}


#--------------------------------------------------------------------------------------------------------------------------------
# III. GENERATE REPEATED SAMPLES (FROM SITE EVIDENCE)
#--------------------------------------------------------------------------------------------------------------------------------

## Sample the evidence at each site nsim times (returns a list of data frames)
## @param evidence        = the possible samples that can be pulled for each site
## @param sampling_method = method of sampling the data (exhaustive, uniform, singleton_ancient, singleton_recent, singleton_random, bracketed)
## @param sampling_effort = the number of samples to take per site (only applies to uniform sampling)
## @param percent_sites   = the percent of sites to apply sampling method to - others sites are exhaustively sampled
## @param nsim            = the number of samples we want to take
get_samples <- function(evidence, sampling_method, sampling_effort, percent_sites, nsim){
  
  # Create an empty list to store our sample data
  sample_data <- vector("list", length = nsim)
  
  # Get sample data
  for (i in 1:nsim){
    sample_data[[i]] <- get_a_sample(evidence, sampling_method, sampling_effort, percent_sites)
    
  }
  
  return(sample_data)
}


#--------------------------------------------------------------------------------------------------------------------------------
# IV. CALIBRATE SAMPLE DATA
#--------------------------------------------------------------------------------------------------------------------------------

## Calibrate the sample data generated using get_samples (returns a list of data frames)
## @param samples    = the sample data (list of data frames)
## @param noramlised = should the calibration curves be normalised?
## @param ncores     = the number of threads to run the calibration across
calibrate_samples <- function(samples, normalised, ncores){
  
  # Create an empty list for the calibrated data (one entry for each sample)
  calibrated_data <- vector("list", length = length(samples))
  
  # For each sample:
  for (s in 1:length(samples)){
    
    ## Create vectors of the ages and their respective errors and ids in the data set
    ages   <- samples[[s]]$age
    errors <- samples[[s]]$error
    
    ## Calibrate these ages using the SHCal20 calibration curves (can calibrate dates up to 55,000 years old)
    calibrated_data[[s]] <- rcarbon::calibrate(x = ages, errors = errors, calCurves = 'shcal20', normalised = normalised, ncores = ncores)
  }
  
  return(calibrated_data)
}


#--------------------------------------------------------------------------------------------------------------------------------
# V. CREATE MULTIPLE SPDS WITH DIFFERENT DATA SETS
#--------------------------------------------------------------------------------------------------------------------------------

## Create a summed probability distribution for each sample
## @param samples         = our list of calibrated sampled data (length = nsim)
## @param timeRange       = the time range we are interested in sampling over
## @param runm            = the running mean to be used when creating the SPD
## @param normalised      = logical for normalising the calibration curves and SPD (TRUE or FALSE)
generate_multiple_spds <- function(calibrated_samples, timeRange, runm, normalised){
  
  # Get calyears
  calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
  
  # Create an empty matrix for the SPDs (one for each sample)
  simulatedSPD <- matrix(NA, nrow = nrow(calyears), ncol = nsim)
  
  # For each sample:
  for (s in 1:length(calibrated_samples)) 
  {
    
    ## Generate the SPD
    tmpSPD <- spd(calibrated_samples[[s]], timeRange = timeRange, runm = runm, spdnormalised = normalised)
    
    ## Save to SPD vector
    simulatedSPD[,s] <- tmpSPD[[2]]$PrDens
    
  }
  
  # Return results
  return(simulatedSPD)
}

#--------------------------------------------------------------------------------------------------------------------------------
# VI. COMPARE SAMPLE SPDS
#--------------------------------------------------------------------------------------------------------------------------------

## Generate and compare the mean SPD
## @param calibrated_samples = our calibrated sample data
## @param timeRange          = the time range we are interested in sampling over
## @param runm               = the running mean to be used when creating the SPD
## @param normalised         = logical for normalising the calibration curves and SPD (TRUE or FALSE)
compare_spds <- function(calibrated_samples, timeRange, runm, normalised){
  
  #-----------------------
  # Get subsampled SPDs
  #-----------------------
  sample_SPDs <- generate_multiple_spds(calibrated_samples = calibrated_samples, timeRange = timeRange, runm = runm, 
                                        normalised = normalised)
  
  #-----------------------
  # Calculate confidence intervals
  #-----------------------
  
  # Combine subsampled spds into list
  subsampled.spds <- list(sample_SPDs)
  
  # Retrieve confidence intervals
  subsampledCIlist <- vector("list",length=length(subsampled.spds))
  for (x in 1:length(subsampled.spds)) {
    subsampledCIlist[[x]] <- cbind(apply(subsampled.spds[[x]], 1, quantile, prob = c(0.025)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.500)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.975)))
  }
  
  # RETURN OUTPUT
  return(list(envelope = subsampledCIlist, raw = subsampled.spds))
  
}


#--------------------------------------------------------------------------------------------------------------------------------
# VII. GENERATE MULTIPLE FREQUENCY DISTRIBUTIONS USING DIFFERENT DATA SETS
#--------------------------------------------------------------------------------------------------------------------------------

## Create a frequency distribution of dates within bins for each sample
## @param data               = our list of uncalibrated data
## @param calibrated_data    = our list of calibrated data
## @param timeRange          = the time range we are interested in sampling over
## @param taphCorrect        = logical specifying whether to taphonomically correct open sites using Williams (2013) (TRUE or FALSE) 
## @param correctForSite     = logical specifying whether to count all radiocarbon dates within a bin or only one per site (TRUE OR FALSE) - note, presently taphCorrect does nothing if this is TRUE
## @param binSize            = size of bin to sort dates into for frequency distribution (200 years used in Williams 2013)
generate_multiple_frequency_dists <- function(data, calibrated_data, timeRange, taphCorrect, correctForSite, binSize){
  
  # Create an empty matrix for the frequency distributions
  simulatedFreq     <- matrix(NA, nrow = (timeRange[1] - timeRange[2])/binSize, ncol = length(data))
  
  # For each sample:
  for (s in 1:length(data)) 
  {
    
    ## Get MedianBP values for each calibrated date
    summary <- summary(calibrated_data[[s]])
    
    ## Get age bins
    data[[s]]$bin <- cut(summary$MedianBP, 
                         breaks = c(timeRange[2], seq(timeRange[2]+binSize, timeRange[1], by = binSize)), 
                         labels = 1:((timeRange[1] - timeRange[2])/binSize)
                         )
    
    ## If we only want to count one occurance of a particular site per bin
    if (correctForSite){
      one_site_per_bin <- unique(data[[s]][c("site", "bin")])
      count.bins <- as.data.frame(table(one_site_per_bin$bin))$Freq
    
    ## Else if we want to count all radiocarbon dates
    } else {
      
      ## Apply taphonomic correction to open sites
      if (taphCorrect){
        open.sites                     <- data[[s]] %>% subset(Open.Closed == "Open")
        count.bins.open                <- as.data.frame(table(open.sites$bin))
        count.bins.open$corrected.freq <- count.bins.open$Freq/(2.107 * 10^7*(median.bin.age + 2754)^(-1.526))
        closed.sites                   <- data[[s]] %>% subset(Open.Closed == "Closed")
        count.bins.closed              <- as.data.frame(as.data.frame(table(closed.sites$bin)))
    
        ## Combine taphonomically corrected open site counts to non-corrected closed site counts to get the corrected frequency distribution in bins
        count.bins                     <- count.bins.open$corrected.freq + count.bins.closed$Freq
     } else {
    
      ## Get uncorrected frequency distribution in bins
      count.bins                     <- as.data.frame(table(data[[s]]$bin))$Freq
      }
    }
    
    ## Save to frequency matrix
    simulatedFreq[,s] <- count.bins
  }
  
  ## Return results
  return(simulatedFreq)
  
}


#--------------------------------------------------------------------------------------------------------------------------------
# VIII. COMPARE SAMPLE FREQUENCY DISTRIBUTIONS TO BASELINE FREQUENCY DISTRIBUTION
#--------------------------------------------------------------------------------------------------------------------------------

## Generate and compare the mean sample frequency distributions to the baseline frequency distribution
## @param evidence            = our full simulated data set
## @param samples             = our sample data
## @param calibrated_evidence = our full calibrated simulated data set
## @param calibrated_samples  = our calibrated sample data
## @param timeRange           = the time range we are interested in sampling over
## @param taphCorrect         = logical specifying whether to taphonomically correct open sites using Williams (2013) (TRUE or FALSE)
## @param correctForSite      = logical specifying whether to count all radiocarbon dates within a bin or only one per site (TRUE OR FALSE) - note, presently taphCorrect does nothing if this is TRUE
## @param binSize             = size of bin to sort dates into for frequency distribution (200 years used in Williams 2013)
compare_frequency_dists <- function(samples, calibrated_samples, timeRange, taphCorrect, correctForSite, 
                                    binSize){
  
  # Get frequency distributions for sub-samples
  sample_freq <- generate_multiple_frequency_dists(data = samples, calibrated_data = calibrated_samples, timeRange = timeRange,
                                                   taphCorrect = taphCorrect, correctForSite = correctForSite, binSize = binSize)

  # Combine subsampled frequency distributions into a list
  sample.freq.dists <- list(sample_freq)
  
  # Standardise the subsampled frequency distributions
  ## (calculate a new standardised variable that expresses the frequencies as a proportion of the total subsample frequency count)
  standardised_sample_freq <- matrix(nrow = nrow(sample_freq), ncol = ncol(sample_freq))
  for (i in 1:ncol(standardised_sample_freq)){
    standardised_sample_freq[,i] = sample_freq[,i]/sum(sample_freq[,i])
  }
  
  # Combine standardised subsampled frequency distributions into a list
  standardised.sample.freq.dists <- list(standardised_sample_freq)
  
  # Retrieve confidence intervals for raw and standardised subsamples
  sample.CI.list <- vector("list",length=length(sample.freq.dists))
  for (x in 1:length(sample.freq.dists)) {
    sample.CI.list[[x]] <- cbind(apply(sample.freq.dists[[x]], 1, quantile, prob = c(0.025)),
                                     apply(sample.freq.dists[[x]], 1, quantile, prob = c(0.500)),
                                     apply(sample.freq.dists[[x]], 1, quantile, prob = c(0.975)))
  }
  standardised.sample.CI.list <- vector("list",length=length(standardised.sample.freq.dists))
  for (x in 1:length(standardised.sample.freq.dists)) {
    standardised.sample.CI.list[[x]] <- cbind(apply(standardised.sample.freq.dists[[x]], 1, quantile, prob = c(0.025)),
                                                  apply(standardised.sample.freq.dists[[x]], 1, quantile, prob = c(0.500)),
                                                  apply(standardised.sample.freq.dists[[x]], 1, quantile, prob = c(0.975)))
  }
  
  # RETURN OUTPUT
  return(list(envelope = sample.CI.list, standardised.envelope = standardised.sample.CI.list,
              raw = sample.freq.dists, standardised.raw = standardised.sample.freq.dists))
}