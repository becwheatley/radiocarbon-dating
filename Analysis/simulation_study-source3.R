#--------------------------------------------------------------------------------------------------------------------------------
# RADIOCARBON DATING PROJECT
# SIMULATION STUDY - INVESTIGATE THE EFFECT OF SAMPLING BIAS ON COMMON ANALYSIS RESULTS
# Source functions
# Code by Rebecca Wheatley
# Last modified 26 November 2021
#--------------------------------------------------------------------------------------------------------------------------------

# MULTIPLE DRAWS OF THE BASELINE DATA SET, USES THE FIRST ONE FOR SUB-SAMPLING

#------------------------------------------------
# I. GET EVIDENCE OF HUMAN OCCUPATION
#------------------------------------------------

## Generate some potential evidence of human occupation at some sites with no biases other than (potentially) pull of the recent
## across the specified time range
## @param timeRange  = the time range we are interested in sampling over
## @param no_sites   = the number of sites we want in our data set
## @param no_samples = the number of samples we want to take per site
## @param pop_trend  = the underlying trend in the population we want to mimic
## @param nsim       = the number of replicates we want to generate
get_available_evidence <- function(timeRange, no_sites, no_samples, pop_trend, nsim)
{
 
  # Create an empty list to store our baseline data replicates
  baseline_data <- vector("list", length = 2)
  
  # Create another two empty lists to store our replicates of "no loss" and "taphonomic loss" baseline data sets
  noloss   <- vector("list", length = nsim)
  taphloss <- vector("list", length = nsim)
  
  # Get baseline data
  for (n in 1:nsim){
 
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
        temp <- extraDistr::rdunif(no_samples*20, min = timeRange[2-3000], max = timeRange[1]+3000)
      
      } else if (pop_trend == "steady growth"){
        samples <- round(EnvStats::rtri(no_samples*100, min = timeRange[2]-3000, max = timeRange[1]+3000, mode = timeRange[2]))
        samples2 <- unlist(lapply(samples, function(x) {if(x >= timeRange[2]-3000){return(x)}}))
        temp <- samples2[1:(no_samples*20)]
    
      } else if (pop_trend == "exponential growth"){
        temp <- round(truncdist::rtrunc(n = no_samples*20, spec = "exp", a = timeRange[2]-3000, b = timeRange[1]+3000, rate = 0.3/500))
    
      } else if (pop_trend == "growth then decline"){
        temp <- round(EnvStats::rtri(no_samples*20, min = timeRange[2]-3000, max = timeRange[1]+3000, mode = (timeRange[1] - timeRange[2])/2))
      
      } else if (pop_trend == "growth decline growth"){
        samples <- round(EnvStats::rtri(no_samples*100, min = timeRange[2] - (timeRange[1] - timeRange[2])/3 - 1000, max = timeRange[2] + (timeRange[1] - timeRange[2])/3 - 1000, mode = timeRange[2]))
        samples2 <- unlist(lapply(samples, function(x) {if(x >= timeRange[2]-1000){return(x)}}))
        temp1 <- samples2[1:(no_samples*10)] 
        temp2 <- round(EnvStats::rtri(no_samples*10, min = timeRange[2] + (timeRange[1] - timeRange[2])/3, max = timeRange[2] + (timeRange[1] - timeRange[2]) + 1000, mode = timeRange[1] - (timeRange[1] - timeRange[2])/3))
        temp <- c(temp1, temp2)
        
      #} else if (pop_trend == "growth then plateau"){
         
      }
      
      # Get rid of samples that are outside the bounds of the region of interest
      samples1 <- unlist(lapply(temp, function(x) {if(x >= timeRange[2] && x <= timeRange[1]){return(x)}}))
      
      # Retain only the number of samples desired (including a buffer for taphonomic loss)
      samples2 <- samples1[1:(no_samples*5)]
      samples <- samples2[order(samples2)]
      
      # Get the occupation history for the site
      occupation_history[s,1] = min(samples)
      occupation_history[s,2] = max(samples)
      
      ## if simulating taphonomic loss, sample from the temp data using weights from an exponential distribution
      ### note that our rate is informed by the results from the taphonomic loss dynamic occupancy model on AustArch
      weights.taphloss <- dexp(samples, rate = 0.04/500)
      data.taphloss1   <- sample(x = samples, size = no_samples, replace = FALSE, prob = weights.taphloss)
      data.taphloss    <- data.taphloss1[order(data.taphloss1)]
      
      ## if not simulating taphonomic loss, sample from the temp data uniformly
      weights.noloss <- dunif(samples, min = timeRange[2], max = timeRange[1])
      data.noloss1   <- sample(x = samples, size = no_samples, replace = FALSE, prob = weights.noloss)
      data.noloss    <- data.noloss1[order(data.noloss1)]
      
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
    
    ## save evidence
    noloss[[n]]   <- evidence.noloss
    taphloss[[n]] <- evidence.taphloss
    
  }
  
  # add our evidence to our big list
  baseline_data[[1]] <- noloss
  baseline_data[[2]] <- taphloss
  
  return(baseline_data)

}

#--------------------------------------------------------------------------------------------------------------------------------
# II. GET THE MEDIAN BASELINE DATA (FROM SITE EVIDENCE SETS)
#--------------------------------------------------------------------------------------------------------------------------------
## Get the median evidence of occupation at some sites with no biases other than (potentially) pull of the recent
## across the specified time range (essentially take the median of the baseline data sets generated in the function above)
## @param data  = the baseline data sets we want to take the median of
get_median_evidence <- function(data)
{
  median.taphloss  <- data.frame(matrix(data = NA, nrow = length(data[[1]][[1]]$age), ncol = 3))
  median.noloss    <- median.taphloss
  
  baseline.noloss         <- matrix(NA, nrow = length(data[[1]][[1]]$age), ncol = length(data[[1]]))
  baseline.taphloss       <- matrix(NA, nrow = length(data[[1]][[2]]$age), ncol = length(data[[2]]))
  baseline.noloss.error   <- matrix(NA, nrow = length(data[[1]][[1]]$error), ncol = length(data[[1]]))
  baseline.taphloss.error <- matrix(NA, nrow = length(data[[1]][[2]]$error), ncol = length(data[[2]]))
  
  for (x in 1:length(data[[1]])){
    baseline.noloss[,x] <- data[[1]][[x]]$age
    baseline.noloss.error[,x] <- data[[1]][[x]]$error
    
    baseline.taphloss[,x] <- data[[2]][[x]]$age
    baseline.taphloss.error[,x] <- data[[2]][[x]]$error
  }
  
  median.noloss[,1]   <- median.taphloss[,1] <- data[[1]][[1]]$site
  median.noloss[,2]   <- apply(baseline.noloss, 1, quantile, prob = c(0.500))
  median.taphloss[,2] <- apply(baseline.taphloss, 1, quantile, prob = c(0.500))
  median.noloss[,3]   <- apply(baseline.noloss.error, 1, quantile, prob = c(0.500))
  median.taphloss[,3] <- apply(baseline.taphloss.error, 1, quantile, prob = c(0.500))
  
  names(median.noloss) <- names(median.taphloss) <- c("site", "age", "error")

  return(list(median.noloss, median.taphloss))
}

#--------------------------------------------------------------------------------------------------------------------------------
# III. GENERATE A SAMPLE (FROM SITE EVIDENCE)
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
  
  # If we are trying to EMULATE THE AUSTARCH data
  } else if (sampling_method == "emulate_AustArch"){
    
    # So that we are not always sampling from the same sites, assign each site a number from 1 to the total number of sites, 
    # in random order, then arrange the data by the site's random number
    site.RN <- sample(1:length(unique(evidence$site)), length(unique(evidence$site)), replace = FALSE)
    evidence$site_RN <- site.RN[rleid(evidence$site)]
    evidence <- evidence %>% 
      arrange(site_RN)
    
    # Create new data frames to store our samples in
    samples1 <- evidence[FALSE,]
    samples2 <- evidence[FALSE,]
    samples34 <- evidence[FALSE,]
    samples5 <- evidence[FALSE,]
    
    # SINGLE SAMPLE: first 50% of sites
    
    ## Set initial integer for loop
    p <- 1
    q <- 1
    s <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first 50% of sites
      if (evidence$site_RN[i] <= (50/100) * length(unique(evidence$site))) {
        
        ## if the site in this row is NOT the same as the site in the next row (i.e. it is the last sample for the site)
        if (evidence$site_RN[i] != evidence$site_RN[i+1]) { 
          
          ## write the line to the new data file (taking the last sample of each site)
          samples1[p,] <- evidence[i,]
          p <- p + 1
        }
        
      ## if we are in the following 20% of sites
      } else if ((evidence$site_RN[i] > (50/100) * length(unique(evidence$site))) &&
        (evidence$site_RN[i] <= (70/100) * length(unique(evidence$site)))){
        
        ## if the site in the current row is not the same as the previous row (i.e. it is the first sample for a new site) OR the 
        ## site in the current row is not the same as the next row (i.e. it is the last sample for the site), add it to the data frame
        if (evidence$site_RN[i] != evidence$site_RN[i-1] || evidence$site_RN[i] != evidence$site_RN[i+1]) {
          
          samples2[q,] <- evidence[i,]
          q <- q + 1
        }
      
      ## if we are in the following 15% of sites  
      } else if ((evidence$site_RN[i] > (70/100) * length(unique(evidence$site))) &&
                 (evidence$site_RN[i] <= (85/100) * length(unique(evidence$site)))){
        
        ## sample exhaustively
        samples5[s,] <- evidence[i,]
        s <- s + 1
        
      ## if we are in the last 15% of sites  
      } else {
        
        # For each site remaining
        for (j in (86/100 * no_sites):no_sites) {
          
          ## get the data for this site only
          just_this_site <- subset(evidence, site_RN == j)
          
          ## take 3 or 4 samples from a uniform distribution
          x <- floor(runif(1, min = 3, max = 5))
            
          ## choose which samples to take from this site
          site_samples <- sample(1:nrow(just_this_site), x, replace = FALSE) ## this will be no_samples long
            
          for (r in 1:length(site_samples))
          {
            samples34 <- rbind(samples34, just_this_site[site_samples[r],])
          }
        }
      }
    }
  }
  
  return(samples)
}


#--------------------------------------------------------------------------------------------------------------------------------
# IV. GENERATE REPEATED SAMPLES (FROM SITE EVIDENCE)
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
# V. CALIBRATE SAMPLE DATA
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
# VI. CREATE MULTIPLE SPDS WITH DIFFERENT DATA SETS
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
  simulatedSPD <- matrix(NA, nrow = nrow(calyears), ncol = nsim+1)
  
  # For each sample:
  for (s in 1:length(calibrated_samples)) 
  {
    
    ## Generate the SPD
    tmpSPD <- spd(calibrated_samples[[s]], timeRange = timeRange, runm = runm, spdnormalised = normalised)
    
    ## Save to SPD vector
    if (s == 1){ simulatedSPD[,s] <- tmpSPD[[2]]$calBP }
    simulatedSPD[,s+1] <- tmpSPD[[2]]$PrDens
    
  }
  
  # Return results
  return(simulatedSPD)
}

#--------------------------------------------------------------------------------------------------------------------------------
# VI. COMPARE SAMPLE SPDS
#--------------------------------------------------------------------------------------------------------------------------------

## Generate and compare the mean SPD
## @param baseline_SPD       = the baseline SPD we want to compare our subsetted SPDs to
## @param calibrated_samples = our calibrated sample data
## @param timeRange          = the time range we are interested in sampling over
## @param runm               = the running mean to be used when creating the SPD
## @param normalised         = logical for normalising the calibration curves and SPD (TRUE or FALSE)
compare_spds <- function(baseline_SPD, calibrated_samples, timeRange, runm, normalised){
  
  #-----------------------
  # Get subsampled SPDs
  #-----------------------
  sample_SPDs <- generate_multiple_spds(calibrated_samples = calibrated_samples, timeRange = timeRange, runm = runm, 
                                        normalised = normalised)
  
  #-----------------------
  # Calculate confidence intervals
  #-----------------------
  
  # Extract calibrated years BP
  sample_calBP <- sample_SPDs[,1]
  
  # Combine subsampled spds into list
  subsampled.spds <- list(sample_SPDs[,2:ncol(sample_SPDs)])
  
  # Retrieve confidence intervals
  subsampledCIlist <- vector("list",length=length(subsampled.spds))
  for (x in 1:length(subsampled.spds)) {
    subsampledCIlist[[x]] <- cbind(apply(subsampled.spds[[x]], 1, quantile, prob = c(0.025)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.500)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.975)))
  }
  
  #--------------------
  # Calculate p-value
  # based on https://rdrr.io/cran/rcarbon/src/R/tests.R
  #--------------------
  
  pValueList <- numeric(length = length(subsampled.spds))
  
  # For each subsample type:
  for (x in 1:length(subsampled.spds)) {
    
    ## Create Vector of Means
    zscoreMean <- apply(subsampled.spds[[x]], 1, mean)
    
    ## Create Vector of SDs
    zscoreSD <- apply(subsampled.spds[[x]], 1, sd)
    
    ## Z-Transform baseline and sub-sampled spds
    tmp.subsampled   <- t(apply(subsampled.spds[[x]], 1, function(p){ return((p - mean(p))/sd(p)) }))
    baseline         <- baseline_SPD$grid$PrDens
    tmp.baseline     <- (baseline - zscoreMean)/zscoreSD
    
    ## Compute CI
    tmp.ci <- t(apply(tmp.subsampled, 1, quantile, prob = c(0.025, 0.975), na.rm = TRUE))
    
    ## Compute expected statistic
    expected.statistic <- abs(apply(tmp.subsampled, 2, function(x, y){ a = x-y; i = which(a<0); return(sum(a[i])) }, y = tmp.ci[,1])) +
      apply(tmp.subsampled, 2, function(x, y){ a = x-y; i = which(a>0); return(sum(a[i])) }, y = tmp.ci[,2])
    
    ## Compute observed statistic
    lower    <- tmp.baseline - tmp.ci[,1]
    indexLow <- which(tmp.baseline < tmp.ci[,1])
    higher   <- tmp.baseline - tmp.ci[,2]
    indexHi  <- which(tmp.baseline > tmp.ci[,2])
    observed.statistic <- sum(abs(lower[indexLow])) + sum(higher[indexHi])
    
    ## Calculate p-value
    pValueList[[x]] <- 1
    if (observed.statistic > 0) {    
      pValueList[[x]] <- c(length(expected.statistic[expected.statistic > observed.statistic]) + 1)/c(nsim + 1)    
    }
  }
  
  # RETURN OUTPUT
  return(list(calBP = sample_calBP, envelope = subsampledCIlist, raw = subsampled.spds, pvalue = pValueList))
  
}


#--------------------------------------------------------------------------------------------------------------------------------
# VII. CALCULATE P VALUES AND TOTAL DISCREPANCY FROM MONTE CARLO SIMULATIONS
#--------------------------------------------------------------------------------------------------------------------------------

## Calculate the global p-value and total discrepancy for monte carlo simulations based on a hypothetical growth model (modelTest 
## results)
## @param model_test_results = a modelTest object comparing the observed curve to Monte Carlo simulations of a hypothetical model
calculate_p_value <- function(model_test_results){
  
  #--------------------------
  # Extract and combine data
  #--------------------------
  
  spd.sims <- model_test_results$sim
  
  colnames(spd.sims) <- paste0("s", 1:model_test_results$nsim)
  
  spd.data <- spd.sims %>% as_tibble %>% 
    mutate(calBP = model_test_results$result$calBP, obs = model_test_results$result$PrDens)
  
  nsim <- model_test_results$nsim
  
  spd.data.long <- spd.data %>% pivot_longer(-calBP, names_to = "rep", values_to = "pd")
  
  #--------------------------
  # Calculate p values
  #--------------------------
  
  spd.data$calBP %>% head
  
  p_sig <- 0.05
  
  dat_sig_p <- 
    spd.data%>% 
    select(-calBP) %>% 
    apply(1, function(row) (row - mean(row))^2) %>% # squared area
    apply(1, cumsum) %>% # cumulative area, i.e., integral
    apply(1, function(row) ((nsim + 2) - rank(row))/(nsim + 1)) %>% # progressive p-value
    {tibble(calBP = spd.data$calBP, p =.["obs",], pd = spd.data$obs)} %>% 
    mutate(sig = p<=p_sig, grp = c(0,diff(sig))) %>% 
    filter(sig) %>% 
    mutate(sig = NULL, grp = cumsum(grp))
  
  p_total <- dat_sig_p %>% slice(n()) %>% pull(p)
  
  plot1 <- 
    dat_sig_p %>% ggplot(aes(-calBP,p)) + geom_line() + ylim(0,0.06) + 
    geom_hline(aes(yintercept = p_sig), lty = "dashed") +
    labs(subtitle = "p-values computed cumulatively from left to right")
  
  #--------------------------
  # Compute total discrepancy
  #--------------------------
  
  discrep <- spd.data.long %>% group_by(calBP) %>% summarise(pd = mean(pd)) %>% 
             mutate(obs = spd.data$obs, 
                    dis = discrepancy_fun(pd,obs), 
                    dis = cumsum(dis)) %>% 
             slice(n()) %>% # take the last value to get the total
             pull(dis)
  
  #--------------------
  # Plot results
  #--------------------
  
  plot2 <- 
    spd.data.long %>% 
    ggplot(aes(-calBP,pd)) +
    geom_line(aes(group = rep), alpha = 0.2, color = blues9[5], size = 0.5) +
    geom_line(col = "grey10", size = 1, lty = "dashed",
              data = ~.x %>% group_by(calBP) %>% summarise(pd = mean(pd))) +
    geom_line(col = blues9[9], size = 0.8,
              data = ~ .x %>% filter(rep == "obs")) +
    geom_line(aes(group = grp), col = "red", data = dat_sig_p, size = 0.9) +
    labs(title = "Summed probabilty distribution",
         subtitle= paste0("Solid line is the emprical curve (red if p < ", p_sig,", blue otherwise)",
                          "\n","p-value for the total data set is ", p_total,
                          "\n", "discrepancy for the total data set is ", discrep))

  #--------------------------
  # Return results
  #--------------------------
  
  return(list(p_sig = p_sig, p_total = p_total, discrepancy = discrep, plot1 = plot1, plot2 = plot2))
  
}

discrepancy_fun <- function(x, y, q = 1, p =2) (abs(x^q - y^q))^p