#--------------------------------------------------------------------------------------------------------------------------------
# RADIOCARBON DATING PROJECT
# SIMULATION STUDY - INVESTIGATE THE EFFECT OF SAMPLING BIAS ON COMMON ANALYSIS RESULTS
# Source functions
# Code by Rebecca Wheatley
# Last modified 30 April 2021
#--------------------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------------------
# I. FUNCTIONS FOR SIMULATING DATA
#--------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------
# GET POTENTIAL EVIDENCE
#------------------------------------------------

## Generate some potential evidence of human occupation at some sites with no biases other than (potentially) pull of the recent
## across the specified time range
## @param timeRange  = the time range we are interested in sampling over
## @param no_sites   = the number of sites we want in our data set
## @param no_samples = the number of samples we want to take per site
## @param pull       = pull of the recent (logical)
get_available_evidence <- function(timeRange, no_sites, no_samples, pull)
{
  
  ## Generate the real occupation period for each site
  ## (year max, year min in cal BP) - at the moment this is just uniform distribution but really it could be more sophisticated, 
  ## e.g. have fewer sites have earlier occupations
  occupation_history <- matrix(data = NA, nrow = no_sites, ncol = 2)
  for (s in 1:no_sites){
    temp <- extraDistr::rdunif(2, min = timeRange[2], max = timeRange[1])
    if (temp[1] < temp[2]) {
      occupation_history[s,1] <- temp[1]
      occupation_history[s,2] <- temp[2]
    } else if (temp[1] > temp[2]) {
      occupation_history[s,1] <- temp[2]
      occupation_history[s,2] <- temp[1]
    }
  }
  
  # For each site, generate some evidence of human occupation that could be sampled from
  evidence <- data.frame(matrix(data = NA, nrow = no_samples * no_sites, ncol = 4))
  names(evidence) <- c("site", "sample", "age", "error")
  for (s in 1:no_sites){
    
    ## if pull of the recent is set to TRUE, draw from a truncated decaying exponential distribution
    if (pull){
      ## what I really want here is to design an exponential distribution that applies across the entire possible date history (e.g. 0 to 55,000 ybp)
      ## and then truncate this single distribution using occupation_history[s,1] and occupation_history[s,2] for the purposes of sampling
      temp1 <- rtrunc(no_samples, spec = "exp", a = occupation_history[s,1], b = occupation_history[s,2])#, rate = 1) ## this code doesn't work
      temp <- round(temp1)
    
    ## if pull of the recent is set to FALSE, draw from a discrete uniform distribution
    } else {
      temp <- rdunif(no_samples, occupation_history[s,1], occupation_history[s,2])
    }
    
    for (i in 1:no_samples)
    {
      evidence[(s-1)*no_samples+i, 1] <- s
      evidence[(s-1)*no_samples+i, 2] <- i
      evidence[(s-1)*no_samples+i, 3] <- temp[i]
    }
  }
  
  ## construct normally distributed errors (mean = 100, sd = 50)
  evidence$error <- round(rtnorm(no_sites * no_samples, mean = 100, sd = 50, a = 0, b = 500))
  
  ## return evidence
  return(evidence)
}

#--------------------------------------------------------------------------------------------------------------------------------
# II. GENERATE A SAMPLE (FROM SITE EVIDENCE)
#--------------------------------------------------------------------------------------------------------------------------------

## Sample the evidence at each site
## @param evidence        = the possible samples that can be pulled for each site
## @param sampling_method = method of sampling the data (exhaustive, uniform, singleton_ancient, singleton_recent, singleton_random, bracketed)
## @param no_samples      = the number of samples to take per site (only applies to uniform sampling)
## @param percent_sites   = the percent of sites to apply sampling method to - others sites are exhaustively sampled
get_a_sample <- function(evidence, sampling_method, no_samples, percent_sites)
{
  # EXHAUSTIVE sampling takes every row of evidence
  if (sampling_method == "exhaustive") {
    samples <- evidence
  
  # UNIFORM sampling takes a number of random samples from a uniform distribution for each site
  } else if (sampling_method == "uniform"){
    
    # Create empty data frame to store the new data
    samples <- evidence[FALSE,]
    
    # For each site in our evidence set...
    for (i in 1:length(unique(evidence$site))) {
      
      ## get the data for this site only
      just_this_site <- subset(evidence, site == i)
      
      ## if we are in the first percent% of sites, take no_samples from a uniform distribution
      if (i <= (percent_sites/100) * length(unique(evidence$site))){
      
        ## choose which samples to take from this site
        site_samples <- rdunif(no_samples, 1, nrow(just_this_site)) ## this will be no_samples long

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
    
    # Create empty data frame to store the new data
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the last row of the data set (i.e. the last sample of the last site)
        if (i == nrow(evidence)) {
          samples[p,] <- evidence[i,]
          p <- p + 1
        
        ## if the site in this row is NOT the same as the site in the next row (i.e. it is the last sample for the site)
        } else if (evidence$site[i] != evidence$site[i+1]) { 
          
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
    
    # Create a new copy of the baseline data set to work with
    random.order <- evidence
    
    # Organise data randomly (to randomise the order of the samples within each site)
    rows <- sample(nrow(random.order))
    random.order <- random.order[rows, ]
    
    # Arrange data by site random number
    random.order <- random.order %>% 
      arrange(site)
    
    # Create empty data frame to store the new data
    samples <- random.order[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(random.order)) {
      
      ## if we are within the first percent% of sites
      if (random.order$site[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the last row of the data set
        if (i == nrow(random.order)) {
          samples[p,] <- evidence[i,]
          p <- p + 1
        
        ## if the site in the current row is not the same as the site in the next row
        } else if (random.order$site[i] != random.order$site[i+1]) { 
          
          ## write the line to the new data file (taking the last sample of each site)
          samples[p,] <- random.order[i,]
          p <- p + 1
        }
        
      ## if we are not within the first percent% of sites
      } else {
        
        ## write the row to the sample data frame
        samples[p,] <- random.order[i,]
        p <- p + 1
      }
    }

    
  # SINGLETON RECENT sampling takes the most recent date only from a certain % of sites (the rest are exhaustively sampled)      
  } else if (sampling_method == "singleton_recent"){
    
    # Create empty data frame to store the new data
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the first row of the data set (i.e. the first sample of the first site)
        if (i == 1) {
          samples[p,] <- evidence[i,]
          p <- p + 1
          
        ## if the site in this row is NOT the same as the site in the previous row (i.e. the first sample for a new site)
        } else if (evidence$site[i] != evidence$site[i-1]) {
          
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
    
    # Create empty data frame to store the new data
    samples <- evidence[FALSE,]
    
    # Set initial integer for loop
    p <- 1
    
    for (i in 1:nrow(evidence)) {
      
      ## if we are within the first percent% of sites
      if (evidence$site[i] <= (percent_sites/100) * length(unique(evidence$site))) {
        
        ## if we are in the very first row (the first sample of our first site) or the last row (the last sample of our last site),
        ## add it to the data frame
        if (i == 1 || i == nrow(evidence)) { 
          samples[p,] <- evidence[i,] 
          p <- p + 1
        
        ## if we are not in the first row of the data frame, and the site in the current row is not the same as the previous row 
        ## (i.e. it is the first sample for a new site) OR the site in the current row is not the same as the next row (i.e. it is the last
        ## sample for the site), add it to the data frame
        } else if (evidence$site[i] != evidence$site[i-1] || evidence$site[i] != evidence$site[i+1]) {
        
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

## Sample the evidence at each site nsim times
## @param evidence        = the possible samples that can be pulled for each site
## @param sampling_method = method of sampling the data (exhaustive, uniform, singleton_ancient, singleton_recent, singleton_random, bracketed)
## @param no_samples      = the number of samples to take per site (only applies to uniform sampling)
## @param percent_sites   = the percent of sites to apply sampling method to - others sites are exhaustively sampled
## @param nsim            = the number of samples we want to take
get_samples <- function(evidence, sampling_method, no_samples, percent_sites, nsim){
  
  # Create an empty list to store our sampled data
  sample_data <- matrix(NA, ncol = nsim) ## work out what nrow is (will depend on sampling_method and percent_sites)
  
  for (i in 1:nsim){
    sample_data[,i] <- get_a_sample(evidence, sampling_method, no_samples, percent_sites)
    
  }
  
  return(sample_data)
}


#--------------------------------------------------------------------------------------------------------------------------------
# II. FUNCTION FOR CREATING MULTIPLE SPDS WITH DIFFERENT DATA SETS
#--------------------------------------------------------------------------------------------------------------------------------

## Simulate data using the functions above, and create a summed probability distribution for each simulated data set
## @param evidence        = our simulated complete data set
## @param timeRange       = the time range we are interested in sampling over
## @param no_samples      = the number of samples we want to take per site (if uniform, other methods have intrinsic sampling rate)
## @param percent_sites   = the percent of sites to use the sampling_method for (remainder are sampled exhaustively)
## @param pull            = pull of the recent (logical)
## @param sampling_method = sampling method to use for the simulated data (exhaustive, uniform, singleton_ancient, singleton_random, singleton_recent, bracketed)
## @param runm            = the running mean to be used when creating the SPD
## @param nsim            = the number of simulations/SPDs to generate
## @param normalised      = logical for normalising the calibration curves and SPD (TRUE or FALSE)
## @param ncores          = number of cores to use in parallel processing when calibrating the dates
generate_multiple_spds <- function(evidence, timeRange, no_sites, no_samples, percent_sites, pull, sampling_method, runm, nsim, normalised, ncores){
  
  ## Get calyears
  calyears <- data.frame(calBP=seq(timeRange[1], timeRange[2],-1))
  
  ## create empty vectors for the spds and frequency plots
  simulatedSPD <- matrix(NA, nrow = nrow(calyears), ncol = nsim)
  
  ## For each "simulation":
  for (s in 1:nsim) 
  {
      
    ## Take samples
    sample <- get_a_sample(evidence, sampling_method, no_samples, percent_sites)
    
    ## Create vectors of the ages and their respective errors and ids in the data set
    ages   <- sample$age
    errors <- sample$error
    
    ## Calibrate these ages using the SHCal20 calibration curves (can calibrate dates up to 55,000 years old)
    calibrated <- rcarbon::calibrate(x = ages, errors = errors, calCurves = 'shcal20', normalised = normalised, ncores = ncores)
      
    ## Create an SPD of the calibrated data
    tmpSPD <- spd(calibrated, timeRange, runm = runm, spdnormalised = normalised)
    
    ## Save to SPD vector
    simulatedSPD[,s] <- tmpSPD[[2]]$PrDens
    
  }
  
  # Return results
  return(simulatedSPD)
}

#--------------------------------------------------------------------------------------------------------------------------------
# III. FUNCTION FOR COMPARING SPDS
#--------------------------------------------------------------------------------------------------------------------------------

## Generate and compare the mean SPDs
## @param evidence        = our full simulated data set
## @param timeRange       = the time range we are interested in sampling over
## @param sampling_method = sampling method to use when subsampling the data (uniform, singleton_ancient, singleton_random, singleton_recent, bracketed)
## @param sampling_effort = the number of samples we want to take per site using the uniform sampling method
## @param percent_sites   = the percent of sites to use the sampling_method for (remainder are sampled exhaustively)
## @param pull            = pull of the recent (logical)
## @param runm            = the running mean to be used when creating the SPD
## @param nsim            = the number of simulations/SPDs to generate
## @param normalised      = logical for normalising the calibration curves and SPD (TRUE or FALSE)
## @param ncores          = number of cores to use in parallel processing when calibrating the dates
compare_spds <- function(evidence, timeRange, sampling_method, sampling_effort, percent_sites, pull, runm, nsim, normalised, ncores){
  
  #-----------------------
  # Get baseline SPD
  #-----------------------
  
  # Create vectors of the ages and their respective errors and ids in the data set
  ages   <- evidence$age
  errors <- evidence$error
  
  # Calibrate these ages using the SHCal20 calibration curves (can calibrate dates up to 55,000 years old)
  calibrated <- rcarbon::calibrate(x = ages, errors = errors, calCurves = 'shcal20', normalised = normalised, ncores = ncores)
  
  # Create an SPD of the calibrated data
  tmpSPD <- spd(calibrated, timeRange, runm = runm, spdnormalised = normalised)
  
  # Save to SPD vector
  baseline.SPD <- tmpSPD[[2]]
  
  #-----------------------
  # Get subsampled SPDs
  #-----------------------

  # % of sites sampled from a uniform distribution
  uniform_SPDs <- generate_multiple_spds(evidence = evidence, timeRange = timeRange, no_samples = no_samples,
                                         percent_sites = percent_sites, pull = pull, sampling_method = "uniform",
                                         runm = runm, nsim = nsim, normalised = normalised, ncores = ncores)
  
  # % of sites sampled using the method specified
  subsampled_SPDs <- generate_multiple_spds(evidence = evidence, timeRange = timeRange, no_samples = no_samples,
                                            percent_sites = percent_sites, pull = pull, sampling_method = sampling_method,
                                            runm = runm, nsim = nsim, normalised = normalised, ncores = ncores)

  #-----------------------
  # Calculate confidence intervals
  #-----------------------
  
  # Combine subsampled spds into list
  subsampled.spds <- list(uniform_SPDs, subsampled_SPDs)
  
  # Retrieve confidence intervals
  subsampledCIlist <- vector("list",length=length(subsampled.spds))
  for (x in 1:length(subsampled.spds)) {
    subsampledCIlist[[x]] <- cbind(apply(subsampled.spds[[x]], 1, quantile, prob = c(0.025)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.500)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.975)))
  }
  
  
  # COMPUTE GLOBAL P-VALUE
  pValueList <- numeric(length = length(subsampled.spds))
  
  # For each subsample type:
  for (x in 1:length(subsampled.spds)) {
    
    ## Create Vector of Means
    zscoreMean <- apply(subsampled.spds[[x]], 1, mean)
    
    ## Create Vector of SDs
    zscoreSD <- apply(subsampled.spds[[x]], 1, sd)
    
    ## Z-Transform baseline and sub-sampled spds
    tmp.subsampled   <- t(apply(subsampled.spds[[x]], 1, function(p){ return((p - mean(p))/sd(p)) }))
    tmp.baseline     <- baseline.SPD$PrDens
    tmp.baseline     <- (tmp.baseline - zscoreMean)/zscoreSD
    
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
  return(list(baseline = baseline.SPD, envelope = subsampledCIlist, raw = subsampled.spds, pValueList = pValueList))
  
}



#--------------------------------------------------------------------------------------------------------------------------------
# III. FUNCTION FOR COMPARING SPDS AND FREQUENCY DISTRIBUTIONS
#--------------------------------------------------------------------------------------------------------------------------------

## Generate and compare the mean SPDs and frequency distributions for (1) the sampling method specified and (2) a uniform sampling method
## against (3) the SPD and frequency distribution for the full evidence
## @param evidence        = our full simulated data set
## @param timeRange       = the time range we are interested in sampling over
## @param sampling_method = sampling method to use when subsampling the data (uniform, singleton_ancient, singleton_random, singleton_recent, bracketed)
## @param sampling_effort = the number of samples we want to take per site using the uniform sampling method
## @param percent_sites   = the percent of sites to use the sampling_method for (remainder are sampled exhaustively)
## @param pull            = pull of the recent (logical)
## @param runm            = the running mean to be used when creating the SPD
## @param nsim            = the number of simulations/SPDs to generate
## @param normalised      = logical for normalising the calibration curves and SPD (TRUE or FALSE)
## @param ncores          = number of cores to use in parallel processing when calibrating the dates
compare_analyses <- function(evidence, timeRange, sampling_method, sampling_effort, percent_sites, pull, runm, nsim, normalised, ncores){
  
  #-----------------------
  # Get baseline SPD
  #-----------------------
  
  # Create vectors of the ages and their respective errors and ids in the data set
  ages   <- evidence$age
  errors <- evidence$error
  
  # Calibrate these ages using the SHCal20 calibration curves (can calibrate dates up to 55,000 years old)
  calibrated <- rcarbon::calibrate(x = ages, errors = errors, calCurves = 'shcal20', normalised = normalised, ncores = ncores)
  
  # Create an SPD of the calibrated data
  tmpSPD <- spd(calibrated, timeRange, runm = runm, spdnormalised = normalised)
  
  # Save to SPD vector
  baseline.SPD <- tmpSPD[[2]]
  
  #-----------------------
  # Get subsampled SPDs
  #-----------------------
  
  # % of sites sampled from a uniform distribution
  uniform_SPDs <- generate_multiple_spds(evidence = evidence, timeRange = timeRange, no_samples = no_samples,
                                         percent_sites = percent_sites, pull = pull, sampling_method = "uniform",
                                         runm = runm, nsim = nsim, normalised = normalised, ncores = ncores)
  
  # % of sites sampled using the method specified
  subsampled_SPDs <- generate_multiple_spds(evidence = evidence, timeRange = timeRange, no_samples = no_samples,
                                            percent_sites = percent_sites, pull = pull, sampling_method = sampling_method,
                                            runm = runm, nsim = nsim, normalised = normalised, ncores = ncores)
  
  #-----------------------
  # Calculate confidence intervals
  #-----------------------
  
  # Combine subsampled spds into list
  subsampled.spds <- list(uniform_SPDs, subsampled_SPDs)
  
  # Retrieve confidence intervals
  subsampledCIlist <- vector("list",length=length(subsampled.spds))
  for (x in 1:length(subsampled.spds)) {
    subsampledCIlist[[x]] <- cbind(apply(subsampled.spds[[x]], 1, quantile, prob = c(0.025)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.500)),
                                   apply(subsampled.spds[[x]], 1, quantile, prob = c(0.975)))
  }
  
  
  # COMPUTE GLOBAL P-VALUE
  pValueList <- numeric(length = length(subsampled.spds))
  
  # For each subsample type:
  for (x in 1:length(subsampled.spds)) {
    
    ## Create Vector of Means
    zscoreMean <- apply(subsampled.spds[[x]], 1, mean)
    
    ## Create Vector of SDs
    zscoreSD <- apply(subsampled.spds[[x]], 1, sd)
    
    ## Z-Transform baseline and sub-sampled spds
    tmp.subsampled   <- t(apply(subsampled.spds[[x]], 1, function(p){ return((p - mean(p))/sd(p)) }))
    tmp.baseline     <- baseline.SPD$PrDens
    tmp.baseline     <- (tmp.baseline - zscoreMean)/zscoreSD
    
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
  return(list(baseline = baseline.SPD, envelope = subsampledCIlist, raw = subsampled.spds, pValueList = pValueList))
  
}

#--------------------------------------------------------------------------------------------------------------------------------
# IV. FUNCTION FOR TAKING REPEATED SUB-SAMPLES AND GENERATING A FREQUENCY DISTRIBUTIONS OF DATES WITHIN BINS FOR EACH SAMPLE
#--------------------------------------------------------------------------------------------------------------------------------

## Generate a frequency distribution of dates within bins for each sample
## @param data        = data set to be sampled over
## @param timeRange   = the time range to generate the plot over; c(MOST ANCIENT DATE, MOST RECENT DATE)
## @param sample      = logical for whether to sub-sample the data or use the full data set
## @param FUN         = function to sub-sample the data
## @param percent    = the random percent of sites to be selected for the sub-sample
## @param nsim        = the number of sub-samples/SPDs to generate
## @param normalised  = logical for normalising the calibration curves (TRUE or FALSE)
## @param taphCorrect = logical specifying whether to taphonomically correct open sites using Williams (2013) (TRUE or FALSE) 
## @param correctForSite = logical specifying whether to count all radiocarbon dates within a bin or only one per site (TRUE OR FALSE) - note, presently taphCorrect does nothing if this is TRUE
## @param binSize     = size of bin to sort dates into for frequency distribution (200 years used in Williams 2013)
## @param ncores      = number of cores to use in parallel processing when calibrating the dates
generate_multiple_frequency_dists <- function(data, timeRange, sample, FUN, percent, nsim, normalised, taphCorrect, correctForSite, binSize, ncores){
  
  ## Create an empty matrix for the frequency distributions
  simulatedFreq     <- matrix(NA, nrow = timeRange[1]/binSize, ncol = nsim)
  
  ## For each "simulation":
  for (s in 1:nsim) 
  {
    
    ## Take a sample, or use the full data set
    if (sample){
      data2 <- FUN(data, percent)
    } else {
      data2 <- data
    }
  
    ## Subset data for terrestrial vs marine (for calibration)
    data.T <- data2 %>% subset(Data.pertinent.for.time.series.analysis.or.calibration == "Terrestrial")
    data.M <- data2 %>% subset(Data.pertinent.for.time.series.analysis.or.calibration == "Marine")
  
    ## Create vectors of the ages and their respective errors and ids in the data set
    ages.T   <- data.T$AGE_NORM
    errors.T <- data.T$ERROR
    ids.T    <- data.T$ADSID
    ages.M   <- data.M$AGE_NORM
    errors.M <- data.M$ERROR
    ids.M    <- data.M$ADSID
    
    ## Calibrate these ages using the SHCal20 and MARINE20 calibration curves (can calibrate dates up to 55,000 years old)
    calibrate.T <- calibrate(x = ages.T, errors = errors.T, ids = ids.T, calCurves = 'shcal20', normalised = normalised, ncores = ncores)
    calibrate.M <- calibrate(x = ages.M, errors = errors.M, ids = ids.M, calCurves = 'marine20', normalised = normalised, ncores = ncores)
      
    ## Combine the two calibrated data sets
    calibrate.all <- combine(calibrate.T, calibrate.M)
    
    ## Get MedianBP values for each calibrated date
    summary <- summary(calibrate.all)
    
    ## Get age bins
    data2$bin <- cut(summary$MedianBP, breaks = c(timeRange[2], seq(binSize, timeRange[1], by = binSize)), labels = 1:(timeRange[1]/binSize))
    
    ## If we only want to count one occurance of a particular site per bin
    if (correctForSite){
      one_site_per_bin <- unique(data2[c("SITE", "bin")])
      count.bins <- as.data.frame(table(one_site_per_bin$bin))$Freq
    
    ## Else if we want to count all radiocarbon dates
    } else {
      
      ## Apply taphonomic correction to open sites
      if (taphCorrect){
        open.sites                     <- data2 %>% subset(Open.or.Closed.Site == "Open")
        count.bins.open                <- as.data.frame(table(open.sites$bin))
        count.bins.open$corrected.freq <- count.bins.open$Freq/(2.107 * 10^7*(median.bin.age + 2754)^(-1.526))
        closed.sites                   <- data2 %>% subset(Open.or.Closed.Site == "Closed")
        count.bins.closed              <- as.data.frame(as.data.frame(table(closed.sites$bin)))
    
        ## Combine taphonomically corrected open site counts to non-corrected closed site counts to get the corrected frequency distribution in bins
        count.bins                     <- count.bins.open$corrected.freq + count.bins.closed$Freq
      } else {
    
        ## Get uncorrected frequency distribution in bins
        count.bins                     <- as.data.frame(table(data2$bin))$Freq
      }
    }
    
    ## Save to frequency matrix
    simulatedFreq[,s] <- count.bins
  }
  
  ## Return results
  return(simulatedFreq)
  
}


#--------------------------------------------------------------------------------------------------------------------------------
# V. FUNCTION FOR COMPARING FREQUENCY DISTRIBUTIONS OF DATES WITHIN BINS
#--------------------------------------------------------------------------------------------------------------------------------

## Generate a frequency distribution of dates within bins (both taphonomically corrected and not)
## @param data        = full data set to be analysed
## @param timeRange   = the time range to generate the plot over; c(MOST ANCIENT DATE, MOST RECENT DATE)
## @param percent     = the random percent of sites to be selected for the sub-sample
## @param runm        = the running mean to be used when creating the SPD
## @param nsim        = the number of sub-samples/SPDs to generate
## @param normalised  = logical for normalising the calibration curves (TRUE or FALSE)
## @param taphCorrect = logical specifying whether to taphonomically correct open sites using Williams (2013) (TRUE or FALSE)
## @param correctForSite = logical specifying whether to count all radiocarbon dates within a bin or only one per site (TRUE OR FALSE) - note, presently taphCorrect does nothing if this is TRUE
## @param binSize     = size of bin to sort dates into for frequency distribution (200 years used in Williams 2013)
## @param ncores      = number of cores to use in parallel processing when calibrating the dates
## Get subsets
### % of sites with 1 randomly pulled sample
compare_frequency_dist_dates_within_bins <- function(data, timeRange, percent, runm, nsim, normalised, taphCorrect, correctForSite, binSize, ncores){
  
  ## Get frequency distribution for baseline data set (nsim = 1, sample = FALSE)
  baseline_freq     <- generate_multiple_frequency_dists(data = data, timeRange = timeRange, sample = FALSE, 
                                                         nsim = 1, normalised = normalised, taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, binSize = binSize, ncores = ncores)
  
  ## Get frequency distributions for sub-samples
  one_random_sample <- generate_multiple_frequency_dists(data = data, timeRange = timeRange, sample = TRUE,
                                                         FUN = get_1_random_sample, percent = percent,
                                                         nsim = nsim, normalised = normalised,
                                                         taphCorrect = taphCorrect, correctForSite = correctForSite, 
                                                         binSize = binSize, ncores = ncores)
  one_ancient_sample <- generate_multiple_frequency_dists(data = data, timeRange = timeRange, sample = TRUE,
                                                          FUN = get_1_ancient_sample, percent = percent,
                                                          nsim = nsim, normalised = normalised,
                                                          taphCorrect = taphCorrect, correctForSite = correctForSite, 
                                                          binSize = binSize, ncores = ncores)
  two_bracketed_samples <- generate_multiple_frequency_dists(data = data, timeRange = timeRange, sample = TRUE,
                                                             FUN = get_2_bracketed_samples, percent = percent,
                                                             nsim = nsim, normalised = normalised,
                                                             taphCorrect = taphCorrect, correctForSite = correctForSite, 
                                                             binSize = binSize, ncores = ncores)
  
  ## Combine subsampled frequency distributions into a list
  subsampled.freq.dists <- list(one_random_sample, one_ancient_sample, two_bracketed_samples)
  
  ## Standardise the baseline and subsampled frequency distributions
  ### (calculate a new standardised variable that expresses the frequencies as a proportion of the total subsample frequency count)
  standardised_baseline_freq <- matrix(nrow = nrow(baseline_freq), ncol = ncol(baseline_freq))
  for (i in 1:ncol(standardised_baseline_freq)){
    standardised_baseline_freq[,i] = baseline_freq[,i]/sum(baseline_freq[,i])
  }
  
  standardised_one_random_sample <- matrix(nrow = nrow(one_random_sample), ncol = ncol(one_random_sample))
  for (i in 1:ncol(standardised_one_random_sample)){
    standardised_one_random_sample[,i] = one_random_sample[,i]/sum(one_random_sample[,i])
  }
  
  standardised_one_ancient_sample <- matrix(nrow = nrow(one_ancient_sample), ncol = ncol(one_ancient_sample))
  for (i in 1:ncol(standardised_one_ancient_sample)){
    standardised_one_ancient_sample[,i] = one_ancient_sample[,i]/sum(one_ancient_sample[,i])
  }
  
  standardised_two_bracketed_samples <- matrix(nrow = nrow(two_bracketed_samples), ncol = ncol(two_bracketed_samples))
  for (i in 1:ncol(standardised_two_bracketed_samples)){
    standardised_two_bracketed_samples[,i] = two_bracketed_samples[,i]/sum(two_bracketed_samples[,i])
  }
  
  ## Combine standardised subsampled frequency distributions into a list
  standardised.subsampled.freq.dists <- list(standardised_one_random_sample, standardised_one_ancient_sample, standardised_two_bracketed_samples)
  
  ## Retrieve confidence intervals for raw and standardised subsamples
  subsampled.CI.list <- vector("list",length=length(subsampled.freq.dists))
  for (x in 1:length(subsampled.freq.dists)) {
    subsampled.CI.list[[x]] <- cbind(apply(subsampled.freq.dists[[x]], 1, quantile, prob = c(0.025)),
                                     apply(subsampled.freq.dists[[x]], 1, quantile, prob = c(0.500)),
                                     apply(subsampled.freq.dists[[x]], 1, quantile, prob = c(0.975)))
  }
  standardised.subsampled.CI.list <- vector("list",length=length(standardised.subsampled.freq.dists))
  for (x in 1:length(standardised.subsampled.freq.dists)) {
    standardised.subsampled.CI.list[[x]] <- cbind(apply(standardised.subsampled.freq.dists[[x]], 1, quantile, prob = c(0.025)),
                                                  apply(standardised.subsampled.freq.dists[[x]], 1, quantile, prob = c(0.500)),
                                                  apply(standardised.subsampled.freq.dists[[x]], 1, quantile, prob = c(0.975)))
  }
  
  # RETURN OUTPUT
  return(list(baseline = baseline_freq, standardised.baseline = standardised_baseline_freq, 
              envelope = subsampled.CI.list, standardised.envelope = standardised.subsampled.CI.list,
              raw = subsampled.freq.dists, standardised.raw = standardised.subsampled.freq.dists))
}



