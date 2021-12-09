#-------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (no population change) - 5 samples/site"
# author: "Rebecca Wheatley"
# date: "2 December 2021"
#-------------------------------------------------------------------------------------------------------------------------------------------

library(rcarbon)
library(ggformula)
library(extraDistr)
library(EnvStats)
library(truncdist)
library(data.table)
library(tidyverse)
library(ggpubr)
theme_set(
  theme_bw() +
    theme(legend.position = "top"))

source("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/simulation_study-source3.R")

#-------------------------------------------------------------------------------------
# I. Simulating no change in population growth over time (without taphonomic loss)
#-------------------------------------------------------------------------------------

# Set the parameters we want for our simulated data:
timeRange        = c(12000, 200)   ## for Holocene sites
no_sites         = 100             ## the number of sites we want in our simulated data set
no_samples       = 5               ## the number of samples we want each site to have
pop_trend        = "no change"     ## the underlying population trend we want to mimic
sampling_effort  = 3               ## the number of samples we want to take using the uniform sampling method
nsim             = 99              ## the number of times we want to replicate each sample

# Get the baseline data set/s:
baseline.data <- get_available_evidence(timeRange = timeRange, no_sites = no_sites, no_samples = no_samples, 
                                        pop_trend = pop_trend, nsim = nsim)

# GET SAMPLES:
sample.uniform.50p           <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "uniform", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 50, 
                                            nsim = nsim)
sample.uniform.75p           <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "uniform", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 75, 
                                            nsim = nsim)
sample.uniform.100p          <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "uniform", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 100, 
                                            nsim = nsim)

sample.singleton.ancient.50p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                             sampling_method = "singleton_ancient", 
                                             sampling_effort = sampling_effort, 
                                             percent_sites = 50, 
                                             nsim = nsim)
sample.singleton.ancient.75p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                             sampling_method = "singleton_ancient", 
                                             sampling_effort = sampling_effort, 
                                             percent_sites = 75, 
                                             nsim = nsim)
sample.singleton.ancient.100p <- get_samples(evidence = baseline.data[[1]][[1]], 
                                             sampling_method = "singleton_ancient", 
                                             sampling_effort = sampling_effort, 
                                             percent_sites = 100, 
                                             nsim = nsim)

sample.singleton.recent.50p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "singleton_recent", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 50, 
                                            nsim = nsim)
sample.singleton.recent.75p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "singleton_recent", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 75, 
                                            nsim = nsim)
sample.singleton.recent.100p <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "singleton_recent", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 100, 
                                            nsim = nsim)

sample.singleton.random.50p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "singleton_random", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 50, 
                                            nsim = nsim)
sample.singleton.random.75p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "singleton_random", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 75, 
                                            nsim = nsim)
sample.singleton.random.100p <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "singleton_random", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 100, 
                                            nsim = nsim)

sample.bracketed.50p         <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "bracketed", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 50, 
                                            nsim = nsim)
sample.bracketed.75p         <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "bracketed", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 75, 
                                            nsim = nsim)
sample.bracketed.100p        <- get_samples(evidence = baseline.data[[1]][[1]], 
                                            sampling_method = "bracketed", 
                                            sampling_effort = sampling_effort, 
                                            percent_sites = 100, 
                                            nsim = nsim)

# Calibrate samples:
normalised      = TRUE ## are calibration curves (and, later, SPDs) normalised?
ncores          = 8    ## the number of threads to use when calibrating the radiocarbon dates

cal.baseline.noloss1 <- rcarbon::calibrate(x = baseline.data[[1]][[1]]$age, 
                                           errors = baseline.data[[1]][[1]]$error, 
                                           calCurves = 'shcal20', 
                                           normalised = normalised)

cal.baseline.noloss <- calibrate_samples(samples = baseline.data[[1]], 
                                         normalised = normalised, 
                                         ncores = ncores)

cal.uniform.50p     <- calibrate_samples(samples = sample.uniform.50p, 
                                         normalised = normalised, 
                                         ncores = ncores)
cal.uniform.75p     <- calibrate_samples(samples = sample.uniform.75p, 
                                         normalised = normalised, 
                                         ncores = ncores)
cal.uniform.100p    <- calibrate_samples(samples = sample.uniform.100p, 
                                         normalised = normalised, 
                                         ncores = ncores)

cal.singleton.ancient.50p  <- calibrate_samples(samples = sample.singleton.ancient.50p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.singleton.ancient.75p  <- calibrate_samples(samples = sample.singleton.ancient.75p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.singleton.ancient.100p <- calibrate_samples(samples = sample.singleton.ancient.100p, 
                                                normalised = normalised, 
                                                ncores = ncores)

cal.singleton.recent.50p   <- calibrate_samples(samples = sample.singleton.recent.50p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.singleton.recent.75p   <- calibrate_samples(samples = sample.singleton.recent.75p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.singleton.recent.100p  <- calibrate_samples(samples = sample.singleton.recent.100p, 
                                                normalised = normalised, 
                                                ncores = ncores)

cal.singleton.random.50p   <- calibrate_samples(samples = sample.singleton.random.50p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.singleton.random.75p   <- calibrate_samples(samples = sample.singleton.random.75p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.singleton.random.100p  <- calibrate_samples(samples = sample.singleton.random.100p, 
                                                normalised = normalised, 
                                                ncores = ncores)

cal.bracketed.50p          <- calibrate_samples(samples = sample.bracketed.50p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.bracketed.75p          <- calibrate_samples(samples = sample.bracketed.75p, 
                                                normalised = normalised, 
                                                ncores = ncores)
cal.bracketed.100p         <- calibrate_samples(samples = sample.bracketed.100p, 
                                                normalised = normalised, 
                                                ncores = ncores)

# Save the workspace
save.image("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/uniform_pop_growth-raw_sample_and_calibrated_data-5s-99sims.RData")

#--------------------------------------------------------------------------------------
# II. COMPARE THE SPD FOR THE BASELINE DATA TO THE MEAN SPD FOR EACH BIASED SUBSAMPLE
#--------------------------------------------------------------------------------------

runm = 100 ## the running mean to use for the SPDs

spd.baseline.noloss1 <- spd(x = cal.baseline.noloss1, timeRange = timeRange, runm = runm, spdnormalised = normalised)

# Get the median and 95% CI for the SPDs for the baseline data sets
spd.baseline.noloss  <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                     calibrated_samples = cal.baseline.noloss, 
                                     subsampling_method = "baseline (replicates)",
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)

# Get the median and 95% CI for the SPDs for the subsampled data sets
spd.uniform.50p      <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                     calibrated_samples = cal.uniform.50p, 
                                     subsampling_method = "uniform (50% sites)", 
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)
spd.uniform.75p      <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                     calibrated_samples = cal.uniform.75p, 
                                     subsampling_method = "uniform (75% sites)",
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)
spd.uniform.100p     <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                     calibrated_samples = cal.uniform.100p,  
                                     subsampling_method = "uniform (100% sites)",
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)

spd.singleton.ancient.50p  <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                           calibrated_samples = cal.singleton.ancient.50p, 
                                           subsampling_method = "singleton ancient (50% sites)", 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.ancient.75p  <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                           calibrated_samples = cal.singleton.ancient.75p, 
                                           subsampling_method = "singleton ancient (75% sites)", 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.ancient.100p <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                           calibrated_samples = cal.singleton.ancient.100p, 
                                           subsampling_method = "singleton ancient (100% sites)",
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)

spd.singleton.recent.50p   <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                           calibrated_samples = cal.singleton.recent.50p,  
                                           subsampling_method = "singleton recent (50% sites)",
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.recent.75p   <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                           calibrated_samples = cal.singleton.recent.75p,   
                                           subsampling_method = "singleton recent (75% sites)",
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.recent.100p  <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                           calibrated_samples = cal.singleton.recent.100p, 
                                           subsampling_method = "singleton recent (100% sites)",
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)

spd.singleton.random.50p  <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                          calibrated_samples = cal.singleton.random.50p, 
                                          subsampling_method = "singleton random (50% sites)",
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.singleton.random.75p  <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                          calibrated_samples = cal.singleton.random.75p,  
                                          subsampling_method = "singleton random (75% sites)", 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.singleton.random.100p <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                          calibrated_samples = cal.singleton.random.100p,   
                                          subsampling_method = "singleton random (100% sites)",
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)

spd.bracketed.50p         <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                          calibrated_samples = cal.bracketed.50p,  
                                          subsampling_method = "bracketed (50% sites)",
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.bracketed.75p         <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                          calibrated_samples = cal.bracketed.75p,  
                                          subsampling_method = "bracketed (75% sites)", 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.bracketed.100p        <- compare_spds(baseline_SPD = spd.baseline.noloss1,
                                          calibrated_samples = cal.bracketed.100p,   
                                          subsampling_method = "bracketed (100% sites)",
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)

# Get total discrepancy scores for subsampled SPDs relative to the baseline SPD

spd.uniform.50p$discrepancy  ## 1.65e-08
spd.uniform.75p$discrepancy  ## 2.21e-08
spd.uniform.100p$discrepancy ## 3.6e-08
spd.uniform.50p$pvalue  ## NA (~1 as mean expected SPD roughly = observed SPD)
spd.uniform.75p$pvalue  ## NA (~1 as mean expected SPD roughly = observed SPD)
spd.uniform.100p$pvalue ## NA (~1 as mean expected SPD roughly = observed SPD)

spd.singleton.ancient.50p$discrepancy  ## 1.38e-06
spd.singleton.ancient.75p$discrepancy  ## 8.62e-05
spd.singleton.ancient.100p$discrepancy ## 0.000188
spd.singleton.ancient.50p$pvalue  ## 0.05
spd.singleton.ancient.75p$pvalue  ## 0.05
spd.singleton.ancient.100p$pvalue ## 0.01 (minimum bounded by #sims)

spd.singleton.recent.50p$discrepancy  ## 5.48e-06
spd.singleton.recent.75p$discrepancy  ## 2.6e-06
spd.singleton.recent.100p$discrepancy ## 0.000155
spd.singleton.recent.50p$pvalue  ## 0.01 (minimum bounded by #sims)
spd.singleton.recent.75p$pvalue  ## 0.01 (minimum bounded by #sims)
spd.singleton.recent.100p$pvalue ## 0.01 (minimum bounded by #sims)

spd.singleton.random.50p$discrepancy  ## 4.41e-08
spd.singleton.random.75p$discrepancy  ## 6.64e-08
spd.singleton.random.100p$discrepancy ## 2.8e-07
spd.singleton.random.50p$pvalue  ## NA (~1 as mean expected SPD roughly = observed SPD)
spd.singleton.random.75p$pvalue  ## NA (~1 as mean expected SPD roughly = observed SPD)
spd.singleton.random.100p$pvalue ## NA (~1 as mean expected SPD roughly = observed SPD)

spd.bracketed.50p$discrepancy  ## 3.71e-06
spd.bracketed.75p$discrepancy  ## 1.5e-05
spd.bracketed.100p$discrepancy ## 6.23e-05
spd.bracketed.50p$pvalue  ## 0.01 (minimum bounded by #sims)
spd.bracketed.75p$pvalue  ## 0.01 (minimum bounded by #sims)
spd.bracketed.100p$pvalue ## 0.01 (minimum bounded by #sims)

#------------------------------------------------------------
# SAVE DATA AND PLOT
#------------------------------------------------------------

# Rearrange data for plotting:
plot.spd.noloss <- data.frame(matrix(NA, nrow = length(spd.baseline.noloss$calBP)*17, ncol = 5))
names(plot.spd.noloss) <- c("sample", "years.ka", "lowerCI", "mean", "upperCI")

for (i in 1:length(spd.baseline.noloss$calBP)){
  plot.spd.noloss[i,1] <- "baseline (replicates)"
  plot.spd.noloss[i,2] <- spd.baseline.noloss$calBP[i]/1000
  plot.spd.noloss[i,3] <- spd.baseline.noloss$envelope[[1]][i,1]
  plot.spd.noloss[i,4] <- spd.baseline.noloss$envelope[[1]][i,2]
  plot.spd.noloss[i,5] <- spd.baseline.noloss$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),1] <- "singleton ancient (50% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),2] <- spd.singleton.ancient.50p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),3] <- spd.singleton.ancient.50p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),4] <- spd.singleton.ancient.50p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),5] <- spd.singleton.ancient.50p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),1] <- "singleton ancient (75% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),2] <- spd.singleton.ancient.75p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),3] <- spd.singleton.ancient.75p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),4] <- spd.singleton.ancient.75p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),5] <- spd.singleton.ancient.75p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),1] <- "singleton ancient (100% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),2] <- spd.singleton.ancient.100p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),3] <- spd.singleton.ancient.100p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),4] <- spd.singleton.ancient.100p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),5] <- spd.singleton.ancient.100p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),1] <- "singleton recent (50% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),2] <- spd.singleton.recent.50p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),3] <- spd.singleton.recent.50p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),4] <- spd.singleton.recent.50p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),5] <- spd.singleton.recent.50p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),1] <- "singleton recent (75% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),2] <- spd.singleton.recent.75p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),3] <- spd.singleton.recent.75p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),4] <- spd.singleton.recent.75p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),5] <- spd.singleton.recent.75p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),1] <- "singleton recent (100% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),2] <- spd.singleton.recent.100p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),3] <- spd.singleton.recent.100p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),4] <- spd.singleton.recent.100p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),5] <- spd.singleton.recent.100p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*7+i),1] <- "singleton random (50% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*7+i),2] <- spd.singleton.random.50p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*7+i),3] <- spd.singleton.random.50p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*7+i),4] <- spd.singleton.random.50p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*7+i),5] <- spd.singleton.random.50p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*8+i),1] <- "singleton random (75% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*8+i),2] <- spd.singleton.random.75p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*8+i),3] <- spd.singleton.random.75p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*8+i),4] <- spd.singleton.random.75p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*8+i),5] <- spd.singleton.random.75p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*9+i),1] <- "singleton random (100% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*9+i),2] <- spd.singleton.random.100p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*9+i),3] <- spd.singleton.random.100p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*9+i),4] <- spd.singleton.random.100p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*9+i),5] <- spd.singleton.random.100p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*10+i),1] <- "bracketed (50% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*10+i),2] <- spd.bracketed.50p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*10+i),3] <- spd.bracketed.50p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*10+i),4] <- spd.bracketed.50p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*10+i),5] <- spd.bracketed.50p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*11+i),1] <- "bracketed (75% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*11+i),2] <- spd.bracketed.75p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*11+i),3] <- spd.bracketed.75p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*11+i),4] <- spd.bracketed.75p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*11+i),5] <- spd.bracketed.75p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*12+i),1] <- "bracketed (100% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*12+i),2] <- spd.bracketed.100p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*12+i),3] <- spd.bracketed.100p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*12+i),4] <- spd.bracketed.100p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*12+i),5] <- spd.bracketed.100p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*13+i),1] <- "uniform (50% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*13+i),2] <- spd.uniform.50p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*13+i),3] <- spd.uniform.50p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*13+i),4] <- spd.uniform.50p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*13+i),5] <- spd.uniform.50p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*14+i),1] <- "uniform (75% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*14+i),2] <- spd.uniform.75p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*14+i),3] <- spd.uniform.75p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*14+i),4] <- spd.uniform.75p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*14+i),5] <- spd.uniform.75p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*15+i),1] <- "uniform (100% sites)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*15+i),2] <- spd.uniform.100p$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*15+i),3] <- spd.uniform.100p$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*15+i),4] <- spd.uniform.100p$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*15+i),5] <- spd.uniform.100p$envelope[[1]][i,3]
  
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),1] <- "baseline"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),2] <- spd.baseline.noloss1$grid$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),3] <- 0
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),4] <- spd.baseline.noloss1$grid$PrDens[i]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),5] <- 0
}

write.csv(plot.spd.noloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/Plot_data-SPD-5samples_99sims.csv")

# Plot (combined):
group.colors <- c("baseline (replicates)" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline" = "black")

# Plot (multi panelled):
uniform <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  ylim(0, 0.0007) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
spd.noloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                        labels = c("A", "B", "C", "D", "E"),
                        ncol = 2, nrow = 3,
                        common.legend = TRUE)
spd.noloss
ggexport(spd.noloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Uniform_population_growth-baseline_SPD_vs_mean_subsamples-5samples_99sims.png",
         width = 1500,
         height = 1000)

#------------------------------------------------------------
# COMPARE SUBSAMPLED SPDS AGAINST THEORETICAL GROWTH MODELS
#------------------------------------------------------------

## Using the first subsample of each grouping
nsim2 = 99

#-------------------
# UNIFORM
#-------------------

# Uniform sampling @ 50% sites
uniform.50p.vs.uniform <- modelTest(cal.uniform.50p[[1]], errors = sample.uniform.50p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                    model = "uniform", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.50p.vs.uniform.p_values <- calculate_p_value(uniform.50p.vs.uniform, "uniform")
uniform.50p.vs.uniform.p_values$plot2
uniform.50p.vs.uniform.p_values$p_total
uniform.50p.vs.uniform.p_values$discrepancy
ggsave(uniform.50p.vs.uniform.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_50p_vs_uniform-5s.png")


# Uniform sampling @ 75% sites
uniform.75p.vs.uniform <- modelTest(cal.uniform.75p[[1]], errors = sample.uniform.75p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                    model = "uniform", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.75p.vs.uniform.p_values <- calculate_p_value(uniform.75p.vs.uniform, "uniform")
uniform.75p.vs.uniform.p_values$plot2
uniform.75p.vs.uniform.p_values$p_total
uniform.75p.vs.uniform.p_values$discrepancy
ggsave(uniform.75p.vs.uniform.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_75p_vs_uniform-5s.png")


# Uniform sampling @ 100% sites
uniform.100p.vs.uniform <- modelTest(cal.uniform.100p[[1]], errors = sample.uniform.100p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                     model = "uniform", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.100p.vs.uniform.p_values <- calculate_p_value(uniform.100p.vs.uniform, "uniform")
uniform.100p.vs.uniform.p_values$plot2
uniform.100p.vs.uniform.p_values$p_total
uniform.100p.vs.uniform.p_values$discrepancy
ggsave(uniform.100p.vs.uniform.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_100p_vs_uniform-5s.png")



#-------------------
# LINEAR
#-------------------

# Uniform sampling @ 50% sites
uniform.50p.vs.linear <- modelTest(cal.uniform.50p[[1]], errors = sample.uniform.50p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                   model = "linear", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.50p.vs.linear.p_values <- calculate_p_value(uniform.50p.vs.linear, "linear")
uniform.50p.vs.linear.p_values$plot2
uniform.50p.vs.linear.p_values$p_total
uniform.50p.vs.linear.p_values$discrepancy
ggsave(uniform.50p.vs.linear.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_50p_vs_linear-5s.png")


# Uniform sampling @ 75% sites
uniform.75p.vs.linear <- modelTest(cal.uniform.75p[[1]], errors = sample.uniform.75p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                   model = "linear", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.75p.vs.linear.p_values <- calculate_p_value(uniform.75p.vs.linear, "linear")
uniform.75p.vs.linear.p_values$plot2
uniform.75p.vs.linear.p_values$p_total
uniform.75p.vs.linear.p_values$discrepancy
ggsave(uniform.75p.vs.linear.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_75p_vs_linear-5s.png")


# Uniform sampling @ 100% sites
uniform.100p.vs.linear <- modelTest(cal.uniform.100p[[1]], errors = sample.uniform.100p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                    model = "linear", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.100p.vs.linear.p_values <- calculate_p_value(uniform.100p.vs.linear, "linear")
uniform.100p.vs.linear.p_values$plot2
uniform.100p.vs.linear.p_values$p_total
uniform.100p.vs.linear.p_values$discrepancy
ggsave(uniform.100p.vs.linear.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_100p_vs_linear-5s.png")


#-------------------
# EXPONENTIAL
#-------------------

# Uniform sampling @ 50% sites
uniform.50p.vs.exponential <- modelTest(cal.uniform.50p[[1]], errors = sample.uniform.50p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                        model = "exponential", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.50p.vs.exponential.p_values <- calculate_p_value(uniform.50p.vs.exponential, "exponential")
uniform.50p.vs.exponential.p_values$plot2
uniform.50p.vs.exponential.p_values$p_total     
uniform.50p.vs.exponential.p_values$discrepancy
ggsave(uniform.50p.vs.exponential.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_50p_vs_exponential-5s.png")


# Uniform sampling @ 75% sites
uniform.75p.vs.exponential <- modelTest(cal.uniform.75p[[1]], errors = sample.uniform.75p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                        model = "exponential", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.75p.vs.exponential.p_values <- calculate_p_value(uniform.75p.vs.exponential, "exponential")
uniform.75p.vs.exponential.p_values$plot2
uniform.75p.vs.exponential.p_values$p_total     
uniform.75p.vs.exponential.p_values$discrepancy
ggsave(uniform.75p.vs.exponential.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_75p_vs_exponential-5s.png")


# Uniform sampling @ 100% sites
uniform.100p.vs.exponential <- modelTest(cal.uniform.100p[[1]], errors = sample.uniform.100p[[1]]$error, nsim = nsim2, timeRange = timeRange, 
                                         model = "exponential", runm = runm, raw = TRUE)
## calculate p-values and discrepancy scores
uniform.100p.vs.exponential.p_values <- calculate_p_value(uniform.100p.vs.exponential, "exponential")
uniform.100p.vs.exponential.p_values$plot2
uniform.100p.vs.exponential.p_values$p_total     
uniform.100p.vs.exponential.p_values$discrepancy
ggsave(uniform.100p.vs.exponential.p_values$plot2, 
       file = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_100p_vs_exponential-5s.png")

