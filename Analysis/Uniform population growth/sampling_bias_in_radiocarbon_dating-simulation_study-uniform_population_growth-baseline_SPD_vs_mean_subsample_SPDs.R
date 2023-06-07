#-------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (UNIFORM population growth)
# subtitle: "Baseline SPD vs mean subsample SPDs"
# author: "Rebecca Wheatley"
# date: "16 May 2023"
#-------------------------------------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required packages
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

# Load required functions
source("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/simulation_study-source3.R")

# Set seed for repeatability
set.seed(1000)

#-------------------------------------------------------------------------------------
# I. SIMULATE SOME BASELINE DATA AND THEN CREATE BIASED SUBSAMPLES
#-------------------------------------------------------------------------------------

# SET PARAMETERS:
timeRange        = c(12000, 1000)  ## for Holocene sites
no_sites         = 100             ## the number of sites we want in our simulated data set
no_samples       = 5               ## the number of samples we want each site to have
pop_trend        = "no change"     ## the underlying population trend we want to mimic
sampling_effort  = 3               ## the number of samples we want to take using the uniform sampling method
nsim             = 5               ## the number of times we want to replicate each sample

# GET BASELINE DATA SET/S:
baseline.data <- get_available_evidence(timeRange = timeRange, no_sites = no_sites, no_samples = no_samples, 
                                        pop_trend = pop_trend, nsim = nsim)

# GENERATE BIASED SUBSAMPLES:
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

# CALIBRATE DATA SETS:
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
save.image("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/uniform_pop_growth-raw_sample_and_calibrated_data-5s-5sims.RData")

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

# SAVE DISCREPANCY SCORES AND TOTAL P-VALUES
pvals <- data.frame(matrix(NA, nrow = 15, ncol = 3))
names(pvals) <- c("sample", "discrepancy", "p-value")

pvals[1,1] <- "uniform.50p"
pvals[1,2] <- spd.uniform.50p$discrepancy
pvals[1,3] <- if(length(spd.uniform.50p$pvalue>0)){ spd.uniform.50p$pvalue }else{ "NA" }
pvals[2,1] <- "uniform.75p"
pvals[2,2] <- spd.uniform.75p$discrepancy
pvals[2,3] <- if(length(spd.uniform.75p$pvalue>0)){ spd.uniform.75p$pvalue }else{ "NA" }
pvals[3,1] <- "uniform.100p"
pvals[3,2] <- spd.uniform.100p$discrepancy
pvals[3,3] <- if(length(spd.uniform.100p$pvalue>0)){ spd.uniform.100p$pvalue }else{ "NA" }

pvals[4,1] <- "singleton.ancient.50p"
pvals[4,2] <- spd.singleton.ancient.50p$discrepancy
pvals[4,3] <- if(length(spd.singleton.ancient.50p$pvalue>0)){ spd.singleton.ancient.50p$pvalue }else{ "NA" }
pvals[5,1] <- "singleton.ancient.75p"
pvals[5,2] <- spd.singleton.ancient.75p$discrepancy
pvals[5,3] <- if(length(spd.singleton.ancient.75p$pvalue>0)){ spd.singleton.ancient.75p$pvalue }else{ "NA" }
pvals[6,1] <- "singleton.ancient.100p"
pvals[6,2] <- spd.singleton.ancient.100p$discrepancy
pvals[6,3] <- if(length(spd.singleton.ancient.100p$pvalue>0)){ spd.singleton.ancient.100p$pvalue }else{ "NA" }

pvals[7,1] <- "singleton.recent.50p"
pvals[7,2] <- spd.singleton.recent.50p$discrepancy
pvals[7,3] <- if(length(spd.singleton.recent.50p$pvalue>0)){ spd.singleton.recent.50p$pvalue }else{ "NA" }
pvals[8,1] <- "singleton.recent.75p"
pvals[8,2] <- spd.singleton.recent.75p$discrepancy
pvals[8,3] <- if(length(spd.singleton.recent.75p$pvalue>0)){ spd.singleton.recent.75p$pvalue }else{ "NA" }
pvals[9,1] <- "singleton.recent.100p"
pvals[9,2] <- spd.singleton.recent.100p$discrepancy
pvals[9,3] <- if(length(spd.singleton.recent.100p$pvalue>0)){ spd.singleton.recent.100p$pvalue }else{ "NA" }

pvals[10,1] <- "singleton.random.50p"
pvals[10,2] <- spd.singleton.random.50p$discrepancy
pvals[10,3] <- if(length(spd.singleton.random.50p$pvalue>0)){ spd.singleton.random.50p$pvalue }else{ "NA" }
pvals[11,1] <- "singleton.random.75p"
pvals[11,2] <- spd.singleton.random.75p$discrepancy
pvals[11,3] <- if(length(spd.singleton.random.75p$pvalue>0)){ spd.singleton.random.75p$pvalue }else{ "NA" }
pvals[12,1] <- "singleton.random.100p"
pvals[12,2] <- spd.singleton.random.100p$discrepancy
pvals[12,3] <- if(length(spd.singleton.random.100p$pvalue>0)){ spd.singleton.random.100p$pvalue }else{ "NA" }

pvals[13,1] <- "bracketed.50p"
pvals[13,2] <- spd.bracketed.50p$discrepancy
pvals[13,3] <- if(length(spd.bracketed.50p$pvalue>0)){ spd.bracketed.50p$pvalue }else{ "NA" }
pvals[14,1] <- "bracketed.75p"
pvals[14,2] <- spd.bracketed.75p$discrepancy
pvals[14,3] <- if(length(spd.bracketed.75p$pvalue>0)){ spd.bracketed.75p$pvalue }else{ "NA" }
pvals[15,1] <- "bracketed.100p"
pvals[15,2] <- spd.bracketed.100p$discrepancy
pvals[14,3] <- if(length(spd.bracketed.100p$pvalue>0)){ spd.bracketed.100p$pvalue }else{ "NA" }

write.csv(pvals, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/uniform_pop_growth-baseline_SPD_vs_mean_subsample_SPDs-5ssamples_5sims.csv")

#------------------------------------------------------------
# III. SAVE DATA AND PLOT
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

write.csv(plot.spd.noloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/Uniform_pop_growth-plot_data-SPD-5samples_5sims.csv")

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
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
singleton.random <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  ylim(0, 0.0007) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
singleton.ancient <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
singleton.recent <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
bracketed <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
baseline <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
spd.noloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                        labels = c("A", "B", "C", "D", "E"),
                        ncol = 2, nrow = 3,
                        common.legend = TRUE)
spd.noloss
ggexport(spd.noloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Uniform_population_growth-baseline_SPD_vs_mean_subsamples-5samples_5sims.png",
         width = 1500,
         height = 1000)

# Save the workspace
save.image("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/uniform_pop_growth-raw_sample_and_calibrated_data_plus_SPDs-5s-5sims.RData")
