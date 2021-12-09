#-------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (no population change) - 5 samples"
# author: "Rebecca Wheatley"
# date: "21 October 2021"
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

source("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/simulation_study-source2.R")

set.seed(1234)

# Set the parameters we want for our simulated data:
timeRange        = c(12000, 200)   ## for Holocene sites
no_sites         = 100             ## the number of sites we want in our simulated data set
no_samples       = 5               ## the number of samples we want each site to have
pop_trend        = "no change"     ## the underlying population trend we want to mimic
#percent_sites    = 50              ## the number of sites we want to apply the sampling method across (the remainder will be exhaustively sampled)
sampling_effort  = 3               ## the number of samples we want to take using the uniform sampling method
nsim             = 20              ## the number of times we want to replicate each sample

#-------------------------------------------------------------------------------------
# I. SIMULATING NO CHANGE IN POPULATION OVER TIME (WITHOUT TAPHONOMIC LOSS)
#-------------------------------------------------------------------------------------

# GET BASELINE DATA SET/S:
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

# CALIBRATE SAMPLES:
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

#-------------------------------------------------
# I.a. Summed probability distributions (SPDs)
#-------------------------------------------------

runm = 100 ## the running mean to use for the SPDs

spd.baseline.noloss1 <- spd(x = cal.baseline.noloss1, timeRange = timeRange, runm = runm, spdnormalised = normalised)

spd.baseline.noloss  <- compare_spds(calibrated_samples = cal.baseline.noloss, 
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)

spd.uniform.50p      <- compare_spds(calibrated_samples = cal.uniform.50p, 
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)
spd.uniform.75p      <- compare_spds(calibrated_samples = cal.uniform.75p, 
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)
spd.uniform.100p     <- compare_spds(calibrated_samples = cal.uniform.100p, 
                                     timeRange = timeRange, 
                                     runm = runm, 
                                     normalised = normalised)

spd.singleton.ancient.50p  <- compare_spds(calibrated_samples = cal.singleton.ancient.50p, 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.ancient.75p  <- compare_spds(calibrated_samples = cal.singleton.ancient.75p, 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.ancient.100p <- compare_spds(calibrated_samples = cal.singleton.ancient.100p, 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)

spd.singleton.recent.50p   <- compare_spds(calibrated_samples = cal.singleton.recent.50p, 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.recent.75p   <- compare_spds(calibrated_samples = cal.singleton.recent.75p, 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)
spd.singleton.recent.100p  <- compare_spds(calibrated_samples = cal.singleton.recent.100p, 
                                           timeRange = timeRange, 
                                           runm = runm, 
                                           normalised = normalised)

spd.singleton.random.50p  <- compare_spds(calibrated_samples = cal.singleton.random.50p, 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.singleton.random.75p  <- compare_spds(calibrated_samples = cal.singleton.random.75p, 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.singleton.random.100p <- compare_spds(calibrated_samples = cal.singleton.random.100p, 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)

spd.bracketed.50p         <- compare_spds(calibrated_samples = cal.bracketed.50p, 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.bracketed.75p         <- compare_spds(calibrated_samples = cal.bracketed.75p, 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)
spd.bracketed.100p        <- compare_spds(calibrated_samples = cal.bracketed.100p, 
                                          timeRange = timeRange, 
                                          runm = runm, 
                                          normalised = normalised)

# Rearrange data for plotting:
plot.spd.noloss <- data.frame(matrix(NA, nrow = length(spd.baseline.noloss$calBP)*17, ncol = 5))
names(plot.spd.noloss) <- c("sample", "years.ka", "lowerCI", "median", "upperCI")

for (i in 1:length(spd.baseline.noloss$calBP)){
  plot.spd.noloss[i,1] <- "baseline (no loss) replicates"
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

  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),1] <- "baseline (no loss)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),2] <- spd.baseline.noloss1$grid$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),3] <- 0
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),4] <- spd.baseline.noloss1$grid$PrDens[i]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*16+i),5] <- 0
}

write.csv(plot.spd.noloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-SPD_noloss-5samples_20sims.csv")

# Plot (combined):
group.colors <- c("baseline (no loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (no loss)" = "black")


# Plot (multi panelled):
uniform <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  ylim(0, 0.0007) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0007) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
spd.noloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                    labels = c("A", "B", "C", "D", "E"),
                    ncol = 2, nrow = 3)
spd.noloss
ggexport(spd.noloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/No_population_change-SPD_noloss-5samples_20sims.png",
         width = 1500,
         height = 1000)

#------------------------------------------------------------
# I.b. Frequency distributions of dates within time bins
#------------------------------------------------------------

taphCorrect = FALSE    ## whether to taphonomically correct the open sites - not currently as haven't specified whether sites are open
correctForSite = FALSE ## per bin: number of sites with at least one date in them, or number of dates within the bin total
binSize = 100          ## the size of the bins to use (years) - use a bin size that actually fits the timeRange specified

# Get median bin age for x axis (same for all plots)
median.bin.age   <- seq(timeRange[2] + binSize/2, timeRange[1] - binSize/2, by = binSize)

# Get frequency distribution for baseline data set
ev2         <- vector("list", length = 1)
ev2[[1]]    <- baseline.data[[1]][[1]]
cal.ev      <- vector("list", length = 1)
cal.ev[[1]] <- cal.baseline.noloss1
baseline.freq.noloss <- generate_multiple_frequency_dists(data = ev2, calibrated_data = cal.ev, timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, correctForSite = correctForSite,
                                                          binSize = binSize)
standardised.baseline.freq.noloss <- matrix(nrow = nrow(baseline.freq.noloss), ncol = ncol(baseline.freq.noloss))
for (i in 1:ncol(standardised.baseline.freq.noloss)){
  standardised.baseline.freq.noloss[,i] = baseline.freq.noloss[,i]/sum(baseline.freq.noloss[,i])
}

# Get sample frequency distributions
fd.baseline.noloss       <- compare_frequency_dists(samples = baseline.data[[1]],
                                                    calibrated_samples = cal.baseline.noloss,
                                                    timeRange = timeRange,
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite,
                                                    binSize = binSize)


fd.uniform.50p           <- compare_frequency_dists(samples = sample.uniform.50p, 
                                                    calibrated_samples = cal.uniform.50p,
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect, 
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.uniform.75p           <- compare_frequency_dists(samples = sample.uniform.75p, 
                                                    calibrated_samples = cal.uniform.75p,
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect, 
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.uniform.100p          <- compare_frequency_dists(samples = sample.uniform.100p, 
                                                    calibrated_samples = cal.uniform.100p,
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect, 
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)

fd.singleton.ancient.50p  <- compare_frequency_dists(samples = sample.singleton.ancient.50p, 
                                                     calibrated_samples = cal.singleton.ancient.50p,
                                                     timeRange = timeRange, 
                                                     taphCorrect = taphCorrect, 
                                                     correctForSite = correctForSite, 
                                                     binSize = binSize)
fd.singleton.ancient.75p  <- compare_frequency_dists(samples = sample.singleton.ancient.75p, 
                                                     calibrated_samples = cal.singleton.ancient.75p,
                                                     timeRange = timeRange, 
                                                     taphCorrect = taphCorrect, 
                                                     correctForSite = correctForSite, 
                                                     binSize = binSize)
fd.singleton.ancient.100p <- compare_frequency_dists(samples = sample.singleton.ancient.100p, 
                                                     calibrated_samples = cal.singleton.ancient.100p,
                                                     timeRange = timeRange, 
                                                     taphCorrect = taphCorrect, 
                                                     correctForSite = correctForSite, 
                                                     binSize = binSize)

fd.singleton.recent.50p  <- compare_frequency_dists(samples = sample.singleton.recent.50p, 
                                                    calibrated_samples = cal.singleton.recent.50p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.singleton.recent.75p  <- compare_frequency_dists(samples = sample.singleton.recent.75p, 
                                                    calibrated_samples = cal.singleton.recent.75p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.singleton.recent.100p <- compare_frequency_dists(samples = sample.singleton.recent.100p, 
                                                    calibrated_samples = cal.singleton.recent.100p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)

fd.singleton.random.50p  <- compare_frequency_dists(samples = sample.singleton.random.50p, 
                                                    calibrated_samples = cal.singleton.random.50p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.singleton.random.75p  <- compare_frequency_dists(samples = sample.singleton.random.75p, 
                                                    calibrated_samples = cal.singleton.random.75p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.singleton.random.100p <- compare_frequency_dists(samples = sample.singleton.random.100p, 
                                                    calibrated_samples = cal.singleton.random.100p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)

fd.bracketed.50p         <- compare_frequency_dists(samples = sample.bracketed.50p, 
                                                    calibrated_samples = cal.bracketed.50p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.bracketed.75p         <- compare_frequency_dists(samples = sample.bracketed.75p, 
                                                    calibrated_samples = cal.bracketed.75p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)
fd.bracketed.100p        <- compare_frequency_dists(samples = sample.bracketed.100p, 
                                                    calibrated_samples = cal.bracketed.100p, 
                                                    timeRange = timeRange, 
                                                    taphCorrect = taphCorrect,
                                                    correctForSite = correctForSite, 
                                                    binSize = binSize)

plot.fd.samples.noloss <- data.frame(matrix(NA, nrow = length(median.bin.age)*17, ncol = 5))
names(plot.fd.samples.noloss) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[i,1] <- "baseline (no loss) replicates"
  plot.fd.samples.noloss[i,2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[i,3] <- fd.baseline.noloss$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[i,4] <- fd.baseline.noloss$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[i,5] <- fd.baseline.noloss$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)+i),1] <- "singleton ancient (50% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)+i),3] <- fd.singleton.ancient.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)+i),4] <- fd.singleton.ancient.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)+i),5] <- fd.singleton.ancient.50p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*2+i),1] <- "singleton ancient (75% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),3] <- fd.singleton.ancient.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),4] <- fd.singleton.ancient.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),5] <- fd.singleton.ancient.75p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*3+i),1] <- "singleton ancient (100% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),3] <- fd.singleton.ancient.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),4] <- fd.singleton.ancient.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),5] <- fd.singleton.ancient.100p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*4+i),1] <- "singleton recent (50% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),3] <- fd.singleton.recent.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),4] <- fd.singleton.recent.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),5] <- fd.singleton.recent.50p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*5+i),1] <- "singleton recent (75% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),3] <- fd.singleton.recent.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),4] <- fd.singleton.recent.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),5] <- fd.singleton.recent.75p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*6+i),1] <- "singleton recent (100% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),3] <- fd.singleton.recent.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),4] <- fd.singleton.recent.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),5] <- fd.singleton.recent.100p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*7+i),1] <- "singleton random (50% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*7+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*7+i),3] <- fd.singleton.random.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*7+i),4] <- fd.singleton.random.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*7+i),5] <- fd.singleton.random.50p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*8+i),1] <- "singleton random (75% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*8+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*8+i),3] <- fd.singleton.random.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*8+i),4] <- fd.singleton.random.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*8+i),5] <- fd.singleton.random.75p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*9+i),1] <- "singleton random (100% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*9+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*9+i),3] <- fd.singleton.random.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*9+i),4] <- fd.singleton.random.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*9+i),5] <- fd.singleton.random.100p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*10+i),1] <- "bracketed (50% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*10+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*10+i),3] <- fd.bracketed.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*10+i),4] <- fd.bracketed.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*10+i),5] <- fd.bracketed.50p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*11+i),1] <- "bracketed (75% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*11+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*11+i),3] <- fd.bracketed.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*11+i),4] <- fd.bracketed.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*11+i),5] <- fd.bracketed.75p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*12+i),1] <- "bracketed (100% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*12+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*12+i),3] <- fd.bracketed.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*12+i),4] <- fd.bracketed.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*12+i),5] <- fd.bracketed.100p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*13+i),1] <- "uniform (50% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*13+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*13+i),3] <- fd.uniform.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*13+i),4] <- fd.uniform.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*13+i),5] <- fd.uniform.50p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*14+i),1] <- "uniform (75% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*14+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*14+i),3] <- fd.uniform.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*14+i),4] <- fd.uniform.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*14+i),5] <- fd.uniform.75p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*15+i),1] <- "uniform (100% sites)"
  plot.fd.samples.noloss[(length(median.bin.age)*15+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*15+i),3] <- fd.uniform.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*15+i),4] <- fd.uniform.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*15+i),5] <- fd.uniform.100p$standardised.envelope[[1]][i,3]

  plot.fd.samples.noloss[(length(median.bin.age)*16+i),1] <- "baseline (no loss)"
  plot.fd.samples.noloss[(length(median.bin.age)*16+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*16+i),3] <- 0
  plot.fd.samples.noloss[(length(median.bin.age)*16+i),4] <- standardised.baseline.freq.noloss[i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*16+i),5] <- 0
}

write.csv(plot.fd.samples.noloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-FD_samples_noloss-5samples_20sims.csv")

# Combined plot
group.colors <- c("baseline (no loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (no loss)" = "black")

# Plot (multi panelled):
uniform <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample,  "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(samples = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
fd.samples.noloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                     labels = c("A", "B", "C", "D", "E"),
                     ncol = 2, nrow = 3)
fd.samples.noloss
ggexport(fd.samples.noloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/No_population_change-FD_samples_noloss-5samples_20sims.png",
         width = 1500,
         height = 1000)

#---------------------------------------------------------------------
# I.c. Frequency distributions of sites with dates within time bins
#---------------------------------------------------------------------

taphCorrect = FALSE    ## whether to taphonomically correct the open sites - not currently as haven't specified whether sites are open
correctForSite = TRUE  ## per bin: number of sites with at least one date in them, or number of dates within the bin total
binSize = 100          ## the size of the bins to use (years) - use a bin size that actually fits the timeRange specified

# Get median bin age for x axis (same for all plots)
median.bin.age   <- seq(timeRange[2] + binSize/2, timeRange[1] - binSize/2, by = binSize)

# Get frequency distribution for baseline data set
ev2         <- vector("list", length = 1)
ev2[[1]]    <- baseline.data[[1]][[1]]
cal.ev      <- vector("list", length = 1)
cal.ev[[1]] <- cal.baseline.noloss1
baseline.freq.noloss.sites <- generate_multiple_frequency_dists(data = ev2, calibrated_data = cal.ev, timeRange = timeRange, 
                                                                taphCorrect = taphCorrect, correctForSite = correctForSite,
                                                                binSize = binSize)
standardised.baseline.freq.noloss.sites <- matrix(nrow = nrow(baseline.freq.noloss.sites), ncol = ncol(baseline.freq.noloss.sites))
for (i in 1:ncol(standardised.baseline.freq.noloss.sites)){
  standardised.baseline.freq.noloss.sites[,i] = baseline.freq.noloss.sites[,i]/sum(baseline.freq.noloss.sites[,i])
}

fd.baseline.sites          <- compare_frequency_dists(samples = baseline.data[[1]],
                                                      calibrated_samples = cal.baseline.noloss,
                                                      timeRange = timeRange,
                                                      taphCorrect = taphCorrect,
                                                      correctForSite = correctForSite,
                                                      binSize = binSize)

fd.uniform.sites.50p       <- compare_frequency_dists(samples = sample.uniform.50p, 
                                                      calibrated_samples = cal.uniform.50p,
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect, 
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)
fd.uniform.sites.75p       <- compare_frequency_dists(samples = sample.uniform.75p, 
                                                      calibrated_samples = cal.uniform.75p,
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect, 
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)
fd.uniform.sites.100p      <- compare_frequency_dists(samples = sample.uniform.100p, 
                                                      calibrated_samples = cal.uniform.100p,
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect, 
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)

fd.singleton.ancient.sites.50p  <- compare_frequency_dists(samples = sample.singleton.ancient.50p, 
                                                           calibrated_samples = cal.singleton.ancient.50p,
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect, 
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.singleton.ancient.sites.75p  <- compare_frequency_dists(samples = sample.singleton.ancient.75p, 
                                                           calibrated_samples = cal.singleton.ancient.75p,
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect, 
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.singleton.ancient.sites.100p <- compare_frequency_dists(samples = sample.singleton.ancient.100p, 
                                                           calibrated_samples = cal.singleton.ancient.100p,
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect, 
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)

fd.singleton.recent.sites.50p   <- compare_frequency_dists(samples = sample.singleton.recent.50p, 
                                                           calibrated_samples = cal.singleton.recent.50p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.singleton.recent.sites.75p   <- compare_frequency_dists(samples = sample.singleton.recent.75p, 
                                                           calibrated_samples = cal.singleton.recent.75p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.singleton.recent.sites.100p  <- compare_frequency_dists(samples = sample.singleton.recent.100p, 
                                                           calibrated_samples = cal.singleton.recent.100p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)

fd.singleton.random.sites.50p   <- compare_frequency_dists(samples = sample.singleton.random.50p, 
                                                           calibrated_samples = cal.singleton.random.50p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.singleton.random.sites.75p   <- compare_frequency_dists(samples = sample.singleton.random.75p, 
                                                           calibrated_samples = cal.singleton.random.75p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.singleton.random.sites.100p  <- compare_frequency_dists(samples = sample.singleton.random.100p, 
                                                           calibrated_samples = cal.singleton.random.100p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)

fd.bracketed.sites.50p          <- compare_frequency_dists(samples = sample.bracketed.50p, 
                                                           calibrated_samples = cal.bracketed.50p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.bracketed.sites.75p          <- compare_frequency_dists(samples = sample.bracketed.75p, 
                                                           calibrated_samples = cal.bracketed.75p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)
fd.bracketed.sites.100p         <- compare_frequency_dists(samples = sample.bracketed.100p, 
                                                           calibrated_samples = cal.bracketed.100p, 
                                                           timeRange = timeRange, 
                                                           taphCorrect = taphCorrect,
                                                           correctForSite = correctForSite, 
                                                           binSize = binSize)


plot.fd.sites.noloss <- data.frame(matrix(NA, nrow = length(median.bin.age)*17, ncol = 5))
names(plot.fd.sites.noloss) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[i,1] <- "baseline (no loss) replicates"
  plot.fd.sites.noloss[i,2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[i,3] <- fd.baseline.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[i,4] <- fd.baseline.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[i,5] <- fd.baseline.sites$standardised.envelope[[1]][i,3]

  plot.fd.sites.noloss[(length(median.bin.age)+i),1] <- "singleton ancient (50% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)+i),3] <- fd.singleton.ancient.sites.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)+i),4] <- fd.singleton.ancient.sites.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)+i),5] <- fd.singleton.ancient.sites.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),1] <- "singleton ancient (75% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),3] <- fd.singleton.ancient.sites.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),4] <- fd.singleton.ancient.sites.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),5] <- fd.singleton.ancient.sites.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),1] <- "singleton ancient (100% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),3] <- fd.singleton.ancient.sites.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),4] <- fd.singleton.ancient.sites.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),5] <- fd.singleton.ancient.sites.100p$standardised.envelope[[1]][i,3]

  plot.fd.sites.noloss[(length(median.bin.age)*4+i),1] <- "singleton recent (50% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),3] <- fd.singleton.recent.sites.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),4] <- fd.singleton.recent.sites.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),5] <- fd.singleton.recent.sites.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),1] <- "singleton recent (75% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),3] <- fd.singleton.recent.sites.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),4] <- fd.singleton.recent.sites.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),5] <- fd.singleton.recent.sites.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),1] <- "singleton recent (100% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),3] <- fd.singleton.recent.sites.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),4] <- fd.singleton.recent.sites.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),5] <- fd.singleton.recent.sites.100p$standardised.envelope[[1]][i,3]

  plot.fd.sites.noloss[(length(median.bin.age)*7+i),1] <- "singleton random (50% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*7+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*7+i),3] <- fd.singleton.random.sites.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*7+i),4] <- fd.singleton.random.sites.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*7+i),5] <- fd.singleton.random.sites.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*8+i),1] <- "singleton random (75% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*8+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*8+i),3] <- fd.singleton.random.sites.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*8+i),4] <- fd.singleton.random.sites.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*8+i),5] <- fd.singleton.random.sites.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*9+i),1] <- "singleton random (100% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*9+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*9+i),3] <- fd.singleton.random.sites.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*9+i),4] <- fd.singleton.random.sites.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*9+i),5] <- fd.singleton.random.sites.100p$standardised.envelope[[1]][i,3]

  plot.fd.sites.noloss[(length(median.bin.age)*10+i),1] <- "bracketed (50% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*10+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*10+i),3] <- fd.bracketed.sites.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*10+i),4] <- fd.bracketed.sites.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*10+i),5] <- fd.bracketed.sites.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*11+i),1] <- "bracketed (75% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*11+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*11+i),3] <- fd.bracketed.sites.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*11+i),4] <- fd.bracketed.sites.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*11+i),5] <- fd.bracketed.sites.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*12+i),1] <- "bracketed (100% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*12+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*12+i),3] <- fd.bracketed.sites.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*12+i),4] <- fd.bracketed.sites.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*12+i),5] <- fd.bracketed.sites.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*13+i),1] <- "uniform (50% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*13+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*13+i),3] <- fd.uniform.sites.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*13+i),4] <- fd.uniform.sites.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*13+i),5] <- fd.uniform.sites.50p$standardised.envelope[[1]][i,3]

  plot.fd.sites.noloss[(length(median.bin.age)*14+i),1] <- "uniform (75% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*14+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*14+i),3] <- fd.uniform.sites.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*14+i),4] <- fd.uniform.sites.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*14+i),5] <- fd.uniform.sites.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.noloss[(length(median.bin.age)*15+i),1] <- "uniform (100% sites)"
  plot.fd.sites.noloss[(length(median.bin.age)*15+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*15+i),3] <- fd.uniform.sites.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*15+i),4] <- fd.uniform.sites.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*15+i),5] <- fd.uniform.sites.100p$standardised.envelope[[1]][i,3]

  plot.fd.sites.noloss[(length(median.bin.age)*16+i),1] <- "baseline (no loss)"
  plot.fd.sites.noloss[(length(median.bin.age)*16+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*16+i),3] <- 0
  plot.fd.sites.noloss[(length(median.bin.age)*16+i),4] <- standardised.baseline.freq.noloss.sites[i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*16+i),5] <- 0
}

write.csv(plot.fd.sites.noloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-FD_sites_noloss-5samples_20sims.csv")

# Plot (combined)
group.colors <- c("baseline (no loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (no loss)" = "black")

# Plot (multi panelled):
uniform <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.07) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
fd.sites.noloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                             labels = c("A", "B", "C", "D", "E"),
                             ncol = 2, nrow = 3)
fd.sites.noloss
ggexport(fd.sites.noloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/No_population_change-FD_sites_noloss-5samples_20sims.png",
         width = 1500,
         height = 1000)

#-------------------------------------------------------------------------------------
# II. SIMULATING NO CHANGE IN POPULATION OVER TIME (WITH TAPHONOMIC LOSS)
#-------------------------------------------------------------------------------------

# Get samples:
sample.uniform.loss.50p           <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "uniform", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 50, 
                                                 nsim = nsim)
sample.uniform.loss.75p           <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "uniform", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 75, 
                                                 nsim = nsim)
sample.uniform.loss.100p          <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "uniform", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 100, 
                                                 nsim = nsim)

sample.singleton.ancient.loss.50p  <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                  sampling_method = "singleton_ancient", 
                                                  sampling_effort = sampling_effort, 
                                                  percent_sites = 50, 
                                                  nsim = nsim)
sample.singleton.ancient.loss.75p  <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                  sampling_method = "singleton_ancient", 
                                                  sampling_effort = sampling_effort, 
                                                  percent_sites = 75, 
                                                  nsim = nsim)
sample.singleton.ancient.loss.100p <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                  sampling_method = "singleton_ancient", 
                                                  sampling_effort = sampling_effort, 
                                                  percent_sites = 100, 
                                                  nsim = nsim)

sample.singleton.recent.loss.50p  <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "singleton_recent", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 50, 
                                                 nsim = nsim)
sample.singleton.recent.loss.75p  <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "singleton_recent", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 75, 
                                                 nsim = nsim)
sample.singleton.recent.loss.100p <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "singleton_recent", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 100, 
                                                 nsim = nsim)

sample.singleton.random.loss.50p  <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "singleton_random", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 50, 
                                                 nsim = nsim)
sample.singleton.random.loss.75p  <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "singleton_random", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 75, 
                                                 nsim = nsim)
sample.singleton.random.loss.100p <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "singleton_random", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 100, 
                                                 nsim = nsim)

sample.bracketed.loss.50p         <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "bracketed", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 50, 
                                                 nsim = nsim)
sample.bracketed.loss.75p         <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "bracketed", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 75, 
                                                 nsim = nsim)
sample.bracketed.loss.100p        <- get_samples(evidence = baseline.data[[2]][[1]], 
                                                 sampling_method = "bracketed", 
                                                 sampling_effort = sampling_effort, 
                                                 percent_sites = 100, 
                                                 nsim = nsim)

# Calibrate samples:
normalised      = TRUE ## are calibration curves (and, later, SPDs) normalised?
ncores          = 6    ## the number of threads to use when calibrating the radiocarbon dates

cal.baseline.taphloss1 <- rcarbon::calibrate(x = baseline.data[[2]][[1]]$age, 
                                             errors = baseline.data[[2]][[1]]$error, 
                                             calCurves = 'shcal20', 
                                             normalised = normalised)

cal.baseline.taphloss  <- calibrate_samples(samples = baseline.data[[2]], 
                                            normalised = normalised, 
                                            ncores = ncores)

cal.uniform.loss.50p           <- calibrate_samples(samples = sample.uniform.loss.50p, 
                                                    normalised = normalised, 
                                                    ncores = ncores)
cal.uniform.loss.75p           <- calibrate_samples(samples = sample.uniform.loss.75p, 
                                                    normalised = normalised, 
                                                    ncores = ncores)
cal.uniform.loss.100p          <- calibrate_samples(samples = sample.uniform.loss.100p, 
                                                    normalised = normalised, 
                                                    ncores = ncores)

cal.singleton.ancient.loss.50p  <- calibrate_samples(samples = sample.singleton.ancient.loss.50p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.singleton.ancient.loss.75p  <- calibrate_samples(samples = sample.singleton.ancient.loss.75p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.singleton.ancient.loss.100p <- calibrate_samples(samples = sample.singleton.ancient.loss.100p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)

cal.singleton.recent.loss.50p   <- calibrate_samples(samples = sample.singleton.recent.loss.50p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.singleton.recent.loss.75p   <- calibrate_samples(samples = sample.singleton.recent.loss.75p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.singleton.recent.loss.100p  <- calibrate_samples(samples = sample.singleton.recent.loss.100p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)

cal.singleton.random.loss.50p   <- calibrate_samples(samples = sample.singleton.random.loss.50p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.singleton.random.loss.75p   <- calibrate_samples(samples = sample.singleton.random.loss.75p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.singleton.random.loss.100p  <- calibrate_samples(samples = sample.singleton.random.loss.100p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)

cal.bracketed.loss.50p          <- calibrate_samples(samples = sample.bracketed.loss.50p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.bracketed.loss.75p          <- calibrate_samples(samples = sample.bracketed.loss.75p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)
cal.bracketed.loss.100p         <- calibrate_samples(samples = sample.bracketed.loss.100p, 
                                                     normalised = normalised, 
                                                     ncores = ncores)

#-----------------------------------------
# II.a. Summed probability distributions
#-----------------------------------------

runm = 100 ## the running mean to use for the SPDs

spd.baseline.taphloss1   <- spd(x = cal.baseline.taphloss1, timeRange = timeRange, runm = runm, spdnormalised = normalised)

spd.baseline.taphloss    <- compare_spds(calibrated_samples = cal.baseline.taphloss, 
                                         timeRange = timeRange, 
                                         runm = runm, 
                                         normalised = normalised)

spd.uniform.loss.50p           <- compare_spds(calibrated_samples = cal.uniform.loss.50p, 
                                               timeRange = timeRange, 
                                               runm = runm, 
                                               normalised = normalised)
spd.uniform.loss.75p           <- compare_spds(calibrated_samples = cal.uniform.loss.75p, 
                                               timeRange = timeRange, 
                                               runm = runm, 
                                               normalised = normalised)
spd.uniform.loss.100p          <- compare_spds(calibrated_samples = cal.uniform.loss.100p, 
                                               timeRange = timeRange, 
                                               runm = runm, 
                                               normalised = normalised)

spd.singleton.ancient.loss.50p  <- compare_spds(calibrated_samples = cal.singleton.ancient.loss.50p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.singleton.ancient.loss.75p  <- compare_spds(calibrated_samples = cal.singleton.ancient.loss.75p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.singleton.ancient.loss.100p <- compare_spds(calibrated_samples = cal.singleton.ancient.loss.100p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)

spd.singleton.recent.loss.50p   <- compare_spds(calibrated_samples = cal.singleton.recent.loss.50p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.singleton.recent.loss.75p   <- compare_spds(calibrated_samples = cal.singleton.recent.loss.75p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.singleton.recent.loss.100p  <- compare_spds(calibrated_samples = cal.singleton.recent.loss.100p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)

spd.singleton.random.loss.50p   <- compare_spds(calibrated_samples = cal.singleton.random.loss.50p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.singleton.random.loss.75p   <- compare_spds(calibrated_samples = cal.singleton.random.loss.75p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.singleton.random.loss.100p  <- compare_spds(calibrated_samples = cal.singleton.random.loss.100p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)

spd.bracketed.loss.50p          <- compare_spds(calibrated_samples = cal.bracketed.loss.50p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.bracketed.loss.75p          <- compare_spds(calibrated_samples = cal.bracketed.loss.75p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)
spd.bracketed.loss.100p         <- compare_spds(calibrated_samples = cal.bracketed.loss.100p, 
                                                timeRange = timeRange, 
                                                runm = runm, 
                                                normalised = normalised)

# Rearrange data for plotting:
plot.spd.taphloss <- data.frame(matrix(NA, nrow = length(spd.baseline.taphloss$calBP)*18, ncol = 5))
names(plot.spd.taphloss) <- c("sample", "years.ka", "lowerCI", "median", "upperCI")

for (i in 1:length(spd.baseline.taphloss$calBP)){
  plot.spd.taphloss[i,1] <- "baseline (taphonomic loss) replicates"
  plot.spd.taphloss[i,2] <- spd.baseline.taphloss$calBP[i]/1000
  plot.spd.taphloss[i,3] <- spd.baseline.taphloss$envelope[[1]][i,1]
  plot.spd.taphloss[i,4] <- spd.baseline.taphloss$envelope[[1]][i,2]
  plot.spd.taphloss[i,5] <- spd.baseline.taphloss$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)+i),1] <- "singleton ancient (50% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)+i),2] <- spd.singleton.ancient.loss.50p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)+i),3] <- spd.singleton.ancient.loss.50p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)+i),4] <- spd.singleton.ancient.loss.50p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)+i),5] <- spd.singleton.ancient.loss.50p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*2+i),1] <- "singleton ancient (75% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*2+i),2] <- spd.singleton.ancient.loss.75p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*2+i),3] <- spd.singleton.ancient.loss.75p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*2+i),4] <- spd.singleton.ancient.loss.75p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*2+i),5] <- spd.singleton.ancient.loss.75p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*3+i),1] <- "singleton ancient (100% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*3+i),2] <- spd.singleton.ancient.loss.100p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*3+i),3] <- spd.singleton.ancient.loss.100p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*3+i),4] <- spd.singleton.ancient.loss.100p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*3+i),5] <- spd.singleton.ancient.loss.100p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*4+i),1] <- "singleton recent (50% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*4+i),2] <- spd.singleton.recent.loss.50p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*4+i),3] <- spd.singleton.recent.loss.50p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*4+i),4] <- spd.singleton.recent.loss.50p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*4+i),5] <- spd.singleton.recent.loss.50p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*5+i),1] <- "singleton recent (75% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*5+i),2] <- spd.singleton.recent.loss.75p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*5+i),3] <- spd.singleton.recent.loss.75p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*5+i),4] <- spd.singleton.recent.loss.75p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*5+i),5] <- spd.singleton.recent.loss.75p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*6+i),1] <- "singleton recent (100% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*6+i),2] <- spd.singleton.recent.loss.100p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*6+i),3] <- spd.singleton.recent.loss.100p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*6+i),4] <- spd.singleton.recent.loss.100p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*6+i),5] <- spd.singleton.recent.loss.100p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*7+i),1] <- "singleton random (50% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*7+i),2] <- spd.singleton.random.loss.50p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*7+i),3] <- spd.singleton.random.loss.50p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*7+i),4] <- spd.singleton.random.loss.50p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*7+i),5] <- spd.singleton.random.loss.50p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*8+i),1] <- "singleton random (75% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*8+i),2] <- spd.singleton.random.loss.75p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*8+i),3] <- spd.singleton.random.loss.75p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*8+i),4] <- spd.singleton.random.loss.75p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*8+i),5] <- spd.singleton.random.loss.75p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*9+i),1] <- "singleton random (100% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*9+i),2] <- spd.singleton.random.loss.100p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*9+i),3] <- spd.singleton.random.loss.100p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*9+i),4] <- spd.singleton.random.loss.100p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*9+i),5] <- spd.singleton.random.loss.100p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*10+i),1] <- "bracketed (50% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*10+i),2] <- spd.bracketed.loss.50p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*10+i),3] <- spd.bracketed.loss.50p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*10+i),4] <- spd.bracketed.loss.50p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*10+i),5] <- spd.bracketed.loss.50p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*11+i),1] <- "bracketed (75% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*11+i),2] <- spd.bracketed.loss.75p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*11+i),3] <- spd.bracketed.loss.75p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*11+i),4] <- spd.bracketed.loss.75p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*11+i),5] <- spd.bracketed.loss.75p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*12+i),1] <- "bracketed (100% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*12+i),2] <- spd.bracketed.loss.100p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*12+i),3] <- spd.bracketed.loss.100p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*12+i),4] <- spd.bracketed.loss.100p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*12+i),5] <- spd.bracketed.loss.100p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*13+i),1] <- "uniform (50% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*13+i),2] <- spd.uniform.loss.50p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*13+i),3] <- spd.uniform.loss.50p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*13+i),4] <- spd.uniform.loss.50p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*13+i),5] <- spd.uniform.loss.50p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*14+i),1] <- "uniform (75% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*14+i),2] <- spd.uniform.loss.75p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*14+i),3] <- spd.uniform.loss.75p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*14+i),4] <- spd.uniform.loss.75p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*14+i),5] <- spd.uniform.loss.75p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*15+i),1] <- "uniform (100% sites)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*15+i),2] <- spd.uniform.loss.100p$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*15+i),3] <- spd.uniform.loss.100p$envelope[[1]][i,1]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*15+i),4] <- spd.uniform.loss.100p$envelope[[1]][i,2]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*15+i),5] <- spd.uniform.loss.100p$envelope[[1]][i,3]
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*16+i),1] <- "baseline (no loss)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*16+i),2] <- spd.baseline.noloss1$grid$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*16+i),3] <- 0
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*16+i),4] <- spd.baseline.noloss1$grid$PrDens[i]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*16+i),5] <- 0
  
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*17+i),1] <- "baseline (taphonomic loss)"
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*17+i),2] <- spd.baseline.taphloss1$grid$calBP[i]/1000
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*17+i),3] <- 0
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*17+i),4] <- spd.baseline.taphloss1$grid$PrDens[i]
  plot.spd.taphloss[(length(spd.baseline.taphloss$calBP)*17+i),5] <- 0
}

write.csv(plot.spd.taphloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-SPD_taphloss-5samples_20sims.csv")

# Plot (multi panelled):
group.colors <- c("baseline (taphonomic loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (taphonomic loss)"    = "black", 
                  "baseline (no loss)"            = "grey48")

uniform <- plot.spd.taphloss[plot.spd.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.00085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.spd.taphloss[plot.spd.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  ylim(0, 0.00085) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.spd.taphloss[plot.spd.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.00085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.spd.taphloss[plot.spd.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.00085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.spd.taphloss[plot.spd.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.00085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.spd.taphloss[plot.spd.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "baseline (taphonomic loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.00085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
spd.taphloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                        labels = c("A", "B", "C", "D", "E"),
                        ncol = 2, nrow = 3)
spd.taphloss
ggexport(spd.taphloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/No_population_change-SPD_taphloss-5samples_20sims.png",
         width = 1500,
         height = 1000)

#-------------------------------------------------------------
# II.a. Frequency distribution of dates in bins (corrected)
#-------------------------------------------------------------

taphCorrect = TRUE     ## whether to taphonomically correct the open sites - not currently as haven't specified whether sites are open
correctForSite = FALSE ## per bin: number of sites with at least one date in them, or number of dates within the bin total
binSize = 100          ## the size of the bins to use (years) - use a bin size that actually fits the timeRange specified

# Get median bin age for x axis (same for all plots)
median.bin.age   <- seq(timeRange[2] + binSize/2, timeRange[1] - binSize/2, by = binSize)

# Get frequency distribution for baseline data set
ev2         <- vector("list", length = 1)
ev2[[1]]    <- baseline.data[[2]][[1]]
cal.ev      <- vector("list", length = 1)
cal.ev[[1]] <- cal.baseline.taphloss1
baseline.freq.taphloss.correct <- generate_multiple_frequency_dists(data = ev2, calibrated_data = cal.ev, timeRange = timeRange, 
                                                                    taphCorrect = taphCorrect, correctForSite = correctForSite,
                                                                    binSize = binSize)
standardised.baseline.freq.taphloss.correct <- matrix(nrow = nrow(baseline.freq.taphloss.correct), 
                                                      ncol = ncol(baseline.freq.taphloss.correct))
for (i in 1:ncol(standardised.baseline.freq.taphloss.correct)){
  standardised.baseline.freq.taphloss.correct[,i] = baseline.freq.taphloss.correct[,i]/sum(baseline.freq.taphloss.correct[,i])
}

# Get sample frequency distributions
fd.baseline.taphloss.correct   <- compare_frequency_dists(samples = baseline.data[[2]],
                                                          calibrated_samples = cal.baseline.taphloss,
                                                          timeRange = timeRange,
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite,
                                                          binSize = binSize)

fd.uniform.loss.correct.50p           <- compare_frequency_dists(samples = sample.uniform.loss.50p, 
                                                         calibrated_samples = cal.uniform.loss.50p,
                                                         timeRange = timeRange, 
                                                         taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, 
                                                         binSize = binSize)
fd.uniform.loss.correct.75p           <- compare_frequency_dists(samples = sample.uniform.loss.75p, 
                                                         calibrated_samples = cal.uniform.loss.75p,
                                                         timeRange = timeRange, 
                                                         taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, 
                                                         binSize = binSize)
fd.uniform.loss.correct.100p          <- compare_frequency_dists(samples = sample.uniform.loss.100p, 
                                                         calibrated_samples = cal.uniform.loss.100p,
                                                         timeRange = timeRange, 
                                                         taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, 
                                                         binSize = binSize)

fd.singleton.ancient.loss.correct.50p  <- compare_frequency_dists(samples = sample.singleton.ancient.loss.50p, 
                                                          calibrated_samples = cal.singleton.ancient.loss.50p,
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, 
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.ancient.loss.correct.75p  <- compare_frequency_dists(samples = sample.singleton.ancient.loss.75p, 
                                                          calibrated_samples = cal.singleton.ancient.loss.75p,
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, 
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.ancient.loss.correct.100p <- compare_frequency_dists(samples = sample.singleton.ancient.loss.100p, 
                                                          calibrated_samples = cal.singleton.ancient.loss.100p,
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, 
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

fd.singleton.recent.loss.correct.50p   <- compare_frequency_dists(samples = sample.singleton.recent.loss.50p, 
                                                          calibrated_samples = cal.singleton.recent.loss.50p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.recent.loss.correct.75p   <- compare_frequency_dists(samples = sample.singleton.recent.loss.75p, 
                                                          calibrated_samples = cal.singleton.recent.loss.75p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.recent.loss.correct.100p  <- compare_frequency_dists(samples = sample.singleton.recent.loss.100p, 
                                                          calibrated_samples = cal.singleton.recent.loss.100p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

fd.singleton.random.loss.correct.50p   <- compare_frequency_dists(samples = sample.singleton.random.loss.50p, 
                                                          calibrated_samples = cal.singleton.random.loss.50p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.random.loss.correct.75p   <- compare_frequency_dists(samples = sample.singleton.random.loss.75p, 
                                                          calibrated_samples = cal.singleton.random.loss.75p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.random.loss.correct.100p  <- compare_frequency_dists(samples = sample.singleton.random.loss.100p, 
                                                          calibrated_samples = cal.singleton.random.loss.100p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

fd.bracketed.loss.correct.50p          <- compare_frequency_dists(samples = sample.bracketed.loss.50p, 
                                                          calibrated_samples = cal.bracketed.loss.50p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.bracketed.loss.correct.75p          <- compare_frequency_dists(samples = sample.bracketed.loss.75p, 
                                                          calibrated_samples = cal.bracketed.loss.75p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.bracketed.loss.correct.100p         <- compare_frequency_dists(samples = sample.bracketed.loss.100p, 
                                                          calibrated_samples = cal.bracketed.loss.100p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

# Organise data for plotting
plot.fd.samples.taphloss.correct <- data.frame(matrix(NA, nrow = length(median.bin.age)*18, ncol = 5))
names(plot.fd.samples.taphloss.correct) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.samples.taphloss.correct[i,1] <- "baseline (taphonomic loss) replicates"
  plot.fd.samples.taphloss.correct[i,2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[i,3] <- fd.baseline.taphloss.correct$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[i,4] <- fd.baseline.taphloss.correct$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[i,5] <- fd.baseline.taphloss.correct$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)+i),1] <- "singleton ancient (50% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)+i),3] <- fd.singleton.ancient.loss.correct.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)+i),4] <- fd.singleton.ancient.loss.correct.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)+i),5] <- fd.singleton.ancient.loss.correct.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*2+i),1] <- "singleton ancient (75% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*2+i),3] <- fd.singleton.ancient.loss.correct.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*2+i),4] <- fd.singleton.ancient.loss.correct.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*2+i),5] <- fd.singleton.ancient.loss.correct.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*3+i),1] <- "singleton ancient (100% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*3+i),3] <- fd.singleton.ancient.loss.correct.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*3+i),4] <- fd.singleton.ancient.loss.correct.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*3+i),5] <- fd.singleton.ancient.loss.correct.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*4+i),1] <- "singleton recent (50% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*4+i),3] <- fd.singleton.recent.loss.correct.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*4+i),4] <- fd.singleton.recent.loss.correct.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*4+i),5] <- fd.singleton.recent.loss.correct.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*5+i),1] <- "singleton recent (75% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*5+i),3] <- fd.singleton.recent.loss.correct.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*5+i),4] <- fd.singleton.recent.loss.correct.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*5+i),5] <- fd.singleton.recent.loss.correct.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*6+i),1] <- "singleton recent (100% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*6+i),3] <- fd.singleton.recent.loss.correct.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*6+i),4] <- fd.singleton.recent.loss.correct.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*6+i),5] <- fd.singleton.recent.loss.correct.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*7+i),1] <- "singleton random (50% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*7+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*7+i),3] <- fd.singleton.random.loss.correct.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*7+i),4] <- fd.singleton.random.loss.correct.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*7+i),5] <- fd.singleton.random.loss.correct.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*8+i),1] <- "singleton random (75% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*8+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*8+i),3] <- fd.singleton.random.loss.correct.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*8+i),4] <- fd.singleton.random.loss.correct.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*8+i),5] <- fd.singleton.random.loss.correct.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*9+i),1] <- "singleton random (100% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*9+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*9+i),3] <- fd.singleton.random.loss.correct.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*9+i),4] <- fd.singleton.random.loss.correct.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*9+i),5] <- fd.singleton.random.loss.correct.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*10+i),1] <- "bracketed (50% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*10+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*10+i),3] <- fd.bracketed.loss.correct.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*10+i),4] <- fd.bracketed.loss.correct.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*10+i),5] <- fd.bracketed.loss.correct.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*11+i),1] <- "bracketed (75% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*11+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*11+i),3] <- fd.bracketed.loss.correct.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*11+i),4] <- fd.bracketed.loss.correct.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*11+i),5] <- fd.bracketed.loss.correct.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*12+i),1] <- "bracketed (100% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*12+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*12+i),3] <- fd.bracketed.loss.correct.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*12+i),4] <- fd.bracketed.loss.correct.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*12+i),5] <- fd.bracketed.loss.correct.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*13+i),1] <- "uniform (50% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*13+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*13+i),3] <- fd.uniform.loss.correct.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*13+i),4] <- fd.uniform.loss.correct.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*13+i),5] <- fd.uniform.loss.correct.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*14+i),1] <- "uniform (75% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*14+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*14+i),3] <- fd.uniform.loss.correct.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*14+i),4] <- fd.uniform.loss.correct.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*14+i),5] <- fd.uniform.loss.correct.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*15+i),1] <- "uniform (100% sites)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*15+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*15+i),3] <- fd.uniform.loss.correct.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*15+i),4] <- fd.uniform.loss.correct.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*15+i),5] <- fd.uniform.loss.correct.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*16+i),1] <- "baseline (no loss)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*16+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*16+i),3] <- 0
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*16+i),4] <- standardised.baseline.freq.noloss[i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*16+i),5] <- 0
  
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*17+i),1] <- "baseline (taphonomic loss)"
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*17+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*17+i),3] <- 0
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*17+i),4] <- standardised.baseline.freq.taphloss.correct[i,1]
  plot.fd.samples.taphloss.correct[(length(median.bin.age)*17+i),5] <- 0
}

write.csv(plot.fd.samples.taphloss.correct, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-FD_samples_taphloss_corrected-5samples_20sims.csv")

# Plot (multi-panelled)
group.colors <- c("baseline (taphonomic loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (taphonomic loss)"    = "black", 
                  "baseline (no loss)"            = "grey48")


uniform <- plot.fd.samples.taphloss.correct[plot.fd.samples.taphloss.correct$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.08) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.samples.taphloss.correct[plot.fd.samples.taphloss.correct$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.08) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.samples.taphloss.correct[plot.fd.samples.taphloss.correct$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.08) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.samples.taphloss.correct[plot.fd.samples.taphloss.correct$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.08) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.samples.taphloss.correct[plot.fd.samples.taphloss.correct$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample,  "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.08) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.fd.samples.taphloss.correct[plot.fd.samples.taphloss.correct$sample %in% c("baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)"),] %>%
  mutate(samples = fct_relevel(sample, "baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.08) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
fd.samples.taphloss.correct <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                                         labels = c("A", "B", "C", "D", "E"),
                                         ncol = 2, nrow = 3)
fd.samples.taphloss.correct
ggexport(fd.samples.taphloss.correct, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/No_population_change-FD_samples_taphloss_corrected-5samples_20sims.png",
         width = 1500,
         height = 1000)

#-------------------------------------------------------------
# II.b. Frequency distribution of dates in bins (uncorrected)
#-------------------------------------------------------------

taphCorrect = FALSE    ## whether to taphonomically correct the open sites - not currently as haven't specified whether sites are open
correctForSite = FALSE ## per bin: number of sites with at least one date in them, or number of dates within the bin total
binSize = 100          ## the size of the bins to use (years) - use a bin size that actually fits the timeRange specified

# Get median bin age for x axis (same for all plots)
median.bin.age   <- seq(timeRange[2] + binSize/2, timeRange[1] - binSize/2, by = binSize)

# Get frequency distribution for baseline data set
ev2         <- vector("list", length = 1)
ev2[[1]]    <- baseline.data[[2]][[1]]
cal.ev      <- vector("list", length = 1)
cal.ev[[1]] <- cal.baseline.taphloss1
baseline.freq.taphloss.uncorrect <- generate_multiple_frequency_dists(data = ev2, calibrated_data = cal.ev, timeRange = timeRange, 
                                                                      taphCorrect = taphCorrect, correctForSite = correctForSite,
                                                                      binSize = binSize)
standardised.baseline.freq.taphloss.uncorrect <- matrix(nrow = nrow(baseline.freq.taphloss.uncorrect), 
                                                        ncol = ncol(baseline.freq.taphloss.uncorrect))
for (i in 1:ncol(standardised.baseline.freq.taphloss.uncorrect)){
  standardised.baseline.freq.taphloss.uncorrect[,i] = baseline.freq.taphloss.uncorrect[,i]/sum(baseline.freq.taphloss.uncorrect[,i])
}

# Get sample frequency distributions
fd.baseline.taphloss.uncorrect   <- compare_frequency_dists(samples = baseline.data[[2]],
                                                            calibrated_samples = cal.baseline.taphloss,
                                                            timeRange = timeRange,
                                                            taphCorrect = taphCorrect,
                                                            correctForSite = correctForSite,
                                                            binSize = binSize)

fd.uniform.loss.uncorrect.50p           <- compare_frequency_dists(samples = sample.uniform.loss.50p, 
                                                                   calibrated_samples = cal.uniform.loss.50p,
                                                                   timeRange = timeRange, 
                                                                   taphCorrect = taphCorrect, 
                                                                   correctForSite = correctForSite, 
                                                                   binSize = binSize)
fd.uniform.loss.uncorrect.75p           <- compare_frequency_dists(samples = sample.uniform.loss.75p, 
                                                                   calibrated_samples = cal.uniform.loss.75p,
                                                                   timeRange = timeRange, 
                                                                   taphCorrect = taphCorrect, 
                                                                   correctForSite = correctForSite, 
                                                                   binSize = binSize)
fd.uniform.loss.uncorrect.100p          <- compare_frequency_dists(samples = sample.uniform.loss.100p, 
                                                                 calibrated_samples = cal.uniform.loss.100p,
                                                                 timeRange = timeRange, 
                                                                 taphCorrect = taphCorrect, 
                                                                 correctForSite = correctForSite, 
                                                                 binSize = binSize)

fd.singleton.ancient.loss.uncorrect.50p  <- compare_frequency_dists(samples = sample.singleton.ancient.loss.50p, 
                                                                  calibrated_samples = cal.singleton.ancient.loss.50p,
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect, 
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.singleton.ancient.loss.uncorrect.75p  <- compare_frequency_dists(samples = sample.singleton.ancient.loss.75p, 
                                                                  calibrated_samples = cal.singleton.ancient.loss.75p,
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect, 
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.singleton.ancient.loss.uncorrect.100p <- compare_frequency_dists(samples = sample.singleton.ancient.loss.100p, 
                                                                  calibrated_samples = cal.singleton.ancient.loss.100p,
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect, 
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)

fd.singleton.recent.loss.uncorrect.50p   <- compare_frequency_dists(samples = sample.singleton.recent.loss.50p, 
                                                                  calibrated_samples = cal.singleton.recent.loss.50p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.singleton.recent.loss.uncorrect.75p   <- compare_frequency_dists(samples = sample.singleton.recent.loss.75p, 
                                                                  calibrated_samples = cal.singleton.recent.loss.75p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.singleton.recent.loss.uncorrect.100p  <- compare_frequency_dists(samples = sample.singleton.recent.loss.100p, 
                                                                  calibrated_samples = cal.singleton.recent.loss.100p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)

fd.singleton.random.loss.uncorrect.50p   <- compare_frequency_dists(samples = sample.singleton.random.loss.50p, 
                                                                  calibrated_samples = cal.singleton.random.loss.50p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.singleton.random.loss.uncorrect.75p   <- compare_frequency_dists(samples = sample.singleton.random.loss.75p, 
                                                                  calibrated_samples = cal.singleton.random.loss.75p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.singleton.random.loss.uncorrect.100p  <- compare_frequency_dists(samples = sample.singleton.random.loss.100p, 
                                                                  calibrated_samples = cal.singleton.random.loss.100p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)

fd.bracketed.loss.uncorrect.50p          <- compare_frequency_dists(samples = sample.bracketed.loss.50p, 
                                                                  calibrated_samples = cal.bracketed.loss.50p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.bracketed.loss.uncorrect.75p          <- compare_frequency_dists(samples = sample.bracketed.loss.75p, 
                                                                  calibrated_samples = cal.bracketed.loss.75p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)
fd.bracketed.loss.uncorrect.100p         <- compare_frequency_dists(samples = sample.bracketed.loss.100p, 
                                                                  calibrated_samples = cal.bracketed.loss.100p, 
                                                                  timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect,
                                                                  correctForSite = correctForSite, 
                                                                  binSize = binSize)

# Organise data for plotting
plot.fd.samples.taphloss.uncorrect <- data.frame(matrix(NA, nrow = length(median.bin.age)*18, ncol = 5))
names(plot.fd.samples.taphloss.uncorrect) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.samples.taphloss.uncorrect[i,1] <- "baseline (taphonomic loss) replicates"
  plot.fd.samples.taphloss.uncorrect[i,2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[i,3] <- fd.baseline.taphloss.uncorrect$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[i,4] <- fd.baseline.taphloss.uncorrect$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[i,5] <- fd.baseline.taphloss.uncorrect$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)+i),1] <- "singleton ancient (50% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)+i),3] <- fd.singleton.ancient.loss.uncorrect.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)+i),4] <- fd.singleton.ancient.loss.uncorrect.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)+i),5] <- fd.singleton.ancient.loss.uncorrect.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*2+i),1] <- "singleton ancient (75% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*2+i),3] <- fd.singleton.ancient.loss.uncorrect.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*2+i),4] <- fd.singleton.ancient.loss.uncorrect.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*2+i),5] <- fd.singleton.ancient.loss.uncorrect.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*3+i),1] <- "singleton ancient (100% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*3+i),3] <- fd.singleton.ancient.loss.uncorrect.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*3+i),4] <- fd.singleton.ancient.loss.uncorrect.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*3+i),5] <- fd.singleton.ancient.loss.uncorrect.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*4+i),1] <- "singleton recent (50% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*4+i),3] <- fd.singleton.recent.loss.uncorrect.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*4+i),4] <- fd.singleton.recent.loss.uncorrect.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*4+i),5] <- fd.singleton.recent.loss.uncorrect.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*5+i),1] <- "singleton recent (75% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*5+i),3] <- fd.singleton.recent.loss.uncorrect.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*5+i),4] <- fd.singleton.recent.loss.uncorrect.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*5+i),5] <- fd.singleton.recent.loss.uncorrect.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*6+i),1] <- "singleton recent (100% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*6+i),3] <- fd.singleton.recent.loss.uncorrect.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*6+i),4] <- fd.singleton.recent.loss.uncorrect.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*6+i),5] <- fd.singleton.recent.loss.uncorrect.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*7+i),1] <- "singleton random (50% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*7+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*7+i),3] <- fd.singleton.random.loss.uncorrect.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*7+i),4] <- fd.singleton.random.loss.uncorrect.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*7+i),5] <- fd.singleton.random.loss.uncorrect.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*8+i),1] <- "singleton random (75% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*8+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*8+i),3] <- fd.singleton.random.loss.uncorrect.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*8+i),4] <- fd.singleton.random.loss.uncorrect.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*8+i),5] <- fd.singleton.random.loss.uncorrect.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*9+i),1] <- "singleton random (100% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*9+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*9+i),3] <- fd.singleton.random.loss.uncorrect.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*9+i),4] <- fd.singleton.random.loss.uncorrect.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*9+i),5] <- fd.singleton.random.loss.uncorrect.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*10+i),1] <- "bracketed (50% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*10+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*10+i),3] <- fd.bracketed.loss.uncorrect.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*10+i),4] <- fd.bracketed.loss.uncorrect.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*10+i),5] <- fd.bracketed.loss.uncorrect.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*11+i),1] <- "bracketed (75% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*11+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*11+i),3] <- fd.bracketed.loss.uncorrect.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*11+i),4] <- fd.bracketed.loss.uncorrect.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*11+i),5] <- fd.bracketed.loss.uncorrect.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*12+i),1] <- "bracketed (100% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*12+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*12+i),3] <- fd.bracketed.loss.uncorrect.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*12+i),4] <- fd.bracketed.loss.uncorrect.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*12+i),5] <- fd.bracketed.loss.uncorrect.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*13+i),1] <- "uniform (50% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*13+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*13+i),3] <- fd.uniform.loss.uncorrect.50p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*13+i),4] <- fd.uniform.loss.uncorrect.50p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*13+i),5] <- fd.uniform.loss.uncorrect.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*14+i),1] <- "uniform (75% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*14+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*14+i),3] <- fd.uniform.loss.uncorrect.75p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*14+i),4] <- fd.uniform.loss.uncorrect.75p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*14+i),5] <- fd.uniform.loss.uncorrect.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*15+i),1] <- "uniform (100% sites)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*15+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*15+i),3] <- fd.uniform.loss.uncorrect.100p$standardised.envelope[[1]][i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*15+i),4] <- fd.uniform.loss.uncorrect.100p$standardised.envelope[[1]][i,2]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*15+i),5] <- fd.uniform.loss.uncorrect.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*16+i),1] <- "baseline (no loss)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*16+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*16+i),3] <- 0
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*16+i),4] <- standardised.baseline.freq.noloss[i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*16+i),5] <- 0
  
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*17+i),1] <- "baseline (taphonomic loss)"
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*17+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*17+i),3] <- 0
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*17+i),4] <- standardised.baseline.freq.taphloss.uncorrect[i,1]
  plot.fd.samples.taphloss.uncorrect[(length(median.bin.age)*17+i),5] <- 0
}

write.csv(plot.fd.samples.taphloss.uncorrect, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-FD_samples_taphloss_uncorrected-5samples_20sims.csv")

# Plot (multi-panelled)
group.colors <- c("baseline (taphonomic loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (taphonomic loss)"    = "black", 
                  "baseline (no loss)"            = "grey48")


uniform <- plot.fd.samples.taphloss.uncorrect[plot.fd.samples.taphloss.uncorrect$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.095) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.samples.taphloss.uncorrect[plot.fd.samples.taphloss.uncorrect$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.095) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.samples.taphloss.uncorrect[plot.fd.samples.taphloss.uncorrect$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.095) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.samples.taphloss.uncorrect[plot.fd.samples.taphloss.uncorrect$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.095) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.samples.taphloss.uncorrect[plot.fd.samples.taphloss.uncorrect$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample,  "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.095) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.fd.samples.taphloss.uncorrect[plot.fd.samples.taphloss.uncorrect$sample %in% c("baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)"),] %>%
  mutate(samples = fct_relevel(sample, "baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  ylim(0, 0.095) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
fd.samples.taphloss.uncorrect <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                                         labels = c("A", "B", "C", "D", "E"),
                                         ncol = 2, nrow = 3)
fd.samples.taphloss.uncorrect
ggexport(fd.samples.taphloss.uncorrect, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figure/No_population_change-FD_samples_taphloss_uncorrected-5samples_20sims.png",
         width = 1500,
         height = 1000)


#-------------------------------------------------------------
# II.c. Frequency distribution of sites in bins
#-------------------------------------------------------------

taphCorrect = FALSE    ## whether to taphonomically correct the open sites - not currently as haven't specified whether sites are open
correctForSite = TRUE ## per bin: number of sites with at least one date in them, or number of dates within the bin total
binSize = 100          ## the size of the bins to use (years) - use a bin size that actually fits the timeRange specified

# Get frequency distribution for baseline data set
ev2         <- vector("list", length = 1)
ev2[[1]]    <- baseline.data[[2]][[1]]
cal.ev      <- vector("list", length = 1)
cal.ev[[1]] <- cal.baseline.taphloss1
baseline.freq.taphloss.sites <- generate_multiple_frequency_dists(data = ev2, calibrated_data = cal.ev, timeRange = timeRange, 
                                                                  taphCorrect = taphCorrect, correctForSite = correctForSite,
                                                                  binSize = binSize)
standardised.baseline.freq.taphloss.sites <- matrix(nrow = nrow(baseline.freq.taphloss.sites), ncol = ncol(baseline.freq.taphloss.sites))
for (i in 1:ncol(standardised.baseline.freq.taphloss.sites)){
  standardised.baseline.freq.taphloss.sites[,i] = baseline.freq.taphloss.sites[,i]/sum(baseline.freq.taphloss.sites[,i])
}

# Get sample frequency distributions
fd.baseline.taphloss   <- compare_frequency_dists(samples = baseline.data[[2]],
                                                  calibrated_samples = cal.baseline.taphloss,
                                                  timeRange = timeRange,
                                                  taphCorrect = taphCorrect,
                                                  correctForSite = correctForSite,
                                                  binSize = binSize)

fd.uniform.loss.50p           <- compare_frequency_dists(samples = sample.uniform.loss.50p, 
                                                         calibrated_samples = cal.uniform.loss.50p,
                                                         timeRange = timeRange, 
                                                         taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, 
                                                         binSize = binSize)
fd.uniform.loss.75p           <- compare_frequency_dists(samples = sample.uniform.loss.75p, 
                                                         calibrated_samples = cal.uniform.loss.75p,
                                                         timeRange = timeRange, 
                                                         taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, 
                                                         binSize = binSize)
fd.uniform.loss.100p          <- compare_frequency_dists(samples = sample.uniform.loss.100p, 
                                                         calibrated_samples = cal.uniform.loss.100p,
                                                         timeRange = timeRange, 
                                                         taphCorrect = taphCorrect, 
                                                         correctForSite = correctForSite, 
                                                         binSize = binSize)

fd.singleton.ancient.loss.50p  <- compare_frequency_dists(samples = sample.singleton.ancient.loss.50p, 
                                                          calibrated_samples = cal.singleton.ancient.loss.50p,
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, 
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.ancient.loss.75p  <- compare_frequency_dists(samples = sample.singleton.ancient.loss.75p, 
                                                          calibrated_samples = cal.singleton.ancient.loss.75p,
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, 
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.ancient.loss.100p <- compare_frequency_dists(samples = sample.singleton.ancient.loss.100p, 
                                                          calibrated_samples = cal.singleton.ancient.loss.100p,
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect, 
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

fd.singleton.recent.loss.50p   <- compare_frequency_dists(samples = sample.singleton.recent.loss.50p, 
                                                          calibrated_samples = cal.singleton.recent.loss.50p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.recent.loss.75p   <- compare_frequency_dists(samples = sample.singleton.recent.loss.75p, 
                                                          calibrated_samples = cal.singleton.recent.loss.75p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.recent.loss.100p  <- compare_frequency_dists(samples = sample.singleton.recent.loss.100p, 
                                                          calibrated_samples = cal.singleton.recent.loss.100p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

fd.singleton.random.loss.50p   <- compare_frequency_dists(samples = sample.singleton.random.loss.50p, 
                                                          calibrated_samples = cal.singleton.random.loss.50p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.random.loss.75p   <- compare_frequency_dists(samples = sample.singleton.random.loss.75p, 
                                                          calibrated_samples = cal.singleton.random.loss.75p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.singleton.random.loss.100p  <- compare_frequency_dists(samples = sample.singleton.random.loss.100p, 
                                                          calibrated_samples = cal.singleton.random.loss.100p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

fd.bracketed.loss.50p          <- compare_frequency_dists(samples = sample.bracketed.loss.50p, 
                                                          calibrated_samples = cal.bracketed.loss.50p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.bracketed.loss.75p          <- compare_frequency_dists(samples = sample.bracketed.loss.75p, 
                                                          calibrated_samples = cal.bracketed.loss.75p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)
fd.bracketed.loss.100p         <- compare_frequency_dists(samples = sample.bracketed.loss.100p, 
                                                          calibrated_samples = cal.bracketed.loss.100p, 
                                                          timeRange = timeRange, 
                                                          taphCorrect = taphCorrect,
                                                          correctForSite = correctForSite, 
                                                          binSize = binSize)

# Organise data for plotting
plot.fd.sites.taphloss <- data.frame(matrix(NA, nrow = length(median.bin.age)*18, ncol = 5))
names(plot.fd.sites.taphloss) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.sites.taphloss[i,1] <- "baseline (taphonomic loss) replicates"
  plot.fd.sites.taphloss[i,2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[i,3] <- fd.baseline.taphloss$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[i,4] <- fd.baseline.taphloss$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[i,5] <- fd.baseline.taphloss$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)+i),1] <- "singleton ancient (50% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)+i),3] <- fd.singleton.ancient.loss.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)+i),4] <- fd.singleton.ancient.loss.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)+i),5] <- fd.singleton.ancient.loss.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*2+i),1] <- "singleton ancient (75% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*2+i),3] <- fd.singleton.ancient.loss.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*2+i),4] <- fd.singleton.ancient.loss.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*2+i),5] <- fd.singleton.ancient.loss.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*3+i),1] <- "singleton ancient (100% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*3+i),3] <- fd.singleton.ancient.loss.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*3+i),4] <- fd.singleton.ancient.loss.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*3+i),5] <- fd.singleton.ancient.loss.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*4+i),1] <- "singleton recent (50% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*4+i),3] <- fd.singleton.recent.loss.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*4+i),4] <- fd.singleton.recent.loss.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*4+i),5] <- fd.singleton.recent.loss.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*5+i),1] <- "singleton recent (75% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*5+i),3] <- fd.singleton.recent.loss.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*5+i),4] <- fd.singleton.recent.loss.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*5+i),5] <- fd.singleton.recent.loss.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*6+i),1] <- "singleton recent (100% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*6+i),3] <- fd.singleton.recent.loss.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*6+i),4] <- fd.singleton.recent.loss.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*6+i),5] <- fd.singleton.recent.loss.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*7+i),1] <- "singleton random (50% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*7+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*7+i),3] <- fd.singleton.random.loss.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*7+i),4] <- fd.singleton.random.loss.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*7+i),5] <- fd.singleton.random.loss.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*8+i),1] <- "singleton random (75% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*8+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*8+i),3] <- fd.singleton.random.loss.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*8+i),4] <- fd.singleton.random.loss.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*8+i),5] <- fd.singleton.random.loss.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*9+i),1] <- "singleton random (100% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*9+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*9+i),3] <- fd.singleton.random.loss.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*9+i),4] <- fd.singleton.random.loss.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*9+i),5] <- fd.singleton.random.loss.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*10+i),1] <- "bracketed (50% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*10+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*10+i),3] <- fd.bracketed.loss.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*10+i),4] <- fd.bracketed.loss.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*10+i),5] <- fd.bracketed.loss.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*11+i),1] <- "bracketed (75% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*11+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*11+i),3] <- fd.bracketed.loss.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*11+i),4] <- fd.bracketed.loss.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*11+i),5] <- fd.bracketed.loss.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*12+i),1] <- "bracketed (100% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*12+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*12+i),3] <- fd.bracketed.loss.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*12+i),4] <- fd.bracketed.loss.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*12+i),5] <- fd.bracketed.loss.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*13+i),1] <- "uniform (50% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*13+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*13+i),3] <- fd.uniform.loss.50p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*13+i),4] <- fd.uniform.loss.50p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*13+i),5] <- fd.uniform.loss.50p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*14+i),1] <- "uniform (75% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*14+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*14+i),3] <- fd.uniform.loss.75p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*14+i),4] <- fd.uniform.loss.75p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*14+i),5] <- fd.uniform.loss.75p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*15+i),1] <- "uniform (100% sites)"
  plot.fd.sites.taphloss[(length(median.bin.age)*15+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*15+i),3] <- fd.uniform.loss.100p$standardised.envelope[[1]][i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*15+i),4] <- fd.uniform.loss.100p$standardised.envelope[[1]][i,2]
  plot.fd.sites.taphloss[(length(median.bin.age)*15+i),5] <- fd.uniform.loss.100p$standardised.envelope[[1]][i,3]
  
  plot.fd.sites.taphloss[(length(median.bin.age)*16+i),1] <- "baseline (no loss)"
  plot.fd.sites.taphloss[(length(median.bin.age)*16+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*16+i),3] <- 0
  plot.fd.sites.taphloss[(length(median.bin.age)*16+i),4] <- standardised.baseline.freq.noloss[i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*16+i),5] <- 0
  
  plot.fd.sites.taphloss[(length(median.bin.age)*17+i),1] <- "baseline (taphonomic loss)"
  plot.fd.sites.taphloss[(length(median.bin.age)*17+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.taphloss[(length(median.bin.age)*17+i),3] <- 0
  plot.fd.sites.taphloss[(length(median.bin.age)*17+i),4] <- standardised.baseline.freq.taphloss.sites[i,1]
  plot.fd.sites.taphloss[(length(median.bin.age)*17+i),5] <- 0
}

write.csv(plot.fd.sites.taphloss, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/No population change/Plot_data-FD_sites_taphloss-5samples_20sims.csv")

# Plot (multi-panelled)
group.colors <- c("baseline (taphonomic loss) replicates" = "grey", 
                  "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                  "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                  "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                  "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                  "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                  "baseline (taphonomic loss)"    = "black", 
                  "baseline (no loss)"            = "grey48")


uniform <- plot.fd.sites.taphloss[plot.fd.sites.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.sites.taphloss[plot.fd.sites.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.sites.taphloss[plot.fd.sites.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.sites.taphloss[plot.fd.sites.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.sites.taphloss[plot.fd.sites.taphloss$sample %in% c("baseline (no loss)", "baseline (taphonomic loss)", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample,  "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
baseline <- plot.fd.sites.taphloss[plot.fd.sites.taphloss$sample %in% c("baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)"),] %>%
  mutate(samples = fct_relevel(sample, "baseline (taphonomic loss) replicates", "baseline (no loss)", "baseline (taphonomic loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  ylim(0, 0.085) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
fd.sites.taphloss <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed, #baseline,
                                           labels = c("A", "B", "C", "D", "E"),
                                           ncol = 2, nrow = 3)
fd.sites.taphloss
ggexport(fd.sites.taphloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/No_population_change-FD_sites_taphloss-5samples_20sims.png",
         width = 1500,
         height = 1000)
