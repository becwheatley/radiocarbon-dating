#-------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (no population change) - 5 samples at 50%"
# author: "Rebecca Wheatley"
# date: "7 October 2021"
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

#-------------------------------------------------------------------------------------
# I. Simulating no change in population growth over time (without taphonomic loss)
#-------------------------------------------------------------------------------------

# Set the parameters we want for our simulated data:
timeRange        = c(12000, 200)   ## for Holocene sites
no_sites         = 100             ## the number of sites we want in our simulated data set
no_samples       = 5               ## the number of samples we want each site to have
pop_trend        = "no change"     ## the underlying population trend we want to mimic
percent_sites    = 50              ## the number of sites we want to apply the sampling method across (the remainder will be exhaustively sampled)
sampling_effort  = 3               ## the number of samples we want to take using the uniform sampling method
nsim             = 20              ## the number of times we want to replicate each sample

# Get the baseline data set/s:
baseline.data <- get_available_evidence(timeRange = timeRange, no_sites = no_sites, no_samples = no_samples, 
                                        pop_trend = pop_trend, nsim = nsim)

# Get samples:
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
                                        percent_sites = percent_sites, 
                                        nsim = nsim)

sample.singleton.random.50p  <- get_samples(evidence = baseline.data[[1]][[1]], 
                                        sampling_method = "singleton_random", 
                                        sampling_effort = sampling_effort, 
                                        percent_sites = percent_sites, 
                                        nsim = nsim)

sample.bracketed.50p         <- get_samples(evidence = baseline.data[[1]][[1]], 
                                        sampling_method = "bracketed", 
                                        sampling_effort = sampling_effort, 
                                        percent_sites = percent_sites, 
                                        nsim = nsim)

# Calibrate samples:
normalised      = TRUE ## are calibration curves (and, later, SPDs) normalised?
ncores          = 6    ## the number of threads to use when calibrating the radiocarbon dates

cal.baseline.noloss1 <- rcarbon::calibrate(x = baseline.data[[1]][[1]]$age, 
                                           errors = baseline.data[[1]][[1]]$error, 
                                           calCurves = 'shcal20', 
                                           normalised = normalised)

cal.baseline.noloss <- calibrate_samples(samples = baseline.data[[1]], 
                                         normalised = normalised, 
                                         ncores = ncores)

cal.uniform           <- calibrate_samples(samples = sample.uniform, 
                                           normalised = normalised, 
                                           ncores = ncores)

cal.singleton.ancient <- calibrate_samples(samples = sample.singleton.ancient, 
                                           normalised = normalised, 
                                           ncores = ncores)

cal.singleton.recent  <- calibrate_samples(samples = sample.singleton.recent, 
                                           normalised = normalised, 
                                           ncores = ncores)

cal.singleton.random  <- calibrate_samples(samples = sample.singleton.random, 
                                           normalised = normalised, 
                                           ncores = ncores)

cal.bracketed         <- calibrate_samples(samples = sample.bracketed, 
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

spd.uniform           <- compare_spds(calibrated_samples = cal.uniform, 
                                      timeRange = timeRange, 
                                      runm = runm, 
                                      normalised = normalised)

spd.singleton.ancient <- compare_spds(calibrated_samples = cal.singleton.ancient, 
                                      timeRange = timeRange, 
                                      runm = runm, 
                                      normalised = normalised)

spd.singleton.recent  <- compare_spds(calibrated_samples = cal.singleton.recent, 
                                      timeRange = timeRange, 
                                      runm = runm, 
                                      normalised = normalised)

spd.singleton.random  <- compare_spds(calibrated_samples = cal.singleton.random, 
                                      timeRange = timeRange, 
                                      runm = runm, 
                                      normalised = normalised)

spd.bracketed         <- compare_spds(calibrated_samples = cal.bracketed, 
                                      timeRange = timeRange, 
                                      runm = runm, 
                                      normalised = normalised)

# Rearrange data for plotting:
plot.spd.noloss <- data.frame(matrix(NA, nrow = length(spd.baseline.noloss$calBP)*7, ncol = 5))
names(plot.spd.noloss) <- c("sample", "years.ka", "lowerCI", "median", "upperCI")

for (i in 1:length(spd.baseline.noloss$calBP)){
  plot.spd.noloss[i,1] <- "baseline (no loss) replicates"
  plot.spd.noloss[i,2] <- spd.baseline.noloss$calBP[i]/1000
  plot.spd.noloss[i,3] <- spd.baseline.noloss$envelope[[1]][i,1]
  plot.spd.noloss[i,4] <- spd.baseline.noloss$envelope[[1]][i,2]
  plot.spd.noloss[i,5] <- spd.baseline.noloss$envelope[[1]][i,3]
}

for (i in 1:length(spd.singleton.ancient$calBP)){
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),1] <- "singleton ancient"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),2] <- spd.singleton.ancient$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),3] <- spd.singleton.ancient$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),4] <- spd.singleton.ancient$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)+i),5] <- spd.singleton.ancient$envelope[[1]][i,3]
}

for (i in 1:length(spd.singleton.recent$calBP)){
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),1] <- "singleton recent"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),2] <- spd.singleton.recent$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),3] <- spd.singleton.recent$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),4] <- spd.singleton.recent$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*2+i),5] <- spd.singleton.recent$envelope[[1]][i,3]
}

for (i in 1:length(spd.singleton.random$calBP)){
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),1] <- "singleton random"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),2] <- spd.singleton.random$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),3] <- spd.singleton.random$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),4] <- spd.singleton.random$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*3+i),5] <- spd.singleton.random$envelope[[1]][i,3]
}

for (i in 1:length(spd.bracketed$calBP)){
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),1] <- "bracketed"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),2] <- spd.bracketed$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),3] <- spd.bracketed$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),4] <- spd.bracketed$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*4+i),5] <- spd.bracketed$envelope[[1]][i,3]
}

for (i in 1:length(spd.uniform$calBP)){
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),1] <- "uniform"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),2] <- spd.uniform$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),3] <- spd.uniform$envelope[[1]][i,1]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),4] <- spd.uniform$envelope[[1]][i,2]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*5+i),5] <- spd.uniform$envelope[[1]][i,3]
}

for (i in 1:length(spd.uniform$calBP)){
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),1] <- "baseline (no loss)"
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),2] <- spd.baseline.noloss1$grid$calBP[i]/1000
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),3] <- 0
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),4] <- spd.baseline.noloss1$grid$PrDens[i]
  plot.spd.noloss[(length(spd.baseline.noloss$calBP)*6+i),5] <- 0
}

# Plot (combined):
group.colors <- c("baseline (no loss) replicates" = "grey", "singleton ancient" = "#FC8D62", "singleton recent" ="#8DA0CB", 
                  "singleton random" = "#E78AC3", "bracketed" = "#A6D854", "uniform" = "#66C2A5", "baseline (no loss)" = "black")

par(mfrow=c(1,1))
p <- plot.spd.noloss %>%
  mutate(sample = fct_relevel(sample, 
                              "baseline (no loss) replicates", "uniform", "bracketed", 
                              "singleton random", "singleton recent", "singleton ancient", 
                              "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ggtitle("Summed probability distribution (no population change)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
p

# Plot (multi panelled):
uniform <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "uniform"),] %>%
  mutate(sample = fct_relevel(sample, "uniform", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "singleton random"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "singleton ancient"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "singleton recent"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.spd.noloss[plot.spd.noloss$sample %in% c("baseline (no loss)", "bracketed"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
figure <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed,
                    labels = c("A", "B", "C", "D", "E"),
                    ncol = 2, nrow = 3)
figure


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
fd.baseline.noloss   <- compare_frequency_dists(samples = baseline.data[[1]],
                                                calibrated_samples = cal.baseline.noloss,
                                                timeRange = timeRange,
                                                taphCorrect = taphCorrect,
                                                correctForSite = correctForSite,
                                                binSize = binSize)


fd.uniform           <- compare_frequency_dists(samples = sample.uniform, 
                                                calibrated_samples = cal.uniform,
                                                timeRange = timeRange, 
                                                taphCorrect = taphCorrect, 
                                                correctForSite = correctForSite, 
                                                binSize = binSize)

fd.singleton.ancient <- compare_frequency_dists(samples = sample.singleton.ancient, 
                                                calibrated_samples = cal.singleton.ancient,
                                                timeRange = timeRange, 
                                                taphCorrect = taphCorrect, 
                                                correctForSite = correctForSite, 
                                                binSize = binSize)

fd.singleton.recent  <- compare_frequency_dists(samples = sample.singleton.recent, 
                                                calibrated_samples = cal.singleton.recent, 
                                                timeRange = timeRange, 
                                                taphCorrect = taphCorrect,
                                                correctForSite = correctForSite, 
                                                binSize = binSize)

fd.singleton.random  <- compare_frequency_dists(samples = sample.singleton.random, 
                                                calibrated_samples = cal.singleton.random, 
                                                timeRange = timeRange, 
                                                taphCorrect = taphCorrect,
                                                correctForSite = correctForSite, 
                                                binSize = binSize)

fd.bracketed         <- compare_frequency_dists(samples = sample.bracketed, 
                                                calibrated_samples = cal.bracketed, 
                                                timeRange = timeRange, 
                                                taphCorrect = taphCorrect,
                                                correctForSite = correctForSite, 
                                                binSize = binSize)

plot.fd.samples.noloss <- data.frame(matrix(NA, nrow = length(median.bin.age)*7, ncol = 5))
names(plot.fd.samples.noloss) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[i,1] <- "baseline (no loss) replicates"
  plot.fd.samples.noloss[i,2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[i,3] <- fd.baseline.noloss$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[i,4] <- fd.baseline.noloss$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[i,5] <- fd.baseline.noloss$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[(length(median.bin.age)+i),1] <- "singleton ancient"
  plot.fd.samples.noloss[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)+i),3] <- fd.singleton.ancient$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)+i),4] <- fd.singleton.ancient$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)+i),5] <- fd.singleton.ancient$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),1] <- "singleton recent"
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),3] <- fd.singleton.recent$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),4] <- fd.singleton.recent$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*2+i),5] <- fd.singleton.recent$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),1] <- "singleton random"
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),3] <- fd.singleton.random$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),4] <- fd.singleton.random$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*3+i),5] <- fd.singleton.random$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),1] <- "bracketed"
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),3] <- fd.bracketed$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),4] <- fd.bracketed$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*4+i),5] <- fd.bracketed$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),1] <- "uniform"
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),3] <- fd.uniform$standardised.envelope[[1]][i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),4] <- fd.uniform$standardised.envelope[[1]][i,2]
  plot.fd.samples.noloss[(length(median.bin.age)*5+i),5] <- fd.uniform$standardised.envelope[[1]][i,3]
}


for (i in 1:length(median.bin.age)){
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),1] <- "baseline (no loss)"
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),3] <- 0
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),4] <- standardised.baseline.freq.noloss[i,1]
  plot.fd.samples.noloss[(length(median.bin.age)*6+i),5] <- 0
}

# Combined plot
group.colors <- c("baseline (no loss) replicates" = "grey", "singleton ancient" = "#FC8D62", "singleton recent" ="#8DA0CB", 
                  "singleton random" = "#E78AC3", "bracketed" = "#A6D854", "uniform" = "#66C2A5", "baseline (no loss)" = "black")

p <- plot.fd.samples.noloss %>%
  mutate(sample = fct_relevel(sample, 
                              "baseline (no loss) replicates", "uniform", "bracketed", 
                              "singleton random", "singleton recent", "singleton ancient", 
                              "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_spline(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5, df = 40) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency") +
  ggtitle("Frequency distribution of dates (no population change)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000)) +
  theme_bw()
p

# Plot (multi panelled):
uniform <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "uniform"),] %>%
  mutate(sample = fct_relevel(sample, "uniform", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "singleton random"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "singleton ancient"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "singleton recent"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.samples.noloss[plot.fd.samples.noloss$sample %in% c("baseline (no loss)", "bracketed"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (dates)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
figure <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed,
                    labels = c("A", "B", "C", "D", "E"),
                    ncol = 2, nrow = 3)
figure


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

fd.uniform.sites           <- compare_frequency_dists(samples = sample.uniform, 
                                                      calibrated_samples = cal.uniform,
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect, 
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)

fd.singleton.ancient.sites <- compare_frequency_dists(samples = sample.singleton.ancient, 
                                                      calibrated_samples = cal.singleton.ancient,
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect, 
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)

fd.singleton.recent.sites  <- compare_frequency_dists(samples = sample.singleton.recent, 
                                                      calibrated_samples = cal.singleton.recent, 
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect,
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)

fd.singleton.random.sites  <- compare_frequency_dists(samples = sample.singleton.random, 
                                                      calibrated_samples = cal.singleton.random, 
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect,
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)

fd.bracketed.sites         <- compare_frequency_dists(samples = sample.bracketed, 
                                                      calibrated_samples = cal.bracketed, 
                                                      timeRange = timeRange, 
                                                      taphCorrect = taphCorrect,
                                                      correctForSite = correctForSite, 
                                                      binSize = binSize)


plot.fd.sites.noloss <- data.frame(matrix(NA, nrow = length(median.bin.age)*7, ncol = 5))
names(plot.fd.sites.noloss) <- c("sample", "median.bin.age", "lowerCI", "median", "upperCI")

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[i,1] <- "baseline (no loss) replicates"
  plot.fd.sites.noloss[i,2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[i,3] <- fd.baseline.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[i,4] <- fd.baseline.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[i,5] <- fd.baseline.sites$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[(length(median.bin.age)+i),1] <- "singleton ancient"
  plot.fd.sites.noloss[(length(median.bin.age)+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)+i),3] <- fd.singleton.ancient.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)+i),4] <- fd.singleton.ancient.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)+i),5] <- fd.singleton.ancient.sites$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),1] <- "singleton recent"
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),3] <- fd.singleton.recent.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),4] <- fd.singleton.recent.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*2+i),5] <- fd.singleton.recent.sites$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),1] <- "singleton random"
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),3] <- fd.singleton.random.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),4] <- fd.singleton.random.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*3+i),5] <- fd.singleton.random.sites$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),1] <- "bracketed"
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),3] <- fd.bracketed.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),4] <- fd.bracketed.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*4+i),5] <- fd.bracketed.sites$standardised.envelope[[1]][i,3]
}

for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),1] <- "uniform"
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),3] <- fd.uniform.sites$standardised.envelope[[1]][i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),4] <- fd.uniform.sites$standardised.envelope[[1]][i,2]
  plot.fd.sites.noloss[(length(median.bin.age)*5+i),5] <- fd.uniform.sites$standardised.envelope[[1]][i,3]
}


for (i in 1:length(median.bin.age)){
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),1] <- "baseline (no loss)"
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),2] <- median.bin.age[i]/1000
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),3] <- 0
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),4] <- standardised.baseline.freq.noloss.sites[i,1]
  plot.fd.sites.noloss[(length(median.bin.age)*6+i),5] <- 0
}

# Plot (combined)
group.colors <- c("baseline (no loss) replicates" = "grey", "singleton ancient" = "#FC8D62", "singleton recent" ="#8DA0CB", 
                  "singleton random" = "#E78AC3", "bracketed" = "#A6D854", "uniform" = "#66C2A5", "baseline (no loss)" = "black")

p <- plot.fd.sites.noloss %>%
  mutate(sample = fct_relevel(sample, 
                              "baseline (no loss) replicates", "uniform", "bracketed", 
                              "singleton random", "singleton recent", "singleton ancient", 
                              "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_spline(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5, df = 40) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency") +
  ggtitle("Frequency distribution of sites (no population change)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000)) +
  theme_bw()
p

# Plot (multi panelled):
uniform <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "uniform"),] %>%
  mutate(sample = fct_relevel(sample, "uniform", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.random <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "singleton random"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.ancient <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "singleton ancient"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
singleton.recent <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "singleton recent"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
bracketed <- plot.fd.sites.noloss[plot.fd.sites.noloss$sample %in% c("baseline (no loss)", "bracketed"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = median.bin.age, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = median.bin.age, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Standardised frequency (sites)") +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
figure <- ggarrange(uniform, singleton.random, singleton.ancient, singleton.recent, bracketed,
                    labels = c("A", "B", "C", "D", "E"),
                    ncol = 2, nrow = 3)
figure
