#-------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study - make plots (5 samples)"
# author: "Rebecca Wheatley"
# date: "21 October 2021"
#-------------------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
theme_set(
  theme_bw() +
    theme(legend.position = "top"))

# Load plot data
setwd("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis")

# No population change over time
no.change.spd.noloss                    <- read.csv("./No population change/Plot_data-SPD_noloss-5samples_20sims.csv", header = T)
no.change.spd.taphloss                  <- read.csv("./No population change/Plot_data-SPD_taphloss-5samples_20sims.csv")
no.change.fd.samples.noloss             <- read.csv("./No population change/Plot_data-FD_samples_noloss-5samples_20sims.csv")
no.change.fd.samples.taphloss.correct   <- read.csv("./No population change/Plot_data-FD_samples_taphloss_corrected-5samples_20sims.csv")
no.change.fd.samples.taphloss.uncorrect <- read.csv("./No population change/Plot_data-FD_samples_taphloss_uncorrected-5samples_20sims.csv")
no.change.fd.sites.noloss               <- read.csv("./No population change/Plot_data-FD_sites_noloss-5samples_20sims.csv")
no.change.fd.sites.taphloss             <- read.csv("./No population change/Plot_data-FD_sites_taphloss-5samples_20sims.csv")

# Steady population growth over time
steady.growth.spd.noloss                    <- read.csv("./Steady population growth/Plot_data-SPD_noloss-5samples_20sims.csv")
steady.growth.spd.taphloss                  <- read.csv("./Steady population growth/Plot_data-SPD_taphloss-5samples_20sims.csv")
steady.growth.fd.samples.noloss             <- read.csv("./Steady population growth/Plot_data-FD_samples_noloss-5samples_20sims.csv")
steady.growth.fd.samples.taphloss.correct   <- read.csv("./Steady population growth/Plot_data-FD_samples_taphloss_corrected-5samples_20sims.csv")
steady.growth.fd.samples.taphloss.uncorrect <- read.csv("./Steady population growth/Plot_data-FD_samples_taphloss_uncorrected-5samples_20sims.csv")
steady.growth.fd.sites.noloss               <- read.csv("./Steady population growth/Plot_data-FD_sites_noloss-5samples_20sims.csv")
steady.growth.fd.sites.taphloss             <- read.csv("./Steady population growth/Plot_data-FD_sites_taphloss-5samples_20sims.csv")

# Exponential population growth over time
exponential.growth.spd.noloss                    <- read.csv("./Exponential population growth/Plot_data-SPD_noloss-5samples_20sims.csv")
exponential.growth.spd.taphloss                  <- read.csv("./Exponential population growth/Plot_data-SPD_taphloss-5samples_20sims.csv")
exponential.growth.fd.samples.noloss             <- read.csv("./Exponential population growth/Plot_data-FD_samples_noloss-5samples_20sims.csv")
exponential.growth.fd.samples.taphloss.correct   <- read.csv("./Exponential population growth/Plot_data-FD_samples_taphloss_corrected-5samples_20sims.csv")
exponential.growth.fd.samples.taphloss.uncorrect <- read.csv("./Exponential population growth/Plot_data-FD_samples_taphloss_uncorrected-5samples_20sims.csv")
exponential.growth.fd.sites.noloss               <- read.csv("./Exponential population growth/Plot_data-FD_sites_noloss-5samples_20sims.csv")
exponential.growth.fd.sites.taphloss             <- read.csv("./Exponential population growth/Plot_data-FD_sites_taphloss-5samples_20sims.csv")

# Initial population growth followed by decline
growth.decline.spd.noloss                    <- read.csv("./Growth then decline/Plot_data-SPD_noloss-5samples_20sims.csv")
growth.decline.spd.taphloss                  <- read.csv("./Growth then decline/Plot_data-SPD_taphloss-5samples_20sims.csv")
growth.decline.fd.samples.noloss             <- read.csv("./Growth then decline/Plot_data-FD_samples_noloss-5samples_20sims.csv")
growth.decline.fd.samples.taphloss.correct   <- read.csv("./Growth then decline/Plot_data-FD_samples_taphloss_corrected-5samples_20sims.csv")
growth.decline.fd.samples.taphloss.uncorrect <- read.csv("./Growth then decline/Plot_data-FD_samples_taphloss_uncorrected-5samples_20sims.csv")
growth.decline.fd.sites.noloss               <- read.csv("./Growth then decline/Plot_data-FD_sites_noloss-5samples_20sims.csv")
growth.decline.fd.sites.taphloss             <- read.csv("./Growth then decline/Plot_data-FD_sites_taphloss-5samples_20sims.csv")


# Plot (SPD no loss for all population trends)
group.colors <- c("baseline (no loss) replicates" = "grey", "baseline (no loss)" = "black")
no.change <- no.change.spd.noloss[no.change.spd.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
steady.growth <- steady.growth.spd.noloss[steady.growth.spd.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
exponential.growth <- exponential.growth.spd.noloss[exponential.growth.spd.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
growth.decline <- growth.decline.spd.noloss[growth.decline.spd.noloss$sample %in% c("baseline (no loss)", "baseline (no loss) replicates"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (no loss) replicates", "baseline (no loss)")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = median, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
spd.noloss <- ggarrange(no.change, steady.growth, exponential.growth, #growth.decline,
                        labels = c("A", "B", "C"),#, "D"),
                        ncol = 2, nrow = 2)
spd.noloss
ggexport(spd.noloss, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Population_trends-SPD_noloss-5samples_20sims.png",
         width = 1500,
         height = 1000)
