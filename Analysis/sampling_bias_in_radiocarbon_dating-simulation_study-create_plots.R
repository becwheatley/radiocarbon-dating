#-------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study - make plots (5 samples)"
# author: "Rebecca Wheatley"
# date: "6 December 2021"
#-------------------------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
theme_set(
  theme_bw() +
    theme(legend.position = "top"))

# Load plot data
setwd("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis")


uniform.spd       <- read.csv("./Uniform population growth/Plot_data-SPD-5samples_99sims.csv", header = T)
linear.spd        <- read.csv("./Linear population growth/Plot_data-SPD-5samples_99sims.csv", header = T)
exponential.spd   <- read.csv("./Exponential population growth/Plot_data-SPD-5samples_99sims.csv", header = T)


# Plot (SPD no loss for all population trends)
group.colors <- c("baseline (replicates)" = "grey", "baseline" = "black")
uniform <- uniform.spd[uniform.spd$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
linear <- linear.spd[linear.spd$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))
exponential <- exponential.spd[exponential.spd$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(12200/1000, 200/1000, by = -6000/1000), limits = c(12200/1000, 200/1000))

spd <- ggarrange(uniform, linear, exponential,
                        labels = c("A", "B", "C"),
                        ncol = 1, nrow = 3,
                 common.legend = TRUE)
spd
ggexport(spd, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Population_trends-SPD-5samples_99sims.png",
         width = 700,
         height = 1200)
