# Clear workspace
rm(list = ls())

# Load R packages required
library(rcarbon)
library(data.table)
library(tidyverse)
library(ggpubr)
theme_set(
  theme_bw() +
    theme(legend.position = "top"))

# Load saved sample and calibrated sample data
load("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/uniform_pop_growth-raw_sample_and_calibrated_data-5s-99sims.RData")

# Loaded stuff
timeRange = c(12000, 1000)
runm = 500
spdnormalised = TRUE
sample <- baseline.data[[1]][[1]]
cal.sample <- cal.baseline.noloss1

# Set number of simulations
nsim = 5

# Generate the observed SPD
observed.spd <- spd(x = cal.sample, timeRange = timeRange, runm = runm, spdnormalised = normalised)
finalSPD<- observed$grid$PrDens

# MONTE CARLO SIMULATIONS
sim <- matrix(NA,nrow=length(observed.spd), ncol=nsim)
fit.time <- seq(timeRange[1],timeRange[2],-1)
pred.time <- fit.time

# Fit desired theoretical growth model to the observed SPD
fit <- NA
if (model=="exponential"){
  fit <- nls(y ~ exp(a + b * x), data=data.frame(x=fit.time, y=finalSPD), start=list(a=a, b=b))
  est <- predict(fit, list(x=pred.time))
  predgrid <- data.frame(calBP=pred.time, PrDens=est)
} else if (model=="uniform"){
  predgrid <- data.frame(calBP=pred.time, PrDens=mean(finalSPD))
} else if (model=="linear"){
  fit <- lm(y ~ x, data=data.frame(x=fit.time, y=finalSPD))
  est <- predict(fit, list(x=pred.time))
  predgrid <- data.frame(calBP=pred.time, PrDens=est)
}

# Prepare the sampling grid for MC simulations
# Prepare Sampling Grid(s)
cragrids = vector("list",length=1)
tmp.grid <- uncalibrate(as.CalGrid(predgrid), calCurves="shcal20", compact=FALSE, verbose=TRUE)
cragrids[[i]] <- tmp.grid
