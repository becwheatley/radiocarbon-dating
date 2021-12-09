#--------------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (uniform population growth) - subsample SPDs vs theoretical growth models"
# author: "Rebecca Wheatley"
# date: "9 December 2021"
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load R packages required
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

# Load saved sample and calibrated sample data
load("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Uniform population growth/uniform_pop_growth-raw_sample_and_calibrated_data-5s-99sims.RData")

# Load source functions
source("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/simulation_study-source3.R")

#----------------------------------------------------------------
# I. COMPARE SUBSAMPLED SPDS AGAINST THEORETICAL GROWTH MODELS
#----------------------------------------------------------------

# Number of Monte-Carlo simulations to run per test
nsim = 499

# The running mean to use for the SPDs
runm = 100

# The number of cores to run the MC simulations over
ncores = 8

# Set up database to store p value and discrepancy scores
pvals <- data.frame(matrix(NA, nrow = (15*length(cal.uniform.50p))*3, ncol = 7))
names(pvals) <- c("sample", "rep", "theoretical growth model", "discrepancy", "p-value", "nsim", "true growth model")

#----------------------
# UNIFORM (50% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
for (i in 1:length(cal.uniform.50p)){
  uniform  <- modelTest(cal.uniform.50p[[i]], errors = sample.uniform.50p[[i]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (50% sites)", replicate = i, theoretical_model = "uniform")
  
  pvals[i,1] <- "uniform (50% sites)"
  pvals[i,2] <- i
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "uniform"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Uniform population growth/Uniform 50p/uniform_pop_growth-uniform_50p_", i, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- (length(cal.uniform.50p)*1)+1
for (i in sl:sl+length(cal.uniform.50p)){
  uniform  <- modelTest(cal.uniform.50p[[i%%(sl-1)]], errors = sample.uniform.50p[[i%%(sl-1)]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (50% sites)", replicate = i%%(sl-1), theoretical_model = "linear")
  
  pvals[i,1] <- "uniform (50% sites)"
  pvals[i,2] <- i%%(sl-1)
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "uniform"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Uniform population growth/Uniform 50p/uniform_pop_growth-uniform_50p_", i%%(sl-1), "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- (length(cal.uniform.50p)*2)+1
for (i in sl:sl+length(cal.uniform.50p)){
  uniform  <- modelTest(cal.uniform.50p[[i%%(sl-1)]], errors = sample.uniform.50p[[i%%(sl-1)]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (50% sites)", replicate = i%%(sl-1), theoretical_model = "exponential")
  
  pvals[i,1] <- "uniform (50% sites)"
  pvals[i,2] <- i%%(sl-1)
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "uniform"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Uniform population growth/Uniform 50p/uniform_pop_growth-uniform_50p_", i%%(sl-1), "_vs_exponential-5s-499sims.png"))
  
}