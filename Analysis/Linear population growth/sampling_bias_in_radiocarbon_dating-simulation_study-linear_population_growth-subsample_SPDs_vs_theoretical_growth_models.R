#--------------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (LINEAR population growth)
# subtitle: "Subsample SPDs vs theoretical growth models"
# author: "Rebecca Wheatley"
# date: "17 April 2023"
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
load("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Linear population growth/linear_pop_growth-raw_sample_and_calibrated_data-5s-99sims.RData")

# Load source functions
source("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/simulation_study-source3.R")

#----------------------------------------------------------------
# COMPARE SUBSAMPLED SPDS AGAINST THEORETICAL GROWTH MODELS
#----------------------------------------------------------------

# Number of Monte-Carlo simulations to run per test
nsim = 499

# The running mean to use for the SPDs
runm = 100

# The number of cores to run the MC simulations over
ncores = 8

# The number of subsamples per biased data set
subnum <- length(cal.uniform.50p)

# Set up database to store p value and discrepancy scores
pvals <- data.frame(matrix(NA, nrow = (15*subnum)*3, ncol = 7))
names(pvals) <- c("sample", "rep", "theoretical growth model", "discrepancy", "p-value", "nsim", "true growth model")

#----------------------
# I. UNIFORM (50% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
for (i in 1:subnum){
  uniform  <- modelTest(cal.uniform.50p[[i]], errors = sample.uniform.50p[[i]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (50% sites)", replicate = i, theoretical_model = "uniform")
  
  pvals[i,1] <- "uniform (50% sites)"
  pvals[i,2] <- i
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 50p/linear_pop_growth-uniform_50p_", i, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*1+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.uniform.50p[[index]], errors = sample.uniform.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (50% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "uniform (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 50p/linear_pop_growth-uniform_50p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*2+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1) } else { index <- subnum }
  
  uniform  <- modelTest(cal.uniform.50p[[index]], errors = sample.uniform.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (50% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "uniform (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 50p/linear_pop_growth-uniform_50p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# II. UNIFORM (75% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*3+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum }
  
  uniform  <- modelTest(cal.uniform.75p[[index]], errors = sample.uniform.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (75% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "uniform (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 75p/linear_pop_growth-uniform_75p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*4+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum }
  
  uniform  <- modelTest(cal.uniform.75p[[index]], errors = sample.uniform.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (75% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "uniform (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 75p/linear_pop_growth-uniform_75p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*5+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.uniform.75p[[index]], errors = sample.uniform.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (75% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "uniform (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 75p/linear_pop_growth-uniform_75p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# III. UNIFORM (100% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*6+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.uniform.100p[[index]], errors = sample.uniform.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (100% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "uniform (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 100p/linear_pop_growth-uniform_100p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*7+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.uniform.100p[[index]], errors = sample.uniform.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (100% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "uniform (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 100p/linear_pop_growth-uniform_100p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*8+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.uniform.100p[[index]], errors = sample.uniform.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "uniform (100% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "uniform (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Uniform 100p/linear_pop_growth-uniform_100p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# IV. SINGLETON ANCIENT (50% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*9+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.50p[[index]], errors = sample.singleton.ancient.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (50% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton ancient (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 50p/linear_pop_growth-singleton_ancient_50p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*10+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.50p[[index]], errors = sample.singleton.ancient.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (50% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton ancient (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 50p/linear_pop_growth-singleton_ancient_50p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*11+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.50p[[index]], errors = sample.singleton.ancient.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (50% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton ancient (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 50p/linear_pop_growth-singleton_ancient_50p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# V. SINGLETON ANCIENT (75% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*12+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.75p[[index]], errors = sample.singleton.ancient.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (75% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton ancient (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 75p/linear_pop_growth-singleton_ancient_75p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*13+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.75p[[index]], errors = sample.singleton.ancient.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (75% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton ancient (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 75p/linear_pop_growth-singleton_ancient_75p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*14+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.75p[[index]], errors = sample.singleton.ancient.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (75% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton ancient (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 75p/linear_pop_growth-singleton_ancient_75p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# VI. SINGLETON ANCIENT (100% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*15+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.100p[[index]], errors = sample.singleton.ancient.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (100% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton ancient (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 100p/linear_pop_growth-singleton_ancient_100p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*16+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.100p[[index]], errors = sample.singleton.ancient.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (100% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton ancient (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 100p/linear_pop_growth-singleton_ancient_100p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*17+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.ancient.100p[[index]], errors = sample.singleton.ancient.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton ancient (100% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton ancient (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton ancient 100p/linear_pop_growth-singleton_ancient_100p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# VII. SINGLETON RECENT (50% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*18+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.50p[[index]], errors = sample.singleton.recent.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (50% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton recent (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 50p/linear_pop_growth-singleton_recent_50p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*19+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.50p[[index]], errors = sample.singleton.recent.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (50% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton recent (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 50p/linear_pop_growth-singleton_recent_50p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*20+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.50p[[index]], errors = sample.singleton.recent.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (50% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton recent (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 50p/linear_pop_growth-singleton_recent_50p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# VIII. SINGLETON RECENT (75% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*21+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.75p[[index]], errors = sample.singleton.recent.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (75% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton recent (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 75p/linear_pop_growth-singleton_recent_75p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*22+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.75p[[index]], errors = sample.singleton.recent.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (75% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton recent (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 75p/linear_pop_growth-singleton_recent_75p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*23+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.75p[[index]], errors = sample.singleton.recent.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (75% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton recent (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 75p/linear_pop_growth-singleton_recent_75p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# IX. SINGLETON RECENT (100% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*24+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.100p[[index]], errors = sample.singleton.recent.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (100% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton recent (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 100p/linear_pop_growth-singleton_recent_100p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*25+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.100p[[index]], errors = sample.singleton.recent.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (100% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton recent (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 100p/linear_pop_growth-singleton_recent_100p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*26+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.recent.100p[[index]], errors = sample.singleton.recent.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton recent (100% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton recent (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton recent 100p/linear_pop_growth-singleton_recent_100p_", index, "_vs_exponential-5s-499sims.png"))
  
}


#----------------------
# X. SINGLETON RANDOM (50% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*27+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.50p[[index]], errors = sample.singleton.random.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (50% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton random (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 50p/linear_pop_growth-singleton_random_50p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*28+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.50p[[index]], errors = sample.singleton.random.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (50% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton random (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 50p/linear_pop_growth-singleton_random_50p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*29+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.50p[[index]], errors = sample.singleton.random.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (50% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton random (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 50p/linear_pop_growth-singleton_random_50p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# XI. SINGLETON RANDOM (75% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*30+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.75p[[index]], errors = sample.singleton.random.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (75% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton random (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 75p/linear_pop_growth-singleton_random_75p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*31+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.75p[[index]], errors = sample.singleton.random.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (75% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton random (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 75p/linear_pop_growth-singleton_random_75p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*32+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.75p[[index]], errors = sample.singleton.random.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (75% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton random (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 75p/linear_pop_growth-singleton_random_75p_", index, "_vs_exponential-5s-499sims.png"))
  
}


#----------------------
# XII. SINGLETON RANDOM (100% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*33+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.100p[[index]], errors = sample.singleton.random.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (100% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "singleton random (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 100p/linear_pop_growth-singleton_random_100p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*34+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.100p[[index]], errors = sample.singleton.random.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (100% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "singleton random (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 100p/linear_pop_growth-singleton_random_100p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*35+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.singleton.random.100p[[index]], errors = sample.singleton.random.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "singleton random (100% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "singleton random (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Singleton random 100p/linear_pop_growth-singleton_random_100p_", index, "_vs_exponential-5s-499sims.png"))
  
}


#----------------------
# XIII.BRACKETED (50% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*36+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.50p[[index]], errors = sample.bracketed.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (50% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "bracketed (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 50p/linear_pop_growth-bracketed_50p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*37+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.50p[[index]], errors = sample.bracketed.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (50% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "bracketed (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 50p/linear_pop_growth-bracketed_50p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*38+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.50p[[index]], errors = sample.bracketed.50p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (50% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "bracketed (50% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 50p/linear_pop_growth-bracketed_50p_", index, "_vs_exponential-5s-499sims.png"))
  
}



#----------------------
# XIV.BRACKETED (75% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*39+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.75p[[index]], errors = sample.bracketed.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (75% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "bracketed (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 75p/linear_pop_growth-bracketed_75p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*40+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.75p[[index]], errors = sample.bracketed.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (75% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "bracketed (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 75p/linear_pop_growth-bracketed_75p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*41+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.75p[[index]], errors = sample.bracketed.75p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (75% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "bracketed (75% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 75p/linear_pop_growth-bracketed_75p_", index, "_vs_exponential-5s-499sims.png"))
  
}

#----------------------
# XV.BRACKETED (100% SITES)
#----------------------

# VS UNIFORM THEORETICAL GROWTH MODEL
sl <- subnum*42+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.100p[[index]], errors = sample.bracketed.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "uniform", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (100% sites)", replicate = index, theoretical_model = "uniform")
  
  pvals[i,1] <- "bracketed (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "uniform"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 100p/linear_pop_growth-bracketed_100p_", index, "_vs_uniform-5s-499sims.png"))
  
}

# VS LINEAR THEORETICAL GROWTH MODEL
sl <- subnum*43+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.100p[[index]], errors = sample.bracketed.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "linear", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (100% sites)", replicate = index, theoretical_model = "linear")
  
  pvals[i,1] <- "bracketed (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "linear"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 100p/linear_pop_growth-bracketed_100p_", index, "_vs_linear-5s-499sims.png"))
  
}

# VS EXPONENTIAL THEORETICAL GROWTH MODEL
sl <- subnum*44+1
for (i in sl:(sl+subnum-1)){
  
  if (i%%(sl-1) > 0) { index <- i%%(sl-1)} else { index <- subnum}
  
  uniform  <- modelTest(cal.bracketed.100p[[index]], errors = sample.bracketed.100p[[index]]$error, nsim = nsim, timeRange = timeRange, 
                        model = "exponential", runm = runm, raw = TRUE, ncores = ncores)
  p_values <- calculate_p_value(uniform, sample_type = "bracketed (100% sites)", replicate = index, theoretical_model = "exponential")
  
  pvals[i,1] <- "bracketed (100% sites)"
  pvals[i,2] <- index
  pvals[i,3] <- "exponential"
  pvals[i,4] <- p_values$discrepancy
  if (length(p_values$p_total) > 0) { pvals[i,5] <- p_values$p_total } else { pvals[i,5] <- "NA" }
  pvals[i,6] <- nsim
  pvals[i,7] <- "linear"
  
  ggsave(p_values$plot2, 
         file = paste0("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/Linear population growth/Bracketed 100p/linear_pop_growth-bracketed_100p_", index, "_vs_exponential-5s-499sims.png"))
  
}


#----------------------
# XVI. SAVE WORKSPACE
#----------------------
write.csv(pvals, "linear_pop_growth-subsamples_vs_theoretical_growth_models-5s_499sims-pvalues.csv")
save.image("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/Linear population growth/linear_pop_growth-subsamples_vs_theoretical_growth_models-5s_499sims.RData")
