#--------------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study (EXPONENTIAL population growth)
# subtitle: "Subsample frequency distributions and SPDs vs theoretical growth models"
# author: "Rebecca Wheatley & Barry Brook"
# date: "22 July 2023"
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list=ls()); options(scipen=999,digits=9)

# Required packages
library("rcarbon")
library("ggplot2")

# Load source functions
source("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Analysis/simulation_study-source3.R")

# Create data frame for saving the number and proportion of correctly identified growth models
correctly_identified <- data.frame(matrix(NA, nrow = 45, ncol = 6))
names(correctly_identified) <- c("sampling_type", "underlying_model", "method", "total_sims", "sims_with_correct_ids", "prop_correct")
correctly_identified[1:45,2] <- "exponential"
correctly_identified[1:45,4] <- 1000 ##nsim

#----------------------------------------
# I. UNIFORM SUBSAMPLING
#----------------------------------------

# 50% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-UNIFORM_50P.RData")
cal_dates          <- sort_calibrated_dates(sample.uniform.50p, cal.uniform.50p) # Sort the calibrated sample dates
uniform.50p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="uniform_50p", true_model="Exponential") # fit and return the best growth model (based on a selected bin width)
uniform.50p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="uniform_50p", true_model="Exponential") # fit and return best growth model based on KDE Gaussian smoothing of the dates
rm(cal.uniform.50p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-UNIFORM_50P.RData")
uniform.50p.spd <- fit_models_spd(spd.uniform.50p, show_plot=FALSE, sampling_type="uniform_50p", true_model="Exponential")
rm(spd.uniform.50p)
rm(spd.baseline.noloss1)
correctly_identified[1:3,1] <- "uniform.50p"
correctly_identified[1,3] <- "fd.binned"
correctly_identified[2,3] <- "fd.kde"
correctly_identified[3,3] <- "spd"
correctly_identified[1,5] <- nrow(subset(uniform.50p.binned, best_model_correct==TRUE))
correctly_identified[2,5] <- nrow(subset(uniform.50p.kde, best_model_correct==TRUE))
correctly_identified[3,5] <- nrow(subset(uniform.50p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 75% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-UNIFORM_75P.RData")
cal_dates          <- sort_calibrated_dates(sample.uniform.75p, cal.uniform.75p)
uniform.75p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="uniform_75p", true_model="Exponential")
uniform.75p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="uniform_75p", true_model="Exponential")
rm(cal.uniform.75p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-UNIFORM_75P.RData")
uniform.75p.spd <- fit_models_spd(spd.uniform.75p, show_plot=FALSE, sampling_type="uniform_75p", true_model="Exponential")
rm(spd.uniform.75p)
rm(spd.baseline.noloss1)
correctly_identified[4:6,1] <- "uniform.75p"
correctly_identified[4,3] <- "fd.binned"
correctly_identified[5,3] <- "fd.kde"
correctly_identified[6,3] <- "spd"
correctly_identified[4,5] <- nrow(subset(uniform.75p.binned, best_model_correct==TRUE))
correctly_identified[5,5] <- nrow(subset(uniform.75p.kde, best_model_correct==TRUE))
correctly_identified[6,5] <- nrow(subset(uniform.75p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 100% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-UNIFORM_100P.RData")
cal_dates          <- sort_calibrated_dates(sample.uniform.100p, cal.uniform.100p)
uniform.100p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="uniform_100p", true_model="Exponential")
uniform.100p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="uniform_100p", true_model="Exponential")
rm(cal.uniform.100p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-UNIFORM_100P.RData")
uniform.100p.spd <- fit_models_spd(spd.uniform.100p, show_plot=FALSE, sampling_type="uniform_100p", true_model="Exponential")
rm(spd.uniform.100p)
rm(spd.baseline.noloss1)
correctly_identified[7:9,1] <- "uniform.100p"
correctly_identified[7,3] <- "fd.binned"
correctly_identified[8,3] <- "fd.kde"
correctly_identified[9,3] <- "spd"
correctly_identified[7,5] <- nrow(subset(uniform.100p.binned, best_model_correct==TRUE))
correctly_identified[8,5] <- nrow(subset(uniform.100p.kde, best_model_correct==TRUE))
correctly_identified[9,5] <- nrow(subset(uniform.100p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

#----------------------------------------
# II. SINGLETON ANCIENT SUBSAMPLING
#----------------------------------------

# 50% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_ANCIENT_50P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.ancient.50p, cal.singleton.ancient.50p)
singleton.ancient.50p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_ancient_50p", true_model="Exponential")
singleton.ancient.50p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_ancient_50p", true_model="Exponential")
rm(cal.singleton.ancient.50p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_ANCIENT_50P.RData")
singleton.ancient.50p.spd <- fit_models_spd(spd.singleton.ancient.50p, show_plot=FALSE, sampling_type="singleton_ancient_50p", true_model="Exponential")
rm(spd.singleton.ancient.50p)
rm(spd.baseline.noloss1)
correctly_identified[10:12,1] <- "singleton.ancient.50p"
correctly_identified[10,3] <- "fd.binned"
correctly_identified[11,3] <- "fd.kde"
correctly_identified[12,3] <- "spd"
correctly_identified[10,5] <- nrow(subset(singleton.ancient.50p.binned, best_model_correct==TRUE))
correctly_identified[11,5] <- nrow(subset(singleton.ancient.50p.kde, best_model_correct==TRUE))
correctly_identified[12,5] <- nrow(subset(singleton.ancient.50p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 75% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_ANCIENT_75P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.ancient.75p, cal.singleton.ancient.75p)
singleton.ancient.75p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_ancient_75p", true_model="Exponential")
singleton.ancient.75p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_ancient_75p", true_model="Exponential")
rm(cal.singleton.ancient.75p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_ANCIENT_75P.RData")
singleton.ancient.75p.spd <- fit_models_spd(spd.singleton.ancient.75p, show_plot=FALSE, sampling_type="singleton_ancient_75p", true_model="Exponential")
rm(spd.singleton.ancient.75p)
rm(spd.baseline.noloss1)
correctly_identified[13:15,1] <- "singleton.ancient.75p"
correctly_identified[13,3] <- "fd.binned"
correctly_identified[14,3] <- "fd.kde"
correctly_identified[15,3] <- "spd"
correctly_identified[13,5] <- nrow(subset(singleton.ancient.75p.binned, best_model_correct==TRUE))
correctly_identified[14,5] <- nrow(subset(singleton.ancient.75p.kde, best_model_correct==TRUE))
correctly_identified[15,5] <- nrow(subset(singleton.ancient.75p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 100% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_ANCIENT_100P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.ancient.100p, cal.singleton.ancient.100p)
singleton.ancient.100p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_ancient_100p", true_model="Exponential")
singleton.ancient.100p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_ancient_100p", true_model="Exponential")
rm(cal.singleton.ancient.100p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_ANCIENT_100P.RData")
singleton.ancient.100p.spd <- fit_models_spd(spd.singleton.ancient.100p, show_plot=FALSE, sampling_type="singleton_ancient_100p", true_model="Exponential")
rm(spd.singleton.ancient.100p)
rm(spd.baseline.noloss1)
correctly_identified[16:18,1] <- "singleton.ancient.100p"
correctly_identified[16,3] <- "fd.binned"
correctly_identified[17,3] <- "fd.kde"
correctly_identified[18,3] <- "spd"
correctly_identified[16,5] <- nrow(subset(singleton.ancient.100p.binned, best_model_correct==TRUE))
correctly_identified[17,5] <- nrow(subset(singleton.ancient.100p.kde, best_model_correct==TRUE))
correctly_identified[18,5] <- nrow(subset(singleton.ancient.100p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

#----------------------------------------
# III. SINGLETON RANDOM SUBSAMPLING
#----------------------------------------

# 50% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_RANDOM_50P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.random.50p, cal.singleton.random.50p)
singleton.random.50p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_random_50p", true_model="Exponential")
singleton.random.50p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_random_50p", true_model="Exponential")
rm(cal.singleton.random.50p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RANDOM_50P.RData")
singleton.random.50p.spd <- fit_models_spd(spd.singleton.random.50p, show_plot=FALSE, sampling_type="singleton_random_50p", true_model="Exponential")
rm(spd.singleton.random.50p)
rm(spd.baseline.noloss1)
correctly_identified[19:21,1] <- "singleton.random.50p"
correctly_identified[19,3] <- "fd.binned"
correctly_identified[20,3] <- "fd.kde"
correctly_identified[21,3] <- "spd"
correctly_identified[19,5] <- nrow(subset(singleton.random.50p.binned, best_model_correct==TRUE))
correctly_identified[20,5] <- nrow(subset(singleton.random.50p.kde, best_model_correct==TRUE))
correctly_identified[21,5] <- nrow(subset(singleton.random.50p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 75% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_RANDOM_75P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.random.75p, cal.singleton.random.75p)
singleton.random.75p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_random_75p", true_model="Exponential")
singleton.random.75p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_random_75p", true_model="Exponential")
rm(cal.singleton.random.75p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RANDOM_75P.RData")
singleton.random.75p.spd <- fit_models_spd(spd.singleton.random.75p, show_plot=FALSE, sampling_type="singleton_random_75p", true_model="Exponential")
rm(spd.singleton.random.75p)
rm(spd.baseline.noloss1)
correctly_identified[22:24,1] <- "singleton.random.75p"
correctly_identified[22,3] <- "fd.binned"
correctly_identified[23,3] <- "fd.kde"
correctly_identified[24,3] <- "spd"
correctly_identified[22,5] <- nrow(subset(singleton.random.75p.binned, best_model_correct==TRUE))
correctly_identified[23,5] <- nrow(subset(singleton.random.75p.kde, best_model_correct==TRUE))
correctly_identified[24,5] <- nrow(subset(singleton.random.75p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 100% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_RANDOM_100P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.random.100p, cal.singleton.random.100p)
singleton.random.100p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_random_100p", true_model="Exponential")
singleton.random.100p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_random_100p", true_model="Exponential")
rm(cal.singleton.random.100p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RANDOM_100P.RData")
singleton.random.100p.spd <- fit_models_spd(spd.singleton.random.100p, show_plot=FALSE, sampling_type="singleton_random_100p", true_model="Exponential")
rm(spd.singleton.random.100p)
rm(spd.baseline.noloss1)
correctly_identified[25:27,1] <- "singleton.random.100p"
correctly_identified[25,3] <- "fd.binned"
correctly_identified[26,3] <- "fd.kde"
correctly_identified[27,3] <- "spd"
correctly_identified[25,5] <- nrow(subset(singleton.random.100p.binned, best_model_correct==TRUE))
correctly_identified[26,5] <- nrow(subset(singleton.random.100p.kde, best_model_correct==TRUE))
correctly_identified[27,5] <- nrow(subset(singleton.random.100p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

#----------------------------------------
# IV. SINGLETON RECENT SUBSAMPLING
#----------------------------------------

# 50% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_RECENT_50P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.recent.50p, cal.singleton.recent.50p)
singleton.recent.50p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_recent_50p", true_model="Exponential")
singleton.recent.50p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_recent_50p", true_model="Exponential")
rm(cal.singleton.recent.50p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RECENT_50P.RData")
singleton.recent.50p.spd <- fit_models_spd(spd.singleton.recent.50p, show_plot=FALSE, sampling_type="singleton_recent_50p", true_model="Exponential")
rm(spd.singleton.recent.50p)
rm(spd.baseline.noloss1)
correctly_identified[28:30,1] <- "singleton.recent.50p"
correctly_identified[28,3] <- "fd.binned"
correctly_identified[29,3] <- "fd.kde"
correctly_identified[30,3] <- "spd"
correctly_identified[28,5] <- nrow(subset(singleton.recent.50p.binned, best_model_correct==TRUE))
correctly_identified[29,5] <- nrow(subset(singleton.recent.50p.kde, best_model_correct==TRUE))
correctly_identified[30,5] <- nrow(subset(singleton.recent.50p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")


# 75% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_RECENT_75P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.recent.75p, cal.singleton.recent.75p)
singleton.recent.75p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_recent_75p", true_model="Exponential")
singleton.recent.75p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_recent_75p", true_model="Exponential")
rm(cal.singleton.recent.75p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RECENT_75P.RData")
singleton.recent.75p.spd <- fit_models_spd(spd.singleton.recent.75p, show_plot=FALSE, sampling_type="singleton_recent_75p", true_model="Exponential")
rm(spd.singleton.recent.75p)
rm(spd.baseline.noloss1)
correctly_identified[31:33,1] <- "singleton.recent.75p"
correctly_identified[31,3] <- "fd.binned"
correctly_identified[32,3] <- "fd.kde"
correctly_identified[33,3] <- "spd"
correctly_identified[31,5] <- nrow(subset(singleton.recent.75p.binned, best_model_correct==TRUE))
correctly_identified[32,5] <- nrow(subset(singleton.recent.75p.kde, best_model_correct==TRUE))
correctly_identified[30,5] <- nrow(subset(singleton.recent.75p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 100% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-SINGLETON_RECENT_100P.RData")
cal_dates                    <- sort_calibrated_dates(sample.singleton.recent.100p, cal.singleton.recent.100p)
singleton.recent.100p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="singleton_recent_100p", true_model="Exponential")
singleton.recent.100p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="singleton_recent_100p", true_model="Exponential")
rm(cal.singleton.recent.100p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RECENT_100P.RData")
singleton.recent.100p.spd <- fit_models_spd(spd.singleton.recent.100p, show_plot=FALSE, sampling_type="singleton_recent_100p", true_model="Exponential")
rm(spd.singleton.recent.100p)
rm(spd.baseline.noloss1)
correctly_identified[34:36,1] <- "singleton.recent.100p"
correctly_identified[34,3] <- "fd.binned"
correctly_identified[35,3] <- "fd.kde"
correctly_identified[36,3] <- "spd"
correctly_identified[34,5] <- nrow(subset(singleton.recent.100p.binned, best_model_correct==TRUE))
correctly_identified[35,5] <- nrow(subset(singleton.recent.100p.kde, best_model_correct==TRUE))
correctly_identified[36,5] <- nrow(subset(singleton.recent.100p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

#----------------------------------------
# V. BRACKETED SUBSAMPLING
#----------------------------------------

# 50% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-BRACKETED_50P.RData")
cal_dates            <- sort_calibrated_dates(sample.bracketed.50p, cal.bracketed.50p)
bracketed.50p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="bracketed_50p", true_model="Exponential")
bracketed.50p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="bracketed_50p", true_model="Exponential")
rm(cal.bracketed.50p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BRACKETED_50P.RData")
bracketed.50p.spd <- fit_models_spd(spd.bracketed.50p, show_plot=FALSE, sampling_type="bracketed_50p", true_model="Exponential")
rm(spd.bracketed.50p)
rm(spd.baseline.noloss1)
correctly_identified[37:39,1] <- "bracketed.50p"
correctly_identified[37,3] <- "fd.binned"
correctly_identified[38,3] <- "fd.kde"
correctly_identified[39,3] <- "spd"
correctly_identified[37,5] <- nrow(subset(bracketed.50p.binned, best_model_correct==TRUE))
correctly_identified[38,5] <- nrow(subset(bracketed.50p.kde, best_model_correct==TRUE))
correctly_identified[39,5] <- nrow(subset(bracketed.50p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 75% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-BRACKETED_75P.RData")
cal_dates            <- sort_calibrated_dates(sample.bracketed.75p, cal.bracketed.75p)
bracketed.75p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="bracketed_75p", true_model="Exponential")
bracketed.75p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="bracketed_75p", true_model="Exponential")
rm(cal.bracketed.75p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BRACKETED_75P.RData")
bracketed.75p.spd <- fit_models_spd(spd.bracketed.75p, show_plot=FALSE, sampling_type="bracketed_75p", true_model="Exponential")
rm(spd.bracketed.75p)
rm(spd.baseline.noloss1)
correctly_identified[40:42,1] <- "bracketed.75p"
correctly_identified[40,3] <- "fd.binned"
correctly_identified[41,3] <- "fd.kde"
correctly_identified[42,3] <- "spd"
correctly_identified[40,5] <- nrow(subset(bracketed.75p.binned, best_model_correct==TRUE))
correctly_identified[41,5] <- nrow(subset(bracketed.75p.kde, best_model_correct==TRUE))
correctly_identified[42,5] <- nrow(subset(bracketed.75p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# 100% sites
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-raw_sample_and_calibrated_data-5s-1000sims-BRACKETED_100P.RData")
cal_dates            <- sort_calibrated_dates(sample.bracketed.100p, cal.bracketed.100p)
bracketed.100p.binned <- fit_models_binned(cal_dates, bin_width=500, show_plot=FALSE, sampling_type="bracketed_100p", true_model="Exponential")
bracketed.100p.kde    <- fit_models_kde(cal_dates, show_plot=FALSE, sampling_type="bracketed_100p", true_model="Exponential")
rm(cal.bracketed.100p)
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BRACKETED_100P.RData")
bracketed.100p.spd <- fit_models_spd(spd.bracketed.100p, show_plot=FALSE, sampling_type="bracketed_100p", true_model="Exponential")
rm(spd.bracketed.100p)
rm(spd.baseline.noloss1)
correctly_identified[43:45,1] <- "bracketed.100p"
correctly_identified[43,3] <- "fd.binned"
correctly_identified[44,3] <- "fd.kde"
correctly_identified[45,3] <- "spd"
correctly_identified[43,5] <- nrow(subset(bracketed.100p.binned, best_model_correct==TRUE))
correctly_identified[44,5] <- nrow(subset(bracketed.100p.kde, best_model_correct==TRUE))
correctly_identified[45,5] <- nrow(subset(bracketed.100p.spd, best_model_correct==TRUE))
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

#----------------------------------------
# VI. COMBINE AND SAVE RESULTS
#----------------------------------------

best.models <- rbind(uniform.50p.binned, uniform.75p.binned, uniform.100p.binned,
                     singleton.ancient.50p.binned, singleton.ancient.75p.binned, singleton.ancient.100p.binned,
                     singleton.recent.50p.binned, singleton.recent.75p.binned, singleton.recent.100p.binned,
                     singleton.random.50p.binned, singleton.random.75p.binned, singleton.random.100p.binned,
                     bracketed.50p.binned, bracketed.75p.binned, bracketed.100p.binned,
                     
                     uniform.50p.kde, uniform.75p.kde, uniform.100p.kde,
                     singleton.ancient.50p.kde, singleton.ancient.75p.kde, singleton.ancient.100p.kde,
                     singleton.random.50p.kde, singleton.random.75p.kde, singleton.random.100p.kde,
                     singleton.recent.50p.kde, singleton.recent.75p.kde, singleton.recent.100p.kde,
                     bracketed.50p.kde, bracketed.75p.kde, bracketed.100p.kde,
                     
                     uniform.50p.spd, uniform.75p.spd, uniform.100p.spd,
                     singleton.ancient.50p.spd, singleton.ancient.75p.spd, singleton.ancient.100p.spd,
                     singleton.random.50p.spd, singleton.random.75p.spd, singleton.random.100p.spd,
                     singleton.recent.50p.spd, singleton.recent.75p.spd, singleton.recent.100p.spd,
                     bracketed.50p.spd, bracketed.75p.spd, bracketed.100p.spd)

write.csv(best.models, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims.csv")

correctly_identified$prop_correct <- correctly_identified$sims_with_correct_ids/correctly_identified$total_sims
write.csv(correctly_identified, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")
