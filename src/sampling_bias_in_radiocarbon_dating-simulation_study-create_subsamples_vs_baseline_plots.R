#--------------------------------------------------------------------------------------------------------------------------------------------------
# title: "Sampling bias in radiocarbon dating project: simulation study
# subtitle: "Make plots"
# author: "Rebecca Wheatley"
# date: "23 July 2023"
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Clear workspace
rm(list=ls()); options(scipen=999,digits=9)

#-----------------------------------------------------------
# BASELINE VS SUBSAMPLE FIGURES
#-----------------------------------------------------------

# LOAD AND SAVE DATA
## Uniform population growth
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-BASELINE_ONLY.RData")
write.csv(spd.baseline.noloss$calBP, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-CALBP.csv")
write.csv(spd.baseline.noloss1$grid, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BASELINE_SINGLE.csv")
write.csv(spd.baseline.noloss$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BASELINE_ONLY.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-UNIFORM_50P.RData")
write.csv(spd.uniform.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-UNIFORM_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-UNIFORM_75P.RData")
write.csv(spd.uniform.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-UNIFORM_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-UNIFORM_100P.RData")
write.csv(spd.uniform.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-UNIFORM_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_RANDOM_50P.RData")
write.csv(spd.singleton.random.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_RANDOM_75P.RData")
write.csv(spd.singleton.random.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_RANDOM_100P.RData")
write.csv(spd.singleton.random.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_ANCIENT_50P.RData")
write.csv(spd.singleton.ancient.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_ANCIENT_75P.RData")
write.csv(spd.singleton.ancient.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_ANCIENT_100P.RData")
write.csv(spd.singleton.ancient.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_RECENT_50P.RData")
write.csv(spd.singleton.recent.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RECENT_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_RECENT_75P.RData")
write.csv(spd.singleton.recent.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RECENT_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-SINGLETON_RECENT_100P.RData")
write.csv(spd.singleton.recent.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RECENT_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-BRACKETED_50P.RData")
write.csv(spd.bracketed.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BRACKETED_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-BRACKETED_75P.RData")
write.csv(spd.bracketed.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BRACKETED_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Uniform population growth/uniform_pop_growth-spds-BRACKETED_100P.RData")
write.csv(spd.bracketed.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BRACKETED_100P.csv")
rm(list=ls())

## Linear population growth
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-BASELINE_ONLY.RData")
write.csv(spd.baseline.noloss$calBP, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-CALBP.csv")
write.csv(spd.baseline.noloss1$grid, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BASELINE_SINGLE.csv")
write.csv(spd.baseline.noloss$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BASELINE_ONLY.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-UNIFORM_50P.RData")
write.csv(spd.uniform.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-UNIFORM_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-UNIFORM_75P.RData")
write.csv(spd.uniform.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-UNIFORM_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-UNIFORM_100P.RData")
write.csv(spd.uniform.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-UNIFORM_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_RANDOM_50P.RData")
write.csv(spd.singleton.random.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_RANDOM_75P.RData")
write.csv(spd.singleton.random.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_RANDOM_100P.RData")
write.csv(spd.singleton.random.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_ANCIENT_50P.RData")
write.csv(spd.singleton.ancient.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_ANCIENT_75P.RData")
write.csv(spd.singleton.ancient.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_ANCIENT_100P.RData")
write.csv(spd.singleton.ancient.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_RECENT_50P.RData")
write.csv(spd.singleton.recent.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RECENT_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_RECENT_75P.RData")
write.csv(spd.singleton.recent.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RECENT_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-SINGLETON_RECENT_100P.RData")
write.csv(spd.singleton.recent.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RECENT_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-BRACKETED_50P.RData")
write.csv(spd.bracketed.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BRACKETED_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-BRACKETED_75P.RData")
write.csv(spd.bracketed.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BRACKETED_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Linear population growth/linear_pop_growth-spds-BRACKETED_100P.RData")
write.csv(spd.bracketed.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BRACKETED_100P.csv")
rm(list=ls())

## Exponential population growth
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BASELINE_ONLY.RData")
write.csv(spd.baseline.noloss$calBP, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-CALBP.csv")
write.csv(spd.baseline.noloss1$grid, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BASELINE_SINGLE.csv")
write.csv(spd.baseline.noloss$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BASELINE_ONLY.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-UNIFORM_50P.RData")
write.csv(spd.uniform.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-UNIFORM_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-UNIFORM_75P.RData")
write.csv(spd.uniform.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-UNIFORM_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-UNIFORM_100P.RData")
write.csv(spd.uniform.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-UNIFORM_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RANDOM_50P.RData")
write.csv(spd.singleton.random.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RANDOM_75P.RData")
write.csv(spd.singleton.random.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RANDOM_100P.RData")
write.csv(spd.singleton.random.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_ANCIENT_50P.RData")
write.csv(spd.singleton.ancient.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_ANCIENT_75P.RData")
write.csv(spd.singleton.ancient.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_ANCIENT_100P.RData")
write.csv(spd.singleton.ancient.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RECENT_50P.RData")
write.csv(spd.singleton.recent.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RECENT_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RECENT_75P.RData")
write.csv(spd.singleton.recent.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RECENT_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-SINGLETON_RECENT_100P.RData")
write.csv(spd.singleton.recent.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RECENT_100P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BRACKETED_50P.RData")
write.csv(spd.bracketed.50p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BRACKETED_50P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BRACKETED_75P.RData")
write.csv(spd.bracketed.75p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BRACKETED_75P.csv")
rm(list=ls())
load("E:/Radiocarbon dating analysis/Analysis/Exponential population growth/exponential_pop_growth-spds-BRACKETED_100P.RData")
write.csv(spd.bracketed.100p$envelope, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BRACKETED_100P.csv")
rm(list=ls())

# LOAD DATA FOR PLOTS
u.calBP <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-CALBP.csv")
u.baseline1 <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BASELINE_SINGLE.csv")
u.baseline <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BASELINE_ONLY.csv")
u.uniform.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-UNIFORM_50P.csv")
u.uniform.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-UNIFORM_75P.csv")
u.uniform.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-UNIFORM_100P.csv")
u.singleton.ancient.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_50P.csv")
u.singleton.ancient.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_75P.csv")
u.singleton.ancient.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_100P.csv")
u.singleton.random.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_50P.csv")
u.singleton.random.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_75P.csv")
u.singleton.random.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_100P.csv")
u.singleton.recent.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RECENT_50P.csv")
u.singleton.recent.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RECENT_75P.csv")
u.singleton.recent.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-SINGLETON_RECENT_100P.csv")
u.bracketed.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BRACKETED_50P.csv")
u.bracketed.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BRACKETED_75P.csv")
u.bracketed.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-BRACKETED_100P.csv")

l.calBP <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-CALBP.csv")
l.baseline1 <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BASELINE_SINGLE.csv")
l.baseline <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BASELINE_ONLY.csv")
l.uniform.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-UNIFORM_50P.csv")
l.uniform.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-UNIFORM_75P.csv")
l.uniform.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-UNIFORM_100P.csv")
l.singleton.ancient.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_50P.csv")
l.singleton.ancient.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_75P.csv")
l.singleton.ancient.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_100P.csv")
l.singleton.random.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_50P.csv")
l.singleton.random.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_75P.csv")
l.singleton.random.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_100P.csv")
l.singleton.recent.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RECENT_50P.csv")
l.singleton.recent.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RECENT_75P.csv")
l.singleton.recent.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-SINGLETON_RECENT_100P.csv")
l.bracketed.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BRACKETED_50P.csv")
l.bracketed.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BRACKETED_75P.csv")
l.bracketed.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-BRACKETED_100P.csv")

e.calBP <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-CALBP.csv")
e.baseline1 <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BASELINE_SINGLE.csv")
e.baseline <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BASELINE_ONLY.csv")
e.uniform.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-UNIFORM_50P.csv")
e.uniform.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-UNIFORM_75P.csv")
e.uniform.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-UNIFORM_100P.csv")
e.singleton.ancient.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_50P.csv")
e.singleton.ancient.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_75P.csv")
e.singleton.ancient.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_ANCIENT_100P.csv")
e.singleton.random.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_50P.csv")
e.singleton.random.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_75P.csv")
e.singleton.random.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RANDOM_100P.csv")
e.singleton.recent.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RECENT_50P.csv")
e.singleton.recent.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RECENT_75P.csv")
e.singleton.recent.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-SINGLETON_RECENT_100P.csv")
e.bracketed.50p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BRACKETED_50P.csv")
e.bracketed.75p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BRACKETED_75P.csv")
e.bracketed.100p <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-BRACKETED_100P.csv")

# REARRANGE DATA FOR PLOTTING:
u.dates <- nrow(u.calBP)
plot.uniform <- data.frame(matrix(NA, nrow = (u.dates*17), ncol = 5))
names(plot.uniform) <- c("sample", "years.ka", "lowerCI", "mean", "upperCI")
for (i in 1:nrow(u.calBP)){
 plot.uniform[i,1] <- "baseline (replicates)"
 plot.uniform[i,2] <- u.calBP[i,2]/1000
 plot.uniform[i,3] <- u.baseline[i,2]
 plot.uniform[i,4] <- u.baseline[i,3]
 plot.uniform[i,5] <- u.baseline[i,4]
   
 plot.uniform[(u.dates+i),1] <- "singleton ancient (50% sites)"
 plot.uniform[(u.dates+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates+i),3] <- u.singleton.ancient.50p[i,2]
 plot.uniform[(u.dates+i),4] <- u.singleton.ancient.50p[i,3]
 plot.uniform[(u.dates+i),5] <- u.singleton.ancient.50p[i,4]

 plot.uniform[(u.dates*2+i),1] <- "singleton ancient (75% sites)"
 plot.uniform[(u.dates*2+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*2+i),3] <- u.singleton.ancient.75p[i,2]
 plot.uniform[(u.dates*2+i),4] <- u.singleton.ancient.75p[i,3]
 plot.uniform[(u.dates*2+i),5] <- u.singleton.ancient.75p[i,4]
   
 plot.uniform[(u.dates*3+i),1] <- "singleton ancient (100% sites)"
 plot.uniform[(u.dates*3+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*3+i),3] <- u.singleton.ancient.100p[i,2]
 plot.uniform[(u.dates*3+i),4] <- u.singleton.ancient.100p[i,3]
 plot.uniform[(u.dates*3+i),5] <- u.singleton.ancient.100p[i,4]
   
 plot.uniform[(u.dates*4+i),1] <- "singleton recent (50% sites)"
 plot.uniform[(u.dates*4+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*4+i),3] <- u.singleton.recent.50p[i,2]
 plot.uniform[(u.dates*4+i),4] <- u.singleton.recent.50p[i,3]
 plot.uniform[(u.dates*4+i),5] <- u.singleton.recent.50p[i,4]
   
 plot.uniform[(u.dates*5+i),1] <- "singleton recent (75% sites)"
 plot.uniform[(u.dates*5+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*5+i),3] <- u.singleton.recent.75p[i,2]
 plot.uniform[(u.dates*5+i),4] <- u.singleton.recent.75p[i,3]
 plot.uniform[(u.dates*5+i),5] <- u.singleton.recent.75p[i,4]
   
 plot.uniform[(u.dates*6+i),1] <- "singleton recent (100% sites)"
 plot.uniform[(u.dates*6+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*6+i),3] <- u.singleton.recent.100p[i,2]
 plot.uniform[(u.dates*6+i),4] <- u.singleton.recent.100p[i,3]
 plot.uniform[(u.dates*6+i),5] <- u.singleton.recent.100p[i,4]
   
 plot.uniform[(u.dates*7+i),1] <- "singleton random (50% sites)"
 plot.uniform[(u.dates*7+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*7+i),3] <- u.singleton.random.50p[i,2]
 plot.uniform[(u.dates*7+i),4] <- u.singleton.random.50p[i,3]
 plot.uniform[(u.dates*7+i),5] <- u.singleton.random.50p[i,4]
   
 plot.uniform[(u.dates*8+i),1] <- "singleton random (75% sites)"
 plot.uniform[(u.dates*8+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*8+i),3] <- u.singleton.random.75p[i,2]
 plot.uniform[(u.dates*8+i),4] <- u.singleton.random.75p[i,3]
 plot.uniform[(u.dates*8+i),5] <- u.singleton.random.75p[i,4]
   
 plot.uniform[(u.dates*9+i),1] <- "singleton random (100% sites)"
 plot.uniform[(u.dates*9+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*9+i),3] <- u.singleton.random.100p[i,2]
 plot.uniform[(u.dates*9+i),4] <- u.singleton.random.100p[i,3]
 plot.uniform[(u.dates*9+i),5] <- u.singleton.random.100p[i,4]
   
 plot.uniform[(u.dates*10+i),1] <- "bracketed (50% sites)"
 plot.uniform[(u.dates*10+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*10+i),3] <- u.bracketed.50p[i,2]
 plot.uniform[(u.dates*10+i),4] <- u.bracketed.50p[i,3]
 plot.uniform[(u.dates*10+i),5] <- u.bracketed.50p[i,4]

 plot.uniform[(u.dates*11+i),1] <- "bracketed (75% sites)"
 plot.uniform[(u.dates*11+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*11+i),3] <- u.bracketed.75p[i,2]
 plot.uniform[(u.dates*11+i),4] <- u.bracketed.75p[i,3]
 plot.uniform[(u.dates*11+i),5] <- u.bracketed.75p[i,4]
   
 plot.uniform[(u.dates*12+i),1] <- "bracketed (100% sites)"
 plot.uniform[(u.dates*12+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*12+i),3] <- u.bracketed.100p[i,2]
 plot.uniform[(u.dates*12+i),4] <- u.bracketed.100p[i,3]
 plot.uniform[(u.dates*12+i),5] <- u.bracketed.100p[i,4]
   
 plot.uniform[(u.dates*13+i),1] <- "uniform (50% sites)"
 plot.uniform[(u.dates*13+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*13+i),3] <- u.uniform.50p[i,2]
 plot.uniform[(u.dates*13+i),4] <- u.uniform.50p[i,3]
 plot.uniform[(u.dates*13+i),5] <- u.uniform.50p[i,4]
 
 plot.uniform[(u.dates*14+i),1] <- "uniform (75% sites)"
 plot.uniform[(u.dates*14+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*14+i),3] <- u.uniform.75p[i,2]
 plot.uniform[(u.dates*14+i),4] <- u.uniform.75p[i,3]
 plot.uniform[(u.dates*14+i),5] <- u.uniform.75p[i,4]
 
 plot.uniform[(u.dates*15+i),1] <- "uniform (100% sites)"
 plot.uniform[(u.dates*15+i),2] <- u.calBP[i,2]/1000
 plot.uniform[(u.dates*15+i),3] <- u.uniform.100p[i,2]
 plot.uniform[(u.dates*15+i),4] <- u.uniform.100p[i,3]
 plot.uniform[(u.dates*15+i),5] <- u.uniform.100p[i,4]
   
 plot.uniform[(u.dates*16+i),1] <- "baseline"
 plot.uniform[(u.dates*16+i),2] <- u.baseline1[i,2]/1000
 plot.uniform[(u.dates*16+i),3] <- 0
 plot.uniform[(u.dates*16+i),4] <- u.baseline1[i,3]
 plot.uniform[(u.dates*16+i),5] <- 0
}
 
write.csv(plot.uniform, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-COMBINED.csv")

l.dates <- nrow(l.calBP)
plot.linear <- data.frame(matrix(NA, nrow = (l.dates*17), ncol = 5))
names(plot.linear) <- c("sample", "years.ka", "lowerCI", "mean", "upperCI")
for (i in 1:nrow(l.calBP)){
  plot.linear[i,1] <- "baseline (replicates)"
  plot.linear[i,2] <- l.calBP[i,2]/1000
  plot.linear[i,3] <- l.baseline[i,2]
  plot.linear[i,4] <- l.baseline[i,3]
  plot.linear[i,5] <- l.baseline[i,4]
  
  plot.linear[(l.dates+i),1] <- "singleton ancient (50% sites)"
  plot.linear[(l.dates+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates+i),3] <- l.singleton.ancient.50p[i,2]
  plot.linear[(l.dates+i),4] <- l.singleton.ancient.50p[i,3]
  plot.linear[(l.dates+i),5] <- l.singleton.ancient.50p[i,4]
  
  plot.linear[(l.dates*2+i),1] <- "singleton ancient (75% sites)"
  plot.linear[(l.dates*2+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*2+i),3] <- l.singleton.ancient.75p[i,2]
  plot.linear[(l.dates*2+i),4] <- l.singleton.ancient.75p[i,3]
  plot.linear[(l.dates*2+i),5] <- l.singleton.ancient.75p[i,4]
  
  plot.linear[(l.dates*3+i),1] <- "singleton ancient (100% sites)"
  plot.linear[(l.dates*3+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*3+i),3] <- l.singleton.ancient.100p[i,2]
  plot.linear[(l.dates*3+i),4] <- l.singleton.ancient.100p[i,3]
  plot.linear[(l.dates*3+i),5] <- l.singleton.ancient.100p[i,4]
  
  plot.linear[(l.dates*4+i),1] <- "singleton recent (50% sites)"
  plot.linear[(l.dates*4+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*4+i),3] <- l.singleton.recent.50p[i,2]
  plot.linear[(l.dates*4+i),4] <- l.singleton.recent.50p[i,3]
  plot.linear[(l.dates*4+i),5] <- l.singleton.recent.50p[i,4]
  
  plot.linear[(l.dates*5+i),1] <- "singleton recent (75% sites)"
  plot.linear[(l.dates*5+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*5+i),3] <- l.singleton.recent.75p[i,2]
  plot.linear[(l.dates*5+i),4] <- l.singleton.recent.75p[i,3]
  plot.linear[(l.dates*5+i),5] <- l.singleton.recent.75p[i,4]
  
  plot.linear[(l.dates*6+i),1] <- "singleton recent (100% sites)"
  plot.linear[(l.dates*6+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*6+i),3] <- l.singleton.recent.100p[i,2]
  plot.linear[(l.dates*6+i),4] <- l.singleton.recent.100p[i,3]
  plot.linear[(l.dates*6+i),5] <- l.singleton.recent.100p[i,4]
  
  plot.linear[(l.dates*7+i),1] <- "singleton random (50% sites)"
  plot.linear[(l.dates*7+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*7+i),3] <- l.singleton.random.50p[i,2]
  plot.linear[(l.dates*7+i),4] <- l.singleton.random.50p[i,3]
  plot.linear[(l.dates*7+i),5] <- l.singleton.random.50p[i,4]
  
  plot.linear[(l.dates*8+i),1] <- "singleton random (75% sites)"
  plot.linear[(l.dates*8+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*8+i),3] <- l.singleton.random.75p[i,2]
  plot.linear[(l.dates*8+i),4] <- l.singleton.random.75p[i,3]
  plot.linear[(l.dates*8+i),5] <- l.singleton.random.75p[i,4]
  
  plot.linear[(l.dates*9+i),1] <- "singleton random (100% sites)"
  plot.linear[(l.dates*9+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*9+i),3] <- l.singleton.random.100p[i,2]
  plot.linear[(l.dates*9+i),4] <- l.singleton.random.100p[i,3]
  plot.linear[(l.dates*9+i),5] <- l.singleton.random.100p[i,4]
  
  plot.linear[(l.dates*10+i),1] <- "bracketed (50% sites)"
  plot.linear[(l.dates*10+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*10+i),3] <- l.bracketed.50p[i,2]
  plot.linear[(l.dates*10+i),4] <- l.bracketed.50p[i,3]
  plot.linear[(l.dates*10+i),5] <- l.bracketed.50p[i,4]
  
  plot.linear[(l.dates*11+i),1] <- "bracketed (75% sites)"
  plot.linear[(l.dates*11+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*11+i),3] <- l.bracketed.75p[i,2]
  plot.linear[(l.dates*11+i),4] <- l.bracketed.75p[i,3]
  plot.linear[(l.dates*11+i),5] <- l.bracketed.75p[i,4]
  
  plot.linear[(l.dates*12+i),1] <- "bracketed (100% sites)"
  plot.linear[(l.dates*12+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*12+i),3] <- l.bracketed.100p[i,2]
  plot.linear[(l.dates*12+i),4] <- l.bracketed.100p[i,3]
  plot.linear[(l.dates*12+i),5] <- l.bracketed.100p[i,4]
  
  plot.linear[(l.dates*13+i),1] <- "uniform (50% sites)"
  plot.linear[(l.dates*13+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*13+i),3] <- l.uniform.50p[i,2]
  plot.linear[(l.dates*13+i),4] <- l.uniform.50p[i,3]
  plot.linear[(l.dates*13+i),5] <- l.uniform.50p[i,4]
  
  plot.linear[(l.dates*14+i),1] <- "uniform (75% sites)"
  plot.linear[(l.dates*14+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*14+i),3] <- l.uniform.75p[i,2]
  plot.linear[(l.dates*14+i),4] <- l.uniform.75p[i,3]
  plot.linear[(l.dates*14+i),5] <- l.uniform.75p[i,4]
  
  plot.linear[(l.dates*15+i),1] <- "uniform (100% sites)"
  plot.linear[(l.dates*15+i),2] <- l.calBP[i,2]/1000
  plot.linear[(l.dates*15+i),3] <- l.uniform.100p[i,2]
  plot.linear[(l.dates*15+i),4] <- l.uniform.100p[i,3]
  plot.linear[(l.dates*15+i),5] <- l.uniform.100p[i,4]
  
  plot.linear[(l.dates*16+i),1] <- "baseline"
  plot.linear[(l.dates*16+i),2] <- l.baseline1[i,2]/1000
  plot.linear[(l.dates*16+i),3] <- 0
  plot.linear[(l.dates*16+i),4] <- l.baseline1[i,3]
  plot.linear[(l.dates*16+i),5] <- 0
}

write.csv(plot.linear, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-COMBINED.csv")

e.dates <- nrow(e.calBP)
plot.exponential <- data.frame(matrix(NA, nrow = (e.dates*17), ncol = 5))
names(plot.exponential) <- c("sample", "years.ka", "lowerCI", "mean", "upperCI")
for (i in 1:nrow(e.calBP)){
  plot.exponential[i,1] <- "baseline (replicates)"
  plot.exponential[i,2] <- e.calBP[i,2]/1000
  plot.exponential[i,3] <- e.baseline[i,2]
  plot.exponential[i,4] <- e.baseline[i,3]
  plot.exponential[i,5] <- e.baseline[i,4]
  
  plot.exponential[(e.dates+i),1] <- "singleton ancient (50% sites)"
  plot.exponential[(e.dates+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates+i),3] <- e.singleton.ancient.50p[i,2]
  plot.exponential[(e.dates+i),4] <- e.singleton.ancient.50p[i,3]
  plot.exponential[(e.dates+i),5] <- e.singleton.ancient.50p[i,4]
  
  plot.exponential[(e.dates*2+i),1] <- "singleton ancient (75% sites)"
  plot.exponential[(e.dates*2+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*2+i),3] <- e.singleton.ancient.75p[i,2]
  plot.exponential[(e.dates*2+i),4] <- e.singleton.ancient.75p[i,3]
  plot.exponential[(e.dates*2+i),5] <- e.singleton.ancient.75p[i,4]
  
  plot.exponential[(e.dates*3+i),1] <- "singleton ancient (100% sites)"
  plot.exponential[(e.dates*3+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*3+i),3] <- e.singleton.ancient.100p[i,2]
  plot.exponential[(e.dates*3+i),4] <- e.singleton.ancient.100p[i,3]
  plot.exponential[(e.dates*3+i),5] <- e.singleton.ancient.100p[i,4]
  
  plot.exponential[(e.dates*4+i),1] <- "singleton recent (50% sites)"
  plot.exponential[(e.dates*4+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*4+i),3] <- e.singleton.recent.50p[i,2]
  plot.exponential[(e.dates*4+i),4] <- e.singleton.recent.50p[i,3]
  plot.exponential[(e.dates*4+i),5] <- e.singleton.recent.50p[i,4]
  
  plot.exponential[(e.dates*5+i),1] <- "singleton recent (75% sites)"
  plot.exponential[(e.dates*5+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*5+i),3] <- e.singleton.recent.75p[i,2]
  plot.exponential[(e.dates*5+i),4] <- e.singleton.recent.75p[i,3]
  plot.exponential[(e.dates*5+i),5] <- e.singleton.recent.75p[i,4]
  
  plot.exponential[(e.dates*6+i),1] <- "singleton recent (100% sites)"
  plot.exponential[(e.dates*6+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*6+i),3] <- e.singleton.recent.100p[i,2]
  plot.exponential[(e.dates*6+i),4] <- e.singleton.recent.100p[i,3]
  plot.exponential[(e.dates*6+i),5] <- e.singleton.recent.100p[i,4]
  
  plot.exponential[(e.dates*7+i),1] <- "singleton random (50% sites)"
  plot.exponential[(e.dates*7+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*7+i),3] <- e.singleton.random.50p[i,2]
  plot.exponential[(e.dates*7+i),4] <- e.singleton.random.50p[i,3]
  plot.exponential[(e.dates*7+i),5] <- e.singleton.random.50p[i,4]
  
  plot.exponential[(e.dates*8+i),1] <- "singleton random (75% sites)"
  plot.exponential[(e.dates*8+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*8+i),3] <- e.singleton.random.75p[i,2]
  plot.exponential[(e.dates*8+i),4] <- e.singleton.random.75p[i,3]
  plot.exponential[(e.dates*8+i),5] <- e.singleton.random.75p[i,4]
  
  plot.exponential[(e.dates*9+i),1] <- "singleton random (100% sites)"
  plot.exponential[(e.dates*9+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*9+i),3] <- e.singleton.random.100p[i,2]
  plot.exponential[(e.dates*9+i),4] <- e.singleton.random.100p[i,3]
  plot.exponential[(e.dates*9+i),5] <- e.singleton.random.100p[i,4]
  
  plot.exponential[(e.dates*10+i),1] <- "bracketed (50% sites)"
  plot.exponential[(e.dates*10+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*10+i),3] <- e.bracketed.50p[i,2]
  plot.exponential[(e.dates*10+i),4] <- e.bracketed.50p[i,3]
  plot.exponential[(e.dates*10+i),5] <- e.bracketed.50p[i,4]
  
  plot.exponential[(e.dates*11+i),1] <- "bracketed (75% sites)"
  plot.exponential[(e.dates*11+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*11+i),3] <- e.bracketed.75p[i,2]
  plot.exponential[(e.dates*11+i),4] <- e.bracketed.75p[i,3]
  plot.exponential[(e.dates*11+i),5] <- e.bracketed.75p[i,4]
  
  plot.exponential[(e.dates*12+i),1] <- "bracketed (100% sites)"
  plot.exponential[(e.dates*12+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*12+i),3] <- e.bracketed.100p[i,2]
  plot.exponential[(e.dates*12+i),4] <- e.bracketed.100p[i,3]
  plot.exponential[(e.dates*12+i),5] <- e.bracketed.100p[i,4]
  
  plot.exponential[(e.dates*13+i),1] <- "uniform (50% sites)"
  plot.exponential[(e.dates*13+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*13+i),3] <- e.uniform.50p[i,2]
  plot.exponential[(e.dates*13+i),4] <- e.uniform.50p[i,3]
  plot.exponential[(e.dates*13+i),5] <- e.uniform.50p[i,4]
  
  plot.exponential[(e.dates*14+i),1] <- "uniform (75% sites)"
  plot.exponential[(e.dates*14+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*14+i),3] <- e.uniform.75p[i,2]
  plot.exponential[(e.dates*14+i),4] <- e.uniform.75p[i,3]
  plot.exponential[(e.dates*14+i),5] <- e.uniform.75p[i,4]
  
  plot.exponential[(e.dates*15+i),1] <- "uniform (100% sites)"
  plot.exponential[(e.dates*15+i),2] <- e.calBP[i,2]/1000
  plot.exponential[(e.dates*15+i),3] <- e.uniform.100p[i,2]
  plot.exponential[(e.dates*15+i),4] <- e.uniform.100p[i,3]
  plot.exponential[(e.dates*15+i),5] <- e.uniform.100p[i,4]
  
  plot.exponential[(e.dates*16+i),1] <- "baseline"
  plot.exponential[(e.dates*16+i),2] <- e.baseline1[i,2]/1000
  plot.exponential[(e.dates*16+i),3] <- 0
  plot.exponential[(e.dates*16+i),4] <- e.baseline1[i,3]
  plot.exponential[(e.dates*16+i),5] <- 0
}

write.csv(plot.exponential, "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-COMBINED.csv")

# LOAD DATA
plot.uniform <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/uniform_pop_growth-1000sims-plot_data-COMBINED.csv")
plot.linear <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/linear_pop_growth-1000sims-plot_data-COMBINED.csv")
plot.exponential <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/Baseline vs subsamples/Plot data/exponential_pop_growth-1000sims-plot_data-COMBINED.csv")

# PLOT
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
theme_set(
  theme_bw() +
    theme(legend.position = "bottom"))
timeRange = c(12000, 1000)       ## for Holocene sites
group.colors <- c("baseline (replicates)" = "grey", 
                 "singleton ancient (50% sites)" = "chocolate1", "singleton ancient (75% sites)" = "chocolate2", "singleton ancient (100% sites)" = "chocolate3",
                 "singleton recent (50% sites)"  = "cadetblue1", "singleton recent (75% sites)"  = "cadetblue2", "singleton recent (100% sites)"  = "cadetblue3",  
                 "singleton random (50% sites)"  = "goldenrod1", "singleton random (75% sites)"  = "goldenrod2", "singleton random (100% sites)"  = "goldenrod3",
                 "bracketed (50% sites)"         = "chartreuse2", "bracketed (75% sites)"        = "chartreuse3", "bracketed (100% sites)"        = "chartreuse4",
                 "uniform (50% sites)"           = "darkorchid1", "uniform (75% sites"           = "darkorchid3", "uniform (100% sites)"          = "darkorchid4",
                 "baseline" = "black")
u.uniform <- plot.uniform[plot.uniform$sample %in% c("baseline", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
   mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline")) %>%
   ggplot() +
   geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
   geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
   scale_fill_manual(values=group.colors) +
   scale_color_manual(values=group.colors) +
   labs(x = "Years (ka)", y = "Probability") +
   ylim(0, 0.0008) +
   scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
u.singleton.random <- plot.uniform[plot.uniform$sample %in% c("baseline", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
   mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline")) %>%
   ggplot() +
   geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
   geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
   scale_fill_manual(values=group.colors) +
   ylim(0, 0.0008) +
   scale_color_manual(values=group.colors) +
   labs(x = "Years (ka)", y = "Probability") +
   scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
u.singleton.ancient <- plot.uniform[plot.uniform$sample %in% c("baseline", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
   mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline")) %>%
   ggplot() +
   geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
   geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
   scale_fill_manual(values=group.colors) +
   scale_color_manual(values=group.colors) +
   labs(x = "Years (ka)", y = "Probability") +
   ylim(0, 0.0008) +
   scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
u.singleton.recent <- plot.uniform[plot.uniform$sample %in% c("baseline", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
   mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline")) %>%
   ggplot() +
   geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
   geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
   scale_fill_manual(values=group.colors) +
   scale_color_manual(values=group.colors) +
   labs(x = "Years (ka)", y = "Probability") +
   ylim(0, 0.0008) +
   scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
u.bracketed <- plot.uniform[plot.uniform$sample %in% c("baseline", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
   mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline")) %>%
   ggplot() +
   geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
   geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
   scale_fill_manual(values=group.colors) +
   scale_color_manual(values=group.colors) +
   labs(x = "Years (ka)", y = "Probability") +
   ylim(0, 0.0008) +
   scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
u.baseline <- plot.uniform[plot.uniform$sample %in% c("baseline", "baseline (replicates)"),] %>%
   mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
   ggplot() +
   geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
   geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
   scale_fill_manual(values=group.colors) +
   scale_color_manual(values=group.colors) +
   labs(x = "Years (ka)", y = "Probability") +
   ylim(0, 0.0003) +
   scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
uniform <- cowplot::plot_grid(u.uniform + theme(axis.text.x = element_blank(),
                                                axis.ticks.x = element_blank(),
                                                axis.title.x = element_blank(),
                                                legend.position = "none"), 
                              u.singleton.random + theme(axis.text.x = element_blank(),
                                                         axis.ticks.x = element_blank(),
                                                         axis.title.x = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.title.y = element_blank(),
                                                         legend.position = "none"),
                              u.singleton.ancient + theme(axis.text.x = element_blank(),
                                                          axis.ticks.x = element_blank(),
                                                          axis.title.x = element_blank(),
                                                          legend.position = "none"),
                              u.singleton.recent + theme(axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.title.y = element_blank(),
                                                         legend.position = "none"),
                              u.bracketed + theme(legend.position = "none"),
                              nrow = 3,
                              ncol = 2,
                              labels = c("A", "B", "C", "D", "E"),
                              align = "hv")
uniform
ggexport(uniform, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Baseline vs subsamples/uniform_population_growth-baseline_SPD_vs_mean_subsamples-5samples_1000sims.png",
          width = 800,
          height = 1000)

l.uniform <- plot.linear[plot.linear$sample %in% c("baseline", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
l.singleton.random <- plot.linear[plot.linear$sample %in% c("baseline", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  ylim(0, 0.0008) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
l.singleton.ancient <- plot.linear[plot.linear$sample %in% c("baseline", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
l.singleton.recent <- plot.linear[plot.linear$sample %in% c("baseline", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
l.bracketed <- plot.linear[plot.linear$sample %in% c("baseline", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
l.baseline <- plot.linear[plot.linear$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0003) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
linear <- cowplot::plot_grid(l.uniform + theme(axis.text.x = element_blank(),
                                                axis.ticks.x = element_blank(),
                                                axis.title.x = element_blank(),
                                               legend.position = "none"), 
                              l.singleton.random + theme(axis.text.x = element_blank(),
                                                         axis.ticks.x = element_blank(),
                                                         axis.title.x = element_blank(),
                                                         axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.title.y = element_blank(),
                                                         legend.position = "none"),
                              l.singleton.ancient + theme(axis.text.x = element_blank(),
                                                          axis.ticks.x = element_blank(),
                                                          axis.title.x = element_blank(),
                                                          legend.position = "none"),
                              l.singleton.recent + theme(axis.text.y = element_blank(),
                                                         axis.ticks.y = element_blank(),
                                                         axis.title.y = element_blank(),
                                                         legend.position = "none"),
                              l.bracketed + theme(legend.position = "none"),
                              nrow = 3,
                              ncol = 2,
                              labels = c("A", "B", "C", "D", "E"),
                              align = "hv")
linear
ggexport(linear, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Baseline vs subsamples/linear_population_growth-baseline_SPD_vs_mean_subsamples-5samples_1000sims.png",
         width = 800,
         height = 1000)


e.uniform <- plot.exponential[plot.exponential$sample %in% c("baseline", "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "uniform (50% sites)", "uniform (75% sites)", "uniform (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.002) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
e.singleton.random <- plot.exponential[plot.exponential$sample %in% c("baseline", "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton random (50% sites)", "singleton random (75% sites)", "singleton random (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  ylim(0, 0.002) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
e.singleton.ancient <- plot.exponential[plot.exponential$sample %in% c("baseline", "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton ancient (50% sites)", "singleton ancient (75% sites)", "singleton ancient (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.002) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
e.singleton.recent <- plot.exponential[plot.exponential$sample %in% c("baseline", "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "singleton recent (50% sites)", "singleton recent (75% sites)", "singleton recent (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.002) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
e.bracketed <- plot.exponential[plot.exponential$sample %in% c("baseline", "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)"),] %>%
  mutate(sample = fct_relevel(sample, "bracketed (50% sites)", "bracketed (75% sites)", "bracketed (100% sites)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.002) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
e.baseline <- plot.exponential[plot.exponential$sample %in% c("baseline", "baseline (replicates)"),] %>%
  mutate(sample = fct_relevel(sample, "baseline (replicates)", "baseline")) %>%
  ggplot() +
  geom_ribbon(aes(x = years.ka, ymin = lowerCI, ymax = upperCI, fill = sample), alpha = 0.2) +
  geom_line(aes(x = years.ka, y = mean, color = sample), lwd = 1.5) +
  scale_fill_manual(values=group.colors) +
  scale_color_manual(values=group.colors) +
  labs(x = "Years (ka)", y = "Probability") +
  ylim(0, 0.0008) +
  scale_x_reverse(breaks = seq(timeRange[1]/1000, timeRange[2]/1000, by = -2000/1000), limits = c(timeRange[1]/1000, timeRange[2]/1000))
exponential <- cowplot::plot_grid(e.uniform + theme(axis.text.x = element_blank(),
                                                    axis.ticks.x = element_blank(),
                                                    axis.title.x = element_blank(),
                                                    legend.position = "none"), 
                                  e.singleton.random + theme(axis.text.x = element_blank(),
                                                             axis.ticks.x = element_blank(),
                                                             axis.title.x = element_blank(),
                                                             axis.text.y = element_blank(),
                                                             axis.ticks.y = element_blank(),
                                                             axis.title.y = element_blank(),
                                                             legend.position = "none"),
                                  e.singleton.ancient + theme(axis.text.x = element_blank(),
                                                              axis.ticks.x = element_blank(),
                                                              axis.title.x = element_blank(),
                                                              legend.position = "none"),
                                  e.singleton.recent + theme(axis.text.y = element_blank(),
                                                             axis.ticks.y = element_blank(),
                                                             axis.title.y = element_blank(),
                                                             legend.position = "none"),
                                  e.bracketed + theme(legend.position = "none"),
                                  nrow = 3,
                                  ncol = 2,
                                  labels = c("A", "B", "C", "D", "E"),
                                  align = "hv")
exponential
ggexport(exponential, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Baseline vs subsamples/exponential_population_growth-baseline_SPD_vs_mean_subsamples-5samples_1000sims.png",
         width = 800,
         height = 1000)


baseline <- cowplot::plot_grid(u.baseline + theme(axis.title.x = element_blank(),
                                                  legend.position="none"), 
                                  l.baseline + theme(#axis.text.y = element_blank(),
                                                     #axis.ticks.y = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     legend.position="none"),
                                  e.baseline + theme(#axis.text.y = element_blank(),
                                                     #axis.ticks.y = element_blank(),
                                                     axis.title.x = element_blank(),
                                                     axis.title.y = element_blank(),
                                                     legend.position="none"),
                                  nrow = 1,
                                  ncol = 3,
                                  labels = c("A", "B", "C"),
                                  align = "hv")
baseline
legend_bottom <- get_legend(u.baseline + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
baseline_legend <- cowplot::plot_grid(baseline, legend_bottom, ncol = 1, rel_heights = c(1, .1))
baseline_legend
ggexport(baseline_legend, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Baseline vs subsamples/baseline_SPD_vs_mean_subsamples-5samples_1000sims.png",
         width = 1000,
         height = 300)

#-----------------------------------------------------------
# SUBSAMPLE VS THEORETICAL GROWTH MODELS BAR CHART
#-----------------------------------------------------------

# Clear workspace
rm(list=ls()); options(scipen=999,digits=9)

# LOAD DATA
uniform <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/uniform_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")
linear <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/linear_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")
exponential <- read.csv("C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Results/exponential_pop_growth-samples_vs_theoretical_growth_models-5s-1000sims-summary.csv")

# SUBSET DATA
uniform.fd.binned <- subset(uniform, method == "fd.binned")
uniform.fd.kde <- subset(uniform, method == "fd.kde")
uniform.spd <- subset(uniform, method == "spd")
linear.fd.binned <- subset(linear, method == "fd.binned")
linear.fd.kde <- subset(linear, method == "fd.kde")
linear.spd <- subset(linear, method == "spd")
exponential.fd.binned <- subset(exponential, method == "fd.binned")
exponential.fd.kde <- subset(exponential, method == "fd.kde")
exponential.spd <- subset(exponential, method == "spd")

# PLOT
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
theme_set(
  theme_bw() +
    theme(legend.position = "bottom"))
# group.colors <- c("singleton.ancient.50p" = "chocolate1", "singleton.ancient.75p" = "chocolate2", "singleton.ancient.100p" = "chocolate3",
#                   "singleton.recent.50p"  = "cadetblue1", "singleton.recent.75p"  = "cadetblue2", "singleton.recent.100p"  = "cadetblue3",  
#                   "singleton.random.50p"  = "goldenrod1", "singleton.random.75p"  = "goldenrod2", "singleton.random.100p"  = "goldenrod3",
#                   "bracketed.50p"         = "chartreuse2", "bracketed.75p"        = "chartreuse3", "bracketed.100p"        = "chartreuse4",
#                   "uniform.50p"           = "darkorchid1", "uniform.75p"           = "darkorchid3", "uniform.100p"          = "darkorchid4"
#                   )
uniform.fd.barplot <-  uniform.fd.binned[uniform.fd.binned$sampling_type %in% 
                                           c("uniform.50p", "uniform.75p", "uniform.100p", 
                                             "singleton.random.50p", "singleton.random.75p", "singleton.random.100p",
                                             "singleton.ancient.50p", "singleton.ancient.75p", "singleton.ancient.100p",
                                             "singleton.recent.50p", "singleton.recent.75p", "singleton.recent.100p",
                                             "bracketed.50p", "bracketed.75p", "bracketed.100p"),] %>%
  mutate(sampling_type = fct_relevel(sampling_type, "uniform.50p", "uniform.75p", "uniform.100p", 
                              "singleton.random.50p", "singleton.random.75p", "singleton.random.100p",
                              "singleton.ancient.50p", "singleton.ancient.75p", "singleton.ancient.100p",
                              "singleton.recent.50p", "singleton.recent.75p", "singleton.recent.100p",
                              "bracketed.50p", "bracketed.75p", "bracketed.100p")) %>%
  ggplot(uniform.fd.binned, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(legend.position="none") +
  ggtitle("Uniform population growth, method = frequency distribution") +
  labs(x = "Subsampling-type", y = "Prop correct")
uniform.kde.barplot <- ggplot(uniform.fd.kde, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(legend.position="none") +
  ggtitle("Uniform population growth, method = frequency distribution (kde)") +
  labs(x = "Subsampling-type", y = "Prop correct")
uniform.spd.barplot <- ggplot(uniform.spd, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") +
  ggtitle("Uniform population growth, method = spd") +
  ylim(0, 1) +
  labs(x = "Subsampling-type", y = "Prop correct")
uniform.barplot <- cowplot::plot_grid(uniform.fd.barplot + theme(axis.text.x = element_blank(),
                                                                 axis.ticks.x = element_blank(),
                                                                 axis.title.x = element_blank() ), 
                                      uniform.kde.barplot + theme(axis.text.x = element_blank(),
                                                                  axis.ticks.x = element_blank(),
                                                                  axis.title.x = element_blank() ),
                                      uniform.spd.barplot,
                                      nrow = 3,
                                      labels = "auto",
                                      align = "h")
ggexport(uniform.barplot, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/uniform_population_growth-subsets_vs_theoretical_growth_models-5samples_1000sims.png",
         width = 1500,
         height = 1000)

linear.fd.barplot <- ggplot(linear.fd.binned, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(legend.position="none") +
  ggtitle("Linear population growth, method = frequency distribution") +
  labs(x = "Subsampling-type", y = "Prop correct")
linear.kde.barplot <- ggplot(linear.fd.kde, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(legend.position="none") +
  ggtitle("Linear population growth, method = frequency distribution (kde)") +
  labs(x = "Subsampling-type", y = "Prop correct")
linear.spd.barplot <- ggplot(linear.spd, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") +
  ggtitle("Linear population growth, method = spd") +
  ylim(0, 1) +
  labs(x = "Subsampling-type", y = "Prop correct")
linear.barplot <- cowplot::plot_grid(linear.fd.barplot + theme(axis.text.x = element_blank(),
                                                                 axis.ticks.x = element_blank(),
                                                                 axis.title.x = element_blank() ), 
                                      linear.kde.barplot + theme(axis.text.x = element_blank(),
                                                                  axis.ticks.x = element_blank(),
                                                                  axis.title.x = element_blank() ),
                                      linear.spd.barplot,
                                      nrow = 3,
                                      labels = "auto",
                                      align = "h")
ggexport(linear.barplot, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/linear_population_growth-subsets_vs_theoretical_growth_models-5samples_1000sims.png",
         width = 1500,
         height = 1000)


exponential.fd.barplot <- ggplot(exponential.fd.binned, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(legend.position="none") +
  ggtitle("Exponential population growth, method = frequency distribution") +
  ylim(0, 1) +
  labs(x = "Subsampling-type", y = "Prop correct")
exponential.kde.barplot <- ggplot(exponential.fd.kde, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(legend.position="none") +
  ggtitle("Exponential population growth, method = frequency distribution (kde)") +
  ylim(0, 1) +
  labs(x = "Subsampling-type", y = "Prop correct")
exponential.spd.barplot <- ggplot(exponential.spd, aes(x = sampling_type, y = prop_correct, color=sampling_type)) +
  geom_bar(stat="identity", fill="white") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="none") +
  ggtitle("Exponential population growth, method = spd") +
  ylim(0, 1) +
  labs(x = "Subsampling-type", y = "Prop correct")
exponential.barplot <- cowplot::plot_grid(exponential.fd.barplot + theme(axis.text.x = element_blank(),
                                                                 axis.ticks.x = element_blank(),
                                                                 axis.title.x = element_blank() ), 
                                      exponential.kde.barplot + theme(axis.text.x = element_blank(),
                                                                  axis.ticks.x = element_blank(),
                                                                  axis.title.x = element_blank() ),
                                      exponential.spd.barplot,
                                      nrow = 3,
                                      labels = "auto",
                                      align = "h")
ggexport(exponential.barplot, filename = "C:/Users/Bec/Work/Projects/Radiocarbon dating/GitHub/Figures/Subsamples vs theoretical growth models/exponential_population_growth-subsets_vs_theoretical_growth_models-5samples_1000sims.png",
         width = 1500,
         height = 1000)
