# radiocarbon-dating
Investigating the effects of possible sampling biases on the results of common analysis techniques for radiocarbon dating of human occupation.
Documentation by Rebecca Wheatley.

Technical details pertaining to the code and where to find the simulated data and fitted models:
* All code written in R v. 4.3.0
* Simulated raw and calibrated data, frequency distributions and summed probability distributions have all been uploaded to the UTas Research Data Portal in the archived project, "RD-SDBIRDS".

Workflow:
* 1. sampling_bias_in_radiocarbon_dating-simulation_study-[trend]_population_growth-baseline_SPD_vs_mean_subsample_SPDs.R files
* * simulate some data
* * create biased subsamples
* * calibrate all the data sets
* * fit and compare the summed probability distribution for the baseline data set to the mean SPD for each biased subsamples
* 2. sampling_bias_in_radiocarbon_dating-[trend]_population_growth-subsample_FDs_vs_theoretical_growth_models.R files
* * fit and compare biased subsample frequency distributions and SPDs to theoretical growth models (uniform, linear, or exponential) and determine the frequency the underlying growth model can be correctly identified
* 3. sampling_bias_in_radiocarbon_dating-simulation_study-create_subsamples_vs_baseline_plots.R file
* * create plots
