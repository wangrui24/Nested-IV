# Nested-IV
This repository contains all code needed to reproduce the simulation and testing results in the paper "Nested Instrumental Variables Analysis: Switcher Average Treatment Effect, Identification, Efficient Estimation and Generalizability". The preprint can be found here: https://arxiv.org/abs/2405.07102



All experiments are implemented in R and designed to be run on a computing cluster. Each R script corresponds to a specific experiment or test, and each `.sh` file is a submission script for running the corresponding R script on a cluster.

## File Structure

| R Script                          | Description                                      | Cluster Script                     |
|----------------------------------|--------------------------------------------------|------------------------------------|
| `cluster1.R`                     | Simulate results for estimation part with continuous outcome        | `run1.sh`                          |
| `cluster2.R`                     | Simulate results for estimation part with binary outcome               | `run2.sh`                          |
| `cluster_test1.R`                | Projection test with linear alternative hypothesis     | `run_test1.sh`                     |
| `cluster_test_nonlinear.R`       | Projection test with nonlinear alternative hypothesis          | `run_test_nonlinear.sh`            |
| `cluster_KS_test.R`              | Nonparametric test with linear alternative hypothesis            | `run_test_KS.sh`                   |
| `cluster_KS_test_nonlinear.R`    | Nonparametric test with nonlinear alternative hypothesis          | `run_test_KS_nonlinear.sh`         |
| `cluster_test1_model_misspecification.R`| Projection test with invalid instrumental variables | `run_test1_misspecification.sh`|
| `cluster_test_time_compare.R`    | Runtime comparison between the projection test and nonparametric test               | (run locally or via custom `.sh`)  |
|`data_clean.R`|clean the PLCO data| (run locally)|
|`real_data_2_biniary`|run the real data analysis| (run locally)|
|`simu_sum`| generate tables and graphs using the simulation output from `cluster1.R` and `cluster2.R` |(run locally)|
|`simu_test_sum`| generate tables and graphs using the simulation output from `cluster_test1.R` and `cluster_KS_test.R` |(run locally)|
|`simu_test_sum_nonlinear`| generate tables and graphs using the simulation output from `cluster_KS_test_nonlinear.R` |(run locally)|
|`simu_test_sum_mis`| generate tables and graphs using the simulation output from `cluster_test1_model_misspecification.R` |(run locally)|


## Additional notes
For replicating the simulation results, please replace the input and output paths in each R script with your own input and output paths.

## A simple illustraing example
The file Simple_example provides a simple synthetic data example illustrating the super-learner–based nonparametric estimator, the projection test, and the nonparametric test. The dataset is generated with a sample size of 1,000. This R script can be run locally and completes within a short amount of time.
