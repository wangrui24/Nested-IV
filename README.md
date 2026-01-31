# Nested-IV
This repository contains all code needed to reproduce the simulation and testing results in the paper "Nested Instrumental Variables Analysis: Switcher Average Treatment Effect, Identification, Efficient Estimation and Generalizability". The preprint can be found here: https://arxiv.org/abs/2405.07102



All experiments are implemented in R and designed to be run on a computing cluster. Each R script corresponds to a specific experiment or test, and each `.sh` file is a submission script for running the corresponding R script on a cluster.

## File Structure and Usage for running on cluster

| R Script                          | Description                                      | Cluster Script                     |
|----------------------------------|--------------------------------------------------|------------------------------------|
| `cluster1.R`                     | Simulation results         | `run1.sh`                          |
| `cluster2.R`                     | Alternative clustering specification             | `run2.sh`                          |
| `cluster_test1.R`                | Hypothesis testing for clustering (linear)       | `run_test1.sh`                     |
| `cluster_test_nonlinear.R`       | Hypothesis testing (nonlinear setting)           | `run_test_nonlinear.sh`            |
| `cluster_KS_test.R`              | Kolmogorov–Smirnov test (linear case)             | `run_test_KS.sh`                   |
| `cluster_KS_test_nonlinear.R`    | Kolmogorov–Smirnov test (nonlinear case)          | `run_test_KS_nonlinear.sh`         |
| `cluster_test_time_compare.R`    | Runtime comparison across methods                | (run locally or via custom `.sh`)  |

