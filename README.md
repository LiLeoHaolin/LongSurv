# Computing Code for the Paper "Design and Analysis of Clinical Trials with Survival Outcome by Incorporating Pre-Randomization Longitudinal Biomarkers"

## Description

This repository contains computing code for the paper "Design and Analysis of Clinical Trials with Survival Outcome by Incorporating Pre-Randomization Longitudinal Biomarkers". 

## Naming Convention 

### Files

We include the following computing code for model fitting based on the proposed robust inference as well as the sample size/power calculation using the illustrative example in the manuscript.

* *data_generation.r* - data generation in simulation studies.
* *Analysis.r* - model fitting based on the proposed robust inference. Note that this will call RcppFunctions.cpp while running the code.
* *power_overall.r* - power calculation for overall treatment effect in the illustrative example.
* *power_interaction.r* - power calculation for treatment-biomarker interaction effect in the illustrative example. 

## References

Li, H., Lu, W., Wang, X., Zeng, D., & Cai, J. (2025). Design and Analysis of Clinical Trials with Survival Outcome by Incorporating Pre-Randomization Longitudinal Biomarkers. Manuscript under review.
