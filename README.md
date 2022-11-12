# Computing Codes for the Paper "Additive Subdistribution Hazards Regression for Competing Risks Data in Case-Cohort Studies"
### Adane F. Wogu, Haolin Li, Shanshan Zhao, Hazel B. Nichols, and Jianwen Cai

## Description

This repository contains computing codes for the paper "Additive Subdistribution Hazards Regression for Competing Risks Data in Case-Cohort Studies". For questions or suggestions about the computing codes, please send an email to haolin@live.unc.edu.

## Folders 

### 1-Data Generation

In this folder, we summarize the computing codes for generating competing risks data following additive subdistribution hazards model in case-cohort studies. The names of the code and the corresponding simulation scenarios in the paper are as follows,

* *data_generation_1.r* - scenario 1; percent censored = 90%; case to noncase ratio = 1:1. 
* *data_generation_2.r* - scenario 1; percent censored = 90%; case to noncase ratio = 1:2. 
* *data_generation_3.r* - scenario 1; percent censored = 90%; case to noncase ratio = 1:3. 
* *data_generation_4.r* - scenario 1; percent censored = 95%; case to noncase ratio = 1:1. 
* *data_generation_5.r* - scenario 1; percent censored = 95%; case to noncase ratio = 1:2. 
* *data_generation_6.r* - scenario 1; percent censored = 95%; case to noncase ratio = 1:3. 
* *data_generation_7.r* - scenario 2; percent censored = 90%; case to noncase ratio = 1:1. 
* *data_generation_8.r* - scenario 2; percent censored = 90%; case to noncase ratio = 1:2. 
* *data_generation_9.r* - scenario 2; percent censored = 90%; case to noncase ratio = 1:3. 
* *data_generation_10.r* - scenario 2; percent censored = 95%; case to noncase ratio = 1:1. 
* *data_generation_11.r* - scenario 2; percent censored = 95%; case to noncase ratio = 1:2. 
* *data_generation_12.r* - scenario 2; percent censored = 95%; case to noncase ratio = 1:3. 

### 2-Analysis

In this folder, we summarize the computing codes analyzing competing risks data from case-cohort studies using the proposed method. The names and descriptions of the files are as follows,

* *1-coefficient.r* - The R code for producing estimates of regression coefficients and standard errors. Note that the "nset" in line 4 and "beta.true" in line 15 needs to be adjusted for each set of simulations.
* *Weighting_rate.r* - The R code for rate of change model (linear regression) using IPW or NRCW. 
* *Weighting_bin.r* - The R code for binary model (logistic regression) using IPW or NRCW. 
* *Weighting_poi.r* - The R code for incidence model (Poisson regression) using IPW or NRCW. 

## References

Cai, J., Zeng, D., Li, H., Butera, N., Baldoni, P., Maitra, P., & Dong, L.(2021+). Comparisons of Statistical Methods for Handling Attrition in a Follow-up Visit with Complex Survey Sampling. Manuscript Submitted for Publication.
