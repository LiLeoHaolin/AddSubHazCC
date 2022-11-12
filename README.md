# Computing Codes for the Paper "Additive Subdistribution Hazards Regression for Competing Risks Data in Case-Cohort Studies"
### Adane F. Wogu, Haolin Li, Shanshan Zhao, Hazel B. Nichols, and Jianwen Cai

## Description

This repository contains computing codes for the paper "Additive Subdistribution Hazards Regression for Competing Risks Data in Case-Cohort Studies". For questions or suggestions about the computing codes, please send an email to haolin@live.unc.edu.

## Folders 

### 1-Data Generation

In this folder, we summarize the computing codes for generating population and sample data sets and creating IPW and NRW adjusted weights. The names and descriptions of the files are as follows,

* *Generate_Population.r* - The R code for generating the population data set. 
* *Generate_Samples.r* - The R code for generating the sample data sets and creating NRCW adjusted weights. 
* *IPW.sas* - The SAS code creating IPW adjusted weights.

### 2-Weighting-Based Approaches

In this folder, we summarize the computing codes for the weighting-based approaches, which include IPW and NRCW. The names and descriptions of the files are as follows,

* *Weighting_diff.r* - The R code for difference model (linear regression) using IPW or NRCW. 
* *Weighting_rate.r* - The R code for rate of change model (linear regression) using IPW or NRCW. 
* *Weighting_bin.r* - The R code for binary model (logistic regression) using IPW or NRCW. 
* *Weighting_poi.r* - The R code for incidence model (Poisson regression) using IPW or NRCW. 

### 3-Multiple Imputation

In this folder, we summarize the computing codes for MI. The names and descriptions of the files are as follows,

* *MI_diff.sas* - The SAS code for difference model (linear regression) using MI. 
* *MI_rate.sas* - The SAS code for rate of change model (linear regression) using MI. 
* *MI_bin.sas* - The SAS code for binary model (logistic regression) using MI. 
* *MI_poi.r* - The R code for incidence model (Poisson regression) using MI.

### 4-Full Information Maximum Likelihood 

In this folder, we summarize the computing codes for FIML. The names and descriptions of the files are as follows,

* *FIML_diff_full.inp* - The Mplus code for the full difference model (linear regression) using FIML.
* *FIML_diff_reduced.inp* - The Mplus code for the reduced difference model (linear regression) using FIML.
* *FIML_rate_full.inp* - The Mplus code for the full rate of change model (linear regression) using FIML.
* *FIML_rate_reduced.inp* - The Mplus code for the reduced rate of change model (linear regression) using FIML.

## References

Cai, J., Zeng, D., Li, H., Butera, N., Baldoni, P., Maitra, P., & Dong, L.(2021+). Comparisons of Statistical Methods for Handling Attrition in a Follow-up Visit with Complex Survey Sampling. Manuscript Submitted for Publication.
