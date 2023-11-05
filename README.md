# Introduction

Individual-level data of recent Alzheimer’s Disease (AD) trials are difficult to obtain. Synthetic/simulated data could be used for preparatory, training or explorative research with low risk of privacy breach. 

We aimed to generate a synthetic version of an original real-world observational dataset, and make our method open-source available. 

# Method

Here we explain the steps in the [R code](generate%20correlated%20data.R)

## Original data

1. Obtain original [real-world data](original_data.csv)
-   _demographic_ (age, sex, education), 
-   _clinical_ (cognition: MMSE and ADAS; function: FAQ; composite cognition/function: CDR, 
ADCOMS) and 
-   _biological_ (genetics: APOE4; cerebrospinal fluid: ABeta, Tau; imaging: PET-SUVR-centiloid) 
-   outcomes at baseline, 6, 12 and/or 18-month follow-up (35 variables), with missing data multiple-imputed to obtain 10 sets of 537 individuals. 
2. Estimate (theoretical) minimum and maximum (all continuous variables) and proportions (all categorical variables). 
3. Rescale to 0-1 range (continuous). 
4. Estimate beta distribution shape parameters (method of moments; continuous). 
5. Transform to cumulative density function (using shape parameters; continuous) and to cumulative probability (categorical). 
6. Convert to a normal distribution. 
7. Estimate variance-covariance matrix. 

## Synthetic data
8. Generate random correlated normal data using Cholesky decomposition of variancecovariance. 
9. Transform to cumulative density function. 
10. Transform to inverse cumulative density function of beta distribution (using beta distribution shape parameters; continuous). 
11. Rescale to original range (using minimum and maximum and proportions from step 2).

![image](https://github.com/ronhandels/synthetic-correlated-data/assets/58787973/f86d9cdb-7bb0-4d26-b512-494ce84f5a92)

## More details
See file [poster syntehtic data ISPOR.pdf](poster%20synthetic%20data%20ISPOR%202023.pdf) for details and supporting figures on an application. 

## Limitations

We believe our method has the following limitations: 
- In case of non-normally distributed data, the variance-covariance matrix is fitted after rescaling to a normal distribution, therewith not representing the variance-covariance of the original data on its original non-normal scale (for example, high values in right skewed data have more impact on the covariance as compared to their covariance after rescaling them to a normal distribution). In case of categorical data, the covariance of the simulated data is based on continuous normal distribution which, after categorization, loses information due to categorization leading to a lower covariance in the categorized data. 
- Our method relies on complete cases. Any conditional missing data in the original data should be handled before creating a synthetic version of the original data. 
- Our method relies on correlation between all variables. It might be limited to represent specific non-linear patterns/associations and interactions. 
- We have simulated 2 outcomes based on the same underlying scale ((CDR global and CDR sum of boxes both variants of CDR). Some combinations of the two outcomes might be inplausible as we did not apply any restriction to their combinations. In addition, CDR-global scores can take values 0, 0.5, 1, 2 or 3. We have rounded them to 0.5 but values other than 0.5 should have been rounded to full numbers. As this is a specific condition of this scale we have not implemented this in our code for reasons to keep our code generic and leave any specific conditions to the user of our code.
- We have not compared our method to alternatives (e.g., R package [synthpop](https://cran.r-project.org/web/packages/synthpop/)). 

# Acknowledgment

## Developers: 
- Ron Handels (Maastricht University, Netherlands)
- Linus Jonsson (Karolinska Institutet, Sweden)
- Lars Lau Raket (Lund University, Sweden)

## Data
For our original purpose we used a real-world dataset named ADNI. For reasons of data protection we created an artificial dataset outside this code loosely based on ADNI data. Random changes have been made on each individual and to mean outcomes making the orginal data available within our code disconnected and non-representative from any original ADNI data. We acknowledge the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu) during development of our method. 
