# Introduction

This R code generates a synthetic version of an original real-world dataset. It was developed with the purpose to generate synthetic trial data in Alzheimer's disease (AD) for sharing data without the risk of privacy violation. 

The following preparations were done: 
- Change working directory. 

The following steps are done: 
(1) Obtain real-world data (and handle any missing data). 
(2) Estimate (theoretical) minimum and maximum (all continuous variables) and proportions (all categorical variables). 
(3) Rescale to 0-1 range (continuous). 
(4) Estimate beta distribution shape parameters (method of moments; continuous). 
(5) Transform to cumulative density function (using shape parameters; continuous) and to cumulative probability (categorical). 
(6) Converted to a normal distribution. 
(7) Estimate variance-covariance matrix. 
(8) Generate random correlated normal data using Cholesky decomposition of variance-covariance. 
(9) Transform to cumulative density function. 
(10) Transform to inverse cumulative density function of beta distribution (using shape parameters; continuous). 
(11) Rescale to original range (using min/max and proportions from step 2). 

# Limitations

We believe our method has the following limitations: 
- The variance-covariance matrix of the original data likely differs from the variance-covariance matrix of the simulated data if the original data is non-normally distributed or of categorical nature. In case of non-normally distributed data, the variance-covariance matrix is fitted after rescaling to a normal distribution, therewith not representing the variance-covariance of the original data on its original non-normal scale (for example, high values in right skewed data have more impact on the covariance as compared to their covariance after rescaling them to a normal distribution). In case of categorical data, the covariance of the simulated data is based on continuous normal distribution which, after categorization, loses information due to categorization leading to a lower covariance in the categorized data.
- Our method relies on complete cases. Any conditional missing data in the original data should be handled before creating a synthetic version of the original data. 
- Our method relies on correlation between all variables. It might be limited to represent specific non-linear patterns/associations and interactions. 
- We have not tested our method for generating synthetic data on a scale/questionnaire item level. When simulating item levels scores, they are not restricted to the limites of their sum-score. 
- We have simulated 2 outcomes (CDR global and CDR sum of boxes) based on the same underlying scale (CDR). Some combinations of the two outcomes might be inplausible as we did not apply any restriction to their combinations. In addition, CDR-global scores can take values 0, 0.5, 1, 2 or 3. We have rounded them to 0.5 but values other than 0.5 should have been rounded to full numbers. As this is a specific condition of this scale we have not implemented this in our code for reasons to keep our code generic and leave any specific conditions to the user of our code. 


# Acknowledgment

## Developers: 
- Ron Handels
- Linus Jonsson
- Lars Lau Raket

## Data
For our original purpose we used a real-world dataset named ADNI. For reasons of data protection we created an artificial dataset outside this code loosely based on ADNI data. Random changes have been made on each individual and to mean outcomes making the orginal data available within our code disconnected and non-representative from any original ADNI data. We acknowledge the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). 
