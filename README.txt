# TO-DO: 

- explain 'user-defined theoretical min/max'
- comment on non-perfect reconstructed correlations (has mainly to do with categorical variables; see 'no longer needed' code in seperate file)
- limitation = correlation on normalized scale is not the same as correlation on original scale
- limitation = only correlations, no functions with non-linear associations or interactions etc. 
- personal remark: if eigenvalues not positie definite check for function to project to nearest positive definite matrix. Possibly variance is too large because of the partly missing data. change names for round1 and round2
- describe we do not simulate on item level
- indicate we are simulating 2 scales using the same items, they are not necessarily create plausible combinations
- add that CDR-global should be rounded to full numbers if above 1 (not done here)
- explain our method uses complete cases only (therewith assumes any missing data are completely at random)
- change working directory and make note on this



# Introduction

This code simulates correlated data. It was developed with the purpose to generate synthetic medical trial data for sharing without the risk of privacy violation. The base case replicates an original dataset. However, part of the code can also relatively easily be used to generate a dataset from user-defined beta distribution shape parameters.  

The main steps are as follows:  
- Open an existing dataset.  
- Rescale variables to 0-1 range and fit beta distribution shape parameters.  
- Convert to normal distribution.  
- Estimate correlation structure.  
- Generate random correlated data using correlation structure.  
- Convert to original distribution.  
- Split data in control and intervention arm.  
- Apply hypothetical treatment effect.  

We believe this method has the following limitations. The variance-covariance matrix of the original data likely differs from the variance-covariance matrix of the simulated data if original data is non-normally distributed or of categorical nature. In case of non-normally distributed data, the variance-covariance matrix is fitted after rescaling to a normal distribution, therewith not representing the variance-covariance of the original data on its original non-normal scale (for example, high values in right skewed data have more impact on the covariance as compared to their covariance after rescaling them to a normal distribution). **DOUBLE CHECK!!!** In case of categorical data, the covariance of the simulated data is based on continuous normal distribution which, after categorization, loses information due to categorization leading to a lower covariance in the categorized data.


# Acknowledgment

## Developers: 
- Ron Handels
- Linus Jonsson
- Lars Lau Raket

## Data
The starting dataset is meant to be an existing real-world dataset. For reasons of data protection we created an artificial dataset outside this code based on ADNI data. Random changes have been made on each individual and to mean outcomes. 
Data used in preparation of this work were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). 
