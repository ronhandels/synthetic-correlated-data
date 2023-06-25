
# TO-DO: 
# 
# explain 'user-defined theoretical min/max'
# comment on non-perfect reconstructed correlations (has mainly to do with categorical variables; see 'no longer needed' code in seperate file)
# limitation = correlation on normalized scale is not the same as correlation on original scale
# limitation = only correlations, no functions with non-linear associations or interactions etc. 
# personal remark: if eigenvalues not positie definite check for function to project to nearest positive definite matrix. Possibly variance is too large because of the partly missing data.
# change names for round1 and round2
# describe we do not simulate on item level
# indicate we are simulating 2 scales using the same items, they are not necessarily create plausible combinations
# add that CDR-global should be rounded to full numbers if above 1 (not done here)
# explain our method uses complete cases only (therewith assumes any missing data are completely at random)
# change working directory and make note on this



######################################## README ########################################

# title: "Generate synthetic trial data (simulate correlated data)"
# author: "Ron Handels"
# date: '2022-11-10'
# developers: 
  # Ron Handels
  # Linus Jonsson
  # Lars Lau Raket


# Introduction
# 
# This code simulates correlated data. It was developed with the purpose to generate synthetic medical trial data for sharing without the risk of privacy violation. The base case replicates an original dataset. However, part of the code can also relatively easily be used to generate a dataset from user-defined beta distribution shape parameters.  
# 
# The main steps are as follows:  
# - Open an existing dataset.  
# - Rescale variables to 0-1 range and fit beta distribution shape parameters.  
# - Convert to normal distribution.  
# - Estimate correlation structure.  
# - Generate random correlated data using correlation structure.  
# - Convert to original distribution.  
# - Split data in control and intervention arm.  
# - Apply hypothetical treatment effect.  
# 
# We believe this method has the following limitations. The variance-covariance matrix of the original data likely differs from the variance-covariance matrix of the simulated data if original data is non-normally distributed or of categorical nature. In case of non-normally distributed data, the variance-covariance matrix is fitted after rescaling to a normal distribution, therewith not representing the variance-covariance of the original data on its original non-normal scale (for example, high values in right skewed data have more impact on the covariance as compared to their covariance after rescaling them to a normal distribution). **DOUBLE CHECK!!!** In case of categorical data, the covariance of the simulated data is based on continuous normal distribution which, after categorization, loses information due to categorization leading to a lower covariance in the categorized data.


# Acknowledgment
# 
# The starting dataset is meant to be an existing real-world dataset. For reasons of data protection we created an artificial dataset outside this code based on ADNI data. Random changes have been made on each individual and to mean outcomes. 
# Data used in preparation of this work were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). 



######################################## PREPARE ########################################

cat("\014") # clear console
rm(list = ls()) # clear environment
setwd("C:/users/Ron/surfdrive/PhD/PAPERS/IPECAD workshop 2023/synthetic trial data/github") # working directory; (USER-INPUT REQUIRED: please change to your prefered location)
#setwd("D:/surfdrive/PhD/PAPERS/IPECAD workshop 2023/synthetic trial data/github") # working directory; (USER-INPUT REQUIRED: please change to your prefered location)



######################################## FUNCTIONS ########################################

# Convert categorical variable to cumulative probabilities (NA not used / ignored).
f.cumsumcat <- function(x) {
  frequency <- table(x, useNA = "no") # vector of frequencies of categorical variable
  probability <- frequency / sum(!is.na(x)) # vector of probabilities of categorical variable
  cumulativesum <- cumsum(probability) # cumulative probability
  out <- rep(x = NA, times = length(!is.na(x))) # create outcome vector with same length as categorical variable x
  for (i in 1:length(cumulativesum)) out[x==as.numeric(names(frequency[i]))] <- cumulativesum[i] # write corresponding cumulative probability to each element of categorical variable x
  return(out)
}


# Estimate beta distribution parameters (method of moments is used: https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance )
f.betamm <- function(mu, var) {
  if (var > (mu*(1-mu))) stop("variance > mu*(1-mu)")
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2 # shape1
  beta <- alpha * (1 / mu - 1) # shape2
  return(params = list(alpha = alpha, beta = beta))
}


# Rescale a vector to 0-1 range. 
f.rescale01 <- function(x, min, max) { (x - min) / (max-min) }


# Back-rescale a vector from 0-1 range. 
f.rescale01_back <- function(x, min, max) {
  x * (max-min) + min
}



######################################## OPEN ORIGINAL DATA ########################################

# open original data to be replicated ( alternatively, data can be generated using user-inputs in the 'm.out' in the next steps)
d_orig <- as.matrix(read.csv(file = "original_data.csv", header = TRUE))

# store column names of continuous variables
v.varscon <- c(
  "AGE","EDUCAT",
  "CDRSB0","CDRSB1","CDRSB2",
  "MMSE0","MMSE1","MMSE2",
  "ABETA0","ABETA1","ABETA2"
)
# define categorical variables
v.varscat <- c(
  "SEX",
  "CDRGLOBAL0","CDRGLOBAL1","CDRGLOBAL2"
)



######################################## DESCRIBE ORIGINAL DATA ########################################

# number of variables
n.varscon <- length(v.varscon)
n.varscat <- length(v.varscat)
n.vars <- n.varscon + n.varscat

# variables statistics (mean, min, max, sd, se), user-defined theoretical min/max, number of decimals for rounding (continuous)
m.out <- matrix(data = NA, nrow = 15, ncol = ncol(d_orig[,v.varscon]), dimnames = list(c("n","mean","min","max","sd","se","95%CIlb","95%CIub","skewness","shape1","shape2","min_th","max_th","round1","round2"),colnames(d_orig[,v.varscon]))) # empty outcome matrix
m.out["n",] <- apply(X = d_orig[,v.varscon], MARGIN = 2, FUN = function(x) sum(!is.na(x))) # estimate mean for each variable in the original data
m.out["mean",] <- apply(X = d_orig[,v.varscon], MARGIN = 2, FUN = mean, na.rm = TRUE)
m.out["min",] <- apply(X = d_orig[,v.varscon], MARGIN = 2, FUN = min, na.rm = TRUE)
m.out["max",] <- apply(X = d_orig[,v.varscon], MARGIN = 2, FUN = max, na.rm = TRUE)
m.out["sd",] <- apply(X = d_orig[,v.varscon], MARGIN = 2, FUN = sd, na.rm = TRUE)
m.out["se",] <- m.out["sd",] / sqrt(m.out["n",])
m.out["95%CIlb",] <- m.out["mean",] - 1.96 * m.out["se",]
m.out["95%CIub",] <- m.out["mean",] + 1.96 * m.out["se",]

# user-defined settings: set theoretical min/max (simulated data uses observed min/max, but one might want to use the scale min/max to allow simulating values outside observed min/max; min/max should be outside observed min/max)
m.out["min_th",c("AGE")] <- 55
m.out["max_th",c("AGE")] <- 90
m.out["min_th",c("EDUCAT")] <- 6
m.out["max_th",c("EDUCAT")] <- 21
m.out["min_th",c("CDRSB0","CDRSB1","CDRSB2")] <- 0
m.out["max_th",c("CDRSB0","CDRSB1","CDRSB2")] <- 18
m.out["min_th",c("MMSE0","MMSE1","MMSE2")] <- 0
m.out["max_th",c("MMSE0","MMSE1","MMSE2")] <- 30
m.out["min_th",c("ABETA0","ABETA1","ABETA2")] <- 199
m.out["max_th",c("ABETA0","ABETA1","ABETA2")] <- 1701

# user-defined settings: set rounding
m.out["round1",c("AGE")] <- 1
m.out["round2",c("AGE")] <- 0
m.out["round1",c("EDUCAT")] <- 1
m.out["round2",c("EDUCAT")] <- 0
m.out["round1",c("CDRSB0","CDRSB1","CDRSB2")] <- 0.5
m.out["round2",c("CDRSB0","CDRSB1","CDRSB2")] <- 0
m.out["round1",c("MMSE0","MMSE1","MMSE2")] <- 1
m.out["round2",c("MMSE0","MMSE1","MMSE2")] <- 0
m.out["round1",c("ABETA0","ABETA1","ABETA2")] <- 1
m.out["round2",c("ABETA0","ABETA1","ABETA2")] <- 0



######################################## RESCALE ORIGINAL DATA AND STORE DISTRIBUTION PARAMETERS ########################################

# rescale variables (continuous)
d_orig_rs <- d_orig # create empty structure with same dimensions and variable names (part 1/2)
d_orig_rs[,] <- NA # create empty structure with same dimensions and variable names (part 2/2); these 2 steps are applied multiple times below
for (i in v.varscon) d_orig_rs[,i] <- f.rescale01(x = d_orig[,i], min = m.out["min_th",i], max = m.out["max_th",i])

# estimate beta distribution parameters (continuous)
for (i in v.varscon) {
  m.out["shape1",i] <- f.betamm(mu = mean(d_orig_rs[,i], na.rm = TRUE), var = var(d_orig_rs[,i], na.rm = TRUE))[["alpha"]] # estimate beta distribution parameter
  m.out["shape2",i] <- f.betamm(mu = mean(d_orig_rs[,i], na.rm = TRUE), var = var(d_orig_rs[,i], na.rm = TRUE))[["beta"]] # estimate beta distribution parameter
}

# convert to probability cumulative density function (CDF)
d_orig_cdf <- d_orig_rs
d_orig_cdf[,] <- NA
for (i in v.varscon) d_orig_cdf[,i] <- pbeta(q = d_orig_rs[,i], shape1 = m.out["shape1",i], shape2 = m.out["shape2",i]) # convert to cumulative density function (continuous)
for (i in v.varscat) d_orig_cdf[,i] <- f.cumsumcat(x = d_orig[,i]) # convert to probability (categorical); alternative: pbinom(q = d_origcat[,"SEX"], size = 1, prob = mean(d_origcat[,"SEX"], na.rm = TRUE))
# replace 0 and 1 to near 0 and near 1
d_orig_cdf2 <- d_orig_cdf
d_orig_cdf2[d_orig_cdf2==1] <- 0.9999 # replace 1 to near-1 (otherwise invalid values when back-transformed)
d_orig_cdf2[d_orig_cdf2==0] <- 0.0001 # replace 0 to near-0 (otherwise invalid values when back-transformed)

# convert to normal distribution
d_orig_norm <- d_orig_cdf2
d_orig_norm[,] <- NA
for (i in 1:n.vars) d_orig_norm[,i] <- qnorm(p = d_orig_cdf2[,i], mean = 0, sd = 1) # convert to estimate on normal distribution

# estimate variance-covariance matrix
m.cor_orig_norm <- cor(d_orig_norm, use = "pairwise.complete.obs")



######################################## GENERATE SYNTHETIC DATA ########################################

# size of simulated data
n.sim <- 10000

# simulate from variance-covariance matrix
m.chol <- chol(m.cor_orig_norm) # Cholesky decomposition (for details, see for example https://towardsdatascience.com/behind-the-models-cholesky-decomposition-b61ef17a65fb )
  # eigen(m.cor_orig_norm)$values # optional check
  # isSymmetric(m.cor_orig_norm) # optional check
set.seed(38567)
m.rnorm <- matrix(data = rnorm(n = n.vars * n.sim), nrow = n.sim, dimnames=list(NULL,colnames(m.chol))) # generate independent random data from normal distribution (other column order)
v.sds <- apply(X = d_orig_norm, MARGIN = 2, FUN = sd, na.rm = TRUE) # estimate standard deviations from original data on normal scale
v.sds[v.varscat] <- 1 # replace SD to 1 for categorical variables (to keep normal distribution with sd=1)
v.means <- colMeans(x = d_orig_norm, na.rm = TRUE) # estimate means from original data on normal scale
v.means[v.varscat] <- 0 # replace mean to 0 for categorical variables (to keep normal distribution with mean=1)
d_sim_norm <- t( t(m.chol) %*% t(m.rnorm) * v.sds + v.means ) # generate correlated data from random normal data, standard deviations and means

# backtransform to probability CDF
d_sim_cdf <- d_sim_norm
d_sim_cdf[,] <- NA
for (i in 1:n.vars) d_sim_cdf[,i] <- pnorm(q = d_sim_norm[,i], mean = 0, sd = 1)

# backtransform to beta distribution (continuous)
d_sim_brs <- d_sim_norm
d_sim_brs[,] <- NA
for (i in v.varscon) d_sim_brs[,i] <- qbeta(p = d_sim_cdf[,i], shape1 = m.out["shape1",i], shape2 = m.out["shape2",i])

# transform to original values
d_sim <- d_sim_brs
d_sim[,] <- NA
# back-rescale from 0-1 range (continuous)
for (i in v.varscon) d_sim[,i] <- f.rescale01_back(x = d_sim_brs[,i], min = m.out["min_th",i], max = m.out["max_th",i])
# transform to categories (categorical)
for (i in v.varscat) {
  proportions <- table(d_orig[,i], useNA = "no")/sum(table(d_orig[,i], useNA = "no"))
  out <- cut(x = d_sim_cdf[,i], breaks = cumsum(c(0, proportions)), include.lowest = TRUE, labels = FALSE)
  values <- as.numeric(labels(table(d_orig[,i], useNA = "no"))[[1]])
  d_sim[,i] <- values[out]
}

# rounding
for (i in v.varscon) d_sim[,i] <- round(x = d_sim[,i] * (1/m.out["round1",i]), digits = m.out["round2",i]) * m.out["round1",i]

# save simulated data
write.csv(x=d_sim, file="synthetic_data.csv", row.names=FALSE)



######################################## COMPARE TO ORIGINAL DATA ########################################

# open data
d_sim <- read.csv(file="synthetic_data.csv")

# plot distribution of original and simulated data
uni01 <- seq(0,1,0.001)
mfcol <- ceiling(sqrt(n.vars))
par(mar = c(2,2,2,1), bg = "white", mfcol = c(mfcol,mfcol), cex = 0.40)
for (i in v.varscon) {
  hist(d_orig[,i], main = i, freq = FALSE, xlim = c(m.out["min_th",i],m.out["max_th",i]))
  lines(density(d_orig[,i], na.rm = TRUE), xlim = c(m.out["min_th",i],m.out["max_th",i]))
  lines(density(f.rescale01_back(x = qbeta(p = uni01, shape1 = m.out["shape1",i], shape2 = m.out["shape2",i]), min = m.out["min_th",i], max = m.out["max_th",i] )), xlim = c(m.out["min_th",i],m.out["max_th",i]), lty = 2)
  lines(density(d_sim[,i]), xlim = c(m.out["min_th",i],m.out["max_th",i]), col = "blue")
}
for (i in v.varscat) {
  hist(d_orig[,i], main = i, freq = FALSE, breaks = 12)
  hist(d_sim[,i], main = i, freq = FALSE, breaks = 12, xaxt = "n", yaxt = "n", col = rgb(0,0,1,0.2), add = TRUE)
}

# simulated data: mean and sd over time
m.out_sim <- matrix(data = NA, nrow = 11, ncol = ncol(d_sim), dimnames = list(c("n","mean","min","max","sd","se","95%CIlb","95%CIub","skewness","shape1","shape2"),colnames(d_sim)))
m.out_sim["n",] <- apply(X = d_sim, MARGIN = 2, FUN = function(x) sum(!is.na(x)))
m.out_sim["mean",] <- apply(X = d_sim, MARGIN = 2, FUN = mean, na.rm = TRUE)
m.out_sim["min",] <- apply(X = d_sim, MARGIN = 2, FUN = min, na.rm = TRUE)
m.out_sim["max",] <- apply(X = d_sim, MARGIN = 2, FUN = max, na.rm = TRUE)
m.out_sim["sd",] <- apply(X = d_sim, MARGIN = 2, FUN = sd, na.rm = TRUE)
m.out_sim["se",] <- m.out_sim["sd",] / sqrt(m.out_sim["n",])
m.out_sim["95%CIlb",] <- m.out_sim["mean",] - 1.96 * m.out_sim["se",]
m.out_sim["95%CIub",] <- m.out_sim["mean",] + 1.96 * m.out_sim["se",]

# plot continuous longitudinal variables
par(mar = c(2,2,2,1), bg = "white", mfcol = c(1,2), cex = 0.50)
time <- c(0,1,2)

plot(x = time, y = m.out["mean",c("CDRSB0","CDRSB1","CDRSB2")], type = "l", main = "CDRSB", xaxp = c(0,2,2), ylim = c(1,4))
arrows(x0 = time, y0 = m.out["95%CIlb",c("CDRSB0","CDRSB1","CDRSB2")], x1 = time, y1 = m.out["95%CIub",c("CDRSB0","CDRSB1","CDRSB2")], length = 0.05, angle = 90, code = 3)
lines(x = time, y = m.out_sim["mean",c("CDRSB0","CDRSB1","CDRSB2")], type = "l", main = "CDRSB", xaxp = c(0,18,3), col = "blue", lty = 2)
arrows(x0 = time, y0 = m.out_sim["95%CIlb",c("CDRSB0","CDRSB1","CDRSB2")], x1 = time, y1 = m.out_sim["95%CIub",c("CDRSB0","CDRSB1","CDRSB2")], length = 0.05, angle = 90, code = 3, col = "blue")

plot(x = time, y = m.out["mean",c("MMSE0","MMSE1","MMSE2")], type = "l", main = "MMSE", xaxp = c(0,2,2), ylim = c(22,28))
arrows(x0 = time, y0 = m.out["95%CIlb",c("MMSE0","MMSE1","MMSE2")], x1 = time, y1 = m.out["95%CIub",c("MMSE0","MMSE1","MMSE2")], length = 0.05, angle = 90, code = 3)
lines(x = time, y = m.out_sim["mean",c("MMSE0","MMSE1","MMSE2")], type = "l", main = "MMSE", xaxp = c(0,18,3), col = "blue", lty = 2)
arrows(x0 = time, y0 = m.out_sim["95%CIlb",c("MMSE0","MMSE1","MMSE2")], x1 = time, y1 = m.out_sim["95%CIub",c("MMSE0","MMSE1","MMSE2")], length = 0.05, angle = 90, code = 3, col = "blue")

# benchmark: correlation matrix (normalized scale)
m.cor_sim_norm <- cor(d_sim_norm)
m.cor_orig_norm <- cor(d_orig, use = "pairwise.complete.obs")
m.cor_dif_norm <- abs(m.cor_sim_norm - m.cor_orig_norm) # absolute difference in correlations
m.cor_dif_norm[lower.tri(x=m.cor_dif_norm, diag=TRUE)] <- NA # empty half of the matrix
par(mar = c(2,2,2,2), bg = "white", mfcol = c(1,1), cex = 0.50) # prepare plot
heatmap(m.cor_dif_norm, Rowv = NA, Colv = NA, revC = TRUE, symm = TRUE, col = rev(heat.colors(9))) # plot heatmap
legend(
  x = "top", 
  legend = c(
    paste(round(min(m.cor_dif_norm, na.rm = TRUE),2), "(min)"), 
    paste(round(median(m.cor_dif_norm, na.rm = TRUE),2), "(median)"), 
    paste(round(max(m.cor_dif_norm, na.rm = TRUE),2), "(max)")), 
  fill = rev(heat.colors(9))[c(1,5,9)], 
  cex = 0.70
)

# some other statistics
mean(m.cor_dif_norm, na.rm=TRUE); sd(m.cor_dif_norm, na.rm=TRUE); min(m.cor_dif_norm, na.rm=TRUE); max(m.cor_dif_norm, na.rm=TRUE)
summary(c(m.cor_dif_norm))
quantile(x=m.cor_dif_norm, probs=c(0.5,0.95), na.rm=TRUE)
hist(m.cor_dif_norm)

# benchmark: correlation matrix (original scale)
m.cor_sim <- cor(d_sim)
m.cor_orig <- cor(d_orig, use = "pairwise.complete.obs")
m.cor_dif <- abs(m.cor_sim - m.cor_orig) # absolute difference in correlations
m.cor_dif[lower.tri(x=m.cor_dif, diag=TRUE)] <- NA # empty half of the matrix
par(mar = c(2,2,2,2), bg = "white", mfcol = c(1,1), cex = 0.50) # prepare plot
heatmap(m.cor_dif, Rowv = NA, Colv = NA, revC = TRUE, symm = TRUE, col = rev(heat.colors(9))) # plot heatmap
legend(
  x = "top", 
  legend = c(
    paste(round(min(m.cor_dif, na.rm = TRUE),2), "(min)"), 
    paste(round(median(m.cor_dif, na.rm = TRUE),2), "(median)"), 
    paste(round(max(m.cor_dif, na.rm = TRUE),2), "(max)")), 
  fill = rev(heat.colors(9))[c(1,5,9)], 
  cex = 0.70
)

