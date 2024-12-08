## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BL22204010)
library(ggplot2)

## -----------------------------------------------------------------------------
 x <- 1:30
x

## -----------------------------------------------------------------------------
gl(2, 6, label=c("Male", "Female"))


## -----------------------------------------------------------------------------
expand.grid(Hight=c(60,80), Width=c(100, 300), sex=c("Male", "Female"))

## -----------------------------------------------------------------------------
ts(1:47, frequency = 12, start = c(1959, 2))

## -----------------------------------------------------------------------------
# Function to generate Rayleigh distributed random samples
generate_rayleigh <- function(n, sigma) {
  # Generate n uniform random variables
  U <- runif(n)
  # Apply the Rayleigh transformation
  X <- sigma * sqrt(-2 * log(U))
  return(X)
}

# Function to check the mode of the generated samples
check_mode <- function(samples, sigma) {
  mode_theoretical <- sigma  # Theoretical mode of Rayleigh distribution
  mode_sample <- density(samples)$x[which.max(density(samples)$y)]  # Sample mode using density estimation
  cat("Theoretical mode: ", mode_theoretical, "\n")
  cat("Sample mode: ", mode_sample, "\n")
}

# Example: Generate Rayleigh samples for different sigma values
set.seed(123)  # For reproducibility
sigma_values <- c(1, 2, 3)  # Different values for sigma
n <- 10000  # Number of samples

# Increase plot margins
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # Set up a 2x2 plotting area with adjusted margins

for (sigma in sigma_values) {
  samples <- generate_rayleigh(n, sigma)
  hist(samples, probability = TRUE, main = paste("Histogram of Rayleigh(", sigma, ")"), 
       xlab = "x", col = "lightblue", breaks = 50)
  check_mode(samples, sigma)
}

# Reset plotting parameters after use
par(mfrow = c(1, 1))  # Reset to single-plot layout


## -----------------------------------------------------------------------------
# Load necessary library
library(ggplot2)

# Function to generate random sample from a mixture of two normal distributions
generate_mixture_sample <- function(n, p1) {
  # Probabilities for components
  p2 <- 1 - p1
  
  # Generate random component choices based on probabilities
  component_choices <- sample(c(1, 2), size = n, replace = TRUE, prob = c(p1, p2))
  
  # Generate samples from N(0, 1) and N(3, 1)
  sample_1 <- rnorm(sum(component_choices == 1), mean = 0, sd = 1)
  sample_2 <- rnorm(sum(component_choices == 2), mean = 3, sd = 1)
  
  # Combine the samples
  mixture_sample <- c(sample_1, sample_2)
  
  return(mixture_sample)
}

# Generate and visualize the mixture sample for p1 = 0.75
set.seed(123)
n <- 1000
p1 <- 0.75
mixture_sample <- generate_mixture_sample(n, p1)

# Plot histogram with density overlay using ggplot2
ggplot(data.frame(x = mixture_sample), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  labs(title = paste("Mixture Distribution (p1 =", p1, ")"), x = "Sample Values", y = "Density") +
  theme_minimal()
# Repeat the process for other values of p1 using base R plots
p1_values <- c(0.25, 0.5, 0.75, 0.9)

# Set up a 2x2 plotting area
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # Adjust margins for better layout

for (p1 in p1_values) {
  mixture_sample <- generate_mixture_sample(n, p1)
  
  # Plot histogram with density overlay using base R
  hist(
    mixture_sample, probability = TRUE, main = paste("p1 =", p1),
    xlab = "Sample Values", col = "lightblue", border = "black", breaks = 30
  )
  lines(density(mixture_sample), col = "red", lwd = 2)
}

# Reset plotting layout to default
par(mfrow = c(1, 1))


## -----------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)

# Function to simulate a Compound Poisson-Gamma process
simulate_compound_poisson_gamma <- function(lambda, alpha, beta, t, n_sim) {
  means <- numeric(n_sim)
  variances <- numeric(n_sim)
  
  # Simulate the process n_sim times
  for (i in 1:n_sim) {
    # Simulate the Poisson process N(t)
    N_t <- rpois(1, lambda * t)  # Poisson process with rate lambda
    
    # Generate jumps from Gamma distribution
    Y <- rgamma(N_t, shape = alpha, rate = beta)
    
    # Compute the sum X(t)
    X_t <- sum(Y)
    
    # Estimate mean and variance
    means[i] <- X_t
    variances[i] <- var(Y) * N_t  # Variance of the sum of N_t gamma variables
  }
  
  # Calculate theoretical mean and variance
  theoretical_mean <- lambda * t * alpha / beta
  theoretical_variance <- lambda * t * (alpha / beta^2 + (alpha / beta)^2)
  
  # Return results
  list(mean_estimate = mean(means),
       variance_estimate = mean(variances),
       theoretical_mean = theoretical_mean,
       theoretical_variance = theoretical_variance)
}

# Set parameters for simulation
lambda <- 2  # Poisson rate
alpha <- 3   # Gamma shape parameter
beta <- 1    # Gamma rate parameter
t <- 10      # Time point to evaluate X(t)
n_sim <- 1000  # Number of simulations

# Run the simulation
results <- simulate_compound_poisson_gamma(lambda, alpha, beta, t, n_sim)

# Display results
cat("Estimated Mean of X(10):", results$mean_estimate, "\n")
cat("Estimated Variance of X(10):", results$variance_estimate, "\n")
cat("Theoretical Mean of X(10):", results$theoretical_mean, "\n")
cat("Theoretical Variance of X(10):", results$theoretical_variance, "\n")


## -----------------------------------------------------------------------------
# Function to estimate the Beta(3, 3) CDF using Monte Carlo simulation
monte_carlo_beta_cdf <- function(x, n_sim) {
  # Generate n_sim samples from Beta(3, 3)
  samples <- rbeta(n_sim, shape1 = 3, shape2 = 3)
  # Estimate the CDF as the proportion of samples <= x
  estimate <- mean(samples <= x)
  return(estimate)
}

# Parameters
n_sim <- 10000  # Number of simulations
x_values <- seq(0.1, 0.9, by = 0.1)  # x values to estimate

# Compute estimates and compare with pbeta
estimates <- sapply(x_values, monte_carlo_beta_cdf, n_sim = n_sim)
pbeta_values <- pbeta(x_values, shape1 = 3, shape2 = 3)

# Results
results <- data.frame(x = x_values, Monte_Carlo_Estimate = estimates, pbeta_Values = pbeta_values)
print(results)


## -----------------------------------------------------------------------------
# Function to generate Rayleigh samples using antithetic variables
generate_rayleigh_antithetic <- function(n, sigma) {
  # Generate n uniform samples
  U <- runif(n)
  
  # Generate antithetic samples
  X1 <- sigma * sqrt(-2 * log(U))
  X2 <- sigma * sqrt(-2 * log(1 - U))  # Antithetic variable
  
  return(list(X1 = X1, X2 = X2))
}

# Parameters
n <- 10000  # Number of samples
sigma <- 1  # Scale parameter for Rayleigh distribution

# Generate Rayleigh samples using antithetic variables
samples <- generate_rayleigh_antithetic(n, sigma)

# Calculate the sum S_antithetic for antithetic variables
S_antithetic <- samples$X1 + samples$X2

# Calculate variances
var_antithetic <- var(S_antithetic)

# Generate independent Rayleigh samples
X1_independent <- sigma * sqrt(-2 * log(runif(n)))
X2_independent <- sigma * sqrt(-2 * log(runif(n)))
S_independent <- X1_independent + X2_independent

# Calculate variance for independent samples
var_independent <- var(S_independent)

# Calculate percent reduction in variance
percent_reduction <- 100 * (var_independent - var_antithetic) / var_independent

# Print results
cat("Variance of X1 + X2 (Independent):", var_independent, "\n")
cat("Variance of X + X2 (Antithetic):", var_antithetic, "\n")
cat("Percent Reduction in Variance:", percent_reduction, "%\n")


## -----------------------------------------------------------------------------
# Load necessary libraries
library(ggplot2)

# Function to measure average sorting time
measure_sort_time <- function(n, simulations = 100) {
  total_time <- 0
  
  for (i in 1:simulations) {
    arr <- sample(1:n)  # Generate random permutation of 1 to n
    start_time <- Sys.time()  # Start time measurement
    sort(arr)  # Sort the array
    total_time <- total_time + as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  }
  
  return(total_time / simulations)  # Return average time
}

# Values of n
n_values <- c(10^4, 2 * 10^4, 4 * 10^4, 6 * 10^4, 8 * 10^4)
average_times <- numeric(length(n_values))

# Run the experiment for each n
for (i in seq_along(n_values)) {
  avg_time <- measure_sort_time(n_values[i])
  average_times[i] <- avg_time
  cat(sprintf("Average sorting time for n=%d: %.4f seconds\n", n_values[i], avg_time))
}

# Calculate t_n = n * log(n)
t_n <- n_values * log(n_values)

# Data frame for regression
results_df <- data.frame(n = n_values, average_time = average_times, t_n = t_n)

# Perform linear regression
model <- lm(average_time ~ t_n, data = results_df)

# Print summary of the regression model
summary(model)

# Plotting the results
ggplot(results_df, aes(x = t_n, y = average_time)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(x = "t_n = n log(n)", y = "Average Sorting Time (seconds)", 
       title = "Average Sorting Time vs. t_n") +
  theme_minimal()

#The plot show regression fitted line


## -----------------------------------------------------------------------------

# Set the random seed for reproducibility
set.seed(42)

# Parameters
num_samples <- 10000  # Number of Monte Carlo simulations
sample_size <- 30     # Size of each sample

# Storage for skewness values
skewness_values <- numeric(num_samples)

# Function to calculate skewness
calculate_skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sd(x)
  skewness <- sum((x - mean_x)^3) / ((n - 1) * sd_x^3)
  return(skewness)
}

# Run Monte Carlo simulation
for (i in 1:num_samples) {
  sample <- rnorm(sample_size)  # Generate normal samples
  skewness <- calculate_skewness(sample)  # Calculate skewness
  skewness_values[i] <- sqrt(abs(skewness))  # Store the square root of the skewness
}

# Estimate quantiles from the Monte Carlo simulation
quantiles_estimated <- quantile(skewness_values, c(0.025, 0.05, 0.95, 0.975))

# Calculate the standard error using the normal approximation
standard_error <- sqrt(6 / sample_size)

# Theoretical quantiles using the normal approximation
quantiles_theoretical <- qnorm(c(0.025, 0.05, 0.95, 0.975), mean = 0, sd = standard_error)

# Display the results
list(
  Estimated_Quantiles = quantiles_estimated,
  Standard_Error = standard_error,
  Theoretical_Quantiles = quantiles_theoretical
)


## -----------------------------------------------------------------------------

# Set the random seed for reproducibility
set.seed(42)

# Parameters
num_samples <- 10000  # Number of simulations
sample_size <- 30     # Sample size for each test
alpha <- 0.05         # Significance level

# Function to calculate the power of a test
calculate_power <- function(p_values, alpha) {
  return(mean(p_values < alpha))
}

# Storage for p-values
p_values_pearson <- numeric(num_samples)
p_values_spearman <- numeric(num_samples)
p_values_kendall <- numeric(num_samples)

# Simulate samples from a bivariate normal distribution
for (i in 1:num_samples) {
  # Generate a bivariate normal sample with correlation 0.5
  x <- rnorm(sample_size)
  y <- 0.5 * x + sqrt(1 - 0.5^2) * rnorm(sample_size)
  
  # Pearson's correlation test
  p_values_pearson[i] <- cor.test(x, y, method = "pearson")$p.value
  
  # Spearman's rank correlation test
  p_values_spearman[i] <- cor.test(x, y, method = "spearman")$p.value
  
  # Kendall's tau test
  p_values_kendall[i] <- cor.test(x, y, method = "kendall")$p.value
}

# Calculate empirical power for each test under bivariate normal distribution
power_pearson <- calculate_power(p_values_pearson, alpha)
power_spearman <- calculate_power(p_values_spearman, alpha)
power_kendall <- calculate_power(p_values_kendall, alpha)

# Display the power of each test
list(
  Power_Pearson = power_pearson,
  Power_Spearman = power_spearman,
  Power_Kendall = power_kendall
)

# Alternative distribution: Bivariate example with a monotonic relationship
p_values_pearson_alt <- numeric(num_samples)
p_values_spearman_alt <- numeric(num_samples)
p_values_kendall_alt <- numeric(num_samples)

for (i in 1:num_samples) {
  # Generate a bivariate distribution with a monotonic but nonlinear relationship
  x <- rnorm(sample_size)
  y <- x^2 + rnorm(sample_size, sd = 0.1)  # Quadratic relationship introduces nonlinearity
  
  # Pearson's correlation test
  p_values_pearson_alt[i] <- cor.test(x, y, method = "pearson")$p.value
  
  # Spearman's rank correlation test
  p_values_spearman_alt[i] <- cor.test(x, y, method = "spearman")$p.value
  
  # Kendall's tau test
  p_values_kendall_alt[i] <- cor.test(x, y, method = "kendall")$p.value
}

# Calculate empirical power for each test under the alternative distribution
power_pearson_alt <- calculate_power(p_values_pearson_alt, alpha)
power_spearman_alt <- calculate_power(p_values_spearman_alt, alpha)
power_kendall_alt <- calculate_power(p_values_kendall_alt, alpha)

# Display the power of each test for the alternative distribution
list(
  Power_Pearson_Alternative = power_pearson_alt,
  Power_Spearman_Alternative = power_spearman_alt,
  Power_Kendall_Alternative = power_kendall_alt
)


## -----------------------------------------------------------------------------
# Data vectors
d1 <- c(-2.961, 0.478, -0.391, -0.869, -0.460, -0.937, 0.779, -1.409, 0.027, -1.569)
d2 <- c(1.608, 1.009, 0.878, 1.600, -0.263, 0.680, 2.280, 2.390, 1.793, 8.091, 1.468)
# Original statistic: mean difference
mean_diff <- mean(d1) - mean(d2)
# Sample standard error of the mean difference
n1 <- length(d1)
n2 <- length(d2)
sample_se <- sqrt(var(d1) / n1 + var(d2) / n2)
# Bootstrap procedure
R <- 10000
boot_diff <- numeric(R)
set.seed(123) 
# Setting a seed for reproducibility
for (i in 1:R) {
boot_sample_d1 <- sample(d1, n1, replace = TRUE)
boot_sample_d2 <- sample(d2, n2, replace = TRUE)
boot_diff[i] <- mean(boot_sample_d1) - mean(boot_sample_d2)
}
# Bootstrap standard error
bootstrap_se <- sd(boot_diff)
# Output the results
cat("Original Mean Difference:", mean_diff, "\n")
cat("Sample Standard Error:", sample_se, "\n")
cat("Bootstrap Standard Error:", bootstrap_se, "\n")
    

## -----------------------------------------------------------------------------
# Load necessary libraries
library(stats)

# Parameters
N <- 1000   # Total number of hypotheses
m <- 10000  # Number of simulation replicates
alpha <- 0.1  # Nominal significance level
num_null <- 950  # Number of null hypotheses
num_alt <- 50   # Number of alternative hypotheses

# Vectors to store FWER, FDR, and TPR for each correction method
fwer_bonf <- numeric(m)
fdr_bonf <- numeric(m)
tpr_bonf <- numeric(m)

fwer_bh <- numeric(m)
fdr_bh <- numeric(m)
tpr_bh <- numeric(m)

# Simulation loop
set.seed(123)  # For reproducibility
for (i in 1:m) {
  # Generate p-values for null hypotheses (uniform distribution)
  pvals_null <- runif(num_null)
  
  # Generate p-values for alternative hypotheses (beta distribution)
  pvals_alt <- rbeta(num_alt, 0.1, 1)
  
  # Combine all p-values
  pvals <- c(pvals_null, pvals_alt)
  
  # Apply Bonferroni correction
  pvals_bonf <- p.adjust(pvals, method = "bonferroni")
  rejections_bonf <- which(pvals_bonf < alpha)
  
  # Apply Benjamini-Hochberg correction
  pvals_bh <- p.adjust(pvals, method = "BH")
  rejections_bh <- which(pvals_bh < alpha)
  
  # Calculate FWER, FDR, and TPR for Bonferroni correction
  num_false_rejections_bonf <- sum(rejections_bonf <= num_null)
  num_true_rejections_bonf <- sum(rejections_bonf > num_null)
  
  fwer_bonf[i] <- ifelse(num_false_rejections_bonf > 0, 1, 0)
  fdr_bonf[i] <- ifelse(length(rejections_bonf) > 0, num_false_rejections_bonf / length(rejections_bonf), 0)
  tpr_bonf[i] <- num_true_rejections_bonf / num_alt
  
  # Calculate FWER, FDR, and TPR for Benjamini-Hochberg correction
  num_false_rejections_bh <- sum(rejections_bh <= num_null)
  num_true_rejections_bh <- sum(rejections_bh > num_null)
  
  fwer_bh[i] <- ifelse(num_false_rejections_bh > 0, 1, 0)
  fdr_bh[i] <- ifelse(length(rejections_bh) > 0, num_false_rejections_bh / length(rejections_bh), 0)
  tpr_bh[i] <- num_true_rejections_bh / num_alt
}

# Calculate average FWER, FDR, and TPR over all simulations for both methods
results <- matrix(c(mean(fwer_bonf), mean(fdr_bonf), mean(tpr_bonf),
                    mean(fwer_bh), mean(fdr_bh), mean(tpr_bh)),
                  nrow = 3, byrow = FALSE)

# Create a table with appropriate row and column names
rownames(results) <- c("FWER", "FDR", "TPR")
colnames(results) <- c("Bonferroni Correction", "B-H Correction")

# Display the results
print("Results:")
print(results)


## -----------------------------------------------------------------------------

# Load the required libraries
library(boot)

# Air-conditioning data
aircondit_data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)

# MLE estimation of lambda (hazard rate) for exponential distribution
mle_lambda <- 1 / mean(aircondit_data)
cat("MLE of lambda:", mle_lambda, "\n")

# Define a function to compute the MLE of lambda from a sample
lambda_mle <- function(data, indices) {
  sample_data <- data[indices]
  return(1 / mean(sample_data))
}

# Bootstrap to estimate the bias and standard error
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data = aircondit_data, statistic = lambda_mle, R = 1000)

# Display bootstrap results
cat("Bootstrap estimate of bias:", mean(bootstrap_results$t) - mle_lambda, "\n")
cat("Bootstrap estimate of standard error:", sd(bootstrap_results$t), "\n")

# Plot the bootstrap distribution
plot(bootstrap_results, main = "Bootstrap Distribution of MLE of Lambda")


## -----------------------------------------------------------------------------

# Load the required libraries
library(boot)

# Air-conditioning data
aircondit_data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)

# Define a function to compute the mean time between failures (1/lambda)
mean_time_failure <- function(data, indices) {
  sample_data <- data[indices]
  return(mean(sample_data))
}

# Perform bootstrap analysis
set.seed(123)  # For reproducibility
bootstrap_results <- boot(data = aircondit_data, statistic = mean_time_failure, R = 1000)

# Compute 95% confidence intervals using different methods
ci_normal <- boot.ci(bootstrap_results, type = "norm")
ci_basic <- boot.ci(bootstrap_results, type = "basic")
ci_percentile <- boot.ci(bootstrap_results, type = "perc")
ci_bca <- boot.ci(bootstrap_results, type = "bca")

# Display the results
cat("95% Bootstrap Confidence Intervals for Mean Time Between Failures (1/lambda):\n")
cat("Standard Normal Method:", ci_normal$normal[2:3], "\n")
cat("Basic Method:", ci_basic$basic[4:5], "\n")
cat("Percentile Method:", ci_percentile$percent[4:5], "\n")
cat("BCa Method:", ci_bca$bca[4:5], "\n")
     

## -----------------------------------------------------------------------------

#Generate or use existing five-dimensional data matrix(example data)
 set.seed(123)
 data_matrix<-matrix(rnorm(100 * 5),ncol=5) #100 observations,5 variables
 #Calculate the sample covariance matrix and eigenvalues for the original data
 cov_matrix<-cov(data_matrix)
 eigenvalues<-eigen(cov_matrix)$values
 eigenvalues<-sort(eigenvalues,decreasing=TRUE)
 #Sample estimate of theta
 theta_hat <-eigenvalues[1]/sum(eigenvalues)
 #Bootstrap Settings
 num_bootstrap<-10000
 theta_bootstrap<-numeric(num_bootstrap)
 #Bootstrap procedure
 for(i in 1:num_bootstrap){
 #Generate a bootstrap sample by resampling rows with replacement
 bootstrap_sample<-data_matrix[sample(1:nrow(data_matrix),replace=TRUE),]
 #Calculate covariance matrix and eigenvalues for the bootstrap sample
 bootstrap_cov<-cov(bootstrap_sample)
 bootstrap_eigenvalues<-eigen(bootstrap_cov)$values
 bootstrap_eigenvalues<-sort(bootstrap_eigenvalues,decreasing=TRUE)
 #Calculate theta for the bootstrap sample
 theta_bootstrap[i]<-bootstrap_eigenvalues[1]/sum(bootstrap_eigenvalues)
 }
 #Bootstrap bias and standard error
 bias_bootstrap<-mean(theta_bootstrap)-theta_hat
 standard_error_bootstrap<-sd(theta_bootstrap)
 #Jackknife procedure
 n<-nrow(data_matrix)
 theta_jackknife<-numeric(n)
 for(i in 1:n){
 #Leave-one-outsample
 jackknife_sample<-data_matrix[-i,]
 #Calculate covariance matrix and eigenvalues for the jackknife sample
 jackknife_cov<-cov(jackknife_sample)
 jackknife_eigenvalues<-eigen(jackknife_cov)$values
 jackknife_eigenvalues<-sort(jackknife_eigenvalues,decreasing=TRUE)
 #Calculate theta for the jackknife sample
 theta_jackknife[i]<-jackknife_eigenvalues[1]/sum(jackknife_eigenvalues)
 }
 #Jackknife bias and standard error
 mean_theta_jackknife<-mean(theta_jackknife)
 bias_jackknife<-(n-1)*(mean_theta_jackknife-theta_hat)
 standard_error_jackknife<-sqrt((n-1)*mean((theta_jackknife-mean_theta_jackknife)^2))
 #Output results
 list(
 theta_hat =theta_hat,
 bootstrap =list(bias=bias_bootstrap,standard_error=standard_error_bootstrap),
 jackknife=list(bias=bias_jackknife,standard_error=standard_error_jackknife)
 )
 

## -----------------------------------------------------------------------------

# Example data (replace with actual data if provided)
 set.seed(123)
 n <- 100
 x <- runif(n, 1, 10)
 y <- 3 + 2 * x- 0.5 * x^2 + 0.05 * x^3 + rnorm(n, 0, 2) # Simulated cubic relationship
 # Prepare data frame
 data <- data.frame(x = x, y = y)
 # Define models
 models <- list(
 linear = function(data) lm(y ~ x, data = data),
 quadratic = function(data) lm(y ~ x + I(x^2), data = data),
 cubic = function(data) lm(y ~ x + I(x^2) + I(x^3), data = data),
 log_linear = function(data) lm(y ~ log(x), data = data)
 )
 # Leave-One-Out Cross-Validation
 loo_cv_mse <- sapply(models, function(model_func) {
 errors <- numeric(n)
 for (i in 1:n) {
 # Leave one out
 train_data <- data[-i, ]
 test_data <- data[i, , drop = FALSE]
 # Fit model and predict left-out observation
 model <- model_func(train_data)
 prediction <- predict(model, newdata = test_data)
 # Calculate squared error for left-out observation
 errors[i] <- (test_data$y- prediction)^2
 }
 # Return mean squared error
 mean(errors)
 })
 # Fit models on entire dataset to get adjusted R^2 values
 adjusted_r_squared <- sapply(models, function(model_func) {
 model <- model_func(data)
 summary(model)$adj.r.squared
 })
 # Find the best model according to LOO-CV MSE and adjusted R^2
 best_model_loo <- names(which.min(loo_cv_mse))
 best_model_adj_r2 <- names(which.max(adjusted_r_squared))
 # Output results
 list(
 loo_cv_mse = loo_cv_mse,
 adjusted_r_squared = adjusted_r_squared,
 best_model_loo = best_model_loo, best_model_adj_r2 = best_model_adj_r2
 )
 

## -----------------------------------------------------------------------------
 # Load the chickwts data and examine its structure
 data(chickwts)
 attach(chickwts)
 attributes(chickwts)
 chickwts$feed
# Define the Cramér-von Mises permutation test function
 Cvm.test_0 <- function(x, y, R = 1000) {
 n <- length(x)
 m <- length(y)
 N <- n + m
 z <- c(x, y)
 # Compute observed Cramér-von Mises statistic
 Fn <- sapply(1:N, function(i) mean(as.integer(z[i] <= x)))
 Gm <- sapply(1:N, function(i) mean(as.integer(z[i] <= y)))
 cv1 <- (n * m) / (n + m)^2 * sum((Fn- Gm)^2)
 # Permutation test
 cv_perm <- replicate(R, {
 d <- sample(1:N)
 # Randomly permute the indices
 z_perm <- z[d]
 x_perm <- z_perm[1:n]
 y_perm <- z_perm[(n + 1):N]
 # Compute the Cramér-von Mises statistic for permuted samples
 Fn_perm <- sapply(1:N, function(i) mean(as.integer(z_perm[i] <= x_perm)))
 Gm_perm <- sapply(1:N, function(i) mean(as.integer(z_perm[i] <= y_perm)))
 (n * m) / (n + m)^2 * sum((Fn_perm- Gm_perm)^2)
 })
 # Calculate p-value
 p_value <- mean(cv_perm >= cv1)
 # Return the observed statistic and p-value
 list(statistic = cv1, p_value = p_value)
 }
 # Extract soybean and linseed weights
 soybean_weights <- chickwts$weight[chickwts$feed == "soybean"]
 linseed_weights <- chickwts$weight[chickwts$feed == "linseed"]
 sunflower_weights <- chickwts$weight[chickwts$feed == "sunflower"]
 linseed_weights <- chickwts$weight[chickwts$feed == "linseed"]
 # Run the test
 result <- Cvm.test_0(soybean_weights, linseed_weights)
 result

result<-Cvm.test_0(sunflower_weights,linseed_weights)
 result

## -----------------------------------------------------------------------------

 #Setup the Spearman rank correlation permutationtest function
 spearman_perm_test<-function(x,y,R=1000){
 #Calculate observed Spearman rank correlation
 observed_corr<-cor(x,y,method="spearman")
 #Initialize a vector to store permutation results
 permuted_corrs<-numeric(R)
 #Generate permutation distribution by shufflingy
 set.seed(123) #for reproducibility
 for(i in 1:R){
 y_permuted<-sample(y)
 permuted_corrs[i]<-cor(x,y_permuted,method="spearman")
 }
 #Calculatep-value for the two-sidedtest
 p_value<-mean(abs(permuted_corrs)>=abs(observed_corr))
 #Return observed correlationand permutationp-value
 list(observed_correlation=observed_corr,perm_p_value=p_value)
 }
 #Generate some example data for x and y
 set.seed(123)
 x <- rnorm(30)
 y <- 0.5 * x + rnorm(30)
 # Run permutation test
 perm_test_result <- spearman_perm_test(x, y)
 # Use cor.test to get the traditional p-value for comparison
 cor_test_result <- cor.test(x, y, method = "spearman")
 # Output both p-values for comparison
 list(
 observed_correlation = perm_test_result$observed_correlation,
 permutation_p_value = perm_test_result$perm_p_value,
 cor_test_p_value = cor_test_result$p.value
 )

## -----------------------------------------------------------------------------


#Metropolis-Hastings Algorithm for Standard Cauchy Distribution}

# Set parameters
set.seed(123)  # For reproducibility
N <- 11000     # Total number of iterations (including burn-in)
burn_in <- 1000
sigma <- 1     # Standard deviation of the proposal distribution

# Initialize
x <- numeric(N)
x[1] <- 0      # Initial value
accept_count <- 0

# Metropolis-Hastings algorithm
for (t in 2:N) {
  # Propose a new value from the proposal distribution (normal distribution)
  x_proposed <- rnorm(1, mean = x[t-1], sd = sigma)
  
  # Compute acceptance probability
  # For the standard Cauchy distribution, f(x) \propto 1 / (1 + x^2)
  f_current <- 1 / (1 + x[t-1]^2)
  f_proposed <- 1 / (1 + x_proposed^2)
  alpha <- min(1, f_proposed / f_current)
  
  # Accept or reject the proposed value
  u <- runif(1)
  if (u <= alpha) {
    x[t] <- x_proposed  # Accept the proposal
    accept_count <- accept_count + 1
  } else {
    x[t] <- x[t-1]      # Reject the proposal
  }
}

# Compute acceptance rate
accept_rate <- accept_count / (N - 1)
cat("Acceptance rate:", accept_rate, "\n") # Discard burn-in samples
x_samples <- x[(burn_in + 1):N]

# Compute deciles of the generated samples
deciles_generated <- quantile(x_samples, probs = seq(0.1, 0.9, by = 0.1))

# Compute theoretical deciles of the standard Cauchy distribution
deciles_theoretical <- qcauchy(seq(0.1, 0.9, by = 0.1))

# Compare the deciles
comparison <- data.frame(
  Decile = seq(0.1, 0.9, by = 0.1),
  Generated = deciles_generated,
  Theoretical = deciles_theoretical
)

print("Comparison of Deciles:") # Plotting the results
par(mfrow = c(1, 2))  # Set up plotting area

# Histogram of the generated samples
hist(x_samples, breaks = 50, probability = TRUE, main = "Histogram of Generated Samples",
     xlab = "Value", xlim = c(-20, 20), col = "lightblue", border = "white")

# Overlay the theoretical standard Cauchy density
curve(dcauchy(x), add = TRUE, col = "red", lwd = 2)

# Q-Q plot to compare the distributions
qqplot(qcauchy(ppoints(length(x_samples))), x_samples,
       main = "Q-Q Plot of Generated Samples vs. Standard Cauchy",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", pch = 20, col = "blue")
abline(0, 1, col = "red", lwd = 2)


## -----------------------------------------------------------------------------
# Set parameters
set.seed(123)    # For reproducibility
N <- 10000       # Total number of iterations
burn_in <- 1000  # Burn-in period
n <- 20          # Parameter of the Binomial distribution
a <- 5           # Beta distribution parameter
b <- 5           # Beta distribution parameter

# Initialize storage for x and y
x <- numeric(N)
y <- numeric(N)

# Initial values
y[1] <- 0.5                      # Starting value for y
x[1] <- rbinom(1, n, y[1])       # Sample initial x from Binomial(n, y[1])

# Gibbs sampler
for (t in 2:N) {
  # Sample x given y
  x[t] <- rbinom(1, n, y[t - 1])
  
  # Sample y given x
  y[t] <- rbeta(1, x[t] + a, n - x[t] + b)
}

# Discard burn-in samples
x_samples <- x[(burn_in + 1):N]
y_samples <- y[(burn_in + 1):N]

# Analyze the samples
mean_x <- mean(x_samples)
var_x <- var(x_samples)
mean_y <- mean(y_samples)
var_y <- var(y_samples)

cat("Estimated mean of x:", mean_x, "\n") 
cat("Estimated variance of x:", var_x, "\n") 
cat("Estimated mean of y:", mean_y, "\n") 
cat("Estimated variance of y:", var_y, "\n") 

# Plotting the results
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # Adjust margins: c(bottom, left, top, right)

# Trace plots
plot(x_samples, type = 'l', main = "Trace Plot of x", xlab = "Iteration", ylab = "x")
plot(y_samples, type = 'l', main = "Trace Plot of y", xlab = "Iteration", ylab = "y")

# Histograms
hist(x_samples, breaks = (n + 1), freq = FALSE, main = "Histogram of x", xlab = "x",
     col = "lightblue", border = "white")
hist(y_samples, breaks = 30, freq = FALSE, main = "Histogram of y", xlab = "y",
     col = "lightgreen", border = "white") 

# Compute the correlation between x and y
cor_xy <- cor(x_samples, y_samples)
cat("Estimated correlation between x and y:", cor_xy, "\n")


## -----------------------------------------------------------------------------



# Load necessary library
library(coda)

# Step 1: Define target density function for standard Cauchy
f <- function(x) {
  return(1 / (pi * (1 + x^2)))
}

# Metropolis-Hastings Algorithm for standard Cauchy
metropolis_hastings <- function(n, sigma, initial) {
  x <- numeric(n)
  x[1] <- initial
  
  for (i in 2:n) {
    # Propose a new value from a normal distribution centered at the current value
    x_new <- rnorm(1, mean = x[i-1], sd = sigma)
    
    # Calculate acceptance ratio
    alpha <- f(x_new) / f(x[i-1])
    
    # Accept or reject the new value
    if (runif(1) < min(1, alpha)) {
      x[i] <- x_new
    } else {
      x[i] <- x[i-1]  # Stay at the current value
    }
  }
  return(x)
}

# Step 2: Run multiple chains and check convergence using Gelman-Rubin
num_chains <- 4
n_samples <- 10000
sigma <- 1
burn_in <- 1000

chains <- list()
for (j in 1:num_chains) {
  chains[[j]] <- metropolis_hastings(n_samples, sigma, rnorm(1))
}

# Convert chains to a matrix for convergence diagnostics
mcmc_chains <- mcmc.list(lapply(chains, mcmc))

# Calculate Gelman-Rubin diagnostic
gelman_diag <- gelman.diag(mcmc_chains, autoburnin = FALSE)
print(gelman_diag) # Step 3: Check if \hat{R} < 1.2
if (all(gelman_diag$psrf[,1] < 1.2)) {
  cat("Chains have converged based on Gelman-Rubin diagnostic.\n")
} else {
  cat("Chains have not yet converged. Increase the number of iterations.\n")
} # Combine chains and discard burn-in samples
combined_samples <- do.call(c, lapply(chains, function(x) x[(burn_in + 1):n_samples]))

# Step 4: Compare deciles with the standard Cauchy distribution
generated_deciles <- quantile(combined_samples, probs = seq(0.1, 0.9, by = 0.1))
theoretical_deciles <- qcauchy(seq(0.1, 0.9, by = 0.1))

# Display comparison of deciles
comparison <- data.frame(
  Decile = seq(0.1, 0.9, by = 0.1),
  Generated = generated_deciles,
  Theoretical = theoretical_deciles
)

print(comparison) # Optional: Plot the generated samples and compare with theoretical deciles
hist(combined_samples, breaks = 50, probability = TRUE, main = "Histogram of Generated Samples vs. Standard Cauchy PDF",
     xlab = "Value")
curve(dcauchy(x), add = TRUE, col = "red", lwd = 2)
points(theoretical_deciles, rep(0, length(theoretical_deciles)), col = "blue", pch = 16)
legend("topright", legend = c("Generated Samples", "Theoretical Deciles"),
       col = c("black", "blue"), pch = c(NA, 16), lty = c(1, NA), lwd = c(2, NA))


## -----------------------------------------------------------------------------


# Load necessary library for the Gelman-Rubin diagnostic
library(coda)

# Parameters
n <- 10          # Number of trials for the binomial
a <- 2           # Alpha parameter for Beta distribution
b <- 3           # Beta parameter for Beta distribution
num_iter <- 10000  # Number of iterations for each Gibbs chain
num_chains <- 4    # Number of chains for the Gelman-Rubin diagnostic
burn_in <- 1000    # Burn-in period to discard initial samples

# Function to run the Gibbs sampler for a single chain
run_gibbs_chain <- function(n, a, b, num_iter) {
  x <- 0          # Initialize x
  y <- 0.5        # Initialize y
  
  x_samples <- numeric(num_iter)
  y_samples <- numeric(num_iter)
  
  for (i in 1:num_iter) {
    # Sample x given y from Binomial(n, y)
    x <- rbinom(1, n, y)
    
    # Sample y given x from Beta(x + a, n - x + b)
    y <- rbeta(1, x + a, n - x + b)
    
    # Store samples
    x_samples[i] <- x
    y_samples[i] <- y
  }
  
  return(list(x_samples = x_samples, y_samples = y_samples))
}

# Run multiple chains and store the results
chains <- vector("list", num_chains)
for (j in 1:num_chains) {
  chains[[j]] <- run_gibbs_chain(n, a, b, num_iter)
}

# Combine the chains into a format suitable for the Gelman-Rubin diagnostic
x_chains <- mcmc.list(lapply(chains, function(chain) mcmc(chain$x_samples)))
y_chains <- mcmc.list(lapply(chains, function(chain) mcmc(chain$y_samples)))

# Calculate Gelman-Rubin diagnostic for both x and y chains
gelman_diag_x <- gelman.diag(x_chains, autoburnin = FALSE)
gelman_diag_y <- gelman.diag(y_chains, autoburnin = FALSE)

cat("Gelman-Rubin Diagnostic for X:\n") # Gelman-Rubin Diagnostic for X:
print(gelman_diag_x) 
cat("\nGelman-Rubin Diagnostic for Y:\n") # Gelman-Rubin Diagnostic for Y:
print(gelman_diag_y)

# Check if both x and y have converged (R_hat < 1.2)
if (all(gelman_diag_x$psrf[, 1] < 1.2) && all(gelman_diag_y$psrf[, 1] < 1.2)) {
  cat("Chains have converged based on Gelman-Rubin diagnostic.\n")
} else {
  cat("Chains have not yet converged. Increase the number of iterations.\n")
}

# Discard burn-in samples and combine samples from all chains
x_samples_combined <- unlist(lapply(chains, function(chain) chain$x_samples[(burn_in + 1):num_iter]))
y_samples_combined <- unlist(lapply(chains, function(chain) chain$y_samples[(burn_in + 1):num_iter]))

# Step 6: Analyze the Samples
# Calculate deciles of the combined samples
x_deciles <- quantile(x_samples_combined, probs = seq(0.1, 0.9, by = 0.1))
y_deciles <- quantile(y_samples_combined, probs = seq(0.1, 0.9, by = 0.1))

# Display deciles
cat("\nDeciles of X Samples:\n") # Deciles of X Samples:
print(x_deciles) 
cat("\nDeciles of Y Samples:\n") # Deciles of Y Samples:
print(y_deciles)

# Optional: Plot the samples 
par(mfrow = c(1, 2))
hist(x_samples_combined, breaks = 30, main = "Histogram of X Samples", xlab = "X")
hist(y_samples_combined, breaks = 30, main = "Histogram of Y Samples", xlab = "Y")



## -----------------------------------------------------------------------------

# Load the necessary library for gamma function
 library(stats)
 # Function to compute the k-th term of the series
 compute_kth_term <- function(a, k, d) {
 # Check if d is a positive integer
 if (d < 1 || (d %% 1) != 0) {
   stop("d must bean integer greater than or equa l to 1.")
 }
 #Calculate the Euclidean norm of vector a
 norm_a<-sqrt(sum(a^2))
 #Compute the k-th term
 term<-((-1)^k /factorial(k)) *(norm_a^(2* k + 2))*
 (gamma((d+ 1)/ 2)/ ((2* k + 1)* (2*k + 2)*
 gamma((k+ 3)/ 2)*gamma(k+(d/ 2)+ 1)))
 return(term)
 }
 #Example usage
 d_value<-2
 k_value<-0
 a_value<-c(1,2)#Corresponding to vector a
 kth_term<-compute_kth_term(a_value,k_value,d_value)
 cat(sprintf("The% d-th term is:%.6f\n",k_value,kth_term))
 ##The 0-th term is:2.500000
 #Function to  compute the sum of the series up to max_k
 compute_series_sum<-function(a,max_k,d){
 #Initializetotal sum
 total_sum<-0
 #Loop through each term from 0 to max_k
 for(k in 0:max_k){
 total_sum<-total_sum +compute_kth_term(a,k,d)
 }
 return(total_sum)
 }
 #Example usage
 max_k_value <-10
 series_sum<-compute_series_sum(a_value,max_k_value,d_value)
 cat(sprintf("The sum of the series up to k=%dis:%.6f\n",max_k_value,series_sum))

## -----------------------------------------------------------------------------

#R code for solving the equation and finding the intersection points of Sk-1(a) and Sk(a)
 library(stats)
 #Function to compute Sk(a)for a given k and a
 evaluate_Sk<-function(k,a){
 if(a^2 >=k +1){
 return(NA) #Avoid invalid values under the square root
 }
 prob<-pt(sqrt((a^2 *k)/ (k+ 1-a^2)),df=k,lower.tail=FALSE)
 return(prob)
 }
 #FunctiontocomputeSk-1(a)foragivenkanda
 evaluate_Sk_minus_1<-function(k,a){
 if(a^2 >=k){
 return(NA) #Avoid invalid values under the square root
 }
 prob<-pt(sqrt((a^2 *(k-1)) / (k-a^2)),df=k-1,lower.tail=FALSE)
 return(prob)
 }
 #Function to find the intersection point A(k) for given k
 find_intersection<-function(k,lower=0.01,upper= sqrt(k) *0.99){
 f<-function(a){
 evaluate_Sk_minus_1(k,a)-evaluate_Sk(k,a)
 }
 result<-tryCatch({
 uniroot(f,lower=lower,upper=upper)$root
 },error=function(e){
 NA #Return NA if the solution cannot be found
 })
 return(result)
 }
 #Calculate intersection points for specified k values
 k_values<-c(4:25,100,500,1000)
 intersection_points<-sapply(k_values,find_intersection)
 #Print the results
 intersection_results<-data.frame(k =k_values,A_k=intersection_points)
 print(intersection_results)
 cat("Intersection points A(k) for k = 4:25, 100, 500, 1000:\n")
  print(intersection_results)
  

## -----------------------------------------------------------------------------

#Given observed data
 observed_data<-c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)
 #Initial estimate of lambda
 lambda_est<-1.0
 tolerance <-1e-6
 max_iterations<-1000
 n<-length(observed_data)
 censored_threshold<-1.0
 #E-Malgorithm
 for(iteration in 1:max_iterations){
 #E-step: Calculate the expected values for censored data
 expected_values<-ifelse(
 observed_data ==censored_threshold,
 censored_threshold +1 / lambda_est,
 observed_data
 )
 #M-step:Update lambda
 new_lambda_est<-n/sum(expected_values)
 #Check for convergence
 if(abs(new_lambda_est-lambda_est)< tolerance){
 break
 }
 lambda_est<-new_lambda_est
 }
 cat(sprintf("Estimated lambda using E-M:%.4f\n",lambda_est))
 ##Estimated lambda using E-M:1.0370
 #Compare with direct MLE
 lambda_mle<-n /sum(observed_data)
 cat(sprintf("MLE estimate ofl ambda:%.4f\n",lambda_mle))

## -----------------------------------------------------------------------------


# List of formulas
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

# List to store model results
lm_results_for <- list()

# Loop through the formulas
for (i in 1:length(formulas)) {
  lm_results_for[[i]] <- lm(formulas[[i]], data = mtcars)
}

# Print summary of models
lapply(lm_results_for, summary)

# Using lapply to apply linear models
lm_results_lapply <- lapply(formulas, function(f) lm(f, data = mtcars))

# Print summary of models
lapply(lm_results_lapply, summary)


## -----------------------------------------------------------------------------


# Function to fit the model mpg ~ disp
fit_model <- function(data) {
  lm(mpg ~ disp, data = data)
}

# Bootstrapping mtcars dataset
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), replace = TRUE)  # Sampling with replacement
  mtcars[rows, ]
})

# List to store model results
lm_results_for <- list()

# Fit the model for each bootstrap replicate using a for loop
for (i in 1:length(bootstraps)) {
  lm_results_for[[i]] <- fit_model(bootstraps[[i]])
}

# Print summaries of models
lapply(lm_results_for, summary)


## -----------------------------------------------------------------------------


# Define the function to extract R^2
rsq <- function(mod) summary(mod)$r.squared

# Assuming lm_results_for is the list of models fitted with a for loop (from previous step)
rsq_results_for <- numeric(length(lm_results_for))

# Extract R^2 for each model
for (i in 1:length(lm_results_for)) {
  rsq_results_for[i] <- rsq(lm_results_for[[i]])
}

# Print R^2 values
print(rsq_results_for)
# Define the function to extract R^2
rsq <- function(mod) summary(mod)$r.squared

# Assuming lm_results_lapply is the list of models fitted with lapply() (from previous step)
rsq_results_lapply <- lapply(lm_results_lapply, rsq)

# Print R^2 values
print(rsq_results_lapply)


## -----------------------------------------------------------------------------


# Simulate 100 trials and perform t-tests
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

# Use sapply() with an anonymous function to extract p-values
p_values <- sapply(trials, function(x) x$p.value)

# Print the p-values
print(p_values)
# Extract p-values using sapply() and [[ directly
p_values_no_anon <- sapply(trials, `[[`, "p.value")

# Print the p-values
print(p_values_no_anon)


## -----------------------------------------------------------------------------

# Custom function that uses Map() and vapply()
parallel_lapply <- function(..., FUN, FUN.VALUE) {
  # Use Map to iterate in parallel over all inputs
  result_list <- Map(FUN, ...)
  
  # Use vapply to ensure the output is stored in a vector (or matrix)
  result <- vapply(result_list, FUN.VALUE = FUN.VALUE, FUN = identity)
  
  return(result)
}

# Example usage:
# Define some example functions and data
x <- 1:5
y <- 6:10
z <- 11:15

# Apply a function (e.g., summing corresponding elements of x, y, z)
parallel_lapply(x, y, z, FUN = function(a, b, c) a + b + c, FUN.VALUE = numeric(1))



## -----------------------------------------------------------------------------

fast_chisq_test <- function(x, y) {
  # Remove missing values
  valid_data <- complete.cases(x, y)
  x <- x[valid_data]
  y <- y[valid_data]
  
  # Calculate observed frequencies (contingency table)
  obs <- table(x, y)
  
  # Calculate row and column totals
  row_totals <- rowSums(obs)
  col_totals <- colSums(obs)
  
  # Calculate expected frequencies under the assumption of independence
  total <- sum(obs)
  expected <- outer(row_totals, col_totals, FUN = "*") / total
  
  # Calculate the chi-square statistic
  chi_square_statistic <- sum((obs - expected)^2 / expected)
  
  return(chi_square_statistic)
}

# Example usage with two numeric vectors
x <- c(1, 2, 1, 2, 1, 2, 3, 3, 2, 3)
y <- c(2, 3, 1, 2, 2, 3, 1, 1, 3, 3)

# Apply the fast chi-square test
chi_square_statistic <- fast_chisq_test(x, y)
print(chi_square_statistic)


## -----------------------------------------------------------------------------

fast_table <- function(x, y) {
  # Ensure both inputs are integer vectors with no missing values
  if (!is.integer(x) || !is.integer(y)) {
    stop("Both x and y must be integer vectors.")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("Both x and y must not contain missing values.")
  }

  # Compute unique pairs (x, y) directly by combining the vectors
  pairs <- cbind(x, y)
  
  # Use a hash-based approach to count the occurrences of each pair
  unique_pairs <- unique(pairs)  # Find unique pairs
  n_pairs <- nrow(unique_pairs)  # Number of unique pairs
  counts <- integer(n_pairs)  # Vector to store counts
  
  # Count each pair's frequency
  for (i in seq_len(n_pairs)) {
    counts[i] <- sum((pairs[,1] == unique_pairs[i,1]) & (pairs[,2] == unique_pairs[i,2]))
  }
  
  # Return the frequency counts as a named vector
  return(setNames(counts, apply(unique_pairs, 1, function(v) paste(v, collapse = "_"))))
}

# Fast chi-square test using optimized table function
fast_chisq_test <- function(x, y) {
  # Check if inputs are integer vectors without missing values
  if (!is.integer(x) || !is.integer(y)) {
    stop("Both inputs must be integer vectors.")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("Both input vectors must not contain missing values.")
  }
  
  # Generate the contingency table using the fast_table() function
  obs <- fast_table(x, y)
  
  # Create the expected table using the standard method
  n <- sum(obs)  # Total number of observations
  row_totals <- sum(obs)  # Row sums (same for one-dimensional table)
  col_totals <- sum(obs)  # Column sums (same for one-dimensional table)
  expected <- row_totals * col_totals / n  # Expected frequency calculation
  
  # Calculate chi-square statistic
  chi_sq_stat <- sum((obs - expected)^2 / expected)
  
  # Return the chi-square statistic
  return(chi_sq_stat)
}

# Example Usage:
x <- as.integer(c(1, 1, 2, 2, 3, 3))  # Make sure the inputs are integer vectors
y <- as.integer(c(1, 2, 1, 2, 1, 2))

fast_chisq_test(x, y)



## -----------------------------------------------------------------------------

# Load necessary library
library(Rcpp)

# Define the Rcpp function using cppFunction
Rcpp::cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List gibbs_sampler(int n, double a, double b, int n_iter) {
  // Initialize variables
  IntegerVector x_chain(n_iter);
  NumericVector y_chain(n_iter);

  // Initial values
  int x = n / 2;  // Arbitrary starting value for x
  double y = 0.5; // Arbitrary starting value for y

  for (int i = 0; i < n_iter; i++) {
    // Sample x | y from Binomial(n, y)
    x = R::rbinom(n, y);
    // Sample y | x from Beta(x + a, n - x + b)
    y = R::rbeta(x + a, n - x + b);
    // Store values in the chains
    x_chain[i] = x;
    y_chain[i] = y;
  }

  // Return the generated chains as a list
  return List::create(
    Named("x_chain") = x_chain,
    Named("y_chain") = y_chain
  );
}
')

# Parameters
n <- 10         # Number of trials in Binomial distribution
a <- 2          # Parameter 'a' for Beta distribution
b <- 2          # Parameter 'b' for Beta distribution
n_iter <- 10000 # Number of iterations for the Gibbs sampler

# Run Gibbs sampler
set.seed(123) # For reproducibility
result <- gibbs_sampler(n, a, b, n_iter)

# Extract chains
x_chain <- result$x_chain
y_chain <- result$y_chain

# Burn-in period (optional)
burn_in <- 1000
x_chain_burned <- x_chain[(burn_in + 1):n_iter]
y_chain_burned <- y_chain[(burn_in + 1):n_iter]

# Plot results
# Trace plot for y
plot(y_chain_burned, type = "l",
     main = "Trace Plot of y",
     xlab = "Iteration",
     ylab = "y")

# Histogram for y
hist(y_chain_burned, breaks = 30,
     main = "Histogram of y",
     xlab = "y",
     probability = TRUE)

# Overlay the theoretical Beta density
curve(dbeta(x, mean(x_chain_burned) + a, (n - mean(x_chain_burned)) + b),
      add = TRUE, col = "red", lwd = 2)

# Summary statistics
cat("Summary of y_chain:\n")
print(summary(y_chain_burned))

# Autocorrelation plot for y
acf(y_chain_burned, main = "Autocorrelation of y Chain")


## -----------------------------------------------------------------------------

 # Load necessary library
library(Rcpp)

# Define the Gibbs sampler function using Rcpp
Rcpp::cppFunction('
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List gibbs_sampler(int n, double a, double b, int n_iter) {
  IntegerVector x_chain(n_iter);
  NumericVector y_chain(n_iter);
  
  // Initial values
  int x = n / 2; // Initial value for x
  double y = 0.5; // Initial value for y
  
  for (int i = 0; i < n_iter; i++) {
    // Sample x | y from Binomial(n, y)
    x = R::rbinom(n, y);
    // Sample y | x from Beta(x + a, n - x + b)
    y = R::rbeta(x + a, n - x + b);
    // Store values in the chains
    x_chain[i] = x;
    y_chain[i] = y;
  }
  
  // Return the generated chains as a list
  return List::create(
    Named("x_chain") = x_chain,
    Named("y_chain") = y_chain
  );
}
')

# Parameters
n <- 10
a <- 2
b <- 2
n_iter <- 10000
burn_in <- 1000

# Run Gibbs sampler
set.seed(123)
result <- gibbs_sampler(n, a, b, n_iter)

# Extract chains and remove burn-in
x_chain <- result$x_chain[(burn_in + 1):n_iter]
y_chain <- result$y_chain[(burn_in + 1):n_iter]

# Generate random numbers using R's built-in functions
set.seed(123)
x_rbinom <- rbinom(length(x_chain), size = n, prob = mean(y_chain))
y_rbeta <- rbeta(length(y_chain), shape1 = mean(x_chain) + a, shape2 = (n - mean(x_chain)) + b)

# QQ plot for x
qqplot(x_rbinom, x_chain, main = "QQ Plot for x", xlab = "rbinom Quantiles", ylab = "Gibbs x Quantiles")
abline(0, 1, col = "red", lwd = 2)

# QQ plot for y
qqplot(y_rbeta, y_chain, main = "QQ Plot for y", xlab = "rbeta Quantiles", ylab = "Gibbs y Quantiles")
abline(0, 1, col = "red", lwd = 2)

 

## -----------------------------------------------------------------------------
 # Load necessary libraries
library(Rcpp)
library(ggplot2)
library(microbenchmark)

# Define the Gibbs sampler function using Rcpp
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix gibbsSamplerRcpp(int n, double a, double b, int N, int burn_in) {
  NumericMatrix samples(N, 2);
  int x = n / 2; // Initial value of x
  double y = 0.5; // Initial value of y
  for (int t = 0; t < N + burn_in; ++t) {
    // Sample x | y
    x = R::rbinom(n, y);
    // Sample y | x
    y = R::rbeta(x + a, n - x + b);
    // Store samples after burn-in period
    if (t >= burn_in) {
      samples(t - burn_in, 0) = x;
      samples(t - burn_in, 1) = y;
    }
  }
  return samples;
}
')

# Gibbs sampler function in R
gibbsSamplerR <- function(n, a, b, N, burn_in) {
  samples <- matrix(0, nrow = N, ncol = 2)
  x <- floor(n / 2) # Initial value of x
  y <- 0.5          # Initial value of y
  for (t in 1:(N + burn_in)) {
    # Sample x | y
    x <- rbinom(1, n, y)
    # Sample y | x
    y <- rbeta(1, x + a, n - x + b)
    # Store samples after burn-in
    if (t > burn_in) {
      samples[t - burn_in, ] <- c(x, y)
    }
  }
  colnames(samples) <- c("x", "y")
  return(samples)
}

# Parameters
n <- 10      # Number of trials
a <- 2       # Beta parameter a
b <- 2       # Beta parameter b
N <- 10000   # Number of samples to generate
burn_in <- 1000 # Burn-in period

set.seed(123) # For reproducibility

# Benchmarking
benchmark_results <- microbenchmark(
  gibbs_Rcpp = gibbsSamplerRcpp(n, a, b, N, burn_in),
  gibbs_R = gibbsSamplerR(n, a, b, N, burn_in),
  times = 10
)

# Print benchmark results
print(benchmark_results)

# Visualize benchmark results
autoplot(benchmark_results)

# Calculate speed-up factor
speedup <- mean(benchmark_results$time[benchmark_results$expr == "gibbs_R"]) /
           mean(benchmark_results$time[benchmark_results$expr == "gibbs_Rcpp"])
speedup_factor <- speedup
cat("Speedup Factor:", speedup_factor, "\n")

 

