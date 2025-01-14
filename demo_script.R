
# Bayes estimation with MCMC ###########################################################################################

# First define the algorithm we'll use to generate the Markov chain ####################################################

set.seed(1234)

# Random-walk Metropolis algorithm
RMW <- function(
    h,                 # possibly unnormalized pdf
    X,                 # domain of h
    steps = 1e3,       # number of steps
    step_size_var = 1, # variance of step size
    use_log = FALSE,
    update_h = NULL
) {
  
  # Find h domain dimension
  domain_dim_length <- length(dim(X))
  if (domain_dim_length > 2) stop("Domain dimension length must be 1 or 2")
  if (domain_dim_length == 0) domain_dim <- 1
  else domain_dim <- ncol(X)
  
  # Find h domain size
  if (domain_dim == 1) domain_size <- length(X)
  else domain_size <- nrow(X)
  
  # Find non-zero initiation point
  x0_found <- FALSE
  sample_range <- 1:domain_size
  while (!x0_found) {
    x0_idx <- sample(sample_range, 1)
    if (domain_dim == 1) x0 <- X[x0_idx]
    else x0 <- X[x0_idx,]
    if (class(h(x0)) != "numeric" || length(h(x0)) != 1) stop("h(X) must be a 1D numeric value")
    if (use_log) good_point <- !is.infinite(h(x0))
    else good_point <- h(x0) > 0
    if (good_point) {
      x0_found <- TRUE
      X_sim <- x0
      if (domain_dim > 1) dim(X_sim) <- c(1, domain_dim)
    } else {
      sample_range <- sample_range[sample_range != x0_idx]
      if (length(sample_range) == 0) stop("No non-zero initiation point found")
    }
  }
  
  # Initialize steps at 0
  step <- 0 
  while (step < steps) {
    
    if (domain_dim == 1) x_current <- X_sim[length(X_sim)]
    else x_current <- X_sim[nrow(X_sim),]
    
    # Generate random step
    x_next <- rnorm(domain_dim, x_current, step_size_var)
    
    h_next <- h(x_next)
    h_current <- h(x_current)
    
    # Calculate acceptance probability
    if (use_log) {
      if (h_next > h_current) r <- 1
      else r <- exp(h_next - h_current)
    } else {
      r <- h(x_next)/h(x_current)
    }
    r <- min(1, r)
    if (is.nan(r) || is.infinite(r)) r <- 0
    
    # Accept or reject step
    if (runif(1) < r) {
      if (domain_dim == 1) X_sim <- c(X_sim, x_next)
      else X_sim <- rbind(X_sim, x_next)
      if (!is.null(update_h)) update_h <- h(x_next)
    }
    
    # Increment step
    step <- step + 1
    
  }
  
  # Return simulated X, which will have pdf h
  return(X_sim)
  
}

# Sanity check: Round-about distribution sampling ######################################################################

# Use RMW to sample from a normal distribution with a given mean and sd
given_mean <- 2
given_sd <- 4
test_sample <- RMW(
  h = function(x) dnorm(x, given_mean, given_sd),
  X = c(0),
  steps = 1e4,
  step_size_var = 0.5
)

# Examine results to ensure this works as expected
hist(test_sample)
abline(v = given_mean, col = "red")
abline(v = given_mean + given_sd, col = "blue")
abline(v = given_mean - given_sd, col = "blue")
cat("\nmean of test sample: ", mean(test_sample), ", should be near: ", given_mean, sep = "")
cat("\nstandard deviation of test sample: ", sd(test_sample), ", should be near: ", given_sd, sep = "")

# Bayesian parameter estimation ########################################################################################
# Now use RMW in conjunction with Bayes' rule to estimate the mean and sd from an observed normal sample

# Define R6 object to conveniently store the prior and calculate updates
posterior_unnormalized_R6 <- R6::R6Class(
  "posterior_unnormalized_R6",
  public = list(
    prior = log(0.05),
    update = function(
    theta, # two-valued vector holding mean and sd
    observed_sample
    ) {
      if (theta[2] <= 0) return(-Inf)
      likelihood <- sum(dnorm(observed_sample, mean = theta[1], sd = theta[2], log = TRUE))
      return(likelihood + self$prior)
    }
  )
)

# Simulate some observations 
given_mean <- 4
given_sd <- 7
observed_sample <- rnorm(1e3, given_mean, given_sd)

# Make a new object to use in estimating/ recovering the sample parameters
posterior_unnormalized <- posterior_unnormalized_R6$new()

# Set the total number of steps to try
n_steps <- 1e4

# Run the MCMC simulation to estimate the parameters 
RMW_Bayes_parameter_sims <- RMW(
  h = function(x) posterior_unnormalized$update(x, observed_sample),
  X = array(runif(20), c(10,2)), # A few random initization values to try
  steps = n_steps,
  step_size_var = 0.05,
  use_log = TRUE,
  update_h = posterior_unnormalized$prior
)

# Clean up results a bit
colnames(RMW_Bayes_parameter_sims) <- c("mean", "sd")

# Count steps taken
n_steps_taken <- nrow(RMW_Bayes_parameter_sims)

# Discard first 25% of steps as burn-in
RMW_Bayes_parameter_sims_end <- RMW_Bayes_parameter_sims[as.integer(n_steps_taken*0.25):n_steps_taken,]

# Examine results 
estimates <- colMeans(RMW_Bayes_parameter_sims_end)
col_quantiles <- apply(RMW_Bayes_parameter_sims_end, 2, quantile, probs = c(0.025, 0.975))
cat("\nExpected mean:", given_mean)
cat("\nExtimated mean: ", estimates[1], ", 95% CI: ", col_quantiles[1,1], " - ", col_quantiles[2,1], sep = "")
cat("\nExpected sd:", given_sd)
cat("\nExtimated sd: ", estimates[2], ", 95% CI: ", col_quantiles[1,2], " - ", col_quantiles[2,2], sep = "")
cat("\nAcceptance rate: ", n_steps_taken/n_steps)
plot(RMW_Bayes_parameter_sims[,1], col = "blue", main = "Trace of estimated mean", type = "l")
abline(h = given_mean, col = "red")
plot(RMW_Bayes_parameter_sims[,2], col = "blue", main = "Trace of estimated sd", type = "l")
abline(h = given_sd, col = "red")





