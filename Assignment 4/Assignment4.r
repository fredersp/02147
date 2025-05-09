library(httpgd)
hgd() 

# Task 1.1.1

# Parameters
a <- 0.9
b <- 1        # bias term (not dynamic input)
sigma1 <- 1
n <- 100
X0 <- 5

# Simulate 5 independent trajectories
set.seed(42)
X_list <- list()

for (j in 1:5) {
  X <- numeric(n)
  X[1] <- X0
  for (t in 2:n) {
    X[t] <- a * X[t - 1] + b + rnorm(1, mean = 0, sd = sigma1)
  }
  X_list[[j]] <- X
}

# Plot
matplot(1:n, do.call(cbind, X_list), type = "l", lty = 1, col = 1:5,
        xlab = "Time", ylab = "State X_t", main = "5 Simulated Realizations of X_t")
legend("topleft", legend = paste("Run", 1:5), col = 1:5, lty = 1)


# Task 1.2

# Use first trajectory from Task 1.1
X <- X_list[[2]]

# Observation noise
sigma2 <- 1

# Simulate noisy observations
Y <- X + rnorm(length(X), mean = 0, sd = sigma2)

# Plot state and observation
plot(X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y)), xlab = "Time", ylab = "Value",
     main = "Latent State $X_t$ and Noisy Observations $Y_t$")
lines(Y, col = "red", lwd = 1.5, lty = 2)
legend("topleft", legend = c("State X_t", "Observation Y_t"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)



# Task 1.3

# The Kalman filter function
myKalmanFilter <- function(
  y,             # Vector of observations y_t
  theta,         # Model parameters for X_{t+1} = a*X_t + b + e_t
  R,             # Measurement noise variance
  x_prior = 0,   # Initial prior mean for X_0
  P_prior = 10   # Initial prior variance for X_0
) {


  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  N <- length(y)

  # Initialize vectors
  x_pred  <- numeric(N)  # Predicted means
  P_pred  <- numeric(N)  # Predicted variances
  x_filt  <- numeric(N)  # Filtered means
  P_filt  <- numeric(N)  # Filtered variances
  innovation     <- numeric(N)  # Innovation = y[t] - x_pred[t]
  innovation_var <- numeric(N)  # Innovation variance

  for (t in seq_len(N)) {
    # Prediction step
    if (t == 1) {
      x_pred[t] <- a * x_prior + b
      P_pred[t] <- a^2 * P_prior + sigma1^2
    } else {
      x_pred[t] <- a * x_filt[t - 1] + b
      P_pred[t] <- a^2 * P_filt[t - 1] + sigma1^2
    }

    # Update step
    innovation[t] <- y[t] - x_pred[t]
    innovation_var[t] <- P_pred[t] + R
    K_t <- P_pred[t] / innovation_var[t]  # Kalman gain: how much to trust observation

    x_filt[t] <- x_pred[t] + K_t * innovation[t]
    P_filt[t] <- (1 - K_t) * P_pred[t]
  }

  return(list(
    x_pred = x_pred,
    P_pred = P_pred,
    x_filt = x_filt,
    P_filt = P_filt,
    innovation = innovation,
    innovation_var = innovation_var
  ))
}

# Define parameters
theta <- c(a, b, sigma1)  # a = 0.9, b = 1, sigma1 = 1 from Task 1.1
R <- sigma2               # Measurement noise variance = 1
X0 <- 5                   # Prior mean (could also use 5, as in simulation)
P0 <- 10                  # Prior variance

# Run the Kalman filter
result <- myKalmanFilter(Y, theta, R, X0, P0)

# Extract predicted state and its variance
x_pred <- result$x_pred
P_pred <- result$P_pred

# Compute 95% confidence interval
ci_upper <- x_pred + 1.96 * sqrt(P_pred)
ci_lower <- x_pred - 1.96 * sqrt(P_pred)

# Time vector
time <- 1:length(Y)

# Plot everything
plot(time, X, type = "l", lwd = 2, col = "black", ylim = range(c(X, Y, ci_upper, ci_lower)),
     ylab = "State / Observation", xlab = "Time", main = "Kalman Filter Output")

# Add observations
lines(time, Y, col = "darkgray", lty = 2)

# Add predicted state
lines(time, x_pred, col = "blue", lwd = 2)

# Add 95% confidence band
polygon(c(time, rev(time)),
        c(ci_upper, rev(ci_lower)),
        col = rgb(0.2, 0.4, 1, alpha = 0.2), border = NA)

# Add legend
legend("topleft",
       legend = c("True state (X)", "Observation (Y)", "Predicted state", "95% CI"),
       col = c("black", "darkgray", "blue", rgb(0.2, 0.4, 1, alpha = 0.5)),
       lty = c(1, 2, 1, NA), lwd = c(2, 1, 2, NA), pch = c(NA, NA, NA, 15),
       pt.cex = 2, bty = "n")


#################### Task 1.4

# myLogLikFun <- function(theta, Y, R, X0, P0) {
#   a <- theta[1]
#   b <- theta[2]
#   sigma1 <- theta[3]
  
#   N <- length(Y)
  
#   # Run Kalman filter
#   kf <- myKalmanFilter(y = Y, theta = theta, R = R, x_prior = X0, P_prior = P0)
  
#   # Compute log-likelihood
#   loglik <- 0
#   for (t in 1:N) {
#     nu_t <- kf$innovation[t]
#     S_t <- kf$innovation_var[t]
#     loglik <- loglik + dnorm(nu_t, mean = 0, sd = sqrt(S_t), log = TRUE)
#   }
  
#   # Return NEGATIVE log-likelihood for minimization
#   return(-loglik)
# }

myLogLikFun <- function(theta, y, R, x_prior = 0, P_prior = 10) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  
  N <- length(y)
  
  # Run Kalman filter
  kf <- myKalmanFilter(y = y, theta = theta, R = R, x_prior = x_prior, P_prior = P_prior)
  
  # Compute log-likelihood
  loglik <- 0
  for (t in 1:N) {
    nu_t <- kf$innovation[t]
    S_t <- kf$innovation_var[t]
    
    if (S_t <= 0 || is.na(S_t)) return(Inf)  # Prevent invalid density calc
    
    loglik <- loglik + dnorm(nu_t, mean = 0, sd = sqrt(S_t), log = TRUE)
  }
  
  return(-loglik)  # NEGATIVE log-likelihood
}


set.seed(42)  # for reproducibility

simulate_and_estimate <- function(a_true, b_true, sigma1_true, n = 100, reps = 100) {
  estimates <- matrix(NA, nrow = reps, ncol = 3)  # Store a, b, sigma1
  colnames(estimates) <- c("a", "b", "sigma1")
  
  for (i in 1:reps) {
    # Simulate latent state
    X <- numeric(n)
    X[1] <- 0
    for (t in 2:n) {
      X[t] <- a_true * X[t - 1] + b_true + rnorm(1, 0, sigma1_true)
    }

    # Simulate noisy observations
    Y <- X + rnorm(n, 0, 1)  # sigma2 = 1

    # Initial guess for optimization
    init_theta <- c(0.5, 0.5, 1)

    # Try to estimate parameters with optim
    result <- tryCatch({
      optim(
        par = init_theta,
        fn = myLogLikFun,
        y = Y,
        R = 1,
        x_prior = 0,
        P_prior = 10,
        method = "L-BFGS-B",
        lower = c(-2, -2, 0.01),
        upper = c(10, 10, 10)
      )
    }, error = function(e) NULL)

    # Save result if optimization succeeded
    if (!is.null(result) && !is.null(result$par)) {
      estimates[i, ] <- result$par
    }
  }

  return(as.data.frame(estimates))
}

res1 <- simulate_and_estimate(1, 0.9, 1, n = 100, reps = 100)
summary(res)

res2 <- simulate_and_estimate(5, 0.9, 1, n = 100, reps = 100)
summary(res2)

res3 <- simulate_and_estimate(1, 0.9, 5, n = 100, reps = 100)
summary(res3)


library(ggplot2)
library(tidyr)
library(dplyr)

# Function to plot estimation results against true values
plot_estimation_results <- function(estimates_df, true_values, title = "Parameter Estimation") {
  
  # Convert to long format
  long <- pivot_longer(estimates_df, cols = c("a", "b", "sigma1"), names_to = "parameter", values_to = "estimate")
  
  # Add true value lookup
  true_df <- data.frame(parameter = names(true_values), true = unname(true_values))
  
  ggplot(long, aes(x = parameter, y = estimate, fill = parameter)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
    geom_point(data = true_df, aes(x = parameter, y = true), color = "red", size = 3, shape = 18) +
    labs(title = title, y = "Estimated Value", x = "Parameter") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
}


plot_estimation_results(res1, c(a = 1, b = 0.9, sigma1 = 1), title = "Estimation Results for a=1, b=0.9, sigma1=1")
plot_estimation_results(res2, c(a = 5, b = 0.9, sigma1 = 1), title = "Estimation Results for a=5, b=0.9, sigma1=1")
plot_estimation_results(res3, c(a = 1, b = 0.9, sigma1 = 5), title = "Estimation Results for a=1, b=0.9, sigma1=5")


# SOF ###

# True parameters
a_true <- 1
b_true <- 0.9
sigma1_true <- 1
sigma2 <- 1  # known
n <- 100
n_sim <- 100

# Storage for parameter estimates
results <- matrix(NA, nrow = n_sim, ncol = 3)
colnames(results) <- c("a_hat", "b_hat", "sigma1_hat")

# Simulation function (you already have this)
simulate_data <- function(n, a, b, sigma1, sigma2) {
  X <- numeric(n)
  Y <- numeric(n)
  X[1] <- 0
  Y[1] <- X[1] + rnorm(1, 0, sigma2)
  for (t in 2:n) {
    X[t] <- a * X[t - 1] + b + rnorm(1, 0, sigma1)
    Y[t] <- X[t] + rnorm(1, 0, sigma2)
  }
  list(X = X, Y = Y)
}


# Run the simulations and estimation
for (i in 1:n_sim) {
  sim <- simulate_data(n, a_true, b_true, sigma1_true, sigma2)
  
  est <- optim(par = c(0.5, 0.5, 0.5), fn = myLogLikFun,
               Y = sim$Y, R = sigma2^2, X0 = 0, P0 = 10,
               method = "L-BFGS-B", lower = c(1/1000, 1/1000, 1/1000))
  print(est$par)  # Print estimated parameters for each simulation
  
  results[i, ] <- est$par
}

# Convert to data frame
results_df <- as.data.frame(results)

# Plot results
library(ggplot2)
library(tidyr)

results_long <- pivot_longer(results_df, cols = everything(), names_to = "parameter", values_to = "estimate")

ggplot(results_long, aes(x = parameter, y = estimate)) +
  geom_boxplot(fill = "lightblue") +
  geom_hline(data = data.frame(parameter = c("a_hat", "b_hat", "sigma1_hat"),
                               true = c(a_true, b_true, sigma1_true)),
             aes(yintercept = true), linetype = "dashed", color = "red") +
  labs(title = "Parameter Estimates from 100 Simulations",
       x = "Parameter", y = "Estimated Value") +
  theme_minimal(base_size = 14)




  