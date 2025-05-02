#library(httpgd)
#hgd() 

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
