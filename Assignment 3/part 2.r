# Set working directory and load data
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 3")
D <- read.csv("datasolar.csv")

# Define model parameters
phi1 <- -0.38
Phi1 <- -0.94
mu <- 5.72

# Compute X_t = log(Y_t) - mu
D$X <- log(D$power) - mu

# Preallocate residual vector
residuals <- rep(NA, nrow(D))

# Loop over valid time steps: t = 13 to 35 → residuals for t+1 = 14 to 36
for (t in 13:(nrow(D) - 1)) {
  X_t1 <- D$X[t + 1]         # X_{t+1}
  X_t  <- D$X[t]             # X_t
  X_tm11 <- D$X[t - 11]      # X_{t-11}
  X_tm12 <- D$X[t - 12]      # X_{t-12}

  # Prediction based on model
  X_pred <- 0.38 * X_t + 0.94 * X_tm11 - 0.3572 * X_tm12

  # Residual
  residuals[t + 1] <- X_t1 - X_pred
}

# Add residuals to data frame
D$residuals <- residuals

# Check number of non-NA residuals
sum(!is.na(D$residuals))  # Should now return 23


# Plot residuals
plot(D$residuals, type = "h", main = "One-step Prediction Residuals", ylab = "Residual", xlab = "Time (month)")


# Plot ACF
acf(na.omit(D$residuals), main = "ACF of Residuals")

# Compare residual variance to theoretical sigma^2 = 0.222
var(na.omit(D$residuals))

# plot predicted values vs actual values
plot(D$power, type = "l", col = "blue", main = "Predicted vs Actual Values", ylab = "Value", xlab = "Time (month)")
lines(D$X + D$residuals, col = "red")
legend("topright", legend = c("Predicted", "Actual"), col = c("blue", "red"), lty = 1)
