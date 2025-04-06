# 2 Predicting monthly solar power

# 2.1
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 3")
# Load the data
# Load the data
D <- read.csv("datasolar.csv")

# Define model parameters
phi1 <- -0.38
Phi1 <- -0.94
mu <- 5.72

# Compute X_t = log(Y_t) - mu
D$X <- log(D$power) - mu

# Preallocate residual vector
residuals <- rep(NA, nrow(D))

# Compute residuals: t from 13 onward (to get t-12)
for (t in D) {
  X_t1 <- D$X[t+1]             # X_t
  X_tm1 <- D$X[t - 1]        # X_{t-1}
  X_tm11 <- D$X[t - 11]     # X_{t-11}
  X_tm12 <- D$X[t - 12]     # X_{t-12}

  # Predicted X_{t+1|t}
  # up to t0 = 12, we have no prediction so we use X_pred <- -phi1 * X_tm1
# and for t > 12 we have X_pred <- -phi1 * X_tm1 -Phi1 * X_tm11 - phi1*Phi1 * X_tm12
if (t <= 12) {
    X_pred <- -phi1 * X_tm1
    } else {
    X_pred <- -phi1 * X_tm1 - Phi1 * X_tm11 - phi1 * Phi1 * X_tm12
    }
  
  # Residual: actual - predicted
  residuals[t + 1] <- X_t1 - X_pred
}



# Add to data frame
D$residuals <- residuals

# Plot residuals
plot(D$residuals, type = "h", main = "One-step Prediction Residuals", ylab = "Residual", xlab = "Time (month)")

# Plot ACF
acf(na.omit(D$residuals), main = "ACF of Residuals")

# Ljung-Box test for white noise
Box.test(na.omit(D$residuals), lag = 10, type = "Ljung-Box")

# Compare residual variance to theoretical sigma^2 = 0.222
var(na.omit(D$residuals))
