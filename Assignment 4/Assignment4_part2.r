library(httpgd)
hgd() 


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

#setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 4")
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 4")


# If needed, load your dataset
data <- read.csv("transformer_data.csv")

# Rename columns for clarity (optional, but helpful)
data <- data %>%
  rename(
    Yt = Y,
    Ta_t = Ta,
    Phi_s_t = S,
    Phi_I_t = I
  )


# Check structure
str(data)

# Reshape for plotting
data_long <- data %>%
  mutate(time = 1:n()) %>%
  pivot_longer(cols = c(Yt, Ta_t, Phi_s_t, Phi_I_t),
               names_to = "variable", values_to = "value")

# Plot all variables over time
ggplot(data_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(title = "Exploratory Time Series Plots",
       x = "Time (hours)",
       y = "Value") +
  theme_minimal()


# 2.2

kf_logLik_dt <- function(par, data) {
  # Parameter unpacking
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  C   <- matrix(par[5], 1, 1)
  Sigma1lt <- matrix(par[6], 1, 1)  # scalar sqrt(sigma1^2)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  Sigma2   <- matrix(par[7]^2, 1, 1)  # scalar variance
  X0       <- matrix(par[8], 1, 1)

  # Data setup
  Y <- as.matrix(data[, "Yt"])             # m x T
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])  # p x T
  Tn <- nrow(Y)

  # Initialization
  n      <- nrow(A)
  x_est  <- X0
  P_est  <- diag(1e1, n)
  logLik <- 0

  for (t in 1:Tn) {
    # PREDICTION
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    # INNOVATION
    y_pred <- C %*% x_pred
    S_t    <- C %*% P_pred %*% t(C) + Sigma2
    innov  <- Y[t, ] - y_pred

    # LOG-LIKELIHOOD
    logLik <- logLik - 0.5 * (log(2 * pi * S_t) + (innov^2) / S_t)

    # UPDATE
    K_t   <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(as.numeric(logLik))
}

# Initial parameters: A, B1, B2, B3, C, sigma1_sqrt, sigma2_sqrt, X0
start_par <- c(0.9, 0.01, 0.01, 0.01, 1, 0.5, 0.5, 20)

# Lower/upper bounds
lower <- c(-1, rep(-10, 3), -10, 1e-6, 1e-6, 10)
upper <- c(  1, rep( 10, 3),  10, 10, 10, 30)


# Optimizer wrapper
estimate_dt <- function(start_par, data, lower=NULL, upper=NULL) {
  negLL <- function(par){ -kf_logLik_dt(par, data) }
  optim(
    par    = start_par, fn = negLL,
    method = "L-BFGS-B",
    lower  = lower, upper = upper,
    control= list(maxit=1000, trace=1)
  )
}
# Estimate
fit <- estimate_dt(start_par, data, lower = lower, upper = upper)

fit$par  # Estimated parameters
fit$value  # Final negative log-likelihood


kalman_predict <- function(par, data) {
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  C   <- matrix(par[5], 1, 1)
  Sigma1lt <- matrix(par[6], 1, 1)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  Sigma2   <- matrix(par[7]^2, 1, 1)
  X0       <- matrix(par[8], 1, 1)

  Y <- as.matrix(data[, "Yt"])
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(Y)

  x_est <- X0
  P_est <- diag(1e1, 1)

  y_pred <- numeric(Tn)

  for (t in 1:Tn) {
    # Prediction
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1
    y_pred[t] <- as.numeric(C %*% x_pred)

    # Update
    S_t    <- C %*% P_pred %*% t(C) + Sigma2
    innov  <- Y[t, ] - C %*% x_pred
    K_t    <- P_pred %*% t(C) %*% solve(S_t)
    x_est  <- x_pred + K_t %*% innov
    P_est  <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(y_pred)
}

library(ggplot2)

# Run prediction
y_hat <- kalman_predict(fit$par, data)

# Add to data
data$Y_hat <- y_hat

# Plot observed vs predicted
ggplot(data, aes(x = time)) +
  geom_line(aes(y = Yt, color = "Observed")) +
  geom_line(aes(y = Y_hat, color = "Predicted")) +
  labs(
    title = "Transformer Temperature: Observed vs Predicted",
    x = "Time (hours)", y = "Temperature (°C)", color = ""
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_minimal(base_size = 14)


# Add predictions and residuals to data
data$residuals <- data$Yt - data$Y_hat

# Plot residuals
ggplot(data, aes(x = time, y = residuals)) +
  geom_line(color = "darkblue") +
  labs(title = "Residuals Over Time", x = "Time (hours)", y = "Residuals") +
  theme_minimal(base_size = 14)

# ACF plot and PACF plot of residuals
par(mfrow = c(1, 2))
acf(data$residuals, main = "ACF of Residuals")
pacf(data$residuals, main = "PACF of Residuals")

# Q-Q plot of residuals
qqnorm(data$residuals)
qqline(data$residuals, col = "red")

# AIC and BIC calculation
logLik_val <- kf_logLik_dt(fit$par, data)
k <- length(fit$par)  # number of estimated parameters
n <- nrow(data)       # number of observations

AIC_val <- -2 * logLik_val + 2 * k
BIC_val <- -2 * logLik_val + log(n) * k

cat("AIC:", AIC_val, "\n")
cat("BIC:", BIC_val, "\n")


# 2.3
A <- matrix(par[1:4], 2, 2)
B <- matrix(par[5:10], 2, 3)
C <- matrix(par[11:12], 1, 2)
Sigma1lt <- matrix(par[13:15], 2, 2)  # lower-triangular
Sigma1 <- Sigma1lt %*% t(Sigma1lt)
Sigma2 <- matrix(par[16]^2, 1, 1)
X0 <- matrix(par[17:18], 2, 1)

kf_logLik_dt <- function(par, data) {
  # Unpack parameters
  A <- matrix(par[1:4], 2, 2)
  B <- matrix(par[5:10], 2, 3)
  C <- matrix(par[11:12], 1, 2)
  
  # Lower-triangular for Sigma1 (2x2)
  Sigma1lt <- matrix(c(par[13], 0, par[14], par[15]), 2, 2)
  Sigma1 <- Sigma1lt %*% t(Sigma1lt)
  
  Sigma2 <- matrix(par[16]^2, 1, 1)  # Ensure positive variance
  X0 <- matrix(par[17:18], 2, 1)

  # Prepare data
  Y <- as.matrix(data$Yt)
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(data)
  
  # Initialize
  x_est <- X0
  P_est <- diag(10, 2)  # prior state covariance
  logLik <- 0

  for (t in 1:Tn) {
    # Prediction
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    # Innovation
    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t] - y_pred

    # Log-likelihood contribution
    logLik <- logLik - 0.5 * (log(2 * pi * S_t) + (innov^2) / S_t)

    # Update
    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(as.numeric(logLik))
}

estimate_dt <- function(start_par, data, lower = NULL, upper = NULL) {
  negLL <- function(par) { -kf_logLik_dt(par, data) }
  
  optim(
    par = start_par,
    fn = negLL,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(maxit = 1000, trace = 1)
  )
}

fit <- estimate_dt(start_par, data, lower = lower, upper = upper)

fit$par  # Estimated parameters
fit$value  # Final negative log-likelihood

start_par <- c(
  # A: 2x2
  0.9, 0.0,
  0.0, 0.9,
  # B: 2x3
  0.01, 0.01, 0.01,
  0.01, 0.01, 0.01,
  # C: 1x2
  1, 0,
  # Sigma1 lower-triangle: 3 values (s1, s2, s3)
  0.1, 0.05, 0.1,
  # Sigma2 sqrt
  0.1,
  # X0: 2 values
  20, 0
)

lower <- c(rep(-1, 12), rep(1e-6, 7))
upper <- c(rep(1, 12), rep(10, 7))

kalman_predict_2D <- function(par, data) {
  # Unpack parameters
  A <- matrix(par[1:4], 2, 2)
  B <- matrix(par[5:10], 2, 3)
  C <- matrix(par[11:12], 1, 2)
  Sigma1lt <- matrix(c(par[13], 0, par[14], par[15]), 2, 2)
  Sigma1 <- Sigma1lt %*% t(Sigma1lt)
  Sigma2 <- matrix(par[16]^2, 1, 1)
  X0 <- matrix(par[17:18], 2, 1)

  Y <- as.matrix(data$Yt)
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(data)

  x_est <- X0
  P_est <- diag(10, 2)
  Y_hat <- numeric(Tn)

  for (t in 1:Tn) {
    # Predict
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1

    # Predict observation
    y_pred <- C %*% x_pred
    Y_hat[t] <- y_pred

    # Innovation
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t] - y_pred

    # Update
    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(Y_hat)
}

# Run predictions using estimated parameters
Y_hat <- kalman_predict_2D(fit$par, data)

# Add to data
data$Y_hat <- Y_hat

# Plot
library(ggplot2)

ggplot(data, aes(x = t)) +
  geom_line(aes(y = Y, color = "Observed")) +
  geom_line(aes(y = Y_hat, color = "Predicted"), linewidth = 1) +
  labs(title = "Transformer Temperature: Observed vs Predicted (2D Kalman Model)",
       x = "Time (hour)", y = "Temperature (°C)", color = "") +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_minimal(base_size = 14)
