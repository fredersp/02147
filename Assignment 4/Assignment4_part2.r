library(httpgd)
hgd() 


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

#setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 4")
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 4")


# If needed, load your dataset
df <- read.csv("transformer_data.csv")

# Rename columns for clarity (optional, but helpful)
df <- df %>%
  rename(
    Yt = Y,
    Ta_t = Ta,
    Phi_s_t = S,
    Phi_I_t = I
  )


# Check structure
str(df)

# Reshape for plotting
df_long <- df %>%
  mutate(time = 1:n()) %>%
  pivot_longer(cols = c(Yt, Ta_t, Phi_s_t, Phi_I_t),
               names_to = "variable", values_to = "value")

# Plot all variables over time
ggplot(df_long, aes(x = time, y = value, color = variable)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(title = "Exploratory Time Series Plots",
       x = "Time (hours)",
       y = "Value") +
  theme_minimal()


# 2.1

kf_logLik_dt <- function(par, df) {
  # Parameter unpacking
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  C   <- matrix(par[5], 1, 1)
  Sigma1lt <- matrix(par[6], 1, 1)  # scalar sqrt(sigma1^2)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  Sigma2   <- matrix(par[7]^2, 1, 1)  # scalar variance
  X0       <- matrix(par[8], 1, 1)

  # Data setup
  Y <- as.matrix(df[, "Y"])             # m x T
  U <- as.matrix(df[, c("Ta", "S", "I")])  # p x T
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
estimate_dt <- function(start_par, df, lower=NULL, upper=NULL) {
  negLL <- function(par){ -kf_logLik_dt(par, df) }
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


kalman_predict <- function(par, df) {
  A   <- matrix(par[1], 1, 1)
  B   <- matrix(par[2:4], 1, 3)
  C   <- matrix(par[5], 1, 1)
  Sigma1lt <- matrix(par[6], 1, 1)
  Sigma1   <- Sigma1lt %*% t(Sigma1lt)
  Sigma2   <- matrix(par[7]^2, 1, 1)
  X0       <- matrix(par[8], 1, 1)

  Y <- as.matrix(df[, "Y"])
  U <- as.matrix(df[, c("Ta", "S", "I")])
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
  geom_line(aes(y = Y, color = "Observed")) +
  geom_line(aes(y = Y_hat, color = "Predicted")) +
  labs(
    title = "Transformer Temperature: Observed vs Predicted",
    x = "Time (hours)", y = "Temperature (°C)", color = ""
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_minimal(base_size = 14)


# Add predictions and residuals to data
data$residuals <- data$Y - data$Y_hat

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
