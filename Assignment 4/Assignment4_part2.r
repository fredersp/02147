library(httpgd)
hgd() 


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

#setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 4")
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 4")
#setwd("C:/Users/frede/OneDrive/Skrivebord/02417 - Time Series Analysis/Assignments/02147/Assignment 4")




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


################### 2.2 ####################

kf_logLik_dt_1D <- function(par, data) {
  A <- matrix(par[1], 1, 1)
  B <- matrix(par[2:4], 1, 3)
  C <- matrix(par[5], 1, 1)
  Sigma1lt <- matrix(par[6], 1, 1)
  Sigma1 <- Sigma1lt %*% t(Sigma1lt)
  Sigma2 <- matrix(par[7]^2, 1, 1)
  X0 <- matrix(par[8], 1, 1)

  Y <- as.matrix(data$Yt)
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(Y)

  x_est <- X0
  P_est <- diag(1e1, 1)
  logLik <- 0

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1
    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t, ] - y_pred

    logLik <- logLik - 0.5 * (log(2 * pi * S_t) + (innov^2) / S_t)
    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(as.numeric(logLik))
}

# 1D Estimation wrapper
estimate_dt_1D <- function(start_par, data, lower = NULL, upper = NULL) {
  negLL <- function(par) { -kf_logLik_dt_1D(par, data) }
  optim(par = start_par, fn = negLL, method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = 1000, trace = 1))
}

# 1D Kalman Prediction
kalman_predict_1D <- function(par, data) {
  A <- matrix(par[1], 1, 1)
  B <- matrix(par[2:4], 1, 3)
  C <- matrix(par[5], 1, 1)
  Sigma1lt <- matrix(par[6], 1, 1)
  Sigma1 <- Sigma1lt %*% t(Sigma1lt)
  Sigma2 <- matrix(par[7]^2, 1, 1)
  X0 <- matrix(par[8], 1, 1)

  Y <- as.matrix(data$Yt)
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(Y)

  x_est <- X0
  P_est <- diag(1e1, 1)
  y_pred <- numeric(Tn)

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1
    y_pred[t] <- as.numeric(C %*% x_pred)
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t, ] - C %*% x_pred
    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(y_pred)
}

# Initial values
start_par_1D <- c(0.9, 0.01, 0.01, 0.01, 1, 0.5, 0.5, 20)
lower_1D <- c(-1, rep(-10, 3), -10, 1e-6, 1e-6, 10)
upper_1D <- c(1, rep(10, 3), 10, 10, 10, 30)

# Fit model
fit_1D <- estimate_dt_1D(start_par_1D, data, lower_1D, upper_1D)
y_hat_1D <- kalman_predict_1D(fit_1D$par, data)
data$Y_hat_1D <- y_hat_1D
data$residuals_1D <- data$Yt - y_hat_1D

# 1D parameters
par_1D <- fit_1D$par
cat("Estimated parameters (1D):\n")
cat("A:", par_1D[1], "\n")
cat("B1 (Outside Air Temperature):", par_1D[2], "\n")
cat("B2 (Horizontal Global Solar Radiation):", par_1D[3], "\n")
cat("B3 (Transformer Load):", par_1D[4], "\n")
cat("C:", par_1D[5], "\n")
cat("Sigma1 sqrt:", par_1D[6], "\n")
cat("Sigma2 sqrt:", par_1D[7], "\n")
cat("X0:", par_1D[8], "\n")
library(ggplot2)
library(patchwork)

# Observed vs predicted
p1 <- ggplot(data, aes(x = time)) +
  geom_line(aes(y = Yt, color = "Observed")) +
  geom_line(aes(y = Y_hat_1D, color = "Predicted")) +
  labs(
    title = "Transformer Temperature: Observed vs Predicted",
    x = "Time (hours)", y = "Temperature (°C)", color = ""
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_minimal(base_size = 14)

# Residuals
p2 <- ggplot(data, aes(x = time, y = residuals_1D)) +
  geom_line(color = "darkblue") +
  labs(title = "Residuals Over Time", x = "Time (hours)", y = "Residuals") +
  theme_minimal(base_size = 14)

# ACF and PACF
acf_data <- acf(data$residuals_1D, plot = FALSE, lag.max = 25)
pacf_data <- pacf(data$residuals_1D, plot = FALSE, lag.max = 25)

acf_df <- data.frame(lag = acf_data$lag, acf = acf_data$acf)
pacf_df <- data.frame(lag = pacf_data$lag, pacf = pacf_data$acf)

p3 <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_segment(aes(xend = lag, yend = 0), color = "red") +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(data$residuals_1D)), linetype = "dashed") +
  labs(title = "ACF of Residuals", x = "Lag", y = "ACF") +
  theme_minimal(base_size = 14)

p4 <- ggplot(pacf_df, aes(x = lag, y = pacf)) +
  geom_segment(aes(xend = lag, yend = 0), color = "blue") +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(data$residuals_1D)), linetype = "dashed") +
  labs(title = "PACF of Residuals", x = "Lag", y = "PACF") +
  theme_minimal(base_size = 14)

# QQ plot
qq <- ggplot(data, aes(sample = residuals_1D)) +
  stat_qq() +
  stat_qq_line(col = "red") +
  labs(title = "Q-Q Plot of Residuals") +
  theme_minimal(base_size = 14)

# Combine all plots
combined <- (p1 / p2 / (p3 | p4 | qq)) & theme_minimal(base_size = 18)

# Show combined figure
combined


# AIC and BIC calculation
logLik_val <- kf_logLik_dt_1D(fit_1D$par, data)
k <- length(fit_1D$par)  # number of estimated parameters
n <- nrow(data)       # number of observations

AIC_val_1D <- -2 * logLik_val + 2 * k
BIC_val_1D <- -2 * logLik_val + log(n) * k

cat("AIC:", AIC_val_1D, "\n")
cat("BIC:", BIC_val_1D, "\n")




########## 2.3 ##########

kf_logLik_dt_2D <- function(par, data) {
  A <- matrix(par[1:4], 2, 2)
  B <- matrix(par[5:10], 2, 3)
  C <- matrix(par[11:12], 1, 2)
  Sigma1lt <- matrix(c(par[13], 0, par[14], par[15]), 2, 2)
  Sigma1 <- Sigma1lt %*% t(Sigma1lt)
  Sigma2 <- matrix(par[16]^2, 1, 1)
  X0 <- matrix(par[17:18], 2, 1)

  Y <- as.matrix(data$Yt)
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(Y)

  x_est <- X0
  P_est <- diag(10, 2)
  logLik <- 0

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1
    y_pred <- C %*% x_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t] - y_pred

    logLik <- logLik - 0.5 * (log(2 * pi * S_t) + (innov^2) / S_t)
    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(as.numeric(logLik))
}

# Estimation
estimate_dt_2D <- function(start_par, data, lower = NULL, upper = NULL) {
  negLL <- function(par) { -kf_logLik_dt_2D(par, data) }
  optim(par = start_par, fn = negLL, method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = 500, trace = 1))
}

# Prediction
kalman_predict_2D <- function(par, data) {
  A <- matrix(par[1:4], 2, 2)
  B <- matrix(par[5:10], 2, 3)
  C <- matrix(par[11:12], 1, 2)
  Sigma1lt <- matrix(c(par[13], 0, par[14], par[15]), 2, 2)
  Sigma1 <- Sigma1lt %*% t(Sigma1lt)
  Sigma2 <- matrix(par[16]^2, 1, 1)
  #X0 <- matrix(c(data$Yt[1], 0), 2, 1)
  X0 <- matrix(par[17:18], 2, 1)

  Y <- as.matrix(data$Yt)
  U <- as.matrix(data[, c("Ta_t", "Phi_s_t", "Phi_I_t")])
  Tn <- nrow(Y)

  x_est <- X0
  P_est <- diag(10, 2)
  Y_hat <- numeric(Tn)

  for (t in 1:Tn) {
    x_pred <- A %*% x_est + B %*% matrix(U[t, ], ncol = 1)
    P_pred <- A %*% P_est %*% t(A) + Sigma1
    y_pred <- C %*% x_pred
    Y_hat[t] <- y_pred
    S_t <- C %*% P_pred %*% t(C) + Sigma2
    innov <- Y[t] - y_pred
    K_t <- P_pred %*% t(C) %*% solve(S_t)
    x_est <- x_pred + K_t %*% innov
    P_est <- P_pred - K_t %*% S_t %*% t(K_t)
  }

  return(Y_hat)
}

# Initial values
init_X0 <- c(data$Yt[1], 0)  # First state roughly aligned with first Yt

start_par_2D <- c(
  0.9, 0.0,  # A matrix
  0.0, 0.9,
  0.01, 0.01, 0.01,  # B row 1
  0.01, 0.01, 0.01,  # B row 2
  1, 0,              # C matrix
  0.1, 0.05, 0.1,    # Sigma1 lower triangle
  0.1,               # Sigma2 (observation)
  init_X0            # <- updated initial state
)

# Initial parameters
# lEG MED VÆRDIERNE HER, LIGE NU ER INTIAL GUESS SAT TIL 10
# start_par_2D <- c(
#   0.9, 0.0,
#   0.0, 0.9,
#   0.01, 0.01, 0.01,
#   0.01, 0.01, 0.01,
#   1, 0,
#   0.1, 0.05, 0.1,
#   0.1,
#   10, 0
# )

lower_2D <- c(rep(-1, 12), rep(1e-6, 4), 10, -10)  # tighter lower bound for X0
upper_2D <- c(rep(1, 12), rep(10, 4), 40, 10)      # tighter upper bound for X0

#lower_2D <- c(rep(-1, 12), rep(1e-6, 7))
#upper_2D <- c(rep(1, 12), rep(10, 7))



#start_par_2D <- start_par_2D[1:16]
#lower_2D <- lower_2D[1:16]
#upper_2D <- upper_2D[1:16]

# Fit 2D model

max_retries <- 5

# KØR DENNE LINJE FØRST OG EFTERFØLGENDE LOOPET NEDENUNDER 
fit_2D <- estimate_dt_2D(start_par_2D, data, lower = lower_2D, upper = upper_2D) 
# DETTE ER FOR AT SIKRE AT OPTIMERINGEN IKKE CHRASHER

for (i in 1:max_retries) {
  if (fit_2D$convergence == 0) break
  fit_2D <- estimate_dt_2D(fit_2D$par, data, lower = lower_2D, upper = upper_2D)
}


y_hat_2D <- kalman_predict_2D(fit_2D$par, data)
data$Y_hat_2D <- y_hat_2D
data$residuals_2D <- data$Yt - y_hat_2D

# Plot observed vs predicted
ggplot(data, aes(x = time)) +
  geom_line(aes(y = Yt, color = "Observed")) +
  geom_line(aes(y = Y_hat_2D, color = "Predicted")) +
  labs(
    title = "Transformer Temperature: Observed vs Predicted",
    x = "Time (hours)", y = "Temperature (°C)", color = ""
  ) +
  scale_color_manual(values = c("Observed" = "black", "Predicted" = "red")) +
  theme_minimal(base_size = 14)


# Plot residuals
ggplot(data, aes(x = time, y = residuals_2D)) +
  geom_line(color = "darkblue") +
  labs(title = "Residuals Over Time", x = "Time (hours)", y = "Residuals") +
  theme_minimal(base_size = 14)

# ACF plot and PACF plot of residuals
par(mfrow = c(1, 2))
acf(data$residuals_2D, main = "ACF of Residuals")
pacf(data$residuals_2D, main = "PACF of Residuals")
# Q-Q plot of residuals
qqnorm(data$residuals_2D) 
qqline(data$residuals_2D, col = "red")

# AIC and BIC calculation
logLik_val_2D <- kf_logLik_dt_2D(fit_2D$par, data)
k_2D <- length(fit_2D$par)  # number of estimated parameters
n_2D <- nrow(data)       # number of observations

AIC_val_2D <- -2 * logLik_val_2D + 2 * k_2D
BIC_val_2D <- -2 * logLik_val_2D + log(n_2D) * k_2D
cat("AIC (2D):", AIC_val_2D, "\n")
cat("BIC (2D):", BIC_val_2D, "\n")



######### 2.4 #########
# Having estimated a 2-dimensional model, you will now explore 
# the meaning of the two latent states.

# Plot the two reconstructed state trajectories over time in a single figure.

# Extract the state estimates
A <- matrix(fit_2D$par[1:4], 2, 2)
B <- matrix(fit_2D$par[5:10], 2, 3)
C <- matrix(fit_2D$par[11:12], 1, 2)

Sigma1lt <- matrix(c(fit_2D$par[13], 0, fit_2D$par[14], fit_2D$par[15]), 2, 2)
Sigma1 <- Sigma1lt %*% t(Sigma1lt)
Sigma2 <- matrix(fit_2D$par[16]^2, 1, 1)
X0 <- matrix(fit_2D$par[17:18], 2, 1)

# Initialize state estimates
x_est <- X0
P_est <- diag(10, 2)
state_estimates <- matrix(NA, nrow = nrow(data), ncol = 2)
# Kalman filter to estimate states
for (t in 1:nrow(data)) {
  U_t <- matrix(as.numeric(data[t, c("Ta_t", "Phi_s_t", "Phi_I_t")]), ncol = 1)

  # Prediction
  x_pred <- A %*% x_est + B %*% U_t
  P_pred <- A %*% P_est %*% t(A) + Sigma1

  # Update
  y_pred <- C %*% x_pred
  S_t <- C %*% P_pred %*% t(C) + Sigma2
  innov <- data$Yt[t] - y_pred

  K_t <- P_pred %*% t(C) %*% solve(S_t)
  x_est <- x_pred + K_t %*% innov
  P_est <- P_pred - K_t %*% S_t %*% t(K_t)

  # Store state estimates
  state_estimates[t, ] <- as.numeric(x_est)
}

# Convert to data frame for plotting
state_estimates_df <- as.data.frame(state_estimates)
colnames(state_estimates_df) <- c("State 1", "State 2")
state_estimates_df$time <- 1:nrow(data)
# Plot the state estimates
ggplot(state_estimates_df, aes(x = time)) +
  geom_line(aes(y = `State 1`, color = "State 1")) +
  geom_line(aes(y = `State 2`, color = "State 2")) +
  labs(title = "Latent State Estimates Over Time",
       x = "Time (hours)", y = "State Value", color = "") +
  scale_color_manual(values = c("State 1" = "blue", "State 2" = "red")) +
  theme_minimal(base_size = 14)
# Den anden state space kunne ligne noget køling
