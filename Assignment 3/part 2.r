# Set working directory and load data
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 3")
D <- read.csv("datasolar.csv")

# Define model parameters
phi1 <- -0.38
Phi1 <- -0.94
mu <- 5.72
sigma <- 0.22

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
  X_pred <- -phi1 * X_t - Phi1 * X_tm11 - phi1*Phi1 * X_tm12

  # Residual
  residuals[t + 1] <- X_t1 - X_pred
}

# Add residuals to data frame
D$residuals <- residuals

sum(!is.na(D$residuals)) 

# Plot residuals
plot(D$residuals, type = "h", main = "One-step Prediction Residuals", ylab = "Residual", xlab = "Time (month)")

# Plot ACF
acf(na.omit(D$residuals), main = "ACF of Residuals")

library(lubridate)


library(ggplot2)
library(patchwork)
library(lubridate)

# Make sure `D$date` exists
D$date <- ymd(paste(D$year, D$month, "01"))

# Base font size
base_size <- 18

# 1. Time series of residuals
p1 <- ggplot(na.omit(D), aes(x = date, y = residuals)) +
  geom_line(color = "darkblue", linewidth = 1) +
  labs(title = "Residuals Over Time", x = "Time", y = "Residual") +
  theme_minimal(base_size = base_size) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. ACF plot
acf_data <- acf(D$residuals, lag.max = 36, plot = FALSE)
acf_df <- data.frame(
  Lag = acf_data$lag,
  ACF = acf_data$acf
)

p2 <- ggplot(acf_df, aes(x = Lag, y = ACF)) +
  geom_segment(aes(xend = Lag, yend = 0), color = "darkred") +
  geom_hline(yintercept = c(-1.96, 1.96)/sqrt(length(na.omit(D$residuals))),
             linetype = "dashed", color = "gray40") +
  labs(title = "ACF of Residuals",
       x = "Lag ", y = "ACF") +
  theme_minimal(base_size = base_size)


# 4. Q-Q plot
qq_data <- na.omit(D$residuals)
qq_df <- data.frame(
  sample = qq_data,
  theoretical = qqnorm(qq_data, plot.it = FALSE)$x
)

p4 <- ggplot(qq_df, aes(sample = sample)) +
  stat_qq(color = "black", size = 1) +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal(base_size = base_size)

# Combine all plots vertically
p1 + p2  + p4


# 2.2 Forecast next 12 months from t = 36 (i.e. months 37 to 48)

# Extend X vector with room for 12 new predictions
X_extended <- c(D$X, rep(NA, 12))

# Forecast step-by-step using previous values (observed or predicted)
for (k in 1:12) {
  t_k <- 36 + k
  
  X_tkm1  <- X_extended[t_k - 1]   # X_{t+k-1}
  X_tkm12 <- X_extended[t_k - 12]  # X_{t+k-12}
  X_tkm13 <- X_extended[t_k - 13]  # X_{t+k-13}
  
  # Forecast using SAR(1)(1)_12 model
  X_extended[t_k] <- 0.38 * X_tkm1 + 0.94 * X_tkm12 - 0.3572 * X_tkm13
}

# Convert forecasted log-values to actual power
Y_forecast <- exp(X_extended[37:48] + mu)

# Make table
forecast_table <- data.frame(
  Month = 36 + 1:12,
  Forecast_MWh = round(Y_forecast, 2)
)

print(forecast_table)

library(ggplot2)

# Combine original and forecasted data into one data frame
D$index <- 1:nrow(D)
D_forecast <- data.frame(
  index = 37:48,
  power = Y_forecast,
  source = "Forecast"
)
D$source <- "Observed"

# Combine both
D_all <- rbind(
  D[, c("index", "power", "source")],
  D_forecast
)

# Plot using ggplot2
ggplot(D_all, aes(x = index, y = power, color = source)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Observed and Forecasted Monthly Solar Power",
    x = "Month Index",
    y = "Power (MWh)",
    color = "Data Type"
  ) +
  scale_color_manual(values = c("Observed" = "blue", "Forecast" = "green")) +
  theme_minimal(base_size = 18)  # Increase font size






#2.3
library(ggplot2)
library(lubridate)

k <- 12 # Number of forecasted months
# Compute sigma_k^2 for each k = 1 to 12
sigma_k2 <- sapply(1:k, function(h) {
  sum(phi1^(2 * (0:(h - 1)))) * sigma^2
})

# Compute standard deviation (for prediction interval)
sigma_k <- sqrt(sigma_k2)

# Create forecast interval
log_forecast_upper <- X_extended[37:48] + qnorm(0.975) * sigma_k
log_forecast_lower <- X_extended[37:48] - qnorm(0.975) * sigma_k

# Transform back to power
power_upper <- exp(log_forecast_upper + mu)
power_lower <- exp(log_forecast_lower + mu)

# Create month labels from original data
start_date <- ymd(paste(D$year[1], D$month[1], "01", sep = "-"))
date_seq <- seq.Date(from = start_date, by = "month", length.out = 48)

# Forecast data frame
D_forecast <- data.frame(
  date = date_seq[37:48],
  power = Y_forecast,
  lower = power_lower,
  upper = power_upper,
  source = "Forecast"
)

# Observed data frame
D_obs <- data.frame(
  date = date_seq[1:36],
  power = D$power,
  lower = NA,
  upper = NA,
  source = "Observed"
)

# Combine
D_all <- rbind(D_obs, D_forecast)

# Plot
ggplot(D_all, aes(x = date, y = power, color = source)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_ribbon(
    data = D_forecast,
    aes(x = date, ymin = lower, ymax = upper),
    fill = "green",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  labs(
    title = "Observed and Forecasted Monthly Solar Generation with Prediction Intervals",
    x = "Month",
    y = "Generation (MWh)",
    color = "Data Type"
  ) +
  scale_color_manual(values = c("Observed" = "blue", "Forecast" = "darkgreen")) +
  theme_minimal(base_size = 18) +  # Increase font size
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
