# Part 3 - Predicting hourly heating (W) dependent on temp and solar radiation

library(ggplot2)
library(httpgd)
hgd()


# Set working directory and load data
setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 3")
#setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 3")

D <- read.csv("box_data_60min.csv")

D$tdate <- as.POSIXct(D$tdate, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

###########################################################
# 3.1
###########################################################

# Plot non lagged data (Ph, Tdelta & Gv)
# Set up plotting layout: 3 rows, 1 column
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))  # Adjust margins if needed

# Plot 1: Ph
plot(D$Ph, type = "l", col = "blue", main = "Heating Power (Ph)", ylab = "W", xlab = "Time (hour)")
grid()

# Plot 2: Tdelta
plot(D$Tdelta, type = "l", col = "red", main = "Temperature Difference (Tdelta)", ylab = "°C", xlab = "Time (hour)")
grid()

# Plot 3: Gv
plot(D$Gv, type = "l", col = "green", main = "Solar Radiation (Gv)", ylab = "W/m²", xlab = "Time (hour)")
grid()

# Notes on the plots:
    # Temp difference and heatning are positive correlated
    # Solar radiation and heating are negative correlated 

###########################################################
# 3.2
###########################################################

# Spit data into traina and test 
teststart <- as.POSIXct("2013-02-06 00:00", tz = "UTC")
Dtrain <- D[D$tdate < teststart, ]
Dtest <- D[D$tdate >= teststart, ]

###########################################################
# 3.3
###########################################################

# Scatter plots (Ph vs Tdelta, Ph vs Gv)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(Dtrain$Tdelta, Dtrain$Ph, col = "red", pch = 16,
     main = "Ph vs Tdelta", xlab = "Tdelta (°C)", ylab = "Ph (W)")

plot(Dtrain$Gv, Dtrain$Ph, col = "green", pch = 16,
     main = "Ph vs Gv", xlab = "Gv (W/m²)", ylab = "Ph (W)")

# Cross
ccf(Dtrain$Tdelta, Dtrain$Ph, lag.max = 30, main = "CCF: Tdelta vs Ph")
ccf(Dtrain$Gv, Dtrain$Ph, lag.max = 30, main = "CCF: Gv vs Ph")

# Notes: the last temp difference for several hours are quite important for the heating power.
# The three lags before radiation is important for the heating power and seasonality 24 hours.



# Autocorrelation
acf(Dtrain$Ph, lag.max = 30, main = "ACF of Ph")
# Notes: seasonality in the data, with a periodicity of 24 hours. Lag 1-5 are important

###########################################################
# 3.4
###########################################################
install.packages("onlineforecast")
library(onlineforecast)
install.packages("patchwork")
library(patchwork)

# Define input and output for one system (e.g., Tdelta → Ph)
x_Tdelta <- Dtrain$Tdelta
x_Gv <- Dtrain$Gv
y <- Dtrain$Ph

# Set maximum lag
p <- 30

# Create lag matrix
X_lags_Tdelta <- lagdf(x_Tdelta, 0:p)
X_lags_Gv <- lagdf(x_Gv, 0:p)

# Fit least squares regression (impulse response)
D_model_Tdelta <- as.data.frame(cbind(y = y, X_lags_Tdelta))
D_model_Gv <- as.data.frame(cbind(y = y, X_lags_Gv))

form <- paste0("y ~ 0 + ", paste0("k", 0:p, collapse = " + "))

fit_Tdelta <- lm(form, data = D_model_Tdelta)
fit_Gv <- lm(form, data = D_model_Gv)


# Extract impulse response coefficients
ir_Tdelta <- data.frame(
  Lag = 0:p,
  IR = coef(fit_Tdelta),
  Variable = "Tdelta"
)

ir_Gv <- data.frame(
  Lag = 0:p,
  IR = coef(fit_Gv),
  Variable = "Gv"
)

# Base text size
base_size <- 30

# Tdelta → Ph
plot_Tdelta <- ggplot(ir_Tdelta, aes(x = Lag, y = IR)) +
  geom_col(fill = "steelblue", width = 0.4) +
  geom_hline(yintercept = 0, color = "gray") +
  labs(title = "Impulse Response: Tdelta → Ph",
       x = "Lag", y = "Estimated Effect") +
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26)
  )

# Gv → Ph
plot_Gv <- ggplot(ir_Gv, aes(x = Lag, y = IR)) +
  geom_col(fill = "seagreen", width = 0.4) +
  geom_hline(yintercept = 0, color = "gray") +
  labs(title = "Impulse Response: Gv → Ph",
       x = "Lag", y = "Estimated Effect") +
  theme_minimal(base_size = base_size) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 26),
    axis.text = element_text(size = 26)
  )


plot_Tdelta + plot_Gv

###########################################################
# 3.5
###########################################################

# Fitting the linear regresseion model
model <- lm(Ph ~ Tdelta + Gv, data = Dtrain)
summary(model)

# One step predictions
Dtrain$Ph_hat <- predict(model)


# residuals
Dtrain$residuals <- Dtrain$Ph - Dtrain$Ph_hat


# Actual vs predicted
ggplot(Dtrain, aes(x = tdate)) +
  geom_line(aes(y = Ph, color = "Actual"), size = 1) +
  geom_line(aes(y = Ph_hat, color = "Predicted"), size = 1) +
  labs(title = "Actual vs Predicted Heating Power (OLS)",
       y = "Ph", x = "Time", color = "Legend") +
  scale_color_manual(values = c("Actual" = "blue", "Predicted" = "red")) +
  theme_minimal(base_size = 14)

# Residuals
ggplot(Dtrain, aes(x = tdate, y = residuals)) +
  geom_line(color = "gray") +
  labs(title = "Residuals from OLS model", y = "Residual", x = "Time") +
  theme_minimal()

# ACF of residuals
acf(Dtrain$residuals, lag.max = 30, main = "ACF of Residuals")
# Notes: The residuals are not white noise, as there is a significant correlation at lag 1, 2, 3, 4 and 23, 24
# could imply the model need some time dependency like an AR or SAR model. 

# ccf plot
ccf(Dtrain$residuals, Dtrain$Tdelta, lag.max = 30, main = "CCF: Residuals vs Tdelta")
ccf(Dtrain$residuals, Dtrain$Gv, lag.max = 30, main = "CCF: Residuals vs Gv")
# Gv is important for the residuals, but Tdelta is not.

###########################################################
# 3.6
###########################################################

# Fit model and predict
fit_arx <- lm(Ph ~ Ph.l1 + Tdelta + Gv, data = Dtrain)
Dtrain$Ph_hat_arx <- predict(fit_arx)
Dtrain$residuals_arx <- Dtrain$Ph - Dtrain$Ph_hat_arx

# actual vs predicted
ggplot(Dtrain, aes(x = tdate)) +
  geom_line(aes(y = Ph, color = "Actual"), size = 1) +
  geom_line(aes(y = Ph_hat_arx, color = "Predicted"), size = 1) +
  labs(title = "Actual vs Predicted Heating Power (ARX)",
       y = "Ph", x = "Time", color = "Legend") +
  scale_color_manual(values = c("Actual" = "blue", "Predicted" = "red")) +
  theme_minimal(base_size = 14)

# Residuals
ggplot(Dtrain, aes(x = tdate, y = residuals_arx)) +
  geom_line(color = "gray") +
  labs(title = "Residuals from ARX model", y = "Residual", x = "Time") +
  theme_minimal()

# ACF of residuals
acf(Dtrain$residuals_arx, lag.max = 30, main = "ACF of Residuals (ARX)")

# ccf plot
ccf(Dtrain$residuals_arx, Dtrain$Tdelta, lag.max = 30, main = "CCF: Residuals vs Tdelta (ARX)")
ccf(Dtrain$residuals_arx, Dtrain$Gv, lag.max = 30, main = "CCF: Residuals vs Gv (ARX)")

###########################################################
# 3.7
###########################################################
fit_arx2 <- lm(Ph ~ Ph.l1 + Ph.l2 + Tdelta + Tdelta.l1 + Gv + Gv.l1, data = Dtrain)
Dtrain$Ph_hat_arx2 <- predict(fit_arx2)
Dtrain$residuals_arx2 <- Dtrain$Ph - Dtrain$Ph_hat_arx2

# actual vs predicted
ggplot(Dtrain, aes(x = tdate)) +
  geom_line(aes(y = Ph, color = "Actual"), size = 1) +
  geom_line(aes(y = Ph_hat_arx2, color = "Predicted"), size = 1) +
  labs(title = "Actual vs Predicted Heating Power (ARX2)",
       y = "Ph", x = "Time", color = "Legend") +
  scale_color_manual(values = c("Actual" = "blue", "Predicted" = "red")) +
  theme_minimal(base_size = 14)

# Residuals
ggplot(Dtrain, aes(x = tdate, y = residuals_arx2)) +
  geom_line(color = "gray") +
  labs(title = "Residuals from ARX2 model", y = "Residual", x = "Time") +
  theme_minimal()

# ACF of residuals
acf(Dtrain$residuals_arx2, lag.max = 30, main = "ACF of Residuals (ARX2)")

#ccf plot
ccf(Dtrain$residuals_arx2, Dtrain$Tdelta, lag.max = 30, main = "CCF: Residuals vs Tdelta (ARX2)")
ccf(Dtrain$residuals_arx2, Dtrain$Gv, lag.max = 30, main = "CCF: Residuals vs Gv (ARX2)")

# Plot BIC and AIC vs. the increasing model order. 
# Set maximum model order
max_order <- 9

# Create data frame to store results
model_metrics <- data.frame(Order = integer(), AIC = numeric(), BIC = numeric())

# Loop through model orders
for (order in 1:max_order) {
  
  # Define column names
  ph_lags <- paste0("Ph.l", 1:order)
  tdelta_lags <- paste0("Tdelta.l", 0:(order - 1))
  gv_lags <- paste0("Gv.l", 0:(order - 1))
  
  # Combine into formula
  predictors <- c(ph_lags, tdelta_lags, gv_lags)
  formula_str <- paste("Ph ~", paste(predictors, collapse = " + "))
  formula <- as.formula(formula_str)
  
  
  # Fit model and store AIC/BIC
  fit <- lm(formula, data = Dtrain)
  model_metrics <- rbind(model_metrics, data.frame(
    Order = order,
    AIC = AIC(fit),
    BIC = BIC(fit)
  ))
}

library(tidyr)
library(ggplot2)

model_metrics_long <- pivot_longer(model_metrics, cols = c("AIC", "BIC"), names_to = "Metric", values_to = "Value")

ggplot(model_metrics_long, aes(x = Order, y = Value, color = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(title = "AIC and BIC vs ARX Model Order",
       x = "Model Order", y = "Information Criterion", color = "Metric") +
  scale_color_manual(values = c("AIC" = "blue", "BIC" = "red")) +
  theme_minimal(base_size = 14)


###########################################################
# 3.8
###########################################################

install.packages("Metrics")
library(Metrics)  # for rmse()

max_order <- 9

# Create data frame to store results
model_metrics <- data.frame(Order = integer(), RMSE = numeric())
all_preds <- data.frame()
# Loop through model orders
for (order in 1:max_order) {
  
  # Define column names
  ph_lags <- paste0("Ph.l", 1:order)
  tdelta_lags <- paste0("Tdelta.l", 0:(order - 1))
  gv_lags <- paste0("Gv.l", 0:(order - 1))
  
  # Combine into formula
  predictors <- c(ph_lags, tdelta_lags, gv_lags)
  formula_str <- paste("Ph ~", paste(predictors, collapse = " + "))
  formula <- as.formula(formula_str)
  
  # Fit model on training set
  fit <- lm(formula, data = Dtrain)
  
  # Predict on test set
  preds <- predict(fit, newdata = Dtest)
  temp_df <- data.frame(
  tdate = Dtest$tdate,
  Actual = Dtest$Ph,
  Predicted = preds,
  Order = paste0("Order ", order)
)

all_preds <- rbind(all_preds, temp_df)
  
  # Compute RMSE
  rmse_val <- rmse(Dtest$Ph, preds)
  
  # Store results
  model_metrics <- rbind(model_metrics, data.frame(
    Order = order,
    RMSE = rmse_val
  ))
}

ggplot(all_preds, aes(x = tdate)) +
  geom_line(aes(y = Actual), color = "black", size = 0.8) +
  geom_line(aes(y = Predicted, color = Order), linewidth = 0.8) +
  labs(title = "Actual vs Predicted Heating Power (Test Set)",
       x = "Time", y = "Ph [W]", color = "Model Order") +
  theme_minimal(base_size = 14)

  ggplot(model_metrics, aes(x = Order, y = RMSE)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_point(size = 2, color = "darkred") +
  labs(title = "RMSE vs ARX Model Order",
       x = "Model Order", y = "RMSE on Test Set") +
  theme_minimal(base_size = 14)