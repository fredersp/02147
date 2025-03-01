### Read training data
# ! Perhaps you need to set the working directory!?
# setwd("/home/pbac/g/course02417/2025/assignment1")

# setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147")
D <- read.csv("DST_BIL54.csv")
str(D)

# See the help
?strftime
D$time <- as.POSIXct(paste0(D$time, "-01"), "%Y-%m-%d", tz = "UTC")
D$time
class(D$time)

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz = "UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]


################################################
# 1 Plot data
################################################


# 1.1 Make an array with train set with the year + month/12
x <- 1900 + as.POSIXlt(Dtrain$time)$year + as.POSIXlt(Dtrain$time)$mon / 12

# plot dtotal with x
plot(x, Dtrain$total, xlab = "Year", ylab = "Total (millions)", main = "Total number of vehicles in Denmark", type = "l", col = "blue", lwd = 2)
grid() # Add grid lines
points(x, Dtrain$total, pch = 16, col = "red") # Add points to the plot


################################################
# 2 Linear Trend model
################################################

# 2.2 Fit a linear model
fit <- lm(Dtrain$total ~ x)
summary(fit)
# intercept and slope
# Plot the estimated mean as a line with the observations as points
plot(x, Dtrain$total, xlab="Year", ylab="Total (millions)", main="Total number of vehicles in Denmark", type="p", col="blue", lwd=2)
grid() # Add grid lines
abline(fit, col="red", lwd=2) # Add the fitted line to the plot


# 2.3 Forecast for next 12 months
xnew <- 1900 + as.POSIXlt(Dtest$time)$year + as.POSIXlt(Dtest$time)$mon / 12
pred <- predict(fit, data.frame(x = xnew))

# prediction interval
pred_int <- predict(fit, data.frame(x = xnew), interval = "prediction")
lwr <- pred_int[, 2]
upr <- pred_int[, 3]


# 2.4 Plot the forecast

library(ggplot2)

# Create a dataframe for the observations (training data)
df_train <- data.frame(Year = x, Total = Dtrain$total, Type = "Observed")
df_test <- data.frame(Year = xnew, Total = Dtest$total, Type = "Test")

# Create a dataframe for the predictions
df_pred <- data.frame(Year = xnew, Total = pred, Type = "Predicted")


# Create the plot
library(ggplot2)

ggplot() +
  # Observed data (training set)
  geom_point(data = df_train, aes(x = Year, y = Total, color = "Observed"), size = 3) +  
  geom_line(data = df_train, aes(x = Year, y = Total, color = "Observed"), linetype = "dashed") +  
  
  # Predicted data
  geom_point(data = df_pred, aes(x = Year, y = Total, color = "Predicted"), size = 3) +  
  geom_line(data = df_pred, aes(x = Year, y = Total, color = "Predicted")) +  
  
  # test data
  geom_point(data = df_test, aes(x = Year, y = Total, color = "Test"), size = 3) +
  # Fitted regression line (now correctly red)
  geom_abline(aes(intercept = coef(fit)[1], slope = coef(fit)[2], color = "Fitted Model"), linewidth = 1) +  

  # Confidence interval for predictions
  geom_ribbon(data = df_pred, aes(x = Year, ymin = lwr, ymax = upr, fill = "Prediction Interval"), alpha = 0.2) +
  
  # Labels
  labs(title = "Total Number of Vehicles in Denmark",
       x = "Year", y = "Total (millions)", color = "Legend", fill = "Legend") +
  
  # Define colors for legend
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "green", "Fitted Model" = "red")) +
  scale_fill_manual(values = c("Prediction Interval" = "red")) +
  
  # Minimal theme
  theme_minimal()


# 2.5 Investigate the residuals of the model. Are the model assumptions fulfilled?

# Residuals
e <- fit$residuals

# Plot residuals
plot(x, e, xlab = "Year", ylab = "Residuals", main = "Residuals of the linear model", type = "p", col = "blue", lwd = 2)
abline(h = 0, col = "red", lwd = 2) # Add a horizontal line at 0
grid() # Add grid lines

# qq plot of residuals:
qqnorm(e)
qqline(e)

# ACF plot of residuals
acf(e)

# Assumptions
# - No correlation between residuals (ACF plot)
# - Constant variance of residuals (plot of residuals)
# - Normally distributed residuals (qq plot of residuals)


################################################
# WLS - Local Linear trend model
################################################

# 3.1. Describe the variance-covariance matrix  for the local model and compare it to the variance-covariance matrix of the corresponding global model

# Variance - covariance for global model, all weights are 1 as we value all observations equally
SIGMA <- diag(n)

# Variance-covariance for local model
# Weights are increasing with time as we value recent observations more
lambda = 0.9
weights <- lambda^((n-1):0)

# Create a diagonal matrix with the weights
SIGMA <- diag(n)*1/weights
W = diag(n)*weights


# 3.2 plot the weights vs. time
plot(x, weights, type = "l", xlab = "Time", ylab = "Weights", main = "Weights for local model", col = "blue", lwd = 2)
grid() # Add grid lines

# 3.3 Sum of all lambda weights
sum(weights)

# Sum of weights in OLS model
# 72 (because we have 72 observations all with weight 1)

# 3.4 Estimate and present theta_hat1 and theta_hat2 for the local model with lambda = 0.9

theta_hat <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Dtrain$total
theta_hat

# 3.5. Make a forecast for the next 12 months - i.e., compute predicted values corresponding to the WLS model with λ = 0.9




Xnew <- cbind(1, xnew)
print(Xnew)
pred_local <- Xnew %*% theta_hat


n <- length(e) # number of samples
p <- 2 # Number of parameters

e_wls <- y - yhat_wls
RSS_wls <- t(e_wls)%*%solve(SIGMA)%*%e_wls
sigma2_wls <- as.numeric(RSS_wls/(n - p))
Vmatrix_pred_local <- sigma2 * (1 + (Xnew %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xnew) )

# prediction interval
y_pred_lwr_wls <- y_pred_wls - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred_local))
y_pred_upr_wls <- y_pred_wls + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred_local))


# Create a dataframe for the predictions
df_pred_local <- data.frame(Year = xnew, Total = pred_local, Type = "Predicted")

# Plot
ggplot(data = df_train, aes(x = Year, y = Total)) +
  geom_point(size = 3) +
  geom_point(data = df_pred, aes(x = Year, y = Total), color = "green") + # Points for predicted values
  geom_point(data = df_pred_local, aes(x = Year, y = Total), color = "purple") + # Points predicted values
  geom_line(data = df_train, aes(x = Year, y = Total), color = "blue", linetype = "dashed") + # Line for observed data
  geom_line(data = df_pred, aes(x = Year, y = Total), color = "green") + # Line for predicted data
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "red") + # line for fit
  geom_abline(intercept = theta_hat[1], slope = theta_hat[2], color = "purple") + # line for local fit
  geom_ribbon(data = df_pred, aes(x = Year, ymin = lwr, ymax = upr), fill = "red", alpha = 0.2) +
  geom_ribbon(data = df_pred_local, aes(x = Year, ymin = y_pred_lwr_wls, ymax = y_pred_upr_wls), fill = "purple", alpha = 0.2) +
  labs(
    title = "Total Number of Vehicles in Denmark",
    x = "Year", y = "Total (millions)"
  ) +
  # scale_color_manual(values = c("Observed" = "blue", "Predicted" = "green")) +
  theme_minimal()

################################################
# RLS - Recursive least squares
################################################

# 4.1 On paper

# 4.2
# Implement the update equations in a for-loop in a computer. Calculate the theta_hat_t up to time t = 3.

p = 2
n = 3

R <- diag(2) * 0.1  # Ensures initial invertibility
Theta <- matrix(0, nrow = n, ncol = p)  # Placeholder for parameter estimates
Theta_0 <- c(0, 0)

for (t in 1:n) {
  xvals <- X[t, ]
  y <- Dtrain$total[t]

  if(t==1){
    R <- R + xvals %*% t(xvals)
    Theta[t, ] <- Theta_0 + solve(R) %*% xvals %*% (y - t(xvals) %*% Theta_0)

  }
  else{

  # Compute R and update with forgetting factor
  R <- R + xvals %*% t(xvals)

  # Update Theta
  Theta[t, ] <- Theta[t-1, ] + solve(R) %*% xvals %*% (y - t(xvals) %*% Theta[t-1, ])
}
}

Theta

# 4.3
# Calculate the estimates of ˆθ_N and compare them to the OLS estimates of θ, are they close? Can
# you find a way to decrease the difference by modifying some initial values and explain why initial
# values are important to get right?

RLS <- function(X, target, p = 2, Theta_0 = c(0, 0)) {

  n = length(X[,1])

  R <- diag(2) * 0.1  # Ensures initial invertibility
  Theta <- matrix(0, nrow = n, ncol = p)  # Placeholder for parameter estimates

  for (t in 1:n) {
    xvals <- X[t, ]
    y <- target[t]

    if(t == 1){
      R <- R + xvals %*% t(xvals)

      Theta[t, ] <- Theta_0 + solve(R) %*% xvals %*% (y - t(xvals) %*% Theta_0)

  }
    else{

    # Compute R and update with forgetting factor
    R <- R + xvals %*% t(xvals)

    # Update Theta
    Theta[t, ] <- Theta[t-1, ] + solve(R) %*% xvals %*% (y - t(xvals) %*% Theta[t-1, ])
}

}

return(Theta)

}

RLS(X, Dtrain$total, p = 2, Theta_0 = c(-110, 0))
coef(fit)

# Note: Man skal være obs på intercept, hvis man vil prøve at ramme noget der minder om den globale OLS


# 4.4
# Implement the RLS algorithm with a forgetting factor λ = 0.99 and λ = 0.7. Calculate the estimates of ˆθ_N for both values of λ.

RLS_lambda <- function(X, target, p = 2, Theta_0 = c(0, 0), lambda = 0.99) {

  n = length(X[,1])

  R <- diag(2) * 0.1  # Ensures initial invertibility
  Theta <- matrix(0, nrow = n, ncol = p)  # Placeholder for parameter estimates

  for (t in 1:n) {
    xvals <- X[t, ]
    y <- target[t]

    if(t == 1){
      R <- lambda * R + xvals %*% t(xvals)

      Theta[t, ] <- Theta_0 + solve(R) %*% xvals %*% (y - t(xvals) %*% Theta_0)

  }
    else{

    # Compute R and update with forgetting factor
    R <- lambda * R + xvals %*% t(xvals)

    # Update Theta
    Theta[t, ] <- Theta[t-1, ] + solve(R) %*% xvals %*% (y - t(xvals) %*% Theta[t-1, ])
}

}

return(Theta)

}

RLS_lambda(X, Dtrain$total, p = 2, Theta_0 = c(-110, 0), lambda = 0.99)
RLS_lambda(X, Dtrain$total, p = 2, Theta_0 = c(-110, 0), lambda = 0.7)

# Compare with WLS estimates

# 4.5
# Make one step-ahead predictions
n <- length(X[, 1])

# Prediction for lambda = 0.7
OneStepPred_1 <- matrix(NA, nrow = n)
Theta_1 <- RLS_lambda(X, Dtrain$total, p = 2, Theta_0 = c(-110, 0), lambda = 0.7)
for (t in 5:(n - 1)) {
  OneStepPred_1[t + 1] <- X[t + 1, ] %*% Theta_1[t, ]
}


# Prediction for lambda = 0.99
OneStepPred_2 <- matrix(NA, nrow = n)
Theta_2 <- RLS_lambda(X, Dtrain$total, p = 2, Theta_0 = c(-110, 0), lambda = 0.99)
for (t in 5:(n - 1)) {
  OneStepPred_2[t + 1] <- X[t + 1, ] %*% Theta_2[t, ]
}


# Plot actual total number of vehicles
plot(x, Dtrain$total, xlab = "Year", ylab = "Total (millions)", 
     main = "Total number of vehicles in Denmark", type = "l", 
     col = "blue", lwd = 2)

# Add grid lines
grid()

# Add predictions
lines(x, as.vector(OneStepPred_1), col = "green", lwd = 2)
lines(x, as.vector(OneStepPred_2), col = "purple", lwd = 2)

# Add a legend
legend("topleft", legend = c("Actual", "Lambda = 0.7", "Lambda = 0.99"),
       col = c("blue", "green", "purple"), lwd = 2, bty = "n")


# Calculate OneStepAhead residuals
OneStepRes_1 <- Dtrain$total[5:(n - 1)] - as.vector(OneStepPred_1[5:(n - 1)])
OneStepRes_2 <- Dtrain$total[5:(n - 1)] - as.vector(OneStepPred_2[5:(n - 1)])

# Plot residuals
plot(x[5:(n - 1)], OneStepRes_1, xlab = "Year", ylab = "Residuals", 
     main = "One-step-ahead residuals for lambda = 0.7", type = "p", 
     col = "blue", lwd = 2)


plot(x[5:(n - 1)], OneStepRes_2, xlab = "Year", ylab = "Residuals",
     main = "One-step-ahead residuals for lambda = 0.99", type = "p",
     col = "blue", lwd = 2)

# Note: Ligner random residuals med forgetting faktor, men der er et mønster i residuals når der ikke er


# 4.6
# Optimize the forgetting for the horizons k = 1, . . . , 12. First calculate the k-step residuals
# then calculate the k-step Root Mean Square Error (RMSE)
# Do this for a sequence of λ values (e.g. 0.5,0.51,. . . ,0.99) and make a plot.
# Comment on: Is there a pattern and how would you choose an optimal value of λ? Would you
# let λ depend on the horizon?

compute_k_step_rmse <- function(X, y, k_values, lambda_values) {
  results <- expand.grid(lambda = lambda_values, k = k_values)
  results$RMSE <- NA

  for (i in seq_along(lambda_values)) {
    lambda <- lambda_values[i]
    
    # Apply RLS
    Theta <- RLS_lambda(X, y, p = 2, Theta_0 = c(-110, 0), lambda = lambda)
    
    for (k in k_values) {
      # k-step-ahead predictions
      n <- length(y)
      OneStepPred <- rep(NA, n)
      
      for (t in (k + 1):n) {
        OneStepPred[t] <- sum(X[t, ] * Theta[t - k, ])  # Forecast k steps ahead
      }
      
      # Compute residuals
      residuals <- y - OneStepPred
      rmse <- sqrt(mean(residuals[(k + 1):n]^2, na.rm = TRUE))
      
      # Store RMSE
      results$RMSE[results$lambda == lambda & results$k == k] <- rmse
    }
  }
  
  return(results)
}

lambda_values <- seq(0.5, 0.99, by = 0.01)
k_values <- 1:12  # Forecast horizons

# Compute RMSE for each (lambda, k)
rmse_results <- compute_k_step_rmse(X, Dtrain$total, k_values, lambda_values)
rmse_results

library(ggplot2)
ggplot(rmse_results, aes(x = lambda, y = RMSE, color = as.factor(k))) +
  geom_line(size = 1) +
  labs(title = "RMSE vs. Forgetting Factor (λ) for Different Horizons",
       x = "Forgetting Factor (λ)", y = "Root Mean Square Error (RMSE)", 
       color = "Forecast Horizon (k)") +
  theme_minimal()

# Note: Jo længere forecast horizon, jo større lambda værdi, altså mindre forgetting factor og omvendt for en kort prediction horizon
# Men mere præcis med kort horizon generelt

# For each k, the RMSE is minimized for a different value of λ, make a list of the optimal λ for each k from 1-12
for(i in k_values){
  print(rmse_results[rmse_results$k == i,][which.min(rmse_results[rmse_results$k == i,]$RMSE),])
}

# Define optimal lambda values for different horizons (k)
opt_lambda <- c(0.5, 0.5, 0.58, 0.59, 0.59, 0.58, 0.74, 0.74, 0.75, 0.99, 0.99, 0.99)

# Number of observations
n <- length(Dtrain$total)

# Initialize a dataframe to store k-step predictions
PredictionMatrix <- data.frame(Year = x, Actual = Dtrain$total)

# Loop through different forecast horizons
for (i in seq_along(k_values)) {  
  k <- k_values[i]  # Forecast horizon
  lambda_k <- opt_lambda[i]  # Optimal lambda for this horizon
  
  # Apply RLS with optimal lambda
  Theta <- RLS_lambda(X, Dtrain$total, p = 2, Theta_0 = c(-110, 0), lambda = lambda_k)
  
  # Initialize predictions for this horizon
  OneStepPred <- rep(NA, n)
  
  # Compute k-step-ahead predictions
  for (t in (k + 1):n) {
    OneStepPred[t] <- sum(X[t, ] * Theta[t - 1, ])  
  }
  
  # Store predictions in dataframe
  PredictionMatrix[[paste0("k_", k, "_step")]] <- OneStepPred
  
  # Print progress
  print(paste("Stored predictions for Horizon k =", k, "with λ =", lambda_k))
}

# Display the stored predictions
PredictionMatrix

# Make one plot for each k-step-ahead prediction with legends and the actual data
library(tidyr)
library(ggplot2)

# Reshape the dataframe for plotting
PredictionMatrix_long <- gather(PredictionMatrix, key = "Horizon", value = "Prediction", -Year, -Actual)

# Plot
ggplot(data = PredictionMatrix_long, aes(x = Year, y = Prediction, color = Horizon)) +
  geom_line(size = 1) +
  geom_line(data = PredictionMatrix_long, aes(x = Year, y = Actual, linetype = "Actual"), color = "black", size = 1) +
  scale_linetype_manual(name = "Legend", values = c("Actual" = "dashed")) +
  labs(title = "k-step-ahead Predictions with Optimal λ for Different Horizons",
       x = "Year", y = "Total (millions)", color = "Forecast Horizon (k)") +
  theme_minimal()

# Plot OLS and WLS
ggplot() +
  # Observed data (training set)
  geom_point(data = df_train, aes(x = Year, y = Total, color = "Observed"), size = 3) +  
  geom_line(data = df_train, aes(x = Year, y = Total, color = "Observed"), linetype = "dashed") +  
  
  # Predicted data
  geom_point(data = df_pred, aes(x = Year, y = Total, color = "Predicted"), size = 3) +  
  geom_line(data = df_pred, aes(x = Year, y = Total, color = "Predicted")) +  
  
  # test data
  geom_point(data = df_test, aes(x = Year, y = Total, color = "Test"), size = 3) +
  # Fitted regression line (now correctly red)
  geom_abline(aes(intercept = coef(fit)[1], slope = coef(fit)[2], color = "Fitted Model"), linewidth = 1) +  

  # Confidence interval for predictions
  geom_ribbon(data = df_pred, aes(x = Year, ymin = lwr, ymax = upr, fill = "Prediction Interval"), alpha = 0.2) +
  
  # Labels
  labs(title = "Total Number of Vehicles in Denmark",
       x = "Year", y = "Total (millions)", color = "Legend", fill = "Legend") +
  
  # Define colors for legend
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "green", "Fitted Model" = "red")) +
  scale_fill_manual(values = c("Prediction Interval" = "red")) +
  
  # Minimal theme
  theme_minimal()




