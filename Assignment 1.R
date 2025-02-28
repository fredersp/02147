### Read training data
#! Perhaps you need to set the working directory!?
#setwd("/home/pbac/g/course02417/2025/assignment1")

#setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147")
D <- read.csv("DST_BIL54.csv")
str(D)

# See the help
?strftime
D$time <- as.POSIXct(paste0(D$time,"-01"), "%Y-%m-%d", tz="UTC")
D$time
class(D$time)

## Year to month for each of them
D$year <- 1900 + as.POSIXlt(D$time)$year + as.POSIXlt(D$time)$mon / 12

## Make the output variable a floating point (i.e.\ decimal number)
D$total <- as.numeric(D$total) / 1E6

## Divide intro train and test set
teststart <- as.POSIXct("2024-01-01", tz="UTC")
Dtrain <- D[D$time < teststart, ]
Dtest <- D[D$time >= teststart, ]


################################################
# 1 Plot data
################################################


# 1.1 Make an array with train set with the year + month/12
x <- 1900 + as.POSIXlt(Dtrain$time)$year + as.POSIXlt(Dtrain$time)$mon / 12

# plot dtotal with x
plot(x, Dtrain$total, xlab="Year", ylab="Total (millions)", main="Total number of vehicles in Denmark", type="l", col="blue", lwd=2)
grid() # Add grid lines
points(x, Dtrain$total, pch=16, col="red") # Add points to the plot


################################################
# 2 Linear Trend model
################################################

# 2.2 Fit a linear model
fit <- lm(Dtrain$total ~ x)
summary(fit)
# Plot the estimated mean as a line with the observations as points
plot(x, Dtrain$total, xlab="Year", ylab="Total (millions)", main="Total number of vehicles in Denmark", type="p", col="blue", lwd=2)
grid() # Add grid lines
abline(fit, col="red", lwd=2) # Add the fitted line to the plot


# 2.3 Forecast for next 12 months
xnew <- 1900 + as.POSIXlt(Dtest$time)$year + as.POSIXlt(Dtest$time)$mon / 12
pred <- predict(fit, data.frame(x=xnew))

# prediction interval
pred_int <- predict(fit, data.frame(x=xnew), interval = "prediction")
lwr = pred_int[,2]
upr = pred_int[,3]


# 2.4 Plot the forecast

library(ggplot2)

# Create a dataframe for the observations (training data)
df_train <- data.frame(Year = x, Total = Dtrain$total, Type = "Observed")

# Create a dataframe for the predictions
df_pred <- data.frame(Year = xnew, Total = pred, Type = "Predicted")

# Plot
ggplot(data = df_train, aes(x=Year, y = Total)) +
  geom_point(size = 3) + geom_point(data = df_pred, aes(x=Year, y =Total), color = "green") +  # Points for both observed and predicted values
  geom_line(data = df_train, aes(x = Year, y = Total), color = "blue", linetype = "dashed") +  # Line for observed data
  geom_line(data = df_pred, aes(x = Year, y = Total), color = "green") +  # Line for predicted data
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "red") + # line for fit
  geom_ribbon(data = df_pred, aes(x = Year, ymin = lwr, ymax = upr), fill = "red", alpha = 0.2) +
  labs(title = "Total Number of Vehicles in Denmark",
       x = "Year", y = "Total (millions)") +
  #scale_color_manual(values = c("Observed" = "blue", "Predicted" = "green")) +
  theme_minimal()

# 2.5 Investigate the residuals of the model. Are the model assumptions fulfilled?

# Residuals
e <- fit$residuals

# Plot residuals
plot(x, e, xlab = "Year", ylab = "Residuals", main = "Residuals of the linear model", type = "p", col="blue", lwd=2)
abline(h=0, col="red", lwd=2) # Add a horizontal line at 0
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
plot(x,weights, type="l", xlab="Time", ylab="Weights", main="Weights for local model", col="blue", lwd=2)
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
ggplot(data = df_train, aes(x=Year, y = Total)) +
  geom_point(size = 3) + geom_point(data = df_pred, aes(x=Year, y =Total), color = "green") +  # Points for predicted values
  geom_point(data = df_pred_local, aes(x=Year, y =Total), color = "purple") +  # Points predicted values
  geom_line(data = df_train, aes(x = Year, y = Total), color = "blue", linetype = "dashed") +  # Line for observed data
  geom_line(data = df_pred, aes(x = Year, y = Total), color = "green") +  # Line for predicted data
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color = "red") + # line for fit
  geom_abline(intercept = theta_hat[1], slope = theta_hat[2], color = "purple") + # line for local fit
  geom_ribbon(data = df_pred, aes(x = Year, ymin = lwr, ymax = upr), fill = "red", alpha = 0.2) +
  geom_ribbon(data = df_pred_local, aes(x = Year, ymin = y_pred_lwr_wls, ymax = y_pred_upr_wls), fill = "purple", alpha = 0.2) +
  labs(title = "Total Number of Vehicles in Denmark",
       x = "Year", y = "Total (millions)") +
  #scale_color_manual(values = c("Observed" = "blue", "Predicted" = "green")) +
  theme_minimal()

################################################
# RLS - Recursive least squares
################################################
# 4.2
# Implement the update equations in a for-loop in a computer. Calculate the ˆθt up to time t = 3.

# Initialize
theta_hat <- matrix(0, nrow = 2, ncol = 3)
P <- diag(2)

# Loop
for (t in 1:3) {
  x_t <- matrix(X[t,], ncol = 1)  # Ensure x_t is a column vector
  y_t <- Dtrain$total[t]
  
  # Update
  K <- P %*% x_t / as.numeric(1 + t(x_t) %*% P %*% x_t)
  theta_hat[,t] <- theta_hat[,t] + K * (y_t - as.numeric(t(x_t) %*% theta_hat[,t]))
  P <- P - K %*% t(x_t) %*% P
}
theta_hat

# Do you think it is intuitive to understand the details in the matrix calculations? If yes, give a short explanaition?
# Mangler lidt foklaring her.

# 4.3
# Calculate the estimates of ˆθ_N and compare them to the OLS estimates of θ, are they close? 
theta_hat <- matrix(0, nrow = 2, ncol = n)
P <- diag(2)
theta_hat[,1] <- c(0, 0)  # Initialize the first value
for (t in 2:n) {
  x_t <- matrix(X[t,], ncol = 1)  # Ensure x_t is a column vector
  y_t <- Dtrain$total[t]
  
  # Update
  K <- P %*% x_t / as.numeric(1 + t(x_t) %*% P %*% x_t)
  theta_hat[,t] <- theta_hat[,t-1] + K * (y_t - as.numeric(t(x_t) %*% theta_hat[,t-1]))
  P <- P - K %*% t(x_t) %*% P
}
theta_hat

# Compare to OLS
coef(fit)

# Can you find a way to decrease the difference by modifying some initial values and explain why initial
# values are important to get right?
# Conclusion: The estimates are close, but not the same, you can choose intial values that is close to OLS
# to minimize the burn in period, but even though it does not come close because we maybe see too much difference
# in the evolution of the data. 

# 4.4
# Now implement RLS with forgetting (you just have to multiply with λ at one position in the Rt update)

# Initialize
theta_hat <- matrix(0, nrow = 2, ncol = n)
P <- diag(2)
theta_hat[,1] <- c(-110, 0.056)  # Initialize the first value
lambda <- 0.9
for (t in 2:n) {
  x_t <- matrix(X[t,], ncol = 1)  # Ensure x_t is a column vector
  y_t <- Dtrain$total[t]
  
  # Update
  K <- P %*% x_t / as.numeric(1 + t(x_t) %*% P %*% x_t)
  theta_hat[,t] <- theta_hat[,t-1] + K * (y_t - as.numeric(t(x_t) %*% theta_hat[,t-1]))
  P <- (1/lambda) * (P - K %*% t(x_t) %*% P)
}
# Parameter estimates
theta_hat


rls_with_lambda <- function(X, y, lambda = 0.7, theta_init = c(0, 0)) {
  # Number of observations and features
  n <- nrow(X)  # Number of time steps
  p <- ncol(X)  # Number of parameters
  
  # Store results for different lambda values
  results <- list()

  # Initialize theta and P
  theta_hat <- matrix(0, nrow = p, ncol = n)  # Each column stores theta at time t
  theta_hat[, 1] <- theta_init  # Set initial theta
  P <- diag(p) * 1000  # Large initial value for P to ensure convergence
    
    # Recursive Least Squares loop
  for (t in 2:n) {
      x_t <- matrix(X[t, ], ncol = 1)  # Ensure x_t is a column vector
      y_t <- y[t]  # Output at time t
      
      # Compute gain (K)
      K <- (P %*% x_t) / as.numeric(1 + t(x_t) %*% P %*% x_t)
      
      # Update theta estimate
      theta_hat[, t] <- theta_hat[, t - 1] + K * (y_t - as.numeric(t(x_t) %*% theta_hat[, t - 1]))
      
      # Update covariance matrix with forgetting factor
      P <- (1 / lambda) * (P - K %*% t(x_t) %*% P)
    }
    
    # Store results
    results[[as.character(lambda)]] <- theta_hat
  
  return(results)  # Return all theta estimates
}

# Run the function
theta_estimates1 <- rls_with_lambda(X, Dtrain$total, lambda = 0.7, theta_init = c(0, 0))
theta_estimates2 <- rls_with_lambda(X, Dtrain$total, lambda = 0.99, theta_init = c(0, 0))

# Plot the theta 1 estimates for each lambda value
plot(x, theta_estimates1$`0.7`[1, ], type = "l", xlab = "Time", ylab = "Theta 1 estimates", main = "Theta estimates for different lambda values", col = "blue", lwd = 2)
lines(x, theta_estimates2$`0.99`[1, ], col = "red", lwd = 2)
legend("topright", legend = c("Lambda = 0.7", "Lambda = 0.99"), col = c("blue", "red"), lty = 1, lwd = 2)
# plot the theta 2 estimates for each lambda value
plot(x, theta_estimates1$`0.7`[2, ], type = "l", xlab = "Time", ylab = "Theta 2 estimates", main = "Theta estimates for different lambda values", col = "blue", lwd = 2)
lines(x, theta_estimates2$`0.99`[2, ], col = "red", lwd = 2)
legend("topright", legend = c("Lambda = 0.7", "Lambda = 0.99"), col = c("blue", "red"), lty = 1, lwd = 2)

# 4.5 make one step ahead predictions
# One-step-ahead prediction using theta_estimates1 (lambda = 0.7)
one_step_ahead_pred <- numeric(n)

for (t in 2:n) {
  x_t <- matrix(X[t, ], ncol = 1)  # Ensure x_t is a column vector
  one_step_ahead_pred[t] <- as.numeric(t(theta_estimates1$`0.7`[, t - 1]) %*% x_t)
}

# One-step ahead prediction for lambda = 0.99
one_step_ahead_pred2 <- numeric(n)

for (t in 2:n) {
  x_t <- matrix(X[t, ], ncol = 1)  # Ensure x_t is a column vector
  one_step_ahead_pred2[t] <- as.numeric(t(theta_estimates2$`0.99`[, t - 1]) %*% x_t)
}

# calculate residuals
residuals1 <- Dtrain$total - one_step_ahead_pred
residuals2 <- Dtrain$total - one_step_ahead_pred2

# Determine the range of residuals for auto-scaling
residuals_range <- range(c(residuals1[5:n], residuals2[5:n]))

# Plot residuals from time 5 to n with auto-scaled y-axis
plot(x[5:n], residuals1[5:n], xlab = "Time", ylab = "Residuals", main = "Residuals of the RLS model", type = "p", col = "blue", lwd = 2, ylim = residuals_range)
points(x[5:n], residuals2[5:n], col = "red", lwd = 2)
abline(h = 0, col = "black", lwd = 2)
legend("topright", legend = c("Lambda = 0.7", "Lambda = 0.99"), col = c("blue", "red"), lty = 1, lwd = 2)
grid() # Add grid lines




# plot predictions with original data
# Create a dataframe for the predictions
df_pred_rls1 <- data.frame(Year = x, Total = one_step_ahead_pred, Type = "Predicted")
df_pred_rls2 <- data.frame(Year = x, Total = one_step_ahead_pred2, Type = "Predicted")

# Plot
ggplot(data = df_train, aes(x=Year, y = Total)) +
  geom_point(size = 3) + geom_point(data = df_pred_rls1, aes(x=Year, y =Total), color = "blue") +  # Points for predicted values
  geom_point(data = df_pred_rls2, aes(x=Year, y =Total), color = "red") +  # Points for predicted values
  geom_line(data = df_train, aes(x = Year, y = Total), color = "blue", linetype = "dashed") +  # Line for observed data
  geom_line(data = df_pred_rls1, aes(x = Year, y = Total), color = "blue") +  # Line for predicted data
  geom_line(data = df_pred_rls2, aes(x = Year, y = Total), color = "red") +  # Line for predicted data
  labs(title = "Total Number of Vehicles in Denmark",
       x = "Year", y = "Total (millions)") +
  #scale_color_manual(values = c("Observed" = "blue", "Predicted" = "green")) +
  theme_minimal()



