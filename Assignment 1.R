### Read training data
#! Perhaps you need to set the working directory!?
#setwd("/home/pbac/g/course02417/2025/assignment1")

setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147")
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

# Variance-covariance matrix of global model

# X is the "design matrix"
X <- cbind(1, x)
print(X)

n <- length(e) # number of samples
p <- 2 # Number of parameters


# calculate sum of squared residuals:
RSS <- sum(e^2)

# calculate sigma^2:
sigma2 <- as.numeric(RSS/(n - p))

diag(sigma2, n, n)

# calculate variance-covariance matrix of _parameters_:
V <- sigma2 * solve(t(X) %*% X)
print(V)

# the variances of the parameters are the values in the diagonal:
diag(V)

# and the standard errors are given by:
sqrt(diag(V))

# Variance-covariance matrix of local model
lambda = 0.9
weights <- lambda^((n-1):0)

SIGMA <- diag(n)
diag(SIGMA) <- 1/weights
W <- diag(weights)

V_local <- sigma2 * solve(t(X) %*% W %*% X)


# Conclusion
# Global variance-covariance matrix: 
V
# Local variance-covariance matrix:
V_local


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

Vmatrix_pred_local <- sigma2 * (1 + (Xnew %*% solve(t(X)%*%solve(SIGMA)%*%X)) %*% t(Xnew) )

# prediction interval
y_pred_lwr_wls <- pred_local - qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred_local))
y_pred_upr_wls <- pred_local + qt(0.975, df=n-1)*sqrt(diag(Vmatrix_pred_local))

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
