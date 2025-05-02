# myKalmanFilter <- function(
#   y,             # Vector of observations y_t
#   theta,       # Model parameters for X_{t+1} = a - b*X_t + c*e_t
#   R,             # Measurement noise variance
#   x_prior = 0,   # Initial prior mean for X_0
#   P_prior = 10   # Initial prior variance for X_0
# ) {
#   a <- theta[1]
#   b <- theta[2]
#   sigma1 <- theta[3]
#   N <- length(y)
#   x_pred  <- numeric(N)  # Predicted means
#   P_pred  <- numeric(N)  # Predicted variances
#   x_filt  <- numeric(N)  # Filtered means
#   P_filt  <- numeric(N)  # Filtered variances
#   innovation     <- numeric(N)  # Pre-fit residuals: y[t] - x_pred[t]
#   innovation_var <- numeric(N)  # Innovation covariance: P_pred[t] + R
  
#   for (t in seq_len(N)) {
#     # the prediction step
#     if (t == 1) {
#       x_pred[t] <- # the mean prediction using the prior
#       P_pred[t] <- # the variance prediction using the prior
#     } else {
#       x_pred[t] <- # the mean prediction using the previous filtered estimate
#       P_pred[t] <- # the variance prediction using the previous filtered estimate
#     }
    
#     # the update step
#     innovation[t] <- # the prediction error
#     innovation_var[t] <- # the prediction error variance
#     K_t <- # the Kalman gain
#     x_filt[t] <- # the filtered estimate
#     P_filt[t] <- # the filtered estimate variance
#   }
  
#   return(list(
#     x_pred = x_pred,
#     P_pred = P_pred,
#     x_filt = x_filt,
#     P_filt = P_filt,
#     innovation = innovation,
#     innovation_var = innovation_var
#   ))
# }


myKalmanFilter <- function(
  y,             # Vector of observations y_t
  theta,         # Model parameters for X_{t+1} = a*X_t + b + e_t
  R,             # Measurement noise variance
  x_prior = 0,   # Initial prior mean for X_0
  P_prior = 10   # Initial prior variance for X_0
) {
  a <- theta[1]
  b <- theta[2]
  sigma1 <- theta[3]
  N <- length(y)

  # Initialize vectors
  x_pred  <- numeric(N)  # Predicted means
  P_pred  <- numeric(N)  # Predicted variances
  x_filt  <- numeric(N)  # Filtered means
  P_filt  <- numeric(N)  # Filtered variances
  innovation     <- numeric(N)  # Innovation = y[t] - x_pred[t]
  innovation_var <- numeric(N)  # Innovation variance

  for (t in seq_len(N)) {
    # Prediction step
    if (t == 1) {
      x_pred[t] <- a * x_prior + b
      P_pred[t] <- a^2 * P_prior + sigma1^2
    } else {
      x_pred[t] <- a * x_filt[t - 1] + b
      P_pred[t] <- a^2 * P_filt[t - 1] + sigma1^2
    }

    # Update step
    innovation[t] <- y[t] - x_pred[t]
    innovation_var[t] <- P_pred[t] + R
    K_t <- P_pred[t] / innovation_var[t]  # Kalman gain

    x_filt[t] <- x_pred[t] + K_t * innovation[t]
    P_filt[t] <- (1 - K_t) * P_pred[t]
  }

  return(list(
    x_pred = x_pred,
    P_pred = P_pred,
    x_filt = x_filt,
    P_filt = P_filt,
    innovation = innovation,
    innovation_var = innovation_var
  ))
}
