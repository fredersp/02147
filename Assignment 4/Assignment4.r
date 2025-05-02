#library(httpgd)
#hgd() 

# Task 1.1.1

# Parameters
a <- 0.9
b <- 1        # bias term (not dynamic input)
sigma1 <- 1
n <- 100
X0 <- 5

# Simulate 5 independent trajectories
set.seed(42)
X_list <- list()

for (j in 1:5) {
  X <- numeric(n)
  X[1] <- X0
  for (t in 2:n) {
    X[t] <- a * X[t - 1] + b + rnorm(1, mean = 0, sd = sigma1)
  }
  X_list[[j]] <- X
}

# Plot
matplot(1:n, do.call(cbind, X_list), type = "l", lty = 1, col = 1:5,
        xlab = "Time", ylab = "State X_t", main = "5 Simulated Realizations of X_t")
legend("topleft", legend = paste("Run", 1:5), col = 1:5, lty = 1)


# Task 1.2

# Use first trajectory from Task 1.1
X <- X_list[[1]]

# Observation noise
sigma2 <- 1

# Simulate noisy observations
Y <- X + rnorm(length(X), mean = 0, sd = sigma2)

# Plot state and observation
plot(X, type = "l", col = "blue", lwd = 2,
     ylim = range(c(X, Y)), xlab = "Time", ylab = "Value",
     main = "Latent State $X_t$ and Noisy Observations $Y_t$")
lines(Y, col = "red", lwd = 1.5, lty = 2)
legend("topleft", legend = c("State X_t", "Observation Y_t"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)


# Task 1.3
