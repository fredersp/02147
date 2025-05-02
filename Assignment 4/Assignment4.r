library(httpgd)
hgd() 
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
