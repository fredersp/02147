library(ggplot2)
# for better ACF extraction (optional)

# Assignment 3: Part 1.1
set.seed(123)
rnorm(5)

phi_1 <- -0.6
phi_2 <- 0.5

# Set up plotting layout for 5 plots in one figure
par(mfrow = c(5, 1), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

# Function to plot each realization
plotrealization <- function(x, i){
  plot(x, ylab = "X", main = paste("Simulation", i))
}

# Generate and plot 5 simulations
for(i in 1:5){
  ts_data <- arima.sim(model = list(ar = c(0.6, -0.5)), n = 200)
  plotrealization(ts_data, i)
}

# Part 1.2: Calculate the theoretical ACF and compare with the empirical ACF
rho_1 <- -phi_1/(phi_2+1)
rho_2 <- -phi_1 * rho_1 - phi_2

# loop to calculate rho for k=2 to 30
rho <- numeric(30)

rho[1] <- rho_1
rho[2] <- rho_2

for (k in 3:30){
  rho[k] <- -phi_1 * rho[k-1] - phi_2 * rho[k-2]
}

rho
# add rho_0 = 1
rho <- c(1, rho)


# Set up a 5x2 layout: each row is one simulation, ACF on left, PACF on right
par(mfrow = c(5, 2), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))

# Plot ACF and PACF for 5 simulations 
for(i in 1:5){
    ts_data <- arima.sim(model = list(ar = c(0.6, -0.5)), n = 200)
    acf(ts_data, lag.max = 30, lwd = 2, main = paste("ACF - Empirical", i))
    plot(0:30, rho, type = "h", lwd = 2, ylim = c(-1, 1),
       xlab = "Lag", ylab = "ACF", main = paste("ACF - Theoretical", i))
  points(0:30, rho, pch = 16)  # add dots at the top of each spike
  abline(h = 0)
}
