# Part 3 - Predicting hourly heating (W) dependent on temp and solar radiation

library(ggplot2)
library(httpgd)
hgd()


# Set working directory and load data
#setwd("C:/Users/sofie/OneDrive/Time Series Analysis/02147/Assignment 3")
setwd("/Users/nicolinesimonesachmann/Documents/DTU/Times Series Analysis/02147/Assignment 3")

D <- read.csv("box_data_60min.csv")

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

# 3.2

# Spit data into traina and test 
teststart <- as.POSIXct("2013-02-06 00:00", tz = "UTC")
Dtrain <- D[D$tdate < teststart, ]
Dtest <- D[D$tdate >= teststart, ]

# 3.3
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


# 3.4
# 3.4
# Estimate impulse response from Tdelta and Gv to Ph (lags 0 to 10)

# # --- Tdelta model ---
# fit_tdelta <- lm(Ph ~ 0 + Tdelta.l0 + Tdelta.l1 + Tdelta.l2 + Tdelta.l3 + Tdelta.l4 +
#                     Tdelta.l5 + Tdelta.l6 + Tdelta.l7 + Tdelta.l8 + Tdelta.l9 + Tdelta.l10,
#                  data = Dtrain)
# impulse_tdelta <- coef(fit_tdelta)

# # --- Gv model ---
# fit_gv <- lm(Ph ~ 0 + Gv.l0 + Gv.l1 + Gv.l2 + Gv.l3 + Gv.l4 +
#                 Gv.l5 + Gv.l6 + Gv.l7 + Gv.l8 + Gv.l9 + Gv.l10,
#              data = Dtrain)
# impulse_gv <- coef(fit_gv)

# # --- Plot both using your impulse audio plotting style ---
# par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
# lags <- 0:10

# # Tdelta
# plot(lags, impulse_tdelta, type = "h", lwd = 3, col = "red",
#      xlab = "Lag", ylab = "Response", main = "Impulse Response: Tdelta → Ph")
# abline(h = 0, lty = 2)

# # Gv
# plot(lags, impulse_gv, type = "h", lwd = 3, col = "darkgreen",
#      xlab = "Lag", ylab = "Response", main = "Impulse Response: Gv → Ph")
# abline(h = 0, lty = 2)




# Impulse response 
x <- Dtrain$Gv
y <- Dtrain$Ph

# ----------------------------------------------------------------
# Plot'em
# Impulse reponse
plot(x, type="h", xlab="", main="Impulse input")
for(i in 1:10){
  if(i == 1){ plot(y, type="n", xlab="", main="Impulse Response") }
  lines(y, col=i)
}

