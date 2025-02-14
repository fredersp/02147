### Read training data
#! Perhaps you need to set the working directory!?
#setwd("/home/pbac/g/course02417/2025/assignment1")
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

# Fit a linear model
fit <- lm(Dtrain$total ~ x)
summary(fit)
plot(x, Dtrain$total, xlab="Year", ylab="Total (millions)", main="Total number of vehicles in Denmark", type="l", col="blue", lwd=2)
abline(fit, col="red", lwd=2)

