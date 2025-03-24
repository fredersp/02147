

# part 1
phi_1 <- -0.7
phi_2 <- -0.2

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

# plot the ACF
plot(rho, type = "h", ylim = c(-1,1), ylab = "rho", xlab = "lag", main = "ACF")

# part 2
# (skal der omvendt fortegn på denne, det ser mere rigtigt ud)
sim1 <- arima.sim(model=list(ar=-0.6), n=2000)
acf(sim1, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="ACF", cex.main=3)  
pacf(sim1, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="PACF", cex.main=3)


sim2 <- arima.sim(model=list(ar=c(rep(0,11),0.9)), n=2000)
acf(sim2, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="ACF", cex.main=3)
pacf(sim2, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="PACF", cex.main=3)


sim3 <- arima.sim(model=list(ar=-0.9, ma=c(rep(0,11),0.7)), n=2000)
acf(sim3, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="ACF", cex.main=3)
pacf(sim3, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="PACF", cex.main=3)

sim4 <- arima.sim(model=list(ar=c(0.6,rep(0,10),0.8, -0.48)), n=2000)
acf(sim4, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="ACF", cex.main=3)
pacf(sim4, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="PACF", cex.main=3)


sim5 <- arima.sim(model=list(ma=c(-0.4,rep(0,10),0.8, 0.32)), n=2000)
acf(sim5, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="ACF", cex.main=3)
pacf(sim5, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="PACF", cex.main=3)


sim6 <- arima.sim(model=list(ar=c(rep(0,11),-0.7),(ma=0.4)), n=2000)
acf(sim6, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="ACF", cex.main=3)
pacf(sim6, lag.max=30, lwd=4, cex.axis=1.75, cex.lab=1.75, main="")
title(main="PACF", cex.main=3)

# Remember to flip the signs when using built in arima sim
# 13 lags when seasonality is same order as AR or MA term
# fx (1,0,0) (1,0,0)_12