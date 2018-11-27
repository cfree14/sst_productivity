
par(mfrow=c(2,3))
sd <- 0.15
plot(arima.sim(model=list(ar=0.1), n=40, sd=sd))
plot(arima.sim(model=list(ar=0.5), n=40, sd=sd))
plot(arima.sim(model=list(ar=0.99), n=40, sd=sd))
plot(arima.sim(model=list(ar=-0.1), n=40, sd=sd))
plot(arima.sim(model=list(ar=-0.5), n=40, sd=sd))
plot(arima.sim(model=list(ar=-0.99), n=40, sd=sd))