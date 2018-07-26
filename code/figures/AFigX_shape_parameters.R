


# Parameters
# 45% B0: p=0.55
# 40% B0: p=0.2
# 37% B0: p=0.01
r <- 0.1
k <- 1
p <- c(0.55, 0.2, 0.01)

# Simulate productivity
x <- seq(0.0, 1, 0.005)
y <- r/p[2]*x*(1-(x/k)^p[2])
b_maxp <- x[which.max(y)]
b_maxp




par(mfrow=c(1,1))
plot(x, y, type="l", bty="n", xlab="Biomass (% of maximum)", ylab="Production")

r*k/4

(r*k)/((p+1)^((p+1)/p))



