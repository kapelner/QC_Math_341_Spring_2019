library(tidyverse)
set.seed(1984)
n = 300

true_lambda_1 = 0.7
true_lambda_2 = 7
true_rho = 0.3


x = c(rpois(n * true_rho, true_lambda_1), rpois(n * (1 - true_rho), true_lambda_2))

ggplot(data.frame(x = x)) + geom_bar(aes(x = x)) + xlab("number of calls in a given day")
#plot


#chains
S = 10000
lambda1s = array(NA, S)
lambda2s = array(NA, S)
rhos = array(NA, S)
Is = matrix(NA, nrow = n, ncol = S)
#start positions
lambda1s[1] = mean(x)
lambda2s[1] = mean(x)
rhos[1] = 0.5
Is[, 1] = 0.5 

for (t in 2 : S){
  lambda1 = lambda1s[t - 1]
  lambda2 = lambda2s[t - 1]
  rho = rhos[t - 1]
  I = Is[, t - 1]
  
  sum_I = sum(I)
  sum_1_min_I = n - sum_I
  lambda1s[t] = rgamma(1, sum(I * x) + 1, sum_I)
  lambda2s[t] = rgamma(1, sum((1 - I) * x) + 1, sum_1_min_I)
  
  for (i in 1 : n){#now draw the Is
    a = exp(-lambda1s[t]) * lambda1s[t]^x[i] * rho
    b = exp(-lambda2s[t]) * lambda2s[t]^x[i] * (1 - rho)
    Is[i, t] = rbinom(1, 1, a / (a + b))
    # cat("a =", a, "b = ", b, "p =", a / (a + b), "I =", Is[i, t], "\n")
  }
  rhos[t] = rbeta(1, 1 + sum_I, 1 + sum_1_min_I)
  # cat("t =", t, "sum_I = ", sum_I, "sum_1_min_I =", sum_1_min_I, "rhos[t] =", rhos[t], "\n")
}



###assess convergence
#plot
B = 20
###
par(mfrow = c(3, 1))
S0 = 50
plot(1 : S0, lambda1s[1 : S0], ylab = "lambda_1")
# abline(h = mean(lambda1s[B : S0]), col = "blue")
# abline(h = true_theta_1, col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, lambda2s[1 : S0], ylab = "lambda_2")
# abline(h = mean(lambda2s[B : S0]), col = "blue")
# abline(h = true_theta_2, col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, rhos[1 : S0], ylab = "rho")
# abline(h = mean(rhos[B : S0]), col = "blue")
# abline(h = sqrt(true_rho), col = "red")
# abline(v = B, col = "grey")
#plot

##assess autocorrelation

par(mfrow = c(3, 1))
Kmax = 10
acf(lambda1s[B : S], xlim = c(0, Kmax), lag.max = Kmax)
acf(lambda2s[B : S], xlim = c(0, Kmax), lag.max = Kmax)
acf(rhos[B : S], xlim = c(0, Kmax), lag.max = Kmax)
THIN = 25
#plot

#burn and thin
lambda1s = lambda1s[B : S]
lambda1s = lambda1s[seq(1, S - B, by = THIN)]
lambda2s = lambda2s[B : S]
lambda2s = lambda2s[seq(1, S - B, by = THIN)]
rhos = rhos[B : S]
rhos = rhos[seq(1, S - B, by = THIN)]
Is = Is[, B : S]
Is = Is[, seq(1, S - B, by = THIN)]

#look at posteriors with post-exp at 95% CI
par(mfrow = c(1, 1))
res = 100

ggplot(data.frame(lambda1s = lambda1s)) + geom_histogram(aes(x = lambda1s), bins = res)
ggplot(data.frame(lambda2s = lambda2s)) + geom_histogram(aes(x = lambda2s), bins = res)
ggplot(data.frame(rhos = rhos)) + geom_histogram(aes(x = rhos), bins = res)
ggplot(data.frame(Is = factor(Is[31, ]))) + geom_bar(aes(x = Is))
ggplot(data.frame(Is = factor(Is[102, ]))) + geom_bar(aes(x = Is))

round(rbind(lambda1s, lambda2s, rhos)[, 1: 20], 3)




# abline(v = mean(lambda1s), col = "blue", lwd = 3)
# abline(v = quantile(lambda1s, 0.025), col = "grey", lwd = 3)
# abline(v = quantile(lambda1s, 0.975), col = "grey", lwd = 3)
# abline(v = true_theta_1, col = "red", lwd = 3)

# hist(lambda2s, br = res)
# abline(v = mean(lambda2s), col = "blue", lwd = 3)
# abline(v = quantile(lambda2s, 0.025), col = "grey", lwd = 3)
# abline(v = quantile(lambda2s, 0.975), col = "grey", lwd = 3)
# abline(v = true_theta_2, col = "red", lwd = 3)
# 
# hist(rhos, br = res)
# abline(v = mean(rhos), col = "blue", lwd = 3)
# abline(v = quantile(rhos, 0.025), col = "grey", lwd = 3)
# abline(v = quantile(rhos, 0.975), col = "grey", lwd = 3)
# abline(v = true_rho, col = "red", lwd = 3)
# #plot