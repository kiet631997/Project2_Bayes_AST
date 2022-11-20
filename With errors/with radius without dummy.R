library(R2jags)
library(mcmcplots)
library(mcmcse)

data <- read.csv('SB_LAW.csv', header = T)


source('functions.R')

lmod <- lm(logL1 ~ logT1 + logR1, data = data)

errT <- sd(data$logT1e)
obsT <- data$logT1
errR <- sd(data$logR1e)
obsR <- data$logR1
erry <- sd(data$logL1e)
obsy <- data$logL1
N <- nrow(data)

jags_data <- list(
  obsR = obsR, 
  obsT = obsT,
  obsy = obsy,
  errT = errT,
  errR = errR,
  erry = erry,
  N = N
)

JAGS_error <- "model{
beta_0 ~ dnorm(0, 1e-3)
beta_1 ~ dnorm(0, 1e-3)
beta_2 ~ dnorm(0, 1e-3)

tau ~ dgamma(0.5, 0.5)

for(i in 1:N)
{
 R[i] ~ dnorm(0, 1e3)
 T[i] ~ dnorm(0, 1e3)
}

for(i in 1:N)
{
  obsy[i] ~ dnorm(y[i],1/erry^2)
  obsR[i] ~ dnorm(R[i],1/errR^2)
  obsT[i] ~ dnorm(T[i],1/errT^2)
  mu[i] <- beta_0 + beta_1*T[i] + beta_2*R[i]
  y[i] ~ dnorm(mu[i], tau)
}
}"

inits <- function(){
  list(
    beta_0 = coef(lmod)[1],
    beta_1 = coef(lmod)[2],
    beta_2 = coef(lmod)[3]
  )
}


params0 <- c('beta_0', 'beta_1', 'beta_2', 'tau')

JAGS_fit <- jags(data = jags_data,
                 inits = inits,
                 parameters = params0,
                 model = textConnection(JAGS_error),
                 n.chains = 1,
                 n.iter = 3e6,
                 n.thin = 1,
                 n.burnin = 0,
                 jags.seed = 5731
)

ref_post_sample <- as.mcmc(JAGS_fit)

apply(ref_post_sample, 2, ess)

acf(ref_post_sample[,1:3])

apply(ref_post_sample, 2, summary)

par(mfrow = c(2,2))

hist_ci(ref_post_sample[,1], name = bquote(beta[0]))
hist_ci(ref_post_sample[,2], name = bquote(beta[1]))
hist_ci(ref_post_sample[,3], name = bquote(beta[2]))
hist_ci(1/ref_post_sample[,5], name = bquote(sigma^2))

par(mfrow = c(1,1))
