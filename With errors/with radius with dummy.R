library(R2jags)
library(mcmcplots)
library(mcmcse)
library(SimTools)

data <- read.csv('SB_LAW.csv', header = T)

source('functions.R')

data$dummy <- ifelse(data$logL1<6*data$logT1 - 21, 'MC', 'RG')
data$dummy <- as.factor(data$dummy)

lmod <- lm(logL1 ~ logT1 + logR1 + dummy, data = data)

errT <- sd(data$logT1e)
obsT <- data$logT1
errR <- sd(data$logR1e)
obsR <- data$logR1
erry <- sd(data$logL1e)
obsy <- data$logL1
dummy <- data$dummy
N <- nrow(data)

jags_data <- list(
  obsR = obsR, 
  obsT = obsT,
  obsy = obsy,
  errT = errT,
  errR = errR,
  erry = erry,
  dummy = dummy,
  N = N
)

JAGS_error <- "model{
beta_0 ~ dnorm(0, 1e-3)
beta_1 ~ dnorm(0, 1e-3)
beta_2 ~ dnorm(0, 1e-3)
beta_3 ~ dnorm(0, 1e-3)

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
  mu[i] <- beta_0 + beta_1*T[i] + beta_2*R[i] + beta_3*dummy[i]
  y[i] ~ dnorm(mu[i], tau)
}
}"

inits <- function(){
  list(
    beta_0 = coef(lmod)[1],
    beta_1 = coef(lmod)[2],
    beta_2 = coef(lmod)[3],
    beta_3 = coef(lmod)[4]
  )
}


params0 <- c('beta_0', 'beta_1', 'beta_2', 'beta_3','tau')

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

head(ref_post_sample)

apply(ref_post_sample, 2, ess)

acf(ref_post_sample[,1:3])

apply(ref_post_sample, 2, summary)

par(mfrow = c(2,2))

hist_ci(ref_post_sample[,1], name = bquote(beta[0]))
hist_ci(ref_post_sample[,2], name = bquote(beta[1]))
hist_ci(ref_post_sample[,3], name = bquote(beta[2]))
hist_ci(ref_post_sample[,4], name = bquote(beta[3]))
par(mfrow = c(1,1))
hist_ci(1/ref_post_sample[,5], name = bquote(sigma^2))
