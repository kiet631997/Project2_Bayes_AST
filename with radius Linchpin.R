rm(list = ls())

library(mvtnorm)
library(modelsummary)
library(ggplot2)
library(plotly)
library(mcmcse)
library(SimTools)
library(MCMCpack)
source('functions.R')

data <- read.csv('SB_LAW.csv', header = T)

plot(logL1 ~ logT1, data = data, pch = 16)


plot(logL1 ~ logT1, data = data, pch = 16)
abline(b = 6, a = -21, col = 'red', lty = 2, lwd = 2)

data$dummy <- ifelse(data$logL1<6*data$logT1 - 21, 'MC', 'RG')
data$dummy <- as.factor(data$dummy)

lmod <- lm(logL1 ~ logT1 + logR1 + dummy, data = data)
b_mod <- coef(lmod)

plot(logL1 ~ logT1, data = data, pch = 16)

modelsummary(lmod)


minESS(p = 5)


sse<-sum(lmod$residuals^2) 

lam_shape<-lmod$df/2

mod_mat<-model.matrix(lmod)

XTXinv<-solve(t(mod_mat) %*% mod_mat)

msim<-1e4

ref_post_sample<-matrix(NA_real_, ncol=5, nrow=msim)

colnames(ref_post_sample) <- c(bquote(beta[0]),bquote(beta[1]),
                               bquote(beta[2]),bquote(beta[2]),
                               bquote(lambda))

for (iter in 1:msim){
  lam<-rgamma(1, shape=lam_shape, scale = 2/sse)
  cmat<-(1/lam)*XTXinv
  ref_post_sample[iter,]<-c(rmvnorm(1,
                                    mean = lmod$coefficients, 
                                    sigma = cmat),lam)
}

ref_post_sample <- as.mcmc(ref_post_sample)
apply(ref_post_sample, 2, ess)

apply(ref_post_sample, 2, summary)

acf(ref_post_sample[,1:3])

acf(ref_post_sample[,3:5])

par(mfrow = c(2,2))

hist_ci(ref_post_sample[,1], name = bquote(beta[0]))
hist_ci(ref_post_sample[,2], name = bquote(beta[1]))
hist_ci(ref_post_sample[,3], name = bquote(beta[2]))
hist_ci(ref_post_sample[,4], name = bquote(beta[3]))

par(mfrow = c(1,1))

plot_ly( x = data$logT1, y = data$logR1, z = data$logL1, color = data$dummy)

b0_post <- ref_post_sample[,1]
b1_post <- ref_post_sample[,2]
b2_post <- ref_post_sample[,3]
b3_post <- ref_post_sample[,4]

df <- data.frame(b0 = as.vector(b0_post),
                 b1 = as.vector(b1_post),
                 b2 = as.vector(b2_post),
                 b3 = as.vector(b3_post))

f02 <- ggplot(data = df, aes(x = b0, y=b2) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() + 
  xlab(bquote(beta[0])) + 
  ylab(bquote(beta[2])) + 
  theme(
    legend.position='none'
  )

f12 <- ggplot(data = df, aes(x = b1, y=b2) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlab(bquote(beta[1])) + 
  ylab(bquote(beta[2])) +
  theme(
    legend.position='none'
  )

f01 <- ggplot(data = df, aes(x = b0, y=b1) ) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlab(bquote(beta[0])) + 
  ylab(bquote(beta[1]))+
  theme(
    legend.position='none'
  )

ggpubr::ggarrange(f02, f12, f01, nrow = 1, ncol = 3)



