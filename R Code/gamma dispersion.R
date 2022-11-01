##########################################
# Overview
#
# How is the dispersion parameter in the gamma
# family estimated? Method of Moments (MOM)?
# Maximum Likelihood Estimation? (MLE)?
#
# The generalized linear model is based on
# Maximum Likelihood Estimation, but sometimes
# MLEs don't have closed forms.
#
# When nuisance parameters don't have closed
# forms and don't affect the final estimates
# of betas, MOMs are used instead of MLEs
# for computational speed and ease of implementation.
# This does not change the weights of the generalized
# linear model, but does change things like the
# likelihood and therefore the likelihood ratio
# tests.
#
# In R, this detail of MOM estimator vs MLE estimator is
# handled under the hood. This is great for just getting
# the job done, but makes recreating the results from
# scratch hard.
#
# Further, many stat courses either completely ignore this
# just mention it in passing, or prove using the MOM
# leads to the same betas as using the MLE.
#
# Questions:
# Can I calculate the MOM estimator of dispersion?
# Can I calculate the MLE estimator of dispersion?
# Can I find a publication that provides data and
# computes the two?
#
# Thankfully, the authors in In All Likelihood give
# a great amount of detail around MOM and MLE for
# the gamma family. This R script follows their examples.
##########################################


library(tidyverse)
library(GlmSimulatoR)

########################
# Make data
########################
dat <- tibble(
  Plot = 1:9,
  Density = c(5, 10, 15, 20, 30, 40, 60, 80, 100),
  Yield = c(122.7, 63, 32.5, 34.5, 31.4, 17.7, 21.9, 21.3, 18.4)
)

# Page 161
dat <- dat %>%
  mutate(
    logDensity = log(Density),
    logDensity = logDensity - mean(logDensity)
  )

# Page 166
model1 <- glm(Yield ~ logDensity, data = dat, family = Gamma(link = "inverse"))
summary(object = model1)$coefficients

########################
# MOM dispersion
########################
# Page 166
summary(model1)$dispersion # Provides estimate of dispersion

########################
# MLE dispersion
########################
MASS::gamma.dispersion(model1)

########################
# Compare weights
########################

# MOM
summary(object = model1)$coefficients

# MLE,
# weights the same, but standard errors change
summary(object = model1, dispersion = MASS::gamma.dispersion(model1))$coefficients

########################
# Compare likelihoods
########################
Y <- matrix(dat$Yield, ncol = 1)
Y_Hat <- matrix(predict(object = model1, newdata = dat, type = "link"), ncol = 1)

neg_ll_gamma_mom <- function(Y, Y_Hat) {
  Y_Hat <- Y_Hat^-1 # Inverse link

  # Method of moments estimates
  n <- nrow(Y)
  p <- 2
  numerator <- (Y - Y_Hat)^2
  denominator <- Y_Hat^2 * (n - p)
  phi <- sum(numerator / denominator)


  ll <- -log(Y) - log(gamma(1 / phi)) + (1 / phi) * (log(Y) - log(Y_Hat) - log(phi)) - Y / (Y_Hat * phi)
  ll <- sum(ll)
  ll <- -1 * ll

  return(ll)
}

-1 * neg_ll_gamma_mom(Y, Y_Hat)
as.numeric(logLik(model1))

# Above uses method of moments to estimate phi
# Could you MLE, but hard to find source that does and
# there is no closed form. Requires numerical routine.
# Multiple software implementations
# use method of moments. Stata and R are an examples.
#
# Results of digging through much of R's docs
# logLik depends on family's AIC. Take that AIC and some algebra
# to get logLik. It turns out the gamma family object uses MLE, but
# but the IRLS inside of glm uses the method of moments.
#
# This dispersion parameter from summary() is different than
# the dispersion parameter in aic function.
#
# Thus the -1 * neg_ll_gamma(Y, Y_Hat) == as.numeric(logLik(model1))
# fails because phi is estimated differently.

# copy and past gamma family's log likelihood logic, but with
# method of moments phi, not MLE. Do the functions match?


##############
# copy and paste of my negative log likelihood
# but with MLE phi
# Does this match logLik?
##############
neg_ll_gamma_mle <- function(Y, Y_Hat, phi) {
  Y_Hat <- Y_Hat^-1 # Inverse link



  ll <- -log(Y) - log(gamma(1 / phi)) + (1 / phi) * (log(Y) - log(Y_Hat) - log(phi)) - Y / (Y_Hat * phi)
  ll <- sum(ll)
  ll <- -1 * ll

  return(ll)
}

# yes. Probably floating point or underflow handling differences
-1 * neg_ll_gamma_mle(Y, Y_Hat, MASS::gamma.dispersion(model1))
as.numeric(logLik(model1))
