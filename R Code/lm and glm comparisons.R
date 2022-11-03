#################################################
# One of the really nice things about glms is
# 99% of lm knowledge directly relates. But
# finding a book that makes these connections
# very clear is hard.
#
# Books and blog posts are easy to find for
# linear regression. Hard to find for glms.
# A lot of glm formulas can be written many
# equivalent ways.
#
# Further some people say deviance when they
# really mean scaled deviance. All this makes
# it very hard to link works with code.
#
# To understand the relation, this script builds
# a lm and glm. Then repeats analysis under
# both models and shows test stats and p-values
# are identical.
#################################################
library(tidyverse)
library(lmtest)

#############
# Make data
#############
set.seed(1)
N <- 50
X1 <- rnorm(N, 0, 1)
X2 <- rnorm(N, 0, 1)
X3 <- rnorm(N, 0, 1)

Y <- 1 + 1 * X1 + 2 * X2 + 3 * X3 + rnorm(N, 0, 10)
Y2 <- factor(1:length(Y))

train <- tibble(Y = Y, Y2 = Y2, X1 = X1, X2 = X2, X3 = X3)
rm(Y, Y2, X1, X2, X3, N)

#############
# Models
#############
model_lm <- lm(Y ~ X1 + X2 + X3, data = train)
model_glm <- glm(Y ~ X1 + X2 + X3, data = train, family = gaussian(link = "identity"))
model_glm_const <- glm(Y ~ 1, data = train, family = gaussian(link = "identity"))
model_glm_saturated <- glm(Y ~ Y2, data = train, family = gaussian(link = "identity"))

#############
# Compare fits
#############
all.equal(model_lm$coefficients, model_glm$coefficients)
all.equal(model_lm$fitted.values, model_glm$fitted.values)
all.equal(model_lm$residuals, model_glm$residuals)

#############
# Linear regressions' SSE is glm's deviance
#############
sum((model_lm$fitted.values - train$Y)^2)
model_glm$deviance

#############
# Linear regressions' SSTO is glm's null.deviance
#############

sum((train$Y - mean(train$Y))^2)
model_glm$null.deviance

#############
# Null deviance is constant model's deviance
#############
model_glm$null.deviance == model_glm_const$deviance
model_glm$df.null == model_glm_const$df.residual

#############
# Saturated model has no error
#############
sum((model_glm_saturated$fitted.values - train$Y)^2)
model_glm_saturated$deviance
model_glm_saturated$df.residual


#############
# Compare likelihood tests
#############
# H0: All weights (excet intercept) are zero.
# HA: Is at least one weight not zero.
lrtest(model_lm)
lrtest(model_glm)

# Calculate test myself
# A lot of different ways of writing the test stat
# General formulation. Usually what is seen in proof books. lr test for contingency tables
# -2 x ln(null dev / resid dev)
# Sometimes negative sign inside ln. (flips ratio)
# 2 * ln(resid dev / null dev)
# Usually what is seen with glm. Algebraically equivalent
# Can be interpreted as "Is the reduction in deviance caused by adding variables significant?"
# 2 x ln(resid dev) -2 x ln(null dev)
#  -2 x ln(null dev) - -2 x ln(resid dev)
#
G <- -2 * (logLik(model_glm_const) - logLik(model_glm))
G <- as.numeric(G)
G_df <- model_glm_const$df.residual - model_glm$df.residual
pchisq(G, G_df, lower.tail = FALSE)

#############
# Sequential analysis via analysis of variance and analysis of deviance
#############
anova(model_lm, test = "F")
anova(model_glm, test = "F")

#############
# Sequential analysis using different tests
#############
anova(model_lm, test = "F")
anova(model_lm, test = "Chisq")

anova(model_glm, test = "F")
anova(model_glm, test = "Chisq") # Changes

#############
# Goodnesss of fit
#############
# H0. Model fits data well.
# HA: The model does not fit the data well.
# Can be either mle or pearson based. For gaussian(link = "identity"), they are the same. Not true for other glms.

# See https://www.amelia.mn/sds291/labs/lab_logistic_regression.html
# for clear application of test to logistic regression.
#
# Often called "likelihood ratio test" to refer to testing model_glm to the perfect model.
# But in general, there are many likelihood ratio tests for glms. The lrtest(model_glm) tests
# model_glm is an improvment over model glm(Y ~ 1) and is also a likelihood test.

lrtest(model_glm_saturated, model_glm)

G <- -2 * (logLik(model_glm) - logLik(model_glm_saturated))
G <- as.numeric(G)
G_df <- model_glm$df.residual - model_glm_saturated$df.residual
pchisq(G, G_df, lower.tail = FALSE)

#############
# residuals and deviance
#############
all(residuals(object = model_glm, type = "deviance") == residuals(object = model_glm, type = "pearson")) # Only true for gaussian family

all(residuals(object = model_glm, type = "deviance") == model_glm$residuals)

sum(residuals(object = model_glm, type = "deviance")^2) == model_glm$deviance
