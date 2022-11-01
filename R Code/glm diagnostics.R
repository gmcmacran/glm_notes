#####################################
# There does not seem to be a standard
# way of assessing diagnostics for glms.
# One author recommends one thing
# The next author recommends something different
# Further one math formula may be called
# something different by different authors.
#
# Even R's documentation is unclear.
#
# This script follows Linear Regressions Analysis's
# diagnostics, but with both lm and glm.
#####################################
library(tidyverse)
library(MPV)

##############
# Get data
##############
data("softdrink")
head(softdrink)

# page 213
nrow(softdrink) == 25
ncol(softdrink) == 3

##############
# Make model
##############
model_lm <- lm(formula = y ~ x1 + x2, data = softdrink)

##############
# diagnostics
##############
diag_lm <- influence.measures(model = model_lm)

# Do results match texkbook? Use page 214
temp <- diag_lm$infmat %>% as_tibble()

temp %>% select(cook.d)
temp %>% select(hat)
temp %>% select(dffit)

##############
# compare equality between lm and glm
##############
model_glm <- glm(formula = y ~ x1 + x2, data = softdrink)
diag_glm <- influence.measures(model = model_glm)

all.equal(diag_glm$infmat, diag_lm$infmat)
all.equal(diag_glm$is.inf, diag_lm$is.inf)
