
# Generalized Linear Model Notes

This repo is a collection of notes for the generalized models. Some
notes are mathematical. Some notes are in code.

## gamma dispersion.R

Both R and Stata’s glm implementation use the method of moments
estimator of dispersion instead of the maximum likelihood estimator.
This R script explores the numerical differences and some of the
statistical consequences of different estimators.

## lm and glm comparisons.R

Many statistical properties of the linear model are applicable to the
generalized linear model. The linear models’ SSE is the generalized
linear model’s deviance. Linear model’s SSTO is the generalized linear
model’s null deviance. etc. Unfortunately, many connections are
obfuscated by different names.

This R script explores connections by performing many calculations with
two models. One trained with lm() and the other glm(family =
gaussian()).

## glm diagnostics.R

The statistical community has not reached a consensus on how to analyze
a generalized linear model’s fit.

- Some textbooks recommend one approach. Other textbooks recommend a
  different approach.
- Some approaches are implemented in R and some are not.
- Some assessment techniques are applicable to the binomial family, but
  not the gamma family and vice versa.
- Some textbooks make a distinction between deviance and scaled
  deviance. Other textbooks don’t.
- Sometimes deviance and scaled deviance are the same so the distinction
  is irrelevant (Poisson and binomial).

These details make it nearly impossible to asses a model’s fit. This
script takes a different strategy. It uses R’s influence.measures on a
linear model. The calculations are checked against examples in Linear
Regressions Analysis by Montgomery, Peck and Vining.

Then a glm model with gaussian family is trained. Equality between R’s
influence.measure’s calculations on the two models is confirmed. This
function can be used across all glm families so I have a simple,
consistent, and easy way of assessing fit across all types of
generalized linear models.

## GLM Exponential Forms.qmd

In deriving the generalized linear model, the fist step is algebraically
manipulating a pdf or pmf to its exponential form. Then defining a set
of functions. Finding resources that clearly state these functions for
all distributions is difficult. Inverse Gaussian is almost always
completely skipped.

This document states all the pieces and is largely based on Generalized
Linear Models and Extensions by Hardin and Hilbe.
