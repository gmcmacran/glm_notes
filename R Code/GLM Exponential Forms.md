
## Overview

When deriving generalized linear models, the probability density or
probability mass function is algebraically manipulated into the
exponential form. This document lists out all the necessary pieces for
this form.

All formulas come out of Generalized Linear Models and Extensions by
Hardin and Hilber. See <https://github.com/gmcmacran/altForm> for
working functions based on these equations.

## Exponential Form

There are other forms, but the exponential form is usually written as

$$ f(y; \theta, \phi) = exp(\frac{y\theta - b(\theta)}{a(\phi)} + c(y,\phi))$$
where a, b, c are functions, $\theta$ is a function of $\mu$, and $\phi$
is a function of the nuisance parameter.

## The Gaussian Family

The usual pdf is

$$ f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} exp(-\frac{(y - \mu)^2}{2\sigma^2}) $$

In exponential form, the pdf is

$$ f(y; \mu, \sigma^2) = exp(\frac{y\mu-\mu^2/2}{\sigma^2} - \frac{y^2}{2\sigma^2} - \frac{1}{2}ln(2\pi\sigma^2)) $$

By inspection:

$$ \theta = \mu $$ $$ \phi = \sigma^2 $$ $$ a(\phi) = \phi $$
$$ b(\theta) = \theta^2/2 $$
$$ c(y, \phi) = - \frac{y^2}{2\phi} -\frac{1}{2}ln(2\pi\phi) $$

## The Gamma Family

This distribution is usually parameter by shape and rate. For the
generalized linear models framework, this pdf is parameterized by $\mu$
and $\phi$

$$ f(y; \mu, \phi) = \frac{1}{y\Gamma(\frac{1}{\phi})} (\frac{y}{\mu\phi})^\frac{1}{\phi} exp(-\frac{y}{\mu\phi})$$

In exponential form, the pdf is

$$ f(y; \mu, \phi) = exp(\frac{\frac{y}{\mu} - (-ln(\mu))}{-\phi} + \frac{1-\phi}{\phi}ln(y) - \frac{ln(\phi)}{\phi} - ln{\Gamma(\frac{1}{\phi})})   $$

By inspection:

$$ \theta = \frac{1}{\mu} $$ $$ \phi = \phi $$ $$ a(\phi) = -\phi $$
$$ b(\theta) = -ln(\frac{1}{\theta}) $$
$$ c(y, \phi) = \frac{1-\phi}{\phi}ln(y) - \frac{ln(\phi)}{\phi} - ln{\Gamma(\frac{1}{\phi})} $$

## The Inverse Gaussian Family

The usual pdf is

$$ f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi y^3 \sigma^2}} exp(-\frac{(y - \mu)^2}{2(\mu\sigma)^2y}) $$

In exponential form, the pdf is

$$ f(y; \mu, \sigma^2) = exp(\frac{\frac{y}{2\mu^2} - \frac{1}{\mu}}{-\sigma^2} + \frac{1}{-2y\sigma^2} - \frac{1}{2}ln(2\pi y^3 \sigma^2)) $$

By inspection:

$$ \theta = \frac{1}{2\mu^2} $$ $$ \phi = \sigma^2 $$
$$ a(\phi) = -\phi $$ $$ b(\theta) = \sqrt{2*\theta} $$
$$ c(y, \phi) = \frac{1}{-2y\phi} - \frac{1}{2}ln(2\pi y^3 \phi)) $$

## The Binomial Family

The usual pdf is

$$ f(y; n, p) = \binom{n}{y} p^y (1-p)^{n-y} $$

In exponential form, the pdf is

$$ f(y; n, p) = exp(y*ln(\frac{p}{1-p}) + n*ln(1-p) + ln\binom{n}{y}) $$

By inspection:

$$ \theta = ln(\frac{p}{1 - p}) $$ $$ \phi = 1 $$ $$ a(\phi) = \phi $$
$$ b(\theta) = n * log(1 + exp(\theta)) $$
$$ c(y, \phi) =ln\binom{n}{y}$$

## The Poisson Family

The usual pdf is

$$ f(y; \mu) = \frac{exp(-\mu) * \mu^y}{y!}  $$

In exponential form, the pdf is

$$ f(y; \mu) = exp(y*ln(\mu) - \mu - ln\Gamma(y+1)) $$

By inspection:

$$ \theta = ln(\mu) $$ $$ \phi = 1 $$ $$ a(\phi) = \phi $$
$$ b(\theta) = exp(\theta) $$ $$ c(y, \phi) = -ln\Gamma(y+1))$$

## The Negative Binomial Family

The usual pdf is

$$ f(y; r, p) = \binom{y + r - 1}{r-1} p^r (1-p)^{y} $$

In exponential form, the pdf is

$$ f(y; r, p) = exp(y*ln(1-p) + r*ln(p) + ln\binom{y + r - 1}{r-1}) $$

By inspection:

$$ \theta = ln(1-p) $$ $$ \phi = 1 $$ $$ a(\phi) = \phi $$
$$ b(\theta) = - r * log(1 - exp(\theta)) $$
$$ c(y, \phi) =ln\binom{y + r - 1}{r-1})$$