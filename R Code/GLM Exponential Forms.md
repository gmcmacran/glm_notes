
## Overview

Finding resources on the mathematical aspects of the generalized linear
model is challenging for a few reasons.

- Many resources talk about some families, but not others. Inverse
  Gaussian and Negative binomial are usually left out.
- Only some formulas are stated while others are not. Common example is
  how to calculate
  ![\phi](https://latex.codecogs.com/svg.latex?%5Cphi "\phi") in terms
  of the pdf’s or pmf’s parameters.
- Some formulas drop constant’s because they do not affect optimization
  in the IRWLS algorithm.
- The starting pdf is assumed known and not stated. For the gamma
  distribution, the starting pdf has a unique parameterization.

The goal of this document is to address these problems by clearly
stating **all** formulas for **all** families in a consistent order. No
proofs. Minimal words. All the necessary pieces.

## Exponential Form and Formula Overview

When deriving generalized linear models, the probability density or
probability mass function is algebraically manipulated into the
exponential form. The idea is to start with a formula for the pdf of a
distribution. For example,

![f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} exp(-\frac{(y - \mu)^2}{2\sigma^2})](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Csigma%5E2%29%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%5Csigma%5E2%7D%7D%20exp%28-%5Cfrac%7B%28y%20-%20%5Cmu%29%5E2%7D%7B2%5Csigma%5E2%7D%29 "f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} exp(-\frac{(y - \mu)^2}{2\sigma^2})")

and change it to an equivalent exponential form

![f(y; \mu, \sigma^2) = exp(\frac{y\mu-\mu^2/2}{\sigma^2} - \frac{y^2}{2\sigma^2} - \frac{1}{2}ln(2\pi\sigma^2))](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Csigma%5E2%29%20%3D%20exp%28%5Cfrac%7By%5Cmu-%5Cmu%5E2%2F2%7D%7B%5Csigma%5E2%7D%20-%20%5Cfrac%7By%5E2%7D%7B2%5Csigma%5E2%7D%20-%20%5Cfrac%7B1%7D%7B2%7Dln%282%5Cpi%5Csigma%5E2%29%29 "f(y; \mu, \sigma^2) = exp(\frac{y\mu-\mu^2/2}{\sigma^2} - \frac{y^2}{2\sigma^2} - \frac{1}{2}ln(2\pi\sigma^2))")

Then define functions a, b, and c so the specific exponential form is a
special case of the general exponential form

![f(y; \theta, \phi) = exp(\frac{y\theta - b(\theta)}{a(\phi)} + c(y,\phi))](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Ctheta%2C%20%5Cphi%29%20%3D%20exp%28%5Cfrac%7By%5Ctheta%20-%20b%28%5Ctheta%29%7D%7Ba%28%5Cphi%29%7D%20%2B%20c%28y%2C%5Cphi%29%29 "f(y; \theta, \phi) = exp(\frac{y\theta - b(\theta)}{a(\phi)} + c(y,\phi))")

The general exponential form is the form optimized by the IRWLS
algorithm. Note f started as a function of
![\mu](https://latex.codecogs.com/svg.latex?%5Cmu "\mu") and
![\sigma^2](https://latex.codecogs.com/svg.latex?%5Csigma%5E2 "\sigma^2")
and ended as a function of
![\theta](https://latex.codecogs.com/svg.latex?%5Ctheta "\theta") and
![\phi](https://latex.codecogs.com/svg.latex?%5Cphi "\phi"). In addition
to specifying a, b, and c, functions to calculate
![\theta](https://latex.codecogs.com/svg.latex?%5Ctheta "\theta") and
![\phi](https://latex.codecogs.com/svg.latex?%5Cphi "\phi") are
required.

## The Gaussian Family

The usual pdf is

![f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} exp(-\frac{(y - \mu)^2}{2\sigma^2})](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Csigma%5E2%29%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%5Csigma%5E2%7D%7D%20exp%28-%5Cfrac%7B%28y%20-%20%5Cmu%29%5E2%7D%7B2%5Csigma%5E2%7D%29 "f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} exp(-\frac{(y - \mu)^2}{2\sigma^2})")

In exponential form, the pdf is

![f(y; \mu, \sigma^2) = exp(\frac{y\mu-\mu^2/2}{\sigma^2} - \frac{y^2}{2\sigma^2} - \frac{1}{2}ln(2\pi\sigma^2))](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Csigma%5E2%29%20%3D%20exp%28%5Cfrac%7By%5Cmu-%5Cmu%5E2%2F2%7D%7B%5Csigma%5E2%7D%20-%20%5Cfrac%7By%5E2%7D%7B2%5Csigma%5E2%7D%20-%20%5Cfrac%7B1%7D%7B2%7Dln%282%5Cpi%5Csigma%5E2%29%29 "f(y; \mu, \sigma^2) = exp(\frac{y\mu-\mu^2/2}{\sigma^2} - \frac{y^2}{2\sigma^2} - \frac{1}{2}ln(2\pi\sigma^2))")

By inspection:

![\theta = \mu](https://latex.codecogs.com/svg.latex?%5Ctheta%20%3D%20%5Cmu "\theta = \mu")

![\phi = \sigma^2](https://latex.codecogs.com/svg.latex?%5Cphi%20%3D%20%5Csigma%5E2 "\phi = \sigma^2")

![a(\phi) = \phi](https://latex.codecogs.com/svg.latex?a%28%5Cphi%29%20%3D%20%5Cphi "a(\phi) = \phi")

![b(\theta) = \theta^2/2](https://latex.codecogs.com/svg.latex?b%28%5Ctheta%29%20%3D%20%5Ctheta%5E2%2F2 "b(\theta) = \theta^2/2")

![c(y, \phi) = - \frac{y^2}{2\phi} -\frac{1}{2}ln(2\pi\phi)](https://latex.codecogs.com/svg.latex?c%28y%2C%20%5Cphi%29%20%3D%20-%20%5Cfrac%7By%5E2%7D%7B2%5Cphi%7D%20-%5Cfrac%7B1%7D%7B2%7Dln%282%5Cpi%5Cphi%29 "c(y, \phi) = - \frac{y^2}{2\phi} -\frac{1}{2}ln(2\pi\phi)")

## The Gamma Family

This distribution is usually parameterized by shape and rate. For the
generalized linear models framework, this pdf is parameterized by
![\mu](https://latex.codecogs.com/svg.latex?%5Cmu "\mu") and
![\phi](https://latex.codecogs.com/svg.latex?%5Cphi "\phi")

![f(y; \mu, \phi) = \frac{1}{y\Gamma(\frac{1}{\phi})} (\frac{y}{\mu\phi})^\frac{1}{\phi} exp(-\frac{y}{\mu\phi})](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Cphi%29%20%3D%20%5Cfrac%7B1%7D%7By%5CGamma%28%5Cfrac%7B1%7D%7B%5Cphi%7D%29%7D%20%28%5Cfrac%7By%7D%7B%5Cmu%5Cphi%7D%29%5E%5Cfrac%7B1%7D%7B%5Cphi%7D%20exp%28-%5Cfrac%7By%7D%7B%5Cmu%5Cphi%7D%29 "f(y; \mu, \phi) = \frac{1}{y\Gamma(\frac{1}{\phi})} (\frac{y}{\mu\phi})^\frac{1}{\phi} exp(-\frac{y}{\mu\phi})")

In exponential form, the pdf is

![f(y; \mu, \phi) = exp(\frac{\frac{y}{\mu} - (-ln(\mu))}{-\phi} + \frac{1-\phi}{\phi}ln(y) - \frac{ln(\phi)}{\phi} - ln{\Gamma(\frac{1}{\phi})})](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Cphi%29%20%3D%20exp%28%5Cfrac%7B%5Cfrac%7By%7D%7B%5Cmu%7D%20-%20%28-ln%28%5Cmu%29%29%7D%7B-%5Cphi%7D%20%2B%20%5Cfrac%7B1-%5Cphi%7D%7B%5Cphi%7Dln%28y%29%20-%20%5Cfrac%7Bln%28%5Cphi%29%7D%7B%5Cphi%7D%20-%20ln%7B%5CGamma%28%5Cfrac%7B1%7D%7B%5Cphi%7D%29%7D%29 "f(y; \mu, \phi) = exp(\frac{\frac{y}{\mu} - (-ln(\mu))}{-\phi} + \frac{1-\phi}{\phi}ln(y) - \frac{ln(\phi)}{\phi} - ln{\Gamma(\frac{1}{\phi})})")

By inspection:

![\theta = \frac{1}{\mu}](https://latex.codecogs.com/svg.latex?%5Ctheta%20%3D%20%5Cfrac%7B1%7D%7B%5Cmu%7D "\theta = \frac{1}{\mu}")

![\phi = \phi](https://latex.codecogs.com/svg.latex?%5Cphi%20%3D%20%5Cphi "\phi = \phi")

![a(\phi) = -\phi](https://latex.codecogs.com/svg.latex?a%28%5Cphi%29%20%3D%20-%5Cphi "a(\phi) = -\phi")

![b(\theta) = -ln(\frac{1}{\theta})](https://latex.codecogs.com/svg.latex?b%28%5Ctheta%29%20%3D%20-ln%28%5Cfrac%7B1%7D%7B%5Ctheta%7D%29 "b(\theta) = -ln(\frac{1}{\theta})")

![c(y, \phi) = \frac{1-\phi}{\phi}ln(y) - \frac{ln(\phi)}{\phi} - ln{\Gamma(\frac{1}{\phi})}](https://latex.codecogs.com/svg.latex?c%28y%2C%20%5Cphi%29%20%3D%20%5Cfrac%7B1-%5Cphi%7D%7B%5Cphi%7Dln%28y%29%20-%20%5Cfrac%7Bln%28%5Cphi%29%7D%7B%5Cphi%7D%20-%20ln%7B%5CGamma%28%5Cfrac%7B1%7D%7B%5Cphi%7D%29%7D "c(y, \phi) = \frac{1-\phi}{\phi}ln(y) - \frac{ln(\phi)}{\phi} - ln{\Gamma(\frac{1}{\phi})}")

## The Inverse Gaussian Family

The usual pdf is

![f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi y^3 \sigma^2}} exp(-\frac{(y - \mu)^2}{2(\mu\sigma)^2y})](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Csigma%5E2%29%20%3D%20%5Cfrac%7B1%7D%7B%5Csqrt%7B2%5Cpi%20y%5E3%20%5Csigma%5E2%7D%7D%20exp%28-%5Cfrac%7B%28y%20-%20%5Cmu%29%5E2%7D%7B2%28%5Cmu%5Csigma%29%5E2y%7D%29 "f(y; \mu, \sigma^2) = \frac{1}{\sqrt{2\pi y^3 \sigma^2}} exp(-\frac{(y - \mu)^2}{2(\mu\sigma)^2y})")

In exponential form, the pdf is

![f(y; \mu, \sigma^2) = exp(\frac{\frac{y}{2\mu^2} - \frac{1}{\mu}}{-\sigma^2} + \frac{1}{-2y\sigma^2} - \frac{1}{2}ln(2\pi y^3 \sigma^2))](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%2C%20%5Csigma%5E2%29%20%3D%20exp%28%5Cfrac%7B%5Cfrac%7By%7D%7B2%5Cmu%5E2%7D%20-%20%5Cfrac%7B1%7D%7B%5Cmu%7D%7D%7B-%5Csigma%5E2%7D%20%2B%20%5Cfrac%7B1%7D%7B-2y%5Csigma%5E2%7D%20-%20%5Cfrac%7B1%7D%7B2%7Dln%282%5Cpi%20y%5E3%20%5Csigma%5E2%29%29 "f(y; \mu, \sigma^2) = exp(\frac{\frac{y}{2\mu^2} - \frac{1}{\mu}}{-\sigma^2} + \frac{1}{-2y\sigma^2} - \frac{1}{2}ln(2\pi y^3 \sigma^2))")

By inspection:

![\theta = \frac{1}{2\mu^2}](https://latex.codecogs.com/svg.latex?%5Ctheta%20%3D%20%5Cfrac%7B1%7D%7B2%5Cmu%5E2%7D "\theta = \frac{1}{2\mu^2}")

![\phi = \sigma^2](https://latex.codecogs.com/svg.latex?%5Cphi%20%3D%20%5Csigma%5E2 "\phi = \sigma^2")

![a(\phi) = -\phi](https://latex.codecogs.com/svg.latex?a%28%5Cphi%29%20%3D%20-%5Cphi "a(\phi) = -\phi")

![b(\theta) = \sqrt{2\*\theta}](https://latex.codecogs.com/svg.latex?b%28%5Ctheta%29%20%3D%20%5Csqrt%7B2%2A%5Ctheta%7D "b(\theta) = \sqrt{2*\theta}")

![c(y, \phi) = \frac{1}{-2y\phi} - \frac{1}{2}ln(2\pi y^3 \phi))](https://latex.codecogs.com/svg.latex?c%28y%2C%20%5Cphi%29%20%3D%20%5Cfrac%7B1%7D%7B-2y%5Cphi%7D%20-%20%5Cfrac%7B1%7D%7B2%7Dln%282%5Cpi%20y%5E3%20%5Cphi%29%29 "c(y, \phi) = \frac{1}{-2y\phi} - \frac{1}{2}ln(2\pi y^3 \phi))")

## The Binomial Family

The usual pdf is

![f(y; n, p) = \binom{n}{y} p^y (1-p)^{n-y}](https://latex.codecogs.com/svg.latex?f%28y%3B%20n%2C%20p%29%20%3D%20%5Cbinom%7Bn%7D%7By%7D%20p%5Ey%20%281-p%29%5E%7Bn-y%7D "f(y; n, p) = \binom{n}{y} p^y (1-p)^{n-y}")

In exponential form, the pdf is

![f(y; n, p) = exp(y\*ln(\frac{p}{1-p}) + n\*ln(1-p) + ln\binom{n}{y})](https://latex.codecogs.com/svg.latex?f%28y%3B%20n%2C%20p%29%20%3D%20exp%28y%2Aln%28%5Cfrac%7Bp%7D%7B1-p%7D%29%20%2B%20n%2Aln%281-p%29%20%2B%20ln%5Cbinom%7Bn%7D%7By%7D%29 "f(y; n, p) = exp(y*ln(\frac{p}{1-p}) + n*ln(1-p) + ln\binom{n}{y})")

By inspection:

![\theta = ln(\frac{p}{1 - p})](https://latex.codecogs.com/svg.latex?%5Ctheta%20%3D%20ln%28%5Cfrac%7Bp%7D%7B1%20-%20p%7D%29 "\theta = ln(\frac{p}{1 - p})")

![\phi = 1](https://latex.codecogs.com/svg.latex?%5Cphi%20%3D%201 "\phi = 1")

![a(\phi) = \phi](https://latex.codecogs.com/svg.latex?a%28%5Cphi%29%20%3D%20%5Cphi "a(\phi) = \phi")

![b(\theta) = n \* ln(1 + exp(\theta))](https://latex.codecogs.com/svg.latex?b%28%5Ctheta%29%20%3D%20n%20%2A%20ln%281%20%2B%20exp%28%5Ctheta%29%29 "b(\theta) = n * ln(1 + exp(\theta))")

![c(y, \phi) =ln\binom{n}{y}](https://latex.codecogs.com/svg.latex?c%28y%2C%20%5Cphi%29%20%3Dln%5Cbinom%7Bn%7D%7By%7D "c(y, \phi) =ln\binom{n}{y}")

## The Poisson Family

The usual pdf is

![f(y; \mu) = \frac{exp(-\mu) \* \mu^y}{y!}](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%29%20%3D%20%5Cfrac%7Bexp%28-%5Cmu%29%20%2A%20%5Cmu%5Ey%7D%7By%21%7D "f(y; \mu) = \frac{exp(-\mu) * \mu^y}{y!}")

In exponential form, the pdf is

![f(y; \mu) = exp(y\*ln(\mu) - \mu - ln\Gamma(y+1))](https://latex.codecogs.com/svg.latex?f%28y%3B%20%5Cmu%29%20%3D%20exp%28y%2Aln%28%5Cmu%29%20-%20%5Cmu%20-%20ln%5CGamma%28y%2B1%29%29 "f(y; \mu) = exp(y*ln(\mu) - \mu - ln\Gamma(y+1))")

By inspection:

![\theta = ln(\mu)](https://latex.codecogs.com/svg.latex?%5Ctheta%20%3D%20ln%28%5Cmu%29 "\theta = ln(\mu)")

![\phi = 1](https://latex.codecogs.com/svg.latex?%5Cphi%20%3D%201 "\phi = 1")

![a(\phi) = \phi](https://latex.codecogs.com/svg.latex?a%28%5Cphi%29%20%3D%20%5Cphi "a(\phi) = \phi")

![b(\theta) = exp(\theta)](https://latex.codecogs.com/svg.latex?b%28%5Ctheta%29%20%3D%20exp%28%5Ctheta%29 "b(\theta) = exp(\theta)")

![c(y, \phi) = -ln\Gamma(y+1))](https://latex.codecogs.com/svg.latex?c%28y%2C%20%5Cphi%29%20%3D%20-ln%5CGamma%28y%2B1%29%29 "c(y, \phi) = -ln\Gamma(y+1))")

## The Negative Binomial Family

The usual pdf is

![f(y; r, p) = \binom{y + r - 1}{r-1} p^r (1-p)^{y}](https://latex.codecogs.com/svg.latex?f%28y%3B%20r%2C%20p%29%20%3D%20%5Cbinom%7By%20%2B%20r%20-%201%7D%7Br-1%7D%20p%5Er%20%281-p%29%5E%7By%7D "f(y; r, p) = \binom{y + r - 1}{r-1} p^r (1-p)^{y}")

In exponential form, the pdf is

![f(y; r, p) = exp(y\*ln(1-p) + r\*ln(p) + ln\binom{y + r - 1}{r-1})](https://latex.codecogs.com/svg.latex?f%28y%3B%20r%2C%20p%29%20%3D%20exp%28y%2Aln%281-p%29%20%2B%20r%2Aln%28p%29%20%2B%20ln%5Cbinom%7By%20%2B%20r%20-%201%7D%7Br-1%7D%29 "f(y; r, p) = exp(y*ln(1-p) + r*ln(p) + ln\binom{y + r - 1}{r-1})")

By inspection:

![\theta = ln(1-p)](https://latex.codecogs.com/svg.latex?%5Ctheta%20%3D%20ln%281-p%29 "\theta = ln(1-p)")

![\phi = 1](https://latex.codecogs.com/svg.latex?%5Cphi%20%3D%201 "\phi = 1")

![a(\phi) = \phi](https://latex.codecogs.com/svg.latex?a%28%5Cphi%29%20%3D%20%5Cphi "a(\phi) = \phi")

![b(\theta) = - r \* ln(1 - exp(\theta))](https://latex.codecogs.com/svg.latex?b%28%5Ctheta%29%20%3D%20-%20r%20%2A%20ln%281%20-%20exp%28%5Ctheta%29%29 "b(\theta) = - r * ln(1 - exp(\theta))")

![c(y, \phi) =ln\binom{y + r - 1}{r-1})](https://latex.codecogs.com/svg.latex?c%28y%2C%20%5Cphi%29%20%3Dln%5Cbinom%7By%20%2B%20r%20-%201%7D%7Br-1%7D%29 "c(y, \phi) =ln\binom{y + r - 1}{r-1})")

## Additional Resources

All formulas come out of Generalized Linear Models and Extensions by
Hardin and Hilber.

See <https://github.com/gmcmacran/altForm> for working functions based
on these equations.
