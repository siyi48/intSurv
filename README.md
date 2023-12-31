
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package description

<!-- badges: start -->
<!-- badges: end -->

The goal of intSurv is to estimate the heterogeneous treatment effect on
the survival outcomes, combining the data from the randomized trial and
observational studies. Under the proportional hazards assumption, the
heterogeneous treatment effect function $\tau(X)$ is defined as the
conditional hazard ratio of the potential survival times.

## Installation

You can install the development version of intSurv from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("siyi48/intSurv")
```

## Models

Let $X\in \mathbb{R}^d$ be the covariates, $A$ be the binary treatment,
$S$ be the data source indicator, and $T$ be the survival time. Also,
denote $\lambda_s(t\mid X, A=a)$ as the hazard function for the actual
survival time in the data source $s$, where $s=1$ indicates the trial
and $s=0$ indicates the observational studies. We posit the Cox models
on each data source as
$$\text{RCT:} \quad \lambda_1(t\mid X, A=a) =\lambda_{1,0}(t)\exp\{c_1(X) + a\tau(X)\},$$
$$\text{OS:} \quad \lambda_0(t\mid X, A=a) =\lambda_{0,0}(t)\exp\{c_0(X) + a\tau(X) + a\lambda(X)\}.$$

Here, $\tau(X)$ is the treatment effect function of interest, which is
further assumed to be a linear structure as
$$\tau(X) = \tau(X;\beta_0) = X^{\text T}\beta_0.$$

We conduct sieve approximation on the nuisance functions
$c_1(X), c_0(X)$, and $\lambda(X)$ using the polynomial basis. The
integrative estimator is solved by minimizing the penalized log partial
likelihood, with the penalty function selected as the adaptive lasso. We
only add the penalization for the covariate information of $\lambda(X)$
(the confounding function), since other nuisance functions do not affect
the efficiency of the integrative estimator.

## Example

The main function `surv.hte` provides two types of estimated parameter
$\beta$ in the treatment effect function $\tau(X)$ as `int` (the
integrative estimator) and `rct` (the trial-only estimator).

In the example data, the true $\tau(X) = 1 + X_1 + X_1^2$. We include
all possible covariates in the nuisance functions
$c_1(X), c_0(X),\lambda(X)$ and set the number of polynomial basis
functions for each nuisance function as $3$ (specified by `q.c1`,
`q.c0`, `q.lambda`) for the continuous covariates $X_1$. For the binary
covariate $X_2$, we input it into `c1.discrete`, `c0.discrete`, and
`lambda.discrete`. For identifiability, we do not include the intercept
term in $c_1(X)$ and $c_0(X)$, but we incorporate the intercept term in
the `lambda.discrete` of $\lambda(X)$.

``` r
library(intSurv)
## basic example code
q.c1 <- 3
q.c0 <- 3
q.lambda <- 3
nfolds <- 5
x1.std <- (dat$x1-mean(dat$x1))/sd(dat$x1)
c1.mat <- matrix(x1.std, ncol = 1)
c1.discrete <- matrix(dat$x2, ncol = 1)
c0.mat <- matrix(x1.std, ncol = 1)
c0.discrete <- matrix(dat$x2, ncol = 1)
tau.mat <- cbind(1, x1.std, x1.std^2)
lambda.mat <-  matrix(x1.std, ncol = 1)
lambda.discrete <- cbind(1, dat$x2)
s <- dat$s
delta <- dat$delta
t <- dat$t
a <- dat$a

res.hte <- surv.hte(tau.mat = tau.mat, c1.mat = c1.mat,
                    c1.discrete = c1.discrete, q.c1 = q.c1,
                    c0.mat = c0.mat, c0.discrete = c0.discrete, q.c0 = q.c0,
                    lambda.mat = lambda.mat,
                    lambda.discrete = lambda.discrete, q.lambda = q.lambda,
                    basis.type = "poly.raw",
                    a = a, s = s, t = t, delta = delta,
                    penalty.type = "adaptive.lasso",
                    nfolds = nfolds, type = c("int", "rct"))
res.hte
#> $beta.est
#> $beta.est$int
#>        ax1        ax2        ax3 
#>  0.4883803 -0.5901431 -0.5324468 
#> 
#> $beta.est$rct
#>        ax1        ax2        ax3 
#>  0.5042487 -0.6106382 -0.5455049 
#> 
#> 
#> $ve.beta
#> $ve.beta$int
#>        ax1        ax2        ax3 
#> 0.02371947 0.02210254 0.01392497 
#> 
#> $ve.beta$rct
#>        ax1        ax2        ax3 
#> 0.02397048 0.02244453 0.01411255 
#> 
#> 
#> $cov.beta
#> $cov.beta$int
#>              [,1]          [,2]          [,3]
#> [1,]  0.023719467 -0.0011182567 -0.0100443157
#> [2,] -0.001118257  0.0221025408 -0.0005352139
#> [3,] -0.010044316 -0.0005352139  0.0139249679
#> 
#> $cov.beta$rct
#>              [,1]          [,2]          [,3]
#> [1,]  0.023970480 -0.0013783538 -0.0102527470
#> [2,] -0.001378354  0.0224445329 -0.0003496444
#> [3,] -0.010252747 -0.0003496444  0.0141125490
#> 
#> 
#> $ate.est
#>         int         rct 
#> -0.07159466 -0.01533104 
#> 
#> $ve.ate
#>         int         rct 
#> 0.008316788 0.005923387 
#> 
#> $psi.est
#>   s0ax11   s0ax12   s0ax13    s0ax2    s0ax3 
#> 1.993929 1.964069 0.000000 1.940985 0.000000
```
