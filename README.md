# CressieReadTest

install packages

```
devtools::install_github("krzysiektr/CressieReadTest")
```

## Discrete distributions

```
set.seed(4101)
u <- rpois(100,2)
```

```
obs(u)
$dist
[1] "pois"

$ob
[1]  8 23 31 25  8  3  1  1

$ex
[1] 11.0803158 24.3766948 26.8143643 19.6638672 10.8151269  4.7586559  1.7448405  0.7461346

chi2.test(u)

	Cressie-Read Goodness-of-Fit Test - Poisson distributions

data:  u
CR = 4.8225, df = 6, lambda = 1, p-value = 0.5668
sample estimates:
       mu    logLik 
   2.2000 -167.7122 
```

```
obs(u, dist = "binom")
$dist
[1] "binom"

$ob
[1]  8 23 31 25  8  3  1  1

$ex
[1]  9.5940663 25.0644983 29.7640917 21.2069153 10.0732848  3.3493672  0.7954747  0.1523017

chi2.test(u, dist = "binom", lambda = 1)

	Cressie-Read Goodness-of-Fit Test - binomial distributions

data:  u
CR = 6.3986, df = 6, lambda = 1, p-value = 0.38
sample estimates:
        size         prob       loglik 
  11.0000000    0.1919192 -166.9101709 
```

## Continuous distributions

```
set.seed(4101)
w <- rnorm(100,2,1)
```

```
chi2.gof(w, dist="norm", lambda=1)

	Cressie-Read Goodness-of-Fit Test - normal distributions

data:  w
CR = 8.42, df = 10, lambda = 1, p-value = 0.5879
sample estimates:
          mu          std       loglik 
   2.0236154    0.9506787 -136.3359381
```

```
nortest::pearson.test(w)

	Pearson chi-square normality test

data:  w
P = 8.42, p-value = 0.5879
```
