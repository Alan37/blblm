# blblm

<!-- badges: start -->
<!-- badges: end -->

## Examples

blblm performs bag of little bootstrap (blb) on lm function:

``` r
library(blblm)
fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
coef(fit)
#> (Intercept)          wt          hp       wt:hp 
#> 48.88428523 -7.88702986 -0.11576659  0.02600976
confint(fit, c("wt", "hp"))
#>           2.5%       97.5%
#> wt -10.7902240 -5.61586271
#> hp  -0.1960903 -0.07049867
sigma(fit)
#> [1] 1.838911
sigma(fit, confidence = TRUE)
#>    sigma      lwr      upr 
#> 1.838911 1.350269 2.276347
predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
#>        1        2 
#> 21.55538 18.80785
predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)
#>        fit      lwr      upr
#> 1 21.55538 20.02457 22.48764
#> 2 18.80785 17.50654 19.71772
```

blblm performs bag of little bootstrap (blb) on glm function:

``` r
library(blblogreg)
fit_blblogreg = blblogreg(Direction ~ Lag1 + Lag2 + Volume, data = ISLR::Smarket, m = 3, B = 100)
coef(fit_blblogreg)
#> (Intercept)        Lag1        Lag2      Volume 
#> -0.19414049 -0.08208848 -0.04627690  0.18828639
confint(fit_blblogreg)
#>             2.5%      97.5%
#>Lag1   -0.2474233 0.07744407
#>Lag2   -0.2306516 0.12280745
#>Volume -0.3550769 0.72407196
sigma(fit_blblogreg)
#> [1] 1.475419
predict(fit_blblogreg, new_data = data.frame("Lag1" = c(0.382, 0.4, 0.5),"Lag2" = c(1.03, 1.4, 0.5), "Volume" = c(1.2, 1.34, 1.1))) 
#>          1           2           3 
#>-0.04721983 -0.03945978 -0.05120815 
```