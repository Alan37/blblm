#' @import purrr
#' @import furrr
#' @import future
#' @import stats
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' Bag of Little Bootstraps Version of Linear Regression
#'
#' blblm performs bag of little bootstrap (blb) on lm function. If users set "parallel = TRUE" in the blblm() function, the plan() function is usually needed for users to specify in front of the blblm() function.
#'
#'
#' @param formula an object of class "formula" (or one that can be converted to that class).
#'
#' @param data a data.frame, environment, or list.
#' @param m the number of splits on the data, m = 10 by default.
#' @param B the number of bootstraps, B = 5000 by default.
#' @param parallel boolean, specify whether to use parallelization or not, parallel = FALSE by default
#'
#' @return a blblm object
#' @export
#' @examples
#'blblm(mpg ~ wt * hp, data = mtcars, m = 10, B = 10000, parallel = TRUE)
#'blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {
  data_list <- split_data(data, m)
  if(parallel){
    estimates<- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B=B)
    )
  }else{
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}




#' split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Print the blblm model (formula) of the fitting model
#'
#' @importFrom utils capture.output
#'
#' @param x a blblm object
#' @param ... other conditions
#'
#' @export
#' @return return the blblm model (formula) of x
#' @method print blblm
#' @examples
#' print(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Compute the residual standard deviation (sigma) of a fitted model
#'
#' @param object a blblm object, the fitted model
#' @param confidence boolean, specify whether to compute the confidence interval
#' @param level double, the significance level 1-alpha, 0.95 by default
#' @param ... other conditions
#'
#' @return return the residual standard deviation (sigma) of the object (the fitted model)
#' @export
#' @method sigma blblm
#'
#' @examples
#'sigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
#'sigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), confidence = TRUE)
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Obtain the coefficients of the fitted model
#'
#' @param object a blblm object, the fitted model
#' @param ... other conditions
#'
#' @return return the of object (the fitted model)
#' @export
#' @method coef blblm
#'
#' @examples
#'coef(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' Obtain the confidence intervals for the parameters of the fitted model
#'
#' @param object a blblm object, the fitted model
#' @param parm a specification of which parameters, e.g a string vector, NULL by default
#' @param level double, the significance level 1-alpha, 0.95 by default
#' @param ... extra conditions
#'
#' @return return the confidence intervals for the param of the object
#' @export
#' @method confint blblm
#'
#' @examples
#' confint(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), c("wt", "hp"))
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' Model predictions for fitted models generated by blblm
#'
#'Acquire the predicted value by providing new data and the fitted model that generated by blblm.
#'
#' @param object a blblm object, the fitted model
#' @param new_data a data.frame, environment, or list
#' @param confidence boolean, specify whether to compute the confidence interval
#' @param level double, the significance level 1-alpha, 0.95 by default
#' @param ... other conditions
#'
#' @return return the predicted value of the object (the fitted model)
#' @export
#' @method predict blblm
#' @examples
#'predict(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), data.frame(wt = c(2.5, 3), hp = c(150, 170)))
#'predict(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}





