#' @aliases NULL
#' @import purrr
#' @import furrr
#' @import tidyverse
#' @import vroom
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @details
#' Logistic Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Bag of Little Bootstraps for Logistic Regression
#' Use Bag of Little Bootstraps to do Logistic Regression
#' @param formula a formula
#' @param data a data frame or a character vector
#' @param m an integer giving the number of groups to split the data in
#' @param B an integer that represents the number of bootstraps
#' @param parallel logical; if TRUE, run plan(multiprocess, workers = 4) before running the blbglm function. In this example, workers = 4 but you can change the value of workers to be the number of workers you want.
#'
#' @return blbglm
#' @export
blbglm <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {
  if (is.character(data)) {
    data_list <- data
    n <- sum(map_dbl(data_list, ~ { length(vroom_lines(., altrep = TRUE, progress = FALSE)) - 1L}))
  }
  else {
    n <- nrow(data)
    data_list <- split_data(data, m)
  }
  if (parallel) {
    mapfunc <- future_map
  }  else {
    mapfunc <- map
  }

  estimates <- mapfunc(
    data_list,
    function(data) {
      if (is.character(data)) {
        data <- vroom(.)
      }
      glm_each_subsample(formula = formula, data = data, n = n, B = B)
    }
  )
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}


#' Splitting Data into Approximately Equal Sizes
#' split data into m parts of approximated equal sizes
#' @param data a data frame or a character vector
#' @param m an integer giving the number of groups to split the data in
#'
#' @return m data frames of approximately equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' Compute the estimates
#' Compute the estimates for the given formula
#' @param formula a formula
#' @param data a data frame
#' @param n numeric
#' @param B numeric
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}


#' Regression estimates for blbglm data set
#' Compute the regression estimates for blbglm data set
#' @param formula a formula
#' @param data a data frame
#' @param n an integer
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs, n)
}


#' Fitting a model
#' Fitting a model using bag of little bootstraps
#' @param formula an object of class called formula
#' @param data a data frame, list, or character vector
#' @param freqs a numeric value
#' @param n a numeric value
#'
#' @return a list of the coefficients and the sigma values
glm1 <- function(formula, data, freqs, n) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- suppressWarnings(glm(formula, data, family = binomial, weights = freqs, maxit = 100))
  while(!fit$converged){
    freqs <- rmultinom(1, n, rep(1, nrow(data)))
    fit <- suppressWarnings(glm(formula, family = binomial, data, weights = freqs, maxit = 100))
  }
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' Compute the coefficients from fit
#' Obtain the coefficients of the parameters of fit
#' @param fit a fitted blbglm model
#'
#' @return numeric value, return the value of the coefficients of the parameters
blbcoef <- function(fit) {
  coef(fit)
}


#' Residual Standard Deviation Sigma
#' Computes the Residual Standard Deviation Sigma
#' @param fit a fitted blbglm model
#'
#' @return a numeric value, the value of the residual standard deviation sigma
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' Print Blbglm Values
#' Print values returned by blbglm model.
#' @param x an R object
#' @param ... additional arguments
#'
#' @export
#' @method print blbglm
print.blbglm <- function(x, ...) {
  cat("blbglm model:", capture.output(x$formula))
  cat("\n")
}


#' Obtain the value of the Standard Deviation Sigma
#' Compute the value of the Standard Deviation Sigma for the Model
#' @param object an object model
#' @param confidence logical, if TRUE you construct the confidence interval of sigma
#' @param level the confidence level
#' @param ... additional arguments
#'
#' @return A numeric value
#' @export
#' @method sigma blbglm
sigma.blbglm <- function(object, confidence = FALSE, level = 0.95, ...) {
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

#' Coefficient of Model
#' Find the coefficients of the model object.
#' @param object an object in which coefficient extraction will be performed on
#' @param ... additional arguments
#'
#' @return the coefficient extract from the object
#' @export
#' @method coef blbglm
coef.blbglm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Model Parameter Confidence Intervals
#' Compute the confidence intervals for a parameter in the model.
#' @param object a fitted model object for which a confidence interval is wanted
#' @param parm the parameters that you want to compute the confidence interval of
#' @param level numeric, the confidence level
#' @param ... additonal arguments
#'
#' @return A matrix or vector that give the lower and upper confidence bounds for the parameter
#' @export
#' @method confint blbglm
confint.blbglm <- function(object, parm = NULL, level = 0.95, ...) {
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

# create a function that calculates the probability using the inverse logit expression
invlogit <- function(x) 1/(1+exp(-x))

#' Prediction Interval for Model Parameters
#' Make prediction intervals for model paramters
#' @param object a model object for which prediction is wanted
#' @param new_data a data frame or character vector
#' @param confidence logical
#' @param level numeric, the prediction level
#' @param ... additonal arguments that impact the predictions
#'
#' @return A matrix or vector with columns, depends on the class of the argument
#' @export
#' @method predict blbglm
predict.blbglm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ invlogit(X %*% .$coef)) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ invlogit(X %*% .$coef)) %>% rowMeans())
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