% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{blblm}
\alias{blblm}
\title{Bag of Little Bootstraps for Linear Regression
Use Bag of Little Bootstraps to do Linear Regression}
\usage{
blblm(formula, data, m = 10, B = 5000, parallel = FALSE)
}
\arguments{
\item{formula}{a formula}

\item{data}{a data frame or a character vector}

\item{m}{an integer giving the number of groups to split the data in}

\item{B}{an integer that represents the number of bootstraps}

\item{parallel}{logical; if TRUE, run plan(multiprocess, workers = 4) before running the blblm function. In this example, workers = 4 but you can change the value of workers to be the number of workers you want.}
}
\value{
blblm
}
\description{
Bag of Little Bootstraps for Linear Regression
Use Bag of Little Bootstraps to do Linear Regression
}
\examples{
fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
coef(fit)
}
