test_that("functions in blbglm file work", {
  iris_sample <- sample(iris, 10000, replace = T)
  fit <- blbglm(Species ~ Petal.Width * Petal.Length, data = iris_sample, m = 3, B = 100, parallel = T)
  expect_s3_class(fit, "blbglm")
  data_m <- split_data(data = iris, m=3)
  expect_equal(length(data_m), 3)
})
