test_that("blblm package works", {
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_s3_class(fit, "blblm")
  data_m <- split_data(data = mtcars, m=3)
  expect_equal(length(data_m), 3)
})
