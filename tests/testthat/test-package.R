test_that("包可以加载", {
  expect_true(requireNamespace("causalR", quietly = TRUE))
})
