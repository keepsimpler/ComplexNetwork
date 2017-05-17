#library(igraph)
library(ComplexNetwork)

context("Jacobians")

#' @title test that the positive-negative interactions should be symmetric
test_that('positive-negative interactions should be symmetric', {
  pn <- function(N, C) {
    adj_type = 'gnm'
    mat_type = 'pn'
    dist_type = 'binary'

    J <- gen_Jacobian(N = N, C = C, adj_type = adj_type, mat_type = mat_type, dist_type = dist_type)
    all(t(J) + J == 0)  # should be symmetric
  }
  
  expect_equal(pn(100, 0.1), TRUE)
})
