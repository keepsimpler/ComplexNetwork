library(ComplexNetwork)

context('Stabilities')

#' @title test localization
test_that('localization of more smooth vector should be lesser', {
  # localization of a is less than that of b
  a = c(1/2, 1/2, 1/2, 1/2)
  #sum(a^2) == 1
  expect_equal(sum(a^2), 1)
  b = c(0, 1/sqrt(2), 1/ sqrt(2),0)
  #sum(b^2) == 1
  expect_equal(sum(b^2), 1)
  #sum(a^4) < sum(b^4)
  expect_lt(sum(a^4), sum(b^4))
})