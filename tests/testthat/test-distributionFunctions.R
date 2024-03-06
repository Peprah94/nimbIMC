devtools::load_all( )

# Unit test for customised Poisson distribution
test_that("Poisson distribution check", {
#density function
checkdpois <- compileNimble(dmypois)
expect_equal(is.numeric(checkdpois(x = c(3, 5, 10), lambda = 2, log = 1)), TRUE) # result should be numeric

#random function
checkrpois <- compileNimble(rmypois)
expect_equal(class(checkrpois(n = 10, lambda = 2, includeLowerBound = 0)), "numeric")

expect_error(checkrpois(n = 10, lambda = 2, c(0,1),includeLowerBound = 1))
})


# Unit test for customised Binomial distribution
test_that("Poisson distribution check", {
  #density function
  checkdbinom <- compileNimble(dmybinom)
  expect_equal(is.numeric(checkdbinom(x = c(1,2,3,1,1), size = 4, prob = 0.5,log = 1)), TRUE) # result should be numeric
  expect_equal(checkdbinom(x = c(1,2,3,1,1), size = 1, prob = 0.5,log = 1), -Inf)
  expect_equal(checkdbinom(x = c(1,2,3,1,1), size = 3, prob = 0,log = 1), -Inf)
  #random function
  checkrbinom <- compileNimble(rmybinom)
  expect_equal(class(checkrbinom(n = 10, size = 1,prob = 0.2)), "integer")
})
