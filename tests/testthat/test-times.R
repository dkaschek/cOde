context("odeC")
library(cOde)
test_that("odeC returns output with times equal to input times", {
  f <- c(
    A = "-k1*A",
    B = "k1*A - k2*B"
  )
  func <- funC(f)
  times <- 0:100
  yini <- c(A = 1, B = 0)
  parms <- c(k1 = 2e-1, k2 = 1e-1)

  out <- odeC(y = yini, times = times, func = func, parms = parms)

  expect_equal(times, out[,1])
})
