context("getSymbols")
test_that("getSymbols works for NULL", {

  #-!Start example code
  #-! library(cOde)
  mysymbols <- getSymbols(c("A*AB+B^2"))
  
  #-!End example code


  # Define your expectations here
  expectation_1 <- c("A", "AB", "B")
  expect_identical(mysymbols, expectation_1)
  expect_null(getSymbols(NULL))
})
  