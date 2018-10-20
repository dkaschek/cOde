myD <- function(eq, x) {
  
  MAX <- function(x, y) y + 0.5*(sign(x-y)+1)*(x-y)
  MIN <- function(x, y) y + 0.5*(sign(y-x)+1)*(x-y)
  
  Deriv::Deriv(eq, "x", cache.exp = FALSE)
  
}


myD("a+MIN(x, a)", "a")

eq <- "a*x^2 + b*x + x^n/(a^n + x^n)"

"x^n/(K^n + x^n) = K^n * (x/K)^n/(1 + (x/K)^n) = K^n*fn(x^n/K^n) with fn(x) = x/(1+x)"
Deriv("log(exp(x))^n", "n")

x <- seq(-10, 10, .1)
plot(x, pmin(x, 3))


plot(x, y)
plot(x, )
