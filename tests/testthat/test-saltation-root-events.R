context("saltation for root-triggered events")

# dx/dt = -k * x^2 * time, with the additive reset x := x + v fired by the
# root xc - x = 0. Both x(t) and every first-order sensitivity have closed
# forms; we use those as ground truth and require agreement to integrator
# precision - no finite differences.

test_that("symbolic sensitivities match analytic solution through root event", {

  skip_if_not_installed("deSolve")

  owd <- setwd(tempdir())
  on.exit(setwd(owd))

  eqns <- c(x = "-k*x^2*time")

  events <- data.frame(
    var    = "x",
    time   = NA,
    value  = "v",
    method = "add",
    root   = "xc - x",
    stringsAsFactors = FALSE
  )

  f.sens <- sensitivitiesSymb(eqns, events = events)

  sens_events <- do.call(rbind, attr(f.sens, "events"))
  user_events <- cbind(events, use_pre_state = FALSE)
  cols <- c("var", "time", "value", "method", "root", "use_pre_state")
  all_events <- rbind(user_events[, cols], sens_events[, cols])

  func <- funC(c(eqns, f.sens), events = all_events,
               modelname = "salttest_analytic")

  pars  <- c(k = 1, v = 1, xc = 0.25)
  yini  <- c(x = 1, attr(f.sens, "yini"))
  # Skip t just around tstar where the analytic piecewise flips - the
  # integrator locates the root to within rtol/atol, not exactly at tstar.
  tstar <- sqrt(2 * (yini["x"] - pars["xc"]) /
                (pars["k"] * yini["x"] * pars["xc"]))
  times <- seq(0, 5, length.out = 400)
  times <- times[abs(times - tstar) > 5e-3]

  out <- odeC(yini, times, func, pars,
              rtol = 1e-10, atol = 1e-12, method = "lsodar")

  # Analytic reference (single-segment piecewise).
  analytic <- function(t, x0, k, v, xc) {
    tstar <- sqrt(2 * (x0 - xc) / (k * x0 * xc))
    pre   <- t < tstar
    D     <- k * t^2 * x0 * xc * (v + xc) - 2 * v * x0 + 2 * xc * (v + xc)

    list(
      x    = ifelse(pre, 2 * x0 / (k * t^2 * x0 + 2),
                         2 * x0 * xc * (v + xc) / D),
      x.x  = ifelse(pre, 4 / (k * t^2 * x0 + 2)^2,
                         4 * xc^2 * (v + xc)^2 / D^2),
      x.k  = ifelse(pre, -2 * t^2 * x0^2 / (k * t^2 * x0 + 2)^2,
                         -2 * t^2 * x0^2 * xc^2 * (v + xc)^2 / D^2),
      x.v  = ifelse(pre, 0 * t,
                         4 * x0^2 * xc^2 / D^2),
      x.xc = ifelse(pre, 0 * t,
                         2 * x0 * (xc * D -
                         (v + xc) * (k * t^2 * x0 * xc^2 + 2 * v * x0 + 2 * xc^2))
                         / D^2)
    )
  }

  ref <- analytic(times, x0 = yini["x"], k = pars["k"],
                  v = pars["v"], xc = pars["xc"])

  tol <- 1e-5
  for (col in c("x", "x.x", "x.k", "x.v", "x.xc")) {
    expect_lt(max(abs(out[, col] - ref[[col]])), tol,
              label = paste0("max|symb - analytic| for ", col))
  }
})
