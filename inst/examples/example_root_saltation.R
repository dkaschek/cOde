\dontrun{

######################################################################
## Saltation across a root-triggered event, validated against the
## closed-form solution.
##
##   dx/dt = -k * x^2 * time
##
##   Event (additive): when  xc - x = 0   -->   x := x + v
##
## Analytic solution has two segments split by
##   t* = sqrt(2 * (x0 - xc) / (k * x0 * xc))
## and closed forms for x, ∂x/∂x0, ∂x/∂k, ∂x/∂v, ∂x/∂xc exist in
## both segments - used below as ground truth (no finite differences).
######################################################################

library(deSolve)

# funC() writes the generated .c / .so into the current working directory;
# run the example in tempdir() so nothing gets left behind next to the Rmd.
owd <- setwd(tempdir())
on.exit(setwd(owd), add = TRUE)

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

# Merge the sensitivity saltation events with the user's event frame.
sens_events <- do.call(rbind, attr(f.sens, "events"))
user_events <- cbind(events, use_pre_state = FALSE)
cols <- c("var", "time", "value", "method", "root", "use_pre_state")
all_events <- rbind(user_events[, cols], sens_events[, cols])

func <- funC(c(eqns, f.sens), events = all_events, modelname = "salt_analytic")

pars  <- c(k = 1, v = 1, xc = 0.25)
yini  <- c(x = 1, attr(f.sens, "yini"))
times <- sort(c(2.449490, seq(0, 5, length.out = 400)))

out <- odeC(yini, times, func, pars,
            rtol = 1e-10, atol = 1e-12, method = "lsodar")

# ---- closed-form reference -------------------------------------------------
analytic <- function(t, x0, k, v, xc) {
  tstar <- sqrt(2 * (x0 - xc) / (k * x0 * xc))
  pre   <- t < tstar

  # shared post-event denominator
  D <- k * t^2 * x0 * xc * (v + xc) - 2 * v * x0 + 2 * xc * (v + xc)

  x_pre  <- 2 * x0 / (k * t^2 * x0 + 2)
  x_post <- 2 * x0 * xc * (v + xc) / D

  x.x_pre   <- 4 / (k * t^2 * x0 + 2)^2
  x.x_post  <- 4 * xc^2 * (v + xc)^2 / D^2

  x.k_pre   <- -2 * t^2 * x0^2 / (k * t^2 * x0 + 2)^2
  x.k_post  <- -2 * t^2 * x0^2 * xc^2 * (v + xc)^2 / D^2

  x.v_pre   <- 0 * t
  x.v_post  <- 4 * x0^2 * xc^2 / D^2

  x.xc_pre  <- 0 * t
  x.xc_post <- 2 * x0 *
    (xc * D - (v + xc) * (k * t^2 * x0 * xc^2 + 2 * v * x0 + 2 * xc^2)) / D^2

  data.frame(
    time = t,
    x    = ifelse(pre, x_pre,    x_post),
    x.x  = ifelse(pre, x.x_pre,  x.x_post),
    x.k  = ifelse(pre, x.k_pre,  x.k_post),
    x.v  = ifelse(pre, x.v_pre,  x.v_post),
    x.xc = ifelse(pre, x.xc_pre, x.xc_post)
  )
}

ref <- analytic(times, x0 = yini["x"], k = pars["k"], v = pars["v"], xc = pars["xc"])

# Compare all columns - should match to ~solver precision.
for (col in c("x", "x.x", "x.k", "x.v", "x.xc")) {
  err <- max(abs(out[, col] - ref[, col]))
  cat(sprintf("%-5s  max|symb - analytic| = %.3e\n", col, err))
}

matplot(times, cbind(out[, "x.x"], ref[, "x.x"]),
        type = "l", lty = c(1, 2), col = c("black", "red"),
        xlab = "time", ylab = "dx/dx0",
        main = "saltation: cOde (solid) vs analytic (dashed)")
abline(v = sqrt(2 * (yini["x"] - pars["xc"]) / (pars["k"] * yini["x"] * pars["xc"])),
       lty = 3, col = "grey40")

}
