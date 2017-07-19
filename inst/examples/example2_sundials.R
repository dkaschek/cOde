\dontrun{

######################################################################
## Sensitivity analysis of ozone formation with Sundials solver
######################################################################


# O2 + O <-> O3
f <- c(
  O3 = " build_O3 * O2 * O - decay_O3 * O3",
  O2 = "-build_O3 * O2 * O + decay_O3 * O3",
  O  = "-build_O3 * O2 * O + decay_O3 * O3"
)

# Generate ODE function
func <- funC(f, solver = "Sundials", modelname = "ozon")

# Change attributes of func to cause the solver evaluate sensitivities
extended <- func
attr(extended, "deriv") <- TRUE
attr(extended, "variables") <- c(
  attr(func, "variables"), 
  attr(func, "variablesSens"))

# Initialize times, states, parameters and forcings
times <- seq(0, 15, by = .1)
yini <- c(O3 = 0, O2 = 3, O = 2)
pars <- c(build_O3 = .1, decay_O3 = .01)


# Solve ODE without sensitivities
out <- odeC(yini, times, func, pars, method = "bdf")
# Solve ODE with sensitivities
out <- odeC(yini, times, extended, pars, method = "bdf")

# Plot solution
par(mfcol=c(2,3))
t <- out[,1]
M1 <- out[,2:4]
M2 <- out[,5:7]
M3 <- out[,8:10]
M4 <- out[,11:13]
M5 <- out[,14:16]
M6 <- out[,17:19]

matplot(t, M1, type="l", lty=1, col=1:3, 
        xlab="time", ylab="value", main="solution")
legend("topright", legend = c("O3", "O2", "O"), lty=1, col=1:3)
matplot(t, M2, type="l", lty=1, col=1:3, 
        xlab="time", ylab="value", main="d/(d O3)")
matplot(t, M3, type="l", lty=1, col=1:3, 
        xlab="time", ylab="value", main="d/(d O2)")
matplot(t, M4, type="l", lty=1, col=1:3, 
        xlab="time", ylab="value", main="d/(d O)")
matplot(t, M5, type="l", lty=1, col=1:3, 
        xlab="time", ylab="value", main="d/(d build_O3)")
matplot(t, M6, type="l", lty=1, col=1:3, 
        xlab="time", ylab="value", main="d/(d decay_O3)")

}
