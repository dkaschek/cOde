\dontrun{

##############################################################################################
## Ozone formation and decay, modified by external forcings
##############################################################################################

library(deSolve)
data(forcData)

# O2 + O <-> O3
f <- c(
  O3 = " (build_O3 + u_build) * O2 * O - (decay_O3 + u_degrade) * O3",
  O2 = "-(build_O3 + u_build) * O2 * O + (decay_O3 + u_degrade) * O3",
  O  = "-(build_O3 + u_build) * O2 * O + (decay_O3 + u_degrade) * O3"
)

# Generate ODE function
forcings <- c("u_build", "u_degrade")
func <- funC(f, forcings = forcings)

# Initialize times, states, parameters and forcings
times <- seq(0, 8, by = .1)
yini <- c(O3 = 0, O2 = 3, O = 2)
pars <- c(build_O3 = 1/6, decay_O3 = 1)
forc <- setForcings(func, forcData)

# Solve ODE
out <- odeC(y = yini, times = times, func = func, parms = pars, forcings = forc)

# Plot solution
par(mfcol=c(1,2))
t1 <- unique(forcData[,2])
M1 <- matrix(forcData[,3], ncol=2)
t2 <- out[,1]
M2 <- out[,-1]

matplot(t1, M1, type="l", lty=1, col=1:2, xlab="time", ylab="value", 
	main="forcings", ylim=c(0, 1.5))
legend("topright", legend = c("u_build", "u_degrade"), lty=1, col=1:2)
matplot(t2, M2, type="l", lty=1, col=1:3, xlab="time", ylab="value", 
	main="response")
legend("topright", legend = c("O3", "O2", "O"), lty=1, col=1:3)

}
