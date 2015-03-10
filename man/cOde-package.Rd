\name{cOde-package}
\alias{cOde-package}
\alias{cOde}
\docType{package}
\title{
Automated C Code Generation for Use with the "deSolve" and "bvpSolve" Packages
}
\description{
Generates all necessary C functions allowing the user to work with the compiled-code interface of ode() and bvptwp(). The implementation supports "forcings" and "events". The package also provides functions to symbolically compute Jacobians, sensitivity equations and adjoint sensitivities being the basis for sensitivity analysis. 
}
\details{
\tabular{ll}{
Package: \tab cOde\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2015-03-05\cr
License: \tab GNU Public License 2 or above\cr
}

The system of ordinary differential equations is defined as a named character vector. This character vector is translated in C syntax and the code output format is matched to the requirements of deSolve and bvpSolve.

}
\author{
Daniel Kaschek <daniel.kaschek@physik.uni-freiburg.de>
}

\keyword{ package }
