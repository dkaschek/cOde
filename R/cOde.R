#' Generate C code for a function and compile it
#' 
#' @param f Named character vector containing the right-hand sides of the ODE.
#'   You may use the key word \code{time} in your equations for non-autonomous
#'   ODEs.
#' @param forcings Character vector with the names of the forcings
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric or character, time point), 
#' "value" (numeric or character, value), "method" (character, either
#' "replace" or "add"). See \link[deSolve]{events}. If "var" and "time" are
#' characters, their values need to be speciefied in the parameter vector
#' when calling \code{\link{odeC}}. An event function is generated and compiled
#' with the ODE.
#' @param fixed  character vector with the names of parameters (initial values
#'   and dynamic) for which no sensitivities are required (will speed up the
#'   integration).
#' @param outputs Named character vector for additional output variables, see
#'   arguments \code{nout} and \code{outnames} of \link[deSolve]{lsode}
#' @param jacobian Character, either "none" (no jacobian is computed), "full"
#'   (full jacobian is computed and written as a function into the C file) or
#'   "inz.lsodes" (only the non-zero elements of the jacobian are determined,
#'   see \link[deSolve]{lsodes})
#' @param rootfunc Named character vector. The root function (see
#'   \link[deSolve]{lsoda}). Besides the variable names (\code{names(f)}) also
#'   other symbols are allowed that are treated like new parameters.
#' @param boundary data.frame with columns name, yini, yend specifying the
#'   boundary condition set-up. NULL if not a boundary value problem
#' @param compile Logical. If FALSE, only the C file is written
#' @param fcontrol Character, either \code{"nospline"} (default, forcings are
#'   handled by deSolve) or \code{"einspline"} (forcings are handled as splines
#'   within the C code based on the einspline library).
#' @param nGridpoints Integer, defining for which time points the ODE is evaluated
#' or the solution is returned: Set \code{-1} to return only the explicitly requested
#' time points (default). If additional time points are introduced through events, they
#' will not be returned. Set \code{>= 0} to introduce additional time points between tmin
#' and tmax where the ODE is evaluated in any case. Additional time points that might be
#' introduced by events will be returned.
#' If splines are used with \code{fcontrol = "einspline"}, \code{nGridpoinnts} also
#' indicates the number of spline nodes.
#' @param includeTimeZero Logical. Include t = 0 in the integration time points if \code{TRUE} 
#' (default). Consequently,
#' integration starts at t = 0 if only positive time points are provided by the user and at
#' tmin, if also negtive time points are provided.
#' @param precision Numeric. Only used when \code{fcontrol = "einspline"}.
#' @param modelname Character. The C file is generated in the working directory
#'   and is named <modelname>.c. If \code{NULL}, a random name starting with
#'   ".f" is chosen, i.e. the file is hidden on a UNIX system.
#' @param verbose Print compiler output to R command line.
#' @param solver Select the solver suite as either \code{deSolve} or
#'   \code{Sundials} (not available any more). Defaults to \code{deSolve}.
#' @details The function replaces variables by arrays \code{y[i]}, etc. and
#'   replaces "^" by pow() in order to have the correct C syntax. The file name
#'   of the C-File is derived from \code{f}. I.e. \code{funC(abc, ...} will
#'   generate a file abc.c in the current directory. Currently, only explicit
#'   ODE specification is supported, i.e. you need to have the right-hand sides
#'   of the ODE.
#'   
#' @return the name of the generated shared object file together with a number
#'   of attributes
#'
#' @examples 
#' \dontrun{
#' # Exponential decay plus constant supply
#' f <- c(x = "-k*x + supply")
#' func <- funC(f, forcings = "supply")
#' 
#' # Example 2: root function
#' f <- c(A = "-k1*A + k2*B", B = "k1*A - k2*B")
#' rootfunc <- c(steadyState = "-k1*A + k2*B - tol")
#' 
#' func <- funC(f, rootfunc = rootfunc, modelname = "test")
#' 
#' yini <- c(A = 1, B = 2)
#' parms <- c(k1 = 1, k2 = 5, tol = 0.1)
#' times <- seq(0, 10, len = 100)
#' 
#' odeC(yini, times, func, parms)
#' }
#' @export
funC <- function(f, forcings = NULL, events = NULL, fixed = NULL, outputs=NULL, 
                 jacobian=c("none", "full", "inz.lsodes", "jacvec.lsodes"), 
                 rootfunc = NULL, boundary = NULL, 
                 compile = TRUE, fcontrol = c("nospline", "einspline"),
                 nGridpoints = -1, includeTimeZero = TRUE, precision = 1e-5, modelname = NULL,
                 verbose = FALSE, solver = c("deSolve", "Sundials")) {
  
  if (!is.null(events[["root"]]) & !is.null(rootfunc))
    stop("Root functions cannot be used in both events and rootfunc argument.")
  
  
  
  f <- unclass(f)
  
  # If f contains line breaks, replace them by ""
  f <- gsub("\n", "", f)
  
  constraints <- NULL # Might be an interesting option in the future
  myattr <- attributes(f)
  equations <- f
  
  if("names"%in%names(myattr)) myattr <- myattr[-which(names(myattr)=="names")]
  
  if(is.null(modelname)) modelname <- paste(c("f", sample(c(letters, 0:9), 8, TRUE)), collapse="")
  dllname <- modelname
  
  ## If boundary conditions are given, sort for leftbc first
  if(!is.null(boundary)) {
    leftbc <- boundary$name[!is.na(boundary$yini)]
    f <- c(f[names(f)%in%leftbc], f[!names(f)%in%leftbc])
  }
  
  ## Check which kind of forcings are used
  fcontrol <- match.arg(fcontrol)
  forcings.t <- paste(forcings, "t", sep=".")
  forc.replace <- forcings
  forc.t.replace <- forcings.t
  if(fcontrol == "einspline" & !is.null(forcings)) {
    forc.replace <- paste0("x[", 1:length(forcings)-1, "]")
    forc.t.replace <- paste0("xdot[", 1:length(forcings.t)-1, "]")
  }
  
     
  ## Analyze f by parser
  variables <- names(f)
  symbols <- getSymbols(c(f, rootfunc, constraints, as.character(events[["value"]]), as.character(events[["time"]]), as.character(events[["root"]])))
  parameters <- symbols[!symbols%in%c(variables, forcings, names(constraints), names(rootfunc), "time")]


  jac <- NULL
  inz <- NULL
  jacobian <- match.arg(jacobian)
  if(jacobian != "none") jac  <- jacobianSymb(f)
  if(jacobian %in% c("inz.lsodes", "jacvec.lsodes")) {
    jac.matrix <- matrix(jac, length(f), length(f))
    inz <- apply(jac.matrix, 2, function(v) which(v != "0"))
    inz <- do.call(rbind, lapply(1:length(inz), function(j) if(length(inz[[j]]) > 0) cbind(i = inz[[j]], j = j)))
  }
  not.zero.jac <- which(jac != "0")
  
  dv <- length(variables)
  dp <- length(parameters)
  if(is.null(forcings)) di <- 0 else di <- length(forcings)
  if(is.null(constraints)) dc <- 0 else dc <- length(constraints)
  if(is.null(outputs)) do <- 0 else do <- length(outputs)
  if(is.null(rootfunc)) dr <- 0 else dr <- length(rootfunc)

  
  
  solver <- match.arg(solver)
  if (solver == "deSolve") {

  ## ------------ deSolve: write code -------------
  
  deSolveSyntax <- function(f) {
    f <- replaceNumbers(f)
    f <- gsub("++", "+", f, fixed = TRUE)
    f <- gsub("+-", "-", f, fixed = TRUE)
    f <- gsub("-+", "-", f, fixed = TRUE)
    f <- gsub("--", "+", f, fixed = TRUE)
    f <- replaceOperation("^", "pow", f)
    f <- replaceOperation("**", "pow", f)
    f <- replaceSymbols(variables, paste0("y[", format(1:length(variables) - 1, trim = TRUE, scientific = FALSE), "]"), f)
    f <- replaceSymbols(names(constraints), paste0("cons[", format(1:dc - 1, trim = TRUE, scientific = FALSE), "]"), f)
    f <- replaceSymbols(forcings, forc.replace, f)
    f <- replaceSymbols(forcings.t, forc.t.replace, f)
    f <- gsub("&", " & ", f, fixed = TRUE)
    return(f)
  }
    
  f <- deSolveSyntax(f)
  
  if(jacobian %in% c("full", "jacvec.lsodes")) 
    jac <- deSolveSyntax(jac)
  
  if(!is.null(constraints)) 
    constraints <- deSolveSyntax(constraints)
  
  if(!is.null(outputs)) 
    outputs <- deSolveSyntax(outputs)
  
  if(!is.null(rootfunc))
    rootfunc <- deSolveSyntax(rootfunc)
  
  
  triggerByRoot <- FALSE
  if(!is.null(events)) {
    
    # Check events if only standard events (no root)
    if (is.null(events[["root"]]) || all(is.na(events[["root"]]))) {
      
      triggerByRoot <- FALSE
      var <- deSolveSyntax(as.character(events[["var"]]))
      time <- deSolveSyntax(as.character(events[["time"]]))
      value <- deSolveSyntax(as.character(events[["value"]]))
      method <- as.character(events[["method"]])
      time.unique <- unique(time)
      eventsfn <- sapply(1:length(time.unique), function(i) {
        t <- time.unique[i]
        paste0("\t if(*t == ", t, " & eventcounter[", i-1, "] == 0) {\n",
               paste(sapply(which(time == t), function(tix) {
                 value[tix] <- gsub("eventcounter__", paste0("eventcounter[", i-1, "]"), value[tix])
                 switch(
                   method[tix],
                   replace = paste0("\t\t", var[tix], " = ", value[tix], ";"),
                   add = paste0("\t\t", var[tix], " = ", var[tix], " + ", value[tix], ";"),
                   multiply = paste0("\t\t", var[tix], " = (", var[tix], ") * (", value[tix], ");")
                 )
               }), collapse = "\n"),
               "\n",
               "\t\t", "eventcounter[", i-1, "] = eventcounter[", i-1, "] + 1.;\n",
               "\t }\n"
        )
      })
      
    } else {
      
      triggerByRoot <- TRUE
      # Generate rootfun
      roots <- apply(
        X = cbind(as.character(events[["time"]]),
                  as.character(events[["root"]])), 
        MARGIN = 1,
        FUN = function(x) {
          if (!is.na(x[1])) {
            out <- paste0("time - ", x[1])
          }
          if (!is.na(x[2])) {
            out <- x[2]
          }
          if (all(is.na(x))) {
            stop("Either time or root must be provided with events data frame")
          }
          return(out)
        }
      )
      names(roots) <- paste0("constr", seq_along(roots))
      rootfunc <- deSolveSyntax(roots)
      rootfunc <- unique(rootfunc)
      dr <- length(rootfunc)
      # Generate conditions for eventfun
      var <- deSolveSyntax(as.character(events[["var"]]))
      condition <- apply(
        X = cbind(as.character(events[["time"]]),
                  as.character(events[["root"]])), 
        MARGIN = 1,
        FUN = function(x) {
          if (!is.na(x[1])) {
            out <- paste0("(time - ", x[1], ") == 0 & eventcounter[(idx)] == 0")
          }
          if (!is.na(x[2])) {
            out <- paste0("fabs(", x[2], ") < 1e-6")
          }
          if (all(is.na(x))) {
            stop("Either time or root must be provided with events data frame")
          }
          out <- deSolveSyntax(out)
          return(out)
        }
      )
      value <- deSolveSyntax(as.character(events[["value"]]))
      method <- as.character(events[["method"]])
      condition.unique <- unique(condition)
      # Generate eventfun
      eventsfn <- sapply(1:length(condition.unique), function(i) {
        t <- condition.unique[i]
        paste0("\t if(", gsub("(idx)", i-1, t, fixed = TRUE), " & eventcounter[", i-1, "] == 0) {\n",
               paste(sapply(which(condition == t), function(tix) {
                 value[tix] <- gsub("eventcounter__", paste0("eventcounter[", i-1, "]"), value[tix])
                 switch(
                   method[tix],
                   replace = paste0("\t\t", var[tix], " = ", value[tix], ";"),
                   add = paste0("\t\t", var[tix], " = ", var[tix], " + ", value[tix], ";"),
                   multiply = paste0("\t\t", var[tix], " = (", var[tix], ") * (", value[tix], ");")
                 )
               }), collapse = "\n"),
               "\n",
               "\t\t", "eventcounter[", i-1, "] = eventcounter[", i-1, "] + 1.;\n",
               "\t }\n"
        )
      })
      
    }
    
    
  }
  
  
  mypath <- system.file(package="cOde")
  splinefile <- paste0("cat ", mypath,"/code/splineCreateEvaluate.c")
  includings <- c("#include <R.h>",
                  "#include <math.h>")
  
  if(fcontrol == "einspline") 
    includings <- c(includings, "#include <einspline/nubspline.h>")
  
  definitions <- paste0("#define ", c(parameters, paste0("y",0:(dv-1),"_0")), " parms[", 0:(dv+dp-1),"]")
  if(!is.null(forcings)) 
    definitions <- c(definitions, paste0("#define ", forcings, " forc[", 0:(di-1),"]"))
  
  
  
  # First unload possible loaded DLL before writing C file (otherwise we have problems on Windows)
  .so <- .Platform$dynlib.ext
  test <- try(dyn.unload(paste0(dllname, .so)), silent = TRUE)
  if (!inherits(test, "try-error")) message(paste("A shared object with name", dllname, "already existed and was overwritten."))
  
  # Write to C file
  filename <- paste0(dllname, ".c")
  sink(filename)
  cat("/** Code auto-generated by cOde", as.character(utils::packageVersion("cOde")), "**/\n")
  cat(paste(includings, "\n"))
  cat("\n")
  cat(paste("static double parms[", dv+dp,"];\n", sep=""))
  cat(paste("static double forc[", di,"];\n", sep=""))
  cat(paste("static double cons[", dc,"];\n", sep=""))
  if (!is.null(events)) cat(paste("static double eventcounter[", length(eventsfn), "];\n", sep = ""))
  cat("static double range[2];\n")
  cat("\n")
  cat(paste("#define nGridpoints",nGridpoints,"\n"))
  cat(paste("#define nSplines", di, "\n"))
  cat(paste("#define precision", precision, "\n"))
  cat("\n")
  cat(paste(definitions, "\n"))
  cat("#define tmin range[0]\n")
  cat("#define tmax range[1]\n")
  cat("\n")
  if(fcontrol == "einspline") {
    splineprogram <- paste(system(splinefile, intern=TRUE), "\n")
    splineprogram <- gsub("createSplines(", paste0(modelname, "_createSplines("), splineprogram, fixed = TRUE)
    splineprogram <- gsub("evaluateSplines(", paste0(modelname, "_evaluateSplines("), splineprogram, fixed = TRUE)
    cat(splineprogram)
  }
  cat("\n")
  cat(paste0("void ", modelname, "_initmod(void (* odeparms)(int *, double *)) {\n"))
  cat(paste("\t int N=", dv+dp,";\n",sep=""))
  cat("\t odeparms(&N, parms);\n")
  if (!is.null(events)) {
    cat("\t for(int i=0; i<", length(eventsfn), "; ++i) eventcounter[i] = 0;\n", sep = "")
  }
  cat("}\n")
  cat("\n")
  cat(paste0("void ", modelname, "_initforc(void (* odeforcs)(int *, double *)) {\n"))
  cat(paste("\t int N=", di,";\n",sep=""))
  cat("\t odeforcs(&N, forc);\n")
  cat("}\n")
  cat("\n")
  
  
  ## Derivative function
  
  cat("/** Derivatives (ODE system) **/\n")
  cat(paste0("void ", modelname, "_derivs (int *n, double *t, double *y, double *ydot, double *RPAR, int *IPAR) {\n"))
  cat("\n")
  if(fcontrol == "einspline") {
    cat("\t double x[nSplines];\n")
    cat("\t double xdot[nSplines];\n")
    cat(paste0("\t ", modelname, "_evaluateSplines(t, x, xdot);\n"))
  }
  cat("\t double time = *t;\n")
  if(!is.null(constraints)) 
    cat(paste("\t cons[", 1:dc-1, "] = ", constraints, ";\n", sep=""))
  
  cat("\n")
  #if(length(reductions)>0) cat(paste("\t double ", reductions, ";\n", sep=""))
  cat(paste("\t ydot[", 0:(dv-1),"] = ", f,";\n", sep=""))
  cat("\n")
  
  # Return forcings and other outputs (only for IVP)
  if(is.null(boundary)) {
    if(di > 0){
      cat(paste0("\t RPAR[", 0:(di-1),"] = ", forc.replace,";\n"))
    }
    if(do > 0){
      cat("\t for(int i= ",di,"; i < ",do+di,"; ++i) RPAR[i] = 0;\n")
      non.zero.outputs <- which(outputs != "0")
      for(i in non.zero.outputs) 
        cat(paste0("\t RPAR[", di + i - 1, "] = ", outputs[i], ";\n")) 
    }
  }
    
  cat("}\n")
  cat("\n")
  
  ## Jacobian of deriv
  if(jacobian == "full") {
    cat("/** Jacobian of the ODE system **/\n")
    cat(paste0("void ", modelname, "_jacobian (int *n, double *t, double *y, double * df, double *RPAR, int *IPAR) {\n"))
    cat("\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat(paste0("\t ", modelname, "_evaluateSplines(t, x, xdot);\n"))
    }
    cat("\n")
    cat("double time = *t;\n")
    cat("\t int i;\n")
    cat("for(i=0; i<(*n) * (*n); i++) df[i] = 0.;\n")
    cat(paste("\t df[", not.zero.jac-1,"] = ", jac[not.zero.jac],";\n", sep=""))
    cat("\n")
    cat("}\n")
    cat("\n")
  }
  
  ## Jacvec of deriv
  if(jacobian == "jacvec.lsodes") {
    vecs <- lapply(1:dv, function(i) matrix(jac, ncol=dv, nrow=dv)[,i])
    not.zero.vec <- lapply(vecs, function(v) which(v != "0"))
    not.zero.columns <- which(sapply(not.zero.vec, function(v) length(v) > 0))
    cat("/** Jacobian vector of the ODE system **/\n")
    cat(paste0("void ", modelname, "_jacvec (int *neq, double *t, double *y, int *j, int *ian, int *jan, double *pdj, double *yout, int *iout) {\n"))
    cat("\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat(paste0("\t ", modelname, "_evaluateSplines(t, x, xdot);\n"))
    }
    cat("\t double time = *t;\n")
    cat("\t int i;\n")
    cat("\t for(i=0; i<*neq; i++) pdj[i] = 0.;\n")
    
    j <- not.zero.columns[1]
    cat(paste("\t if(*j ==", j, ") {\n"))
    cat(paste("\t pdj[", not.zero.vec[[j]]-1,"] = ", vecs[[j]][not.zero.vec[[j]]],";\n", sep=""))
    cat("\t }\n")
    for(j in not.zero.columns[-1]) {
      cat(paste("\t else if(*j ==", j, ") {\n"))
      cat(paste("\t pdj[", not.zero.vec[[j]]-1,"] = ", vecs[[j]][not.zero.vec[[j]]],";\n", sep=""))
      cat("\t }\n")
    }
    
    cat("\n")
    cat("}\n")
    cat("\n")
    
  }
  
  if (!is.null(events)) {
    
    cat("/** Event function **/\n")
    cat(paste0("void ", modelname, "_myevent(int *n, double *t, double *y) {\n"))
    cat("\n")
    cat("\t double time = *t;\n")
    cat("\n")
    cat(eventsfn, sep = "\n")
    cat("\n")
    cat("}\n")
    
    
  }
  
  if(!is.null(rootfunc)) {
    
    cat("/** Root function **/\n")
    cat(paste0("void ", modelname, "_myroot(int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip ) {\n"))
    cat("\n")
    cat("\t double time = *t;\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat(paste0("\t ", modelname, "_evaluateSplines(t, x, xdot);\n"))
    }
    cat(paste("\t gout[", 0:(dr-1),"] = ", rootfunc,";\n", sep=""))
    cat("\n")
    cat("}\n")
    
  }
  
  if(!is.null(boundary)) {
    
    ## Check length of boundary conditions
    nbc <- length(which(!is.na(c(boundary$yini, boundary$yend))))
    if(nbc != dv) {
      sink()
      warning("Number of boundary conditions not correct\n")
      return()
    }
    
    boundary <- boundary[match(variables, boundary$name),]
    
    leftbc <- which(!is.na(boundary$yini))
    rightbc <- which(!is.na(boundary$yend))
    myorder <- c(leftbc, rightbc)
    
    ## Boundary Condition (for compatibility with bvpSolve)
    
    cat("/** Boundary Conditions **/\n")
    cat(paste0("void ", modelname, "_gsub(int *i, int *n, double *z, double *g, double *RPAR, int *IPAR) {\n"))
    cat("\n")
    cat(paste("\t if (*i==", 1,") *g=z[", myorder[1]-1, "]-y", 0, "_0;\n", sep=""))
    if(dv>1) cat(paste("\t else if (*i==", 2:dv,") *g=z[", myorder[-1]-1, "]-y", 2:dv-1, "_0;\n", sep=""))
    cat("\n")
    cat("}\n")
    cat("\n")
    
    ## Jacobian of Boundary Condition (for compatibility with bvpSolve)
    
    cat("/** Jacobian of the Boundary Conditions **/\n")
    cat(paste0("void ", modelname, "_dgsub(int *i, int *n, double *z, double *dg, double *RPAR, int *IPAR) {\n"))
    cat("\n")
    cat("\t int j;\n")
    cat("\t for (j = 0; j< *n; j++) dg[j] = 0;\n")
    
    cat(paste("\t if (*i==", 1,") dg[", myorder[1]-1, "] = 1.;\n", sep=""))
    if(dv>1) cat(paste("\t else if (*i==", 2:dv,") dg[", myorder[-1]-1, "]=1.;\n", sep=""))
    cat("\n")
    cat("}\n")
  }
  
  
  
  sink()
  

  
  } else if (solver == "Sundials") {
    
    stop("Sundials support has been removed. If you were an active user of the Sundials implementation, please get in touch.")
    
  }
  
  
  ## ---- Compile and load shared object file ----
  if (compile) compileAndLoad(filename, dllname, fcontrol, verbose)
  
  
  ## ----------- function return -----------

  f <- dllname
  attributes(f) <- c(attributes(f), myattr)

  attr(f, "equations") <- equations
  attr(f, "variables") <- variables
  attr(f, "parameters") <- parameters
  attr(f, "forcings") <- forcings
  attr(f, "events") <- events
  attr(f, "triggerByRoot") <- triggerByRoot
  attr(f, "outputs") <- outputs
  attr(f, "jacobian") <- jacobian
  attr(f, "inz") <- inz
  attr(f, "boundary") <- boundary
  attr(f, "rootfunc") <- rootfunc
  attr(f, "nGridpoints") <- nGridpoints
  attr(f, "includeTimeZero") <- includeTimeZero
  attr(f, "fcontrol") <- fcontrol
  attr(f, "solver" ) <- solver
  attr(f, "modelname") <- modelname
  
  return(f)
  
}



#' Compile and load shared object implementing the ODE system.
#' 
#' @param filename Full file name of the source file.
#' @param dllname Base name for source and dll file.
#' @param fcontrol Interpolation method for forcings.
#' @param verbose Print compiler output or not.
#'   
#' @author Daniel Kaschek, \email{daniel.kaschek@@physik.uni-freiburg.de}
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
compileAndLoad <- function(filename, dllname, fcontrol, verbose) {
  if (fcontrol == "nospline") 
    shlibOut <- system(paste0(R.home(component = "bin"), "/R CMD SHLIB ", filename), intern = TRUE)
  if (fcontrol == "einspline")
    shlibOut <- system(paste0(R.home(component = "bin"), "/R CMD SHLIB ", filename, " -leinspline"), intern = TRUE)
  if (verbose) {
    cat(shlibOut)
  }
  
  .so <- .Platform$dynlib.ext
  mydllname <- paste0(dllname, .so)
  soExists <- file.exists(mydllname)
  if (soExists) {
    try(dyn.unload(mydllname), silent = TRUE)
    dyn.load(mydllname)
  }
}




#' Generate interpolation spline for the forcings and write into list of matrices
#' 
#' @param func result from funC()
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric)
#' @return list of matrices with time points and values assigned to the forcings interface of deSolve
#' @details Splines are generated for each name in forcings and both, function value and first
#' derivative are evaluated at the time points of the data frame.
#' @examples
#' \dontrun{
#' f <- c(x = "-k*x + a - b")
#' func <- funC(f, forcings = c("a", "b"))
#' forcData <- rbind(
#'   data.frame(name = "a", time = c(0, 1, 10), value = c(0, 5, 2)),
#'   data.frame(name = "b", time = c(0, 5, 10), value = c(1, 3, 6)))
#' forc <- setForcings(func, forcData) 
#' }
#' @export
#' @importFrom stats splinefun
setForcings <- function(func, forcings) {
  
  #loadDLL(func)
  
  inputs <- attr(func, "forcings")
  fcontrol <- attr(func, "fcontrol")
  nGridpoints <- attr(func, "nGridpoints")
  if (nGridpoints < 2 & fcontrol == "einspline")
    stop("When using spline interpolation, nGridpoints (the number of spline nodes) must be >= 2.")
  
  trange <- range(forcings$time)
  tspan <- seq(trange[1], trange[2], len=max(2, nGridpoints))
  
  times <- NULL
  values <- NULL
  
  
  out <- do.call(c, lapply(inputs, function(i) {
    
    t <- forcings[forcings$name == i, "time"]
    x <- forcings[forcings$name == i, "value"]
    
    if(fcontrol == "nospline") {
      mat <- list(cbind(t, x))
      names(mat) <- i
    }
    if(fcontrol == "einspline") {
      myfun <- splinefun(t, x)
      out <- myfun(tspan)
      mat <- list(cbind(tspan, out))
      names(mat) <- i
    }
 
    return(mat)
    
  }))
  
  
  if(fcontrol == "einspline") {
   
    times <- do.call(c, lapply(out, function(o) o[,1]))
    values <- do.call(c, lapply(out, function(o) o[,2]))
    cfunc <- getNativeSymbolInfo(paste0(func, "_createSplines"))
    #.C(cfunc, as.double(times), as.double(values))
    do.call(".C", list(cfunc, as.double(times), as.double(values)))
    
  }
  
  return(out)
  
}


#' Interface to ode()
#' 
#' @param y named vector of type numeric. Initial values for the integration
#' @param times vector of type numeric. Integration times. If \code{includeTimeZero}
#' is \code{TRUE} (see \link{funC}), the times vector is augmented by t = 0. If
#' \code{nGridpoints} (see \link{funC}) was set >= 0, uniformly distributed time points
#' between the first and last time point are introduced and the solution is returned
#' for these time points, too. Any additional time points that are introduced during
#' integration (e.g. event time points) are returned unless nGridpoints = -1 (the default).
#' 
#' @param func return value from funC()
#' @param parms named vector of type numeric.
#' @param ... further arguments going to \code{ode()}
#' @details See deSolve-package for a full description of possible arguments
#' @return matrix with times and states
#' @example inst/examples/example1.R
#' @export
odeC <- function(y, times, func, parms, ...) {
  
  # Sundials::cvodes
  if (attr(func, "solver") == "Sundials") {
    
    stop("Sundials is not supported any more. If you were an active user of the Sundials implementation, please get in touch.")
    
  }
  
  ## deSolve ----------------------------------------------------------------------
  nGridpoints <- attr(func, "nGridpoints")
  includeTimeZero <- attr(func, "includeTimeZero")
  modelname <- attr(func, "modelname")
  triggerByRoot <- attr(func, "triggerByRoot")
  
  # Evaluate initial values and parameters
  y <- y[attr(func, "variables")]
  parms <- parms[attr(func, "parameters")]
  parms <- c(parms, rep(0, length(y)))
  
  # First handle events to get additional time points
  events <- attr(func, "events")
  eventlist <- NULL
  if (!is.null(events)) {
    
    eventfunc <- paste(func, "myevent", sep = "_")
    time.expr <- parse(text = paste0("c(", paste(as.character(events[["time"]]), collapse = ", "), ")"))
    # evaluate event parameters and times
    eventtime <- with(as.list(parms), eval(time.expr))
    eventtime <- eventtime[!is.na(eventtime)]
    
    if (triggerByRoot) {
      eventlist <- list(func = eventfunc, root = TRUE)
      triggertimes <- unlist(lapply(eventtime, function(mytime) seq(mytime - 1e-5, mytime + 1e-5, 1e-6)))
    } else {
      eventlist <- list(func = eventfunc, time = sort(unique(eventtime)))
      triggertimes <- NULL
    }
    
  }
  
  
  times.outer <- times.inner <- times
  if (includeTimeZero)
    times.inner <- union(times.inner, 0)
  if (nGridpoints > 0)
    times.inner <- union(
      times.inner, 
      seq(min(times.inner), max(times.inner), len=nGridpoints)
    )
  if (!is.null(events))
   times.inner <- union(times.inner, triggertimes)
  
  times.inner <- sort(unique(times.inner))
  which.times <- match(times.outer, times.inner)
  
  yout <- c(attr(func, "forcings"), names(attr(func, "outputs")))
  
  arglist <- list(y = y, times = times.inner, func = paste0(func, "_derivs"), parms = parms, dllname = modelname, initfunc = paste0(func, "_initmod"))
  
  
  if (attr(func, "jacobian") == "full")
    arglist <- c(arglist, list(jacfunc = paste0(func, "_jacobian")))
    
  if (attr(func, "jacobian") == "inz.lsodes") {
    inz <- attr(func, "inz")
    lrw <- 20 + 3*dim(inz)[1] + 20*length(y)
    arglist <- c(arglist, list(sparsetype = "sparseusr", inz = inz, lrw = lrw))
  }
  
  if (attr(func, "jacobian") == "jacvec.lsodes") {
    inz <- attr(func, "inz")
    arglist <- c(arglist, list(sparsetype = "sparseusr", jacvec = paste0(func, "_jacvec"), inz = inz))
  }
    
  if (!is.null(attr(func, "forcings")) & attr(func, "fcontrol") == "nospline") 
    arglist <- c(arglist, list(initforc = paste0(func, "_initforc")))
 
  if(!is.null(attr(func, "rootfunc")))
    arglist <- c(arglist, list(rootfunc = paste0(func, "_myroot"), nroot = length(attr(func, "rootfunc"))))
   
  if (!is.null(yout)) {
    arglist <- c(arglist, list(nout = length(yout), outnames = yout))
  }
   
  
  # Replace arguments and add new ones
  moreargs <- list(...)
  if(length(moreargs) > 0) {
    if (any(names(moreargs)=="forcings") & attr(func, "fcontrol") == "einspline") {
      moreargs <- moreargs[-which(names(moreargs) == "forcings")]
    }
  }
  
  i <- match(names(moreargs), names(arglist))
  is.overlap <- which(!is.na(i))
  is.new <- which(is.na(i))
  arglist[i[is.overlap]] <- moreargs[is.overlap]
  arglist <- c(arglist, moreargs[is.new])
  
  # Merge events and event time points
  arglist[["events"]] <- c(arglist[["events"]], eventlist)
   
  
  #loadDLL(func)
  
  out <- do.call(deSolve::ode, arglist)
  
  if (nGridpoints == -1) {
    # Return only requested time points
    out <- out[match(times.outer, out[, 1]), , drop = FALSE]
  } else {
    # Return all time points between tmin and tmax (additional time points after
    # tmax might be introduced by the integrator)
    out <- out[out[, 1] >= min(times.inner) & out[, 1] <= max(times.inner), ]  
  }
  
  return(out)
  
  
}


#' Interface to bvptwp()
#' 
#' 
#' @param yini named vector of type numeric. Initial values to be overwritten.
#' @param x vector of type numeric. Integration times
#' @param func return value from funC() with a boundary argument. 
#' @param yend named vector of type numeric. End values to be overwritten.
#' @param parms named vector of type numeric. The dynamic parameters.
#' @param xguess vector of type numeric, the x values
#' @param yguess matrix with as many rows as variables and columns as x values
#' @param ... further arguments going to \code{bvptwp()}
#' @details See bvpSolve-package for a full description of possible arguments
#' @return matrix with times and states
#' @example inst/examples/example4.R
#' @export
bvptwpC <- function(yini=NULL, x, func, yend=NULL, parms, xguess=NULL, yguess=NULL,  ...) {
  
  #loadDLL(func)
  
  dynpar <- parms[attr(func, "parameters")]
  boundary <- attr(func, "boundary")
  leftbc <- boundary$name[!is.na(boundary$yini)]
  rightbc <- boundary$name[!is.na(boundary$yend)]
  
  ## Fill yini/yend with values from func. If yini/yend are given,
  ## set their values.
  bini <- boundary$yini
  names(bini) <- boundary$name
  bini <- bini[!is.na(bini)]
  
  bend <- boundary$yend
  names(bend) <- boundary$name
  bend <- bend[!is.na(bend)]
  
  if(!is.null(yini)) bini[names(yini)] <- yini
  if(!is.null(yend)) bend[names(yend)] <- yend
  if(!is.null(attr(func, "forcings")) & attr(func, "fcontrol") == "nospline") 
    initforc <- paste0(func, "_initforc")
  else
    initforc <- NULL
  
  
  
  posbound <- c(rep(min(x), length(bini)), rep(max(x), length(bend)))
  
  
  statepars <- c(bini, bend)
  newparms <- c(dynpar, statepars)
  
  moreargs <- list(...)
  if(any(names(moreargs)=="forcings") & attr(func, "fcontrol") == "einspline") 
    moreargs <- moreargs[-which(names(moreargs)=="forcings")]

  out <- do.call(bvpSolve::bvptwp, c(moreargs, list(
    x = x, parms = newparms, xguess = xguess, yguess = yguess, posbound=posbound,
    func = paste0(func, "_derivs"), jacfunc = paste0(func, "_jacobian"), bound = paste0(func, "_gsub"), jacbound = paste0(func, "_dgsub"), 
    initfunc = paste0(func, "_initmod"), initforc = initforc,
    dllname = func,
    ncomp = length(statepars)
  )))
  

  colnames(out) <- c("x", attr(func, "variables"))
  
  return(out)
  
  
}

