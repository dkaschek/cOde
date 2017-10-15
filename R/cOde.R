#' Generate C code for a function and compile it
#' 
#' @param f Named character vector containing the right-hand sides of the ODE.
#'   You may use the key word \code{time} in your equations for non-autonomous
#'   ODEs.
#' @param forcings Character vector with the names of the forcings
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
#' @param nGridpoints Integer, defining the number of grid points between tmin
#'   and tmax where the ODE is computed in any case. Indicates also the number
#'   of spline nodes if \code{fcontrol = "einspline"}.
#' @param precision Numeric. Only used when \code{fcontrol = "einspline"}.
#' @param modelname Character. The C file is generated in the working directory
#'   and is named <modelname>.c. If \code{NULL}, a random name starting with
#'   ".f" is chosen, i.e. the file is hidden on a UNIX system.
#' @param verbose Print compiler output to R command line.
#' @param solver Select the solver suite as either \code{deSolve} or
#'   \code{Sundials}. Defaults to \code{deSolve}.
#' @details The function replaces variables by arrays \code{y[i]}, etc. and
#'   replaces "^" by pow() in order to have the correct C syntax. The file name
#'   of the C-File is derived from \code{f}. I.e. \code{funC(abc, ...} will
#'   generate a file abc.c in the current directory. Currently, only explicit
#'   ODE specification is supported, i.e. you need to have the right-hand sides
#'   of the ODE.
#'   
#' @return the name of the generated shared object file together with a number
#'   of attributes
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
funC <- function(f, forcings = NULL, fixed = NULL, outputs=NULL, 
                 jacobian=c("none", "full", "inz.lsodes", "jacvec.lsodes"), 
                 rootfunc = NULL, boundary = NULL, 
                 compile = TRUE, fcontrol = c("nospline", "einspline"),
                 nGridpoints = 500, precision = 1e-5, modelname = NULL,
                 verbose = FALSE, solver = c("deSolve", "Sundials")) {
  
  f <- unclass(f)
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
  symbols <- getSymbols(c(f, rootfunc, constraints))
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
    f <- replaceOperation("^", "pow", f)
    f <- replaceOperation("**", "pow", f)
    f <- replaceSymbols(variables, paste0("y[", 1:length(variables) - 1, "]"), f)
    f <- replaceSymbols(names(constraints), paste0("cons[", 1:dc - 1, "]"), f)
    f <- replaceSymbols(forcings, forc.replace, f)
    f <- replaceSymbols(forcings.t, forc.t.replace, f)
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
  cat("/** Code auto-generated by cOde", as.character(packageVersion("cOde")), "**/\n")
  cat(paste(includings, "\n"))
  cat("\n")
  cat(paste("static double parms[", dv+dp,"];\n", sep=""))
  cat(paste("static double forc[", di,"];\n", sep=""))
  cat(paste("static double cons[", dc,"];\n", sep=""))
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
  
  if(!is.null(rootfunc)) {
    
    cat("/** Root function **/\n")
    cat(paste0("void ", modelname, "_myroot(int *neq, double *t, double *y, int *ng, double *gout, double *out, int *ip ) {\n"))
    cat("\n")
    cat("\n")
    if(fcontrol == "einspline") {
      cat("\t double x[nSplines];\n")
      cat("\t double xdot[nSplines];\n")
      cat(paste0("\t ", modelname, "_evaluateSplines(t, x, xdot);\n"))
    }
    cat("\t double time = *t;\n")
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
    ## ---- Sundials::cvodes: write code ----
    # Construct sensitivity equations
    fSens <- sensitivitiesSymb(f, 
                               states = setdiff(variables, fixed), 
                               parameters = setdiff(parameters, fixed), 
                               inputs = forcings,
                               reduce = FALSE)
    variablesSens <- names(fSens)

    
    #### Assemble source code
    dllname <- paste0(modelname, "_sdcv")
    program <- c(sundialsIncludes(),
                 sundialsOde(f, variables, parameters, dllname),
                 sundialsSensOde(fSens, variables, variablesSens, parameters, dllname))
    
    if (jacobian != "none") {
      program <- c(program,
                   sundialsJac(f, variables, parameters, dllname))
    }
    
    # First unload possible loaded DLL before writing C file (otherwise we have problems on Windows)
    .so <- .Platform$dynlib.ext
    test <- try(dyn.unload(paste0(dllname, .so)), silent = TRUE)
    if (!inherits(test, "try-error")) message(paste("A shared object with name", dllname, "already existed and was overwritten."))
    
    
    ## Write source code
    filename <- paste0(dllname, ".cpp")
    cppFile <- file(filename, open = "w")
    writeLines(c(program,
                 "}" # Matches extern "C" {
    ), cppFile)
    close(cppFile)
    
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
  attr(f, "outputs") <- outputs
  attr(f, "jacobian") <- jacobian
  attr(f, "inz") <- inz
  attr(f, "boundary") <- boundary
  attr(f, "rootfunc") <- rootfunc
  attr(f, "nGridpoints") <- nGridpoints
  attr(f, "fcontrol") <- fcontrol
  attr(f, "solver" ) <- solver
  attr(f, "modelname") <- modelname
  
  if (solver == "Sundials") {
    attr(f, "equationsSens") <- fSens
    attr(f, "variablesSens") <- variablesSens
    #attr(f, "adrDynamics") <- getNativeSymbolInfo("dynamics", PACKAGE = dllname)$address
    #attr(f, "adrSensitivities") <- getNativeSymbolInfo("sensitivities", PACKAGE = dllname)$address
    #attr(f, "adrDynamicsJac") <- if (jacobian != "none") getNativeSymbolInfo("dynamicsJac", PACKAGE = dllname)$address else NULL
    # Attribute deriv is used in odeC to decide if sensitivities are to be
    # integrated or not. As f is stored in func by odemodel(), deriv is set
    # to false. odemodel() stores this f also in extended with deriv = TRUE.
    attr(f, "deriv") <- FALSE 
  }
  
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



#' Sundials::cvodes: Replace powers and symbols to get correct C++ syntax for
#' cvodes.
#' 
#' @param f Named character vector.
#' @param variables Variables appearing on the ODE, \code{names(f)}.
#' @param parameters Parameters appearing in the ODE.
#' @param varSym Symbol of the variables to appear in the source code.
#' @param parSym Symbol of the parameter to appear in the source code.
#' @param forcings Inhomogenity of the ODE. Currently not used.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
cvodesSyntax <- function(f, variables, parameters, varSym = "y", parSym = "p", forcings = NULL) {
  f <- replaceOperation("^", "pow", f)
  f <- replaceOperation("**", "pow", f)
  f <- replaceSymbols(variables, paste0(varSym, "[", 1:length(variables) - 1, "]"), f)
  f <- replaceSymbols(parameters, paste0(parSym, "[", 1:length(parameters) - 1, "]"), f)
  return(f)
}



#' Infos and includes for sundials::cvodes C++ code: ODE, Jacobian.
#' 
#' @return C++ source code as a character vector.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @importFrom utils packageVersion
sundialsIncludes <- function() {
  ## General preparations
  #### User info
  prolog <- paste("/**",
                  paste0("** Code auto-generated by cOde ", as.character(packageVersion("cOde")), "."),
                  "** Do not edit by hand.",
                  "**/", sep = "\n")
  
  
  #### Includes and function signature
  includes <- paste("#include <array>",
                    "#include <vector>",
                    "#include <cmath>",
                    "using std::array;",
                    "using std::vector;", sep = "\n")
  
  
  return(c(prolog, includes, "\n\n", 'extern "C" {\n'))
}



#' Implement the ODE system for sundials::cvodes.
#' 
#' @param f Named character vector containing the right-hand sides of the ODE.
#'   You may use the key word.
#' @param variables Variables appearing on the ODE, \code{names(f)}.
#' @param parameters Parameters appearing in the ODE.
#' @param modelname Base name of the dll.
#' @return C++ source code as a character vector.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
sundialsOde <- function(f, variables, parameters, modelname) {
  ## Header
  odeHead <- paste("/** Derivatives **/",
                    paste0("vector<double> ", modelname, "_dynamics(const double& time, const vector<double>& y,"),
                    "                               const vector<double>& p, const vector<double>& f) {",
                    "    vector<double> ydot(y.size());", sep = "\n")
  
  
  ## ODE system
  f <- cvodesSyntax(f, variables, parameters)
  dv <- length(variables)

  odeSystem <- paste0("    ydot[", 0:(dv - 1),"] = ", f,";")
  
  
  ## Return statement
  ret <- paste("    array<vector<double>, 2> output = {{ydot, ydot}};",
                  "    return output;",
                  "}", sep = "\n")
  ret <- paste("    return ydot;",
               "}", sep = "\n")
  
  
  #### Return entire program
  return(c(odeHead, odeSystem, ret, "\n"))
}



#' Implement sensitivities for an ODE system, sundials::cvodes.
#' 
#' @param f Named character vector containing the right-hand sides of the ODE. 
#'   You may use the key word.
#' @param variablesOde Variables appearing on the ODE of the dynamic system.
#' @param variablesSens Variables appearing on the ODE of the sensitivities of
#'   the dynamic system.
#' @param parameters Parameters appearing in the ODE.
#' @param modelname Base name of the dll.
#' @return C++ source code as a character vector.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
sundialsSensOde <- function(f, variablesOde, variablesSens, parameters, modelname) {
  ## Header
  odeHead <- paste("/** Derivatives of sensitivities **/",
                   paste0("vector<double> ", modelname, "_sensitivities (const double& time,"),
                   "                              const vector<double>& y, const vector<double>& yS,",
                   "                              const vector<double>& p) {",
                   "    vector<double> ySdot(y.size() * (y.size() + p.size()));", sep = "\n")
  
  
  ## Sensitivity ODE system
  f <- cvodesSyntax(f, variablesOde, parameters)
  f <- cvodesSyntax(f, variablesSens, parameters, "yS")
  hasSens.idx <- c(attr(f, "hasSensStatesIdx"), attr(f, "hasSensParametersIdx"))
  hasSens.position <- (1:length(hasSens.idx))[hasSens.idx]
  sensOde <- paste0("    ySdot[", hasSens.position - 1,"] = ", f,";")
  
  
  ## Return statement
  ret <- paste("    return ySdot;",
               "}", sep = "\n")
  
  
  #### Return entire program
  return(c(odeHead, sensOde, ret, "\n"))
}



#' Implement the Jacobian of the ODE system for sundials::cvodes.
#' 
#' @param f Named character vector containing the right-hand sides of the ODE.
#'   You may use the key word.
#' @param variables Variables appearing on the ODE, \code{names(f)}.
#' @param parameters Parameters appearing in the ODE.
#' @param modelname Base name of the dll.
#' @return C++ source code as a character vector.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
sundialsJac <- function(f, variables, parameters, modelname) {
  ## Header
  jacHead <- paste("/** Jacobian **/",
                   paste0("vector<double> ", modelname, "_dynamicsJac(const double& time, const std::vector<double>& y, "),
                   "                        const std::vector<double>& p,",
                   "                        const std::vector<double>& f) {",
                   "    vector<double> yJac(y.size()*y.size());", sep = "\n")
  
  
  ## Jacobian of the ODE system
  # On the cvodes side, the jacobian is filled up column by column which is in
  # line with the format of jac
  jac  <- jacobianSymb(f)
  jac <- cvodesSyntax(jac, variables, parameters)
  not.zero.jac <- which(jac != "0")
  jacobian <- paste0("    yJac[", not.zero.jac - 1,"] = ", jac[not.zero.jac],";")
  
  
  ## Return statement
  ret <- paste("    return yJac;",
               "}", sep = "\n")
  
  
  #### Return entire program
  return(c(jacHead, jacobian, ret, "\n"))
}



#' Generate interpolation spline for the forcings and write into list of matrices
#' 
#' @param func result from funC()
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric)
#' @return list of matrices with time points and values assigned to the forcings interface of deSolve
#' @details Splines are generated for each name in forcings and both, function value and first
#' derivative are evaluated at the time points of the data frame.
#' @examples
#' f <- c(x = "-k*x + a - b")
#' func <- funC(f, forcings = c("a", "b"))
#' forcData <- rbind(
#'   data.frame(name = "a", time = c(0, 1, 10), value = c(0, 5, 2)),
#'   data.frame(name = "b", time = c(0, 5, 10), value = c(1, 3, 6)))
#' forc <- setForcings(func, forcData) 
#' @export
#' @importFrom stats splinefun
setForcings <- function(func, forcings) {
  
  #loadDLL(func)
  
  inputs <- attr(func, "forcings")
  fcontrol <- attr(func, "fcontrol")
  nGridpoints <- attr(func, "nGridpoints")
  trange <- range(forcings$time)
  tspan <- seq(trange[1], trange[2], len=nGridpoints)
  
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
#' @param times vector of type numeric. Integration times
#' @param func return value from funC()
#' @param parms named vector of type numeric.
#' @param ... further arguments going to \code{ode()}
#' @details See deSolve-package for a full description of possible arguments
#' @return matrix with times and states
#' @example inst/examples/example1.R
#' @export
#' @import Rcpp
#' @useDynLib cOde
odeC <- function(y, times, func, parms, ...) {
  
  # Sundials::cvodes
  if (attr(func, "solver") == "Sundials") {
    
    
    if (any(is.na(y))) {
      stop("At least one state is not provided with an initial value.")
    }

    if (any(is.na(parms))) {
      stop("At least one parameter value is unspecified.")
    }
    

    # Extract cvodes settings from ...
    varargs <- list(...)
    varnames <- names(varargs)
    
    
    # Handle events
    if (all(varnames != "events")) {
      events <- data.frame()
    }
    else {
      events <- varargs$events$data
      # On the C++ side, we can't work with species names, only with their
      # position in the systems of odes.
      events$var <- match(events$var, names(y)) - 1
      methods <- c("replace", "add", "multiply")
      events$method <- match(events$method, methods)
    }
    
    
    settings <- list(rtol = 1e-6,
                     atol = 1e-6,
                     maxsteps = 1e3,
                     maxord = 5,
                     hini = 0,
                     hmin = 0,
                     hmax = 100,
                     maxerr = 5,
                     maxnonlin = 10,
                     maxconvfail = 10,
                     method = "bdf",
                     jacobian = if (attr(func, "jacobian") == "none") FALSE
                                else TRUE,
                     minimum = -1e-4,
                     positive = 1,
                     which_states = length(attr(func, "equations")),
                     which_observed = 0,
                     stability = TRUE,
                     sensitivities = if (attr(func, "deriv")) TRUE else FALSE
                     )

    userSettings <- intersect(names(settings), varnames)
    settings[userSettings] <- varargs[userSettings]
    
    
    
    
    # Check if a model is present
    adrDynamics <- try(getNativeSymbolInfo(name = paste0(func, "_dynamics")), silent = TRUE)
    if (inherits(adrDynamics, "try-error")) {
      stop("No model provided. Create one by calling odemodel().")
    }
    
    
    # Check consistency of user setting and user provided jacobian.
    if (settings["jacobian"] == TRUE ) {  
      adrDynamicsJac <- try(getNativeSymbolInfo(name = paste0(func, "_dynamicsJac")), silent = TRUE)
      if (inherits(adrDynamicsJac, "try-error")) {
        stop("You said to provide the jacobian, but you did not.")
      }
    } else {
      adrDynamicsJac <- NULL
    }
    
    
    # Check consistency of user setting and user provided sensitivities.
    if (settings["sensitivities"] == TRUE ) {  
      adrSensitivities <- try(getNativeSymbolInfo(name = paste0(func, "_sensitivities")), silent = TRUE)
      if (inherits(adrSensitivities, "try-error")) {
        stop("You said to provide sensitivities, but you did not.")
      }
    } else {
      adrSensitivities <- NULL
    }
    
    
    # Setup sensitivities
    # This is a pain, as deSolve and Sundials have different approaches.
    # While the parameter y of odeC seems to provied initals for sensitivities,
    # for Sundials they are not usefull, as only initials for sensitivities are
    # reported which are not reduced by reduceSensitivities. Therefore, we have
    # to construct from scratch.
    # initSens must not be NULL when passed to wrap_cvodes
    nStates <- length(attr(func, "equations"))
    nPars <- length(parms)
    initSens <- c(diag(nStates), matrix(0, nStates, nPars))
    
    
    #cat("times:", times, "\n")
    #cat("states:", y[names(attr(func, "equations"))], "\n")
    #cat("parms:", parms[attr(func, "parameters")], "\n")

    
    # Call integrator
    out = wrap_cvodes(times = times,
                      states_ = y[names(attr(func, "equations"))],
                      parameters_ = parms[attr(func, "parameters")],
                      initSens_ = initSens,
                      events_ = events,
                      settings = settings,
                      model_ = as.list(adrDynamics)$address,
                      jacobian_ = as.list(adrDynamicsJac)$address,
                      sens_ = as.list(adrSensitivities)$address)


    # Prepare return
    whichStates <- settings[["which_states"]]
    whichObserved <- settings[["which_observed"]]
    outnames <- c("time",
                  attr(func, "variables")[0:whichStates],
                  make.unique(rep("obs", whichObserved), sep = ""))
    if (settings["sensitivities"] == TRUE) {
      outnames <- c(outnames, attr(func, "variables")[-(1:whichStates + whichObserved)])
    }
    colnames(out) <- outnames

    return(out)
  }
  
  # deSolve
  nGridpoints <- attr(func, "nGridpoints")
  modelname <- attr(func, "modelname")
  times.inner <- seq(min(c(times, 0)), max(times), len=nGridpoints)
  times.inner <- sort(unique(c(times, times.inner)))
  which.times <- match(times, times.inner)
  yout <- c(attr(func, "forcings"), names(attr(func, "outputs")))
  
  y <- y[attr(func, "variables")]
  parms <- parms[attr(func, "parameters")]
  parms <- c(parms, rep(0, length(y)))
  
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
  if(any(names(moreargs)=="forcings") & attr(func, "fcontrol") == "einspline") 
    moreargs <- moreargs[-which(names(moreargs) == "forcings")]
  
  i <- match(names(moreargs), names(arglist))
  is.overlap <- which(!is.na(i))
  is.new <- which(is.na(i))
  arglist[i[is.overlap]] <- moreargs[is.overlap]
  arglist <- c(arglist, moreargs[is.new])
  
  #loadDLL(func)
  
  out <- do.call(deSolve::ode, arglist)
  out.index <- unique(c(which.times[which.times <= nrow(out)], nrow(out)))
  out <- matrix(out[out.index, ], nrow = length(out.index), dimnames = list(NULL, colnames(out)))
  
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

  print(initforc)
   
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

