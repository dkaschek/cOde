
#' Compute sensitivity equations of a function symbolically
#' 
#' @param f named vector of type character, the functions
#' @param states Character vector. Sensitivities are computed with respect to initial
#' values of these states
#' @param parameters Character vector. Sensitivities are computed with respect to initial
#' values of these parameters
#' @param inputs Character vector. Input functions or forcings. They are excluded from
#' the computation of sensitivities.
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric or character, time point), 
#' "value" (numeric or character, value), "method" (character, either
#' "replace" or "add"). See \link[deSolve]{events}.
#' Within \code{sensitivitiesSymb()} a \code{data.frame} of additional events is generated to 
#' reset the sensitivities appropriately, depending on the event method. 
#' @param reduce Logical. Attempts to determine vanishing sensitivities, removes their
#' equations and replaces their right-hand side occurences by 0.
#' @details The sensitivity equations are ODEs that are derived from the original ODE f.
#' They describe the sensitivity of the solution curve with respect to parameters like 
#' initial values and other parameters contained in f. These equtions are also useful
#' for parameter estimation by the maximum-likelihood method. For consistency with the
#' time-continuous setting provided by \link{adjointSymb}, the returned equations contain
#' attributes for the chisquare functional and its gradient.
#' @return Named vector of type character with the sensitivity equations. Furthermore,
#' attributes "chi" (the integrand of the chisquare functional), "grad" (the integrand
#' of the gradient of the chisquare functional), "forcings" (Character vector of the 
#' additional forcings being necessare to compute \code{chi} and \code{grad}) and "yini" (
#' The initial values of the sensitivity equations) are returned.
#' @example inst/examples/example2.R
#' @example inst/examples/example3.R
#' @export
sensitivitiesSymb <- function(f, states = names(f), parameters = NULL, inputs = NULL, events = NULL, reduce = FALSE) {
  
  # If f contains line breaks, replace them by ""
  f <- gsub("\n", "", f)
  
  variables <- names(f)
  states <- states[!states%in%inputs]
  
  if (is.null(parameters)) {
    pars <- getSymbols(c(f,
                         as.character(events[["value"]]),
                         as.character(events[["time"]]),
                         as.character(events[["root"]])),
                       exclude = c(variables, inputs, "time"))
  } else {
    pars <- parameters[!parameters%in%inputs]
  }
  
  if (length(states) == 0 & length(pars) == 0)
    stop("Attempt to compute sensitivities although both states and parameters had length 0.")
  
  Dyf <- jacobianSymb(f, variables)
  Dpf <- jacobianSymb(f, pars)
  
  df <- length(f)
  dv <- length(variables)
  ds <- length(states)
  dp <- length(pars)
  
  # generate sensitivity variable names and names with zero entries then
  # write sensitivity equations in matrix form
  Dy0y <- Dpy <- NULL
  sensParVariablesY0 <- sensParVariablesP <- NULL
  
  if (ds > 0) {
    mygridY0 <- expand.grid.alt(variables, states)
    sensParVariablesY0 <- apply(mygridY0, 1, paste, collapse = ".")
    Dy0y <- matrix(sensParVariablesY0, ncol = ds, nrow = dv)
  }
  
  if (dp > 0) {
    mygridP <- expand.grid.alt(variables, pars)
    sensParVariablesP <- apply(mygridP, 1, paste, collapse = ".")
    Dpy <- matrix(sensParVariablesP, ncol = dp, nrow = dv)
  }
  
  
  gl <- NULL
  if (!is.null(Dy0y)) {
    gl <- c(gl, as.vector(prodSymb(matrix(Dyf, ncol = dv), Dy0y)))
  }
  if (!is.null(Dpy)) {
    gl <- c(gl, as.vector(sumSymb(prodSymb(matrix(Dyf, ncol = dv), Dpy), matrix(Dpf, nrow = dv))))
  }
  
  newfun <- gl
  newvariables.grid <- expand.grid.alt(variables, c(states, pars))
  newvariables <- apply(newvariables.grid, 1, paste, collapse=".")
  names(newfun) <- newvariables
  
  
  events.addon <- eventframe <- NULL
  if (!is.null(events)) {
    
    if (is.null(events[["root"]])) events[["root"]] <- NA
    
    events.addon <- lapply(1:nrow(events), function(i) {
      
      # A data frame representing the i'th event
      myevent <- events[i,]
      
      # Get the info from the event
      xone <- as.character(myevent[["var"]])
      tau <- as.character(myevent[["time"]])
      xi <- as.character(myevent[["value"]])
      root <- as.character(myevent[["root"]])
      xk <- setdiff(variables, xone)
      rootstate <- intersect(getSymbols(root), variables)[1]
      norootstate <- setdiff(variables, rootstate)
      rootpar <- intersect(getSymbols(root), pars)[1]
      
      # Derive some quantities needed for all event methods
      Ji1 <- jacobianSymb(f, variables = xone)
      J11 <- Ji1[paste0(xone, ".", xone)]
      Jk1 <- Ji1[paste0(xk, ".", xone)]
      f1 <- f[xone]
      fk <- f[xk]
      fr <- f[rootstate]
      
      # Collect ODE parameters 
      odepars <- c(states, pars)
      
      
      # Generate additional events for **replacement event**
      if (myevent[["method"]] == "replace") {
        
        
        d <- list()
        
        vars <- intersect(paste0(xone, ".", odepars), newvariables)
        if (length(vars) > 0 & is.na(root)) {
          d[[1]] <- data.frame(
            var = vars,
            time = tau,
            value = 0,
            root = root,
            method = "replace"
          )
        }
        
        # Resetting of tau sensitivities for multiple roots
        vars <- intersect(paste0(c(xone, xk), ".", tau), newvariables)
        if (length(vars) > 0 & !is.na(root)) {
          d[[2]] <- data.frame(
            var = vars,
            time = tau,
            value = 0,
            root = root,
            method = "replace"
          )
        }
        
        
        vars <- intersect(paste0(xone, ".", tau), newvariables)
        if (length(vars) > 0) {
          d[[3]] <- data.frame(
            var = vars,
            time = tau,
            value = paste0("(", J11, ")*(", xone, "- (", xi, ")) - (", f1, ")"),
            root = root,
            method = "add"
          )
        }
          
        vars <- intersect(paste0(xk, ".", tau), newvariables)
        if (length(vars) > 0) {
          d[[4]] <- data.frame(
            var = vars,
            time = tau,
            value = paste0("(", Jk1, ")*(", xone, "- (", xi, "))"),
            root = root,
            method = "add"
          )
        }
        
        vars <- intersect(paste0(xone, ".", xi), newvariables)
        if (length(vars) > 0) {
          d[[5]] <- data.frame(
            var = vars,
            time = tau,
            value = 1,
            root = root,
            method = "add"
          )
        }
        
        
        vars <- intersect(paste0(c(xone, xk), ".", rootpar), newvariables)
        if (length(vars) > 0) {
          d[[6]] <- do.call(rbind, mapply(function(mystate, mypar) {
            data.frame(
              var = paste0(mystate, ".", mypar),
              time = tau,
              value = paste0("(eventcounter__ + 1)*", mystate, ".", tau , "/(", fr, ")"),
              root = root,
              method = "replace"
            )
          }, 
          mystate = sapply(vars, function(v) strsplit(v, ".", fixed = TRUE)[[1]][1]),
          mypar =  sapply(vars, function(v) strsplit(v, ".", fixed = TRUE)[[1]][2]),
          SIMPLIFY = FALSE))
          d[[6]] <- d[[6]][order(d[[6]][["method"]], decreasing = TRUE), ]
        }
   
        excludedPars <- NULL
        if (!is.na(root)) excludedPars <- c(rootpar, tau, xi)
        vars <- outer(c(xone, xk), setdiff(odepars, excludedPars), function(x, y) paste(x, y, sep = "."))
        vars <- intersect(vars, newvariables)
        # sort first non-root states then root state
        vars <- c(vars[!grepl(paste0("^", rootstate, "\\."), vars)], vars[grepl(paste0("^", rootstate, "\\."), vars)])
        if (length(vars) > 0 & !is.na(root)) {
          d[[7]] <- do.call(rbind, mapply(function(mystate, mypar) {
            data.frame(
              var = paste0(mystate, ".", mypar),
              time = tau,
              value = paste0(mystate, ".", tau, "*(-", rootstate, ".", mypar, ")/(", fr, ")"), 
              root = root,
              method = ifelse(mystate == xone, "replace", "add")
            )
          }, 
          mystate = sapply(vars, function(v) strsplit(v, ".", fixed = TRUE)[[1]][1]),
          mypar =  sapply(vars, function(v) strsplit(v, ".", fixed = TRUE)[[1]][2]),
          SIMPLIFY = FALSE))
          #d[[7]] <- d[[7]][order(d[[7]][["method"]], decreasing = TRUE), ]
        }
        
        
        d[["stringsAsFactors"]] <- FALSE
        
        return(do.call(rbind, d))
        
        
        
      }
      
      # Generate additional events for **additive event**
      else if (myevent[["method"]] == "add") {
        
        d <- list()
        
        vars <- intersect(paste0(xone, ".", odepars), newvariables)
        if (length(vars) > 0) {
          d[[1]] <- data.frame(
            var = vars,
            time = tau,
            value = 0,
            root = root,
            method = "add"
          )
        }
        
        vars <- intersect(paste0(xone, ".", tau), newvariables)
        if (length(vars) > 0) {
          d[[2]] <- data.frame(
            var = vars,
            time = tau,
            value = paste0("-(", J11, ") * (", xi, ")"),
            root = root,
            method = "add"
          )
        }
        
        vars <- intersect(paste0(xk, ".", tau), newvariables)
        if (length(vars) > 0) {
          d[[3]] <- data.frame(
            var = vars,
            time = tau,
            value = paste0("-(", Jk1, ") * (", xi, ")"),
            root = root,
            method = "add"
          )
        }
        
        vars <- intersect(paste0(xone, ".", xi), newvariables)
        if (length(vars) > 0) {
          d[[4]] <- data.frame(
            var = vars,
            time = tau,
            value = 1,
            root = root,
            method = "add"
          )
        }
        
        
        d[["stringsAsFactors"]] <- FALSE
        
        return(do.call(rbind, d))
        
        
        
      }
      
      # generate additional events for **multiplicative event**
      else if (myevent[["method"]] == "multiply") {
        
        d <- list()
        
        vars <- intersect(paste0(xone, ".", odepars), newvariables)
        if (length(vars) > 0) {
          d[[1]] <- data.frame(
            var = vars,
            time = tau,
            value = xi,
            root = root,
            method = "multiply" 
          )
        }
        
        vars <- intersect(paste0(xone, ".", tau), newvariables)
        if (length(vars) > 0) {
          d[[2]] <- data.frame(
            var = vars,
            time = tau,
            value = paste0("(1 - (", xi, "))*((", J11, ")*(", xone, ") - (", f1, "))"),
            root = root,
            method = "add"
          )
        }
        
        vars <- intersect(paste0(xk, ".", tau), newvariables)
        if (length(vars) > 0) {
          d[[3]] <- data.frame(
            var = vars,
            time = tau,
            value = paste0("(1 - (", xi, "))*((", Jk1, ")*(", xone, "))"),
            root = root,
            method = "add"
          )
        }
        
        vars <- intersect(paste0(xone, ".", xi), newvariables)
        if (length(vars) > 0) {
          d[[4]] <- data.frame(
            var = vars,
            time = tau,
            value = xone,
            root = root,
            method = "add"
          )
        }
        
        d[["stringsAsFactors"]] <- FALSE
        
        
        
        return(do.call(rbind, d))
        
        
        
      }
      
      else {
        
        stop("Event method must be either 'replace', 'add' or 'multiply'.")
        
      }
      
    })
    
    # Overwrite events
    events <- events.addon
    
    # Add rownames that allow to trace back records in eventframe to the list of events
    for (i in seq_along(events)) {
      rownames(events[[i]]) <- paste(i, seq_along(events[[i]][[1]]), sep = "_")
    }
    
    # Make sure all columns are characters
    eventframe <- do.call(rbind, events)
    for (i in seq_along(eventframe)) {
      eventframe[[i]] <- as.character(eventframe[[i]])
    }
    
    
    # Get (numeric) values in eventframe value column
    symbols <- getSymbols(eventframe$value)
    symbols.vals <- structure(stats::rnorm(length(symbols)), names = symbols)
    values.char <- paste0("c(", paste(eventframe$value, collapse = ", "), ")")
    values.vals <- with(as.list(symbols.vals), eval(parse(text = values.char)))
    
    # Get sensitivities which are initialized by 0 and do not change due to additional events
    is_ini_zero <- sapply(strsplit(eventframe$var, ".", fixed = TRUE), function(x) x[1] != x[2])
    var_reset_to_zero <- sapply(unique(eventframe$var[is_ini_zero]), function(myvar) {
      all(values.vals[eventframe$var == myvar] == 0)
    })
    
    # Determine records which can be removed from eventframe either because they are neutral or the variable can completely be removed
    is_neutral <- (eventframe$method == "add" & values.vals == 0)
    
    # Reduce eventframe
    # eventframe <- eventframe[!is_neutral,]
    
    # Reduce events (by name matching)
    events <- lapply(events, function(e) e[rownames(e) %in% rownames(eventframe), ])
    
    # Check if any root-like events
    # events <- lapply(events, function(e){
    #   if (all(is.na(eventframe[["root"]]))) {
    #     e <- e[-match("root", names(e))]
    #   }
    #   return(e)
    # })
    
    # No factors in events
    events <- lapply(events, function(e) {
      e <- lapply(e, as.character)
      e <- as.data.frame(e, stringsAsFactors = FALSE)
      return(e)
    })
    
    
  }
 
  
  # Reduce the sensitivities
  #vanishing <- c(sensParVariablesY0[!(sensParVariablesY0 %in% eventframe[["var"]][eventframe[["value"]] != "0"])], 
  #               sensParVariablesP[Dpf == "0" & !(sensParVariablesP %in% eventframe[["var"]][eventframe[["value"]] != "0"])])
  vanishing <- setdiff(c(sensParVariablesY0, sensParVariablesP[Dpf == "0"]), eventframe[["var"]])
  
  if(reduce) {
    newfun <- reduceSensitivities(newfun, vanishing)
    is.zero.sens <- names(newfun) %in% attr(newfun,"is.zero")
  } else {
    is.zero.sens <- rep(FALSE, length(newfun))
  }
  events <- lapply(events, function(myevents) {
    myevents[!as.character(myevents[["var"]]) %in% names(newfun)[is.zero.sens], ]
  })
  newfun <- newfun[!is.zero.sens]
  output.reduction <- structure(rep(0, length(which(is.zero.sens))), names = newvariables[is.zero.sens])
  
  # Append initial values
  initials <- rep(0, length(newfun))
  names(initials) <- newvariables[!is.zero.sens]
  ones <- which(apply(newvariables.grid, 1, function(row) row[1] == row[2]))
  initials[newvariables[ones]] <- 1
  
  
  # Construct index vector for states and parameters indicating non-vanishing
  # sensitivity equations.
  # States
  hasSens.stateNames <- intersect(sensParVariablesY0, names(initials))
  hasSens.StatesIdx <- rep(FALSE, length(sensParVariablesY0))
  names(hasSens.StatesIdx) <- sensParVariablesY0
  hasSens.StatesIdx[hasSens.stateNames] <- TRUE
  
  # Parameters
  hasSens.parameterNames <- intersect(sensParVariablesP, names(initials))
  hasSens.ParameterIdx <- rep(FALSE, length(sensParVariablesP))
  names(hasSens.ParameterIdx) <- sensParVariablesP
  hasSens.ParameterIdx[hasSens.parameterNames] <- TRUE
  
  
  # Compute wrss
  pars <- c(pars, states)
  
  statesD <- paste0(states, "D")
  weightsD <- paste0("weight", statesD)
  
  res <- paste0(weightsD,"*(", states, "-", statesD, ")")
  sqres <- paste0(res, "^2")
  chi <- c(chi = paste0(sqres, collapse =" + "))
    
  sensitivities <- lapply(pars, function(p) paste0(states, ".", p))
  names(sensitivities) <- pars
    
  grad <- sapply(pars, function(p) paste0(paste("2*", res, "*", sensitivities[[p]]), collapse=" + "))
  names(grad) <- paste("chi", pars, sep=".")
  grad <- replaceSymbols(newvariables[is.zero.sens], "0", grad)
  
  attr(newfun, "chi") <- chi
  attr(newfun, "grad") <- grad
  attr(newfun, "outputs") <- output.reduction
  attr(newfun, "forcings") <- c(statesD, weightsD)
  attr(newfun, "yini") <- initials
  attr(newfun, "events") <- events
  attr(newfun, "hasSensStatesIdx") <- hasSens.StatesIdx
  attr(newfun, "hasSensParametersIdx") <- hasSens.ParameterIdx
    
  
  return(newfun)
  
  
}




#' Compute adjoint equations of a function symbolically
#' 
#' @param f Named vector of type character, the functions
#' @param states Character vector of the ODE states for which observations are available
#' @param inputs Character vector of the "variable" input states, i.e. time-dependent parameters
#' (in contrast to the forcings).
#' @param parameters Character vector of the parameters
#' @details The adjoint equations are computed with respect to the functional 
#' \deqn{(x, u)\mapsto \int_0^T \|x(t)-x^D(t)\|^2 + \|u(t) - u^D(t)\|^2 dt,}{(x, u) -> int( ||x(t) - xD(t)||^2 + ||u(t) - uD(t)||^2, dt),} 
#' where x are the states being constrained
#' by the ODE, u are the inputs and xD and uD indicate the trajectories to be best
#' possibly approached. When the ODE is linear with respect to u, the attribute \code{inputs}
#' of the returned equations can be used to replace all occurences of u by the corresponding
#' character in the attribute. This guarantees that the input course is optimal with
#' respect to the above function.
#' @return Named vector of type character with the adjoint equations. The vector has attributes
#' "chi" (integrand of the chisquare functional), "grad" (integrand of the gradient of the chisquare functional),
#' "forcings" (character vector of the forcings necessary for integration of the adjoint equations) and
#' "inputs" (the input expressed as a function of the adjoint variables).
#' @example inst/examples/example5.R
#' @export
adjointSymb <- function(f, states=names(f), parameters = NULL, inputs=NULL) {
  
  n <- length(f)
  
  adjNames <- paste0("adj", names(f))
  
  ## Compute adjoint sensitivities  
  jac <- matrix(jacobianSymb(f), n)
  
  negadj <- as.vector(prodSymb(t(jac), matrix(adjNames, ncol=1)))
  adj <- paste0("-(", negadj, ")")
  names(adj) <- adjNames
  
  
  negres <- paste0("(", states, "D - ", states, ")") 
  wres <- paste(negres, paste0("weight", states, "D"), sep="*") 
  
  
  adj[paste0("adj", states)] <- paste(adj[paste0("adj", states)], wres, sep="+")
  
  ## Input equations
  u <- NULL
  if(!is.null(inputs)) {
    
    jac <- matrix(jacobianSymb(f, inputs), n )
    u <- as.vector(
      sumSymb(paste0(inputs, "D"), 
              matrix(paste0("-(",  
                            as.vector(
                              prodSymb(t(jac), matrix(adjNames, ncol=1))),
                            ")*eps/weight", inputs, "D"),
                     ncol=1)
      )
    )
    u <-  paste0("(", u, ")")
    names(u) <- inputs
    
  }
  
  ## Forcings required by the BVP solver
  forcings <- paste0(c(states, inputs), "D")
  forcings <- c(forcings, paste0("weight", forcings))
  
  
  
  ## Time-continous log-likelihood
  res <- paste0("(", c(states, inputs), "-", c(states, inputs), "D)")
  sres <- paste0(res, "^2")
  wsres <- paste(paste0("weight", c(states, inputs), "D"),sres,  sep="*") 
  
  chi <- paste(wsres, collapse = " + ")
  names(chi) <- "chi"
  attr(adj, "chi") <- chi
  
  ## Adjoint "gradient wrt parameters" equations
  if(is.null(parameters)) {
    symbols <- getSymbols(f)
    parameters <- symbols[!(symbols%in%inputs) & !(symbols%in%names(f))]
  }
  if(length(parameters)>0) {
    jac <- matrix(jacobianSymb(f, parameters), n)
    grad <- as.vector(prodSymb(matrix(adjNames, nrow=1), jac))
    gradP <- paste0("2*(", grad, ")")
    names(gradP) <- paste("chi", parameters, sep=".")
    attr(adj, "grad") <- gradP  
  }
  
  
  
  attr(adj, "forcings") <- forcings
  attr(adj, "inputs") <- u
  
  
  
  
  return(adj)
  
}


#'@title Forcings data.frame
#'@name forcData
#'@docType data
#'@description Forcings data.frame
#' 
NULL

#'@title Time-course data of O, O2 and O3
#'@name oxygenData
#'@docType data 
#'@description Forcings data.frame
NULL
