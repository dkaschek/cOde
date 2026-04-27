
#' Compute sensitivity equations of a function symbolically
#' 
#' @param f named vector of type character, the functions
#' @param states Character vector. Sensitivities are computed with respect to initial
#' values of these states
#' @param parameters Character vector. Sensitivities are computed with respect to initial
#' values of these parameters
#' @param inputs Character vector. Input functions or forcings. They are excluded from
#' the computation of sensitivities.
#' @param events data.frame of events with columns \code{var} (character,
#' the name of the state to be affected), \code{time} (numeric or character,
#' event time), \code{value} (numeric or character, value used by the reset),
#' \code{method} (\code{"replace"}, \code{"add"} or \code{"multiply"}) and
#' optional \code{root} (character, a root expression r(x, t, p) = 0).
#' Each row must set exactly one of \code{time} and \code{root} to \code{NA}:
#' a timed event uses \code{root = NA} and a fixed or symbolic time; a
#' root-triggered event uses \code{time = NA} and a non-empty \code{root}.
#' See \link[deSolve]{events}. Within \code{sensitivitiesSymb()} additional
#' events are generated to apply the saltation correction to the forward
#' sensitivities at each event firing.
#' @param reduce Logical. Attempts to determine vanishing sensitivities, removes their
#' equations and replaces their right-hand side occurences by 0.
#' @details The sensitivity equations are ODEs that are derived from the original ODE f.
#' They describe the sensitivity of the solution curve with respect to parameters like
#' initial values and other parameters contained in f. These equtions are also useful
#' for parameter estimation by the maximum-likelihood method. For consistency with the
#' time-continuous setting provided by \link{adjointSymb}, the returned equations contain
#' attributes for the chisquare functional and its gradient.
#' At each event the saltation matrix correction
#' \eqn{S_i^+ = \sum_j (\partial g_i/\partial x_j) S_{j,p} + \partial g_i/\partial p + \Delta_i \cdot d\tau/dp}
#' with
#' \eqn{\Delta_i = \sum_j (\partial g_i/\partial x_j) f_j^- - f_i^+}
#' and, for root-triggered events,
#' \eqn{d\tau/dp = -(\partial r/\partial p + \sum_k (\partial r/\partial x_k) S_{k,p}) /
#'                 (\partial r/\partial t + \sum_k (\partial r/\partial x_k) f_k^-)}
#' is emitted as a reset event with \code{use_pre_state = TRUE} so that
#' \code{\link{funC}} wires the saltation RHS to a pre-event snapshot of the
#' state vector.
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
    event_chars <- c(as.character(events[["value"]]),
                     as.character(events[["time"]]),
                     as.character(events[["root"]]))
    event_chars <- event_chars[!is.na(event_chars) & nzchar(event_chars)]
    pars <- getSymbols(c(f, event_chars),
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

    if (is.null(events[["root"]])) events[["root"]] <- NA_character_

    # -- Saltation via direct chain-rule --
    #
    # For an event at time τ with reset x_i -> g_i(x, p):
    #
    #   S_{i,p}^+ = Σ_j (∂g_i/∂x_j) · S_{j,p}^-
    #             + ∂g_i/∂p
    #             + Δ_i · dτ/dp
    #
    #   Δ_i  = Σ_j (∂g_i/∂x_j) · f_j^-  -  f_i^+
    #
    # Root-triggered (root r(x,t,p) = 0, `time` may be NA):
    #   dτ/dp = -(∂r/∂p + Σ_k (∂r/∂x_k) · S_{k,p}^-) / (∂r/∂t + Σ_k (∂r/∂x_k) · f_k^-)
    #
    # Time-triggered (root NA, `time` = τ(p)):
    #   dτ/dp = ∂(τ-expr)/∂p
    #
    # All RHS references S_{j,p}^- and x_j^- are pre-event values — funC emits
    # a y_pre[] snapshot once at the top of myevent, and events emitted here
    # carry use_pre_state = TRUE so their value expressions are rewritten to
    # read from y_pre[].

    is_zero <- function(x) {
      is.null(x) || length(x) == 0L || any(is.na(x)) || !nzchar(x) || identical(x, "0")
    }
    # Whitespace-insensitive string equality. replaceSymbols() re-parses its
    # input and strips whitespace, so f_pre and f_post can differ only by
    # spaces even when semantically identical.
    same_expr <- function(a, b) {
      identical(gsub(" ", "", a, fixed = TRUE),
                gsub(" ", "", b, fixed = TRUE))
    }
    sum_terms <- function(terms) {
      terms <- terms[!vapply(terms, is_zero, logical(1))]
      if (length(terms) == 0L) return("0")
      paste(terms, collapse = " + ")
    }
    prod_term <- function(a, b) {
      if (is_zero(a) || is_zero(b)) return("0")
      if (identical(a, "1")) return(b)
      if (identical(b, "1")) return(a)
      paste0("(", a, ") * (", b, ")")
    }
    deriv1 <- function(expr, sym) {
      if (is_zero(expr) || is.null(sym) || is.na(sym) || !nzchar(sym)) return("0")
      e <- try(parse(text = expr), silent = TRUE)
      if (inherits(e, "try-error")) return("0")
      d <- try(stats::D(e[[1]], sym), silent = TRUE)
      if (inherits(d, "try-error")) return("0")
      gsub(" ", "", paste(deparse(d), collapse = ""), fixed = TRUE)
    }

    # Sensitivity parameters: state-initial-values first, then true parameters.
    # Expressions (xi, r, tau_expr) do not depend on state-initial-value
    # parameters symbolically — their contribution enters only through
    # Σ_k (∂r/∂x_k) · S_{k, p_init} in dτ/dp.
    sens_pars <- c(states, pars)
    is_init <- stats::setNames(sens_pars %in% states, sens_pars)

    events.addon <- lapply(seq_len(nrow(events)), function(i) {

      myevent  <- events[i, ]
      xone     <- as.character(myevent[["var"]])
      tau_expr <- as.character(myevent[["time"]])
      xi       <- as.character(myevent[["value"]])
      root     <- as.character(myevent[["root"]])
      method   <- as.character(myevent[["method"]])

      if (!method %in% c("replace", "add", "multiply"))
        stop("Event method must be 'replace', 'add' or 'multiply'. Got: ", method)
      time_na <- is.na(tau_expr) || !nzchar(tau_expr)
      root_na <- is.na(root)     || !nzchar(root)
      if (time_na && root_na)
        stop("Event row ", i, ": exactly one of 'time' and 'root' must be set; both are NA.")
      if (!time_na && !root_na)
        stop("Event row ", i, ": exactly one of 'time' and 'root' must be set; both are given (", tau_expr, " / ", root, ").")

      xk <- setdiff(variables, xone)

      # Pre- and post-event ODE right-hand side.
      f_pre  <- f
      xone_post <- switch(
        method,
        replace  = paste0("(", xi, ")"),
        add      = paste0("(", xone, " + (", xi, "))"),
        multiply = paste0("(", xone, " * (", xi, "))")
      )
      f_post <- replaceSymbols(xone, xone_post, f)
      names(f_post) <- names(f)

      # ∂xi/∂x_j  and  ∂xi/∂p
      dxi_dx <- stats::setNames(vapply(variables, function(j) deriv1(xi, j), character(1)),
                         variables)
      dxi_dp <- stats::setNames(vapply(sens_pars, function(p) {
        if (is_init[[p]]) "0" else deriv1(xi, p)
      }, character(1)), sens_pars)

      # Reset-map jacobians for the reset state xone (non-reset states have
      # g_k = x_k, so ∂g_k/∂x_j = δ_kj and ∂g_k/∂p = 0 — handled below).
      dg_dx <- switch(
        method,
        replace  = dxi_dx,
        add      = {
          out <- dxi_dx
          out[[xone]] <- sum_terms(c("1", dxi_dx[[xone]]))
          out
        },
        multiply = {
          out <- stats::setNames(character(length(variables)), variables)
          for (j in variables) {
            if (j == xone) {
              out[[j]] <- sum_terms(c(xi, prod_term(xone, dxi_dx[[xone]])))
            } else {
              out[[j]] <- prod_term(xone, dxi_dx[[j]])
            }
          }
          out
        }
      )
      dg_dp <- switch(
        method,
        replace  = dxi_dp,
        add      = dxi_dp,
        multiply = stats::setNames(vapply(sens_pars, function(p) prod_term(xone, dxi_dp[[p]]),
                                   character(1)), sens_pars)
      )

      # Δ_i = Σ_j (∂g_i/∂x_j) · f_j^-  -  f_i^+
      salt <- stats::setNames(character(length(variables)), variables)
      salt[[xone]] <- sum_terms(c(
        sum_terms(vapply(variables, function(j) prod_term(dg_dx[[j]], f_pre[[j]]),
                         character(1))),
        paste0("-(", f_post[[xone]], ")")
      ))
      for (i_ in xk) {
        a <- f_pre[[i_]]; b <- f_post[[i_]]
        salt[[i_]] <- if (same_expr(a, b)) "0" else sum_terms(c(a, paste0("-(", b, ")")))
      }

      # dτ/dp
      if (!is.na(root)) {
        dr_dx <- stats::setNames(vapply(variables, function(k) deriv1(root, k), character(1)),
                          variables)
        dr_dt <- deriv1(root, "time")
        denom <- sum_terms(c(
          dr_dt,
          vapply(variables, function(k) prod_term(dr_dx[[k]], f_pre[[k]]), character(1))
        ))
        dtau_dp <- stats::setNames(vapply(sens_pars, function(p) {
          dr_dp <- if (is_init[[p]]) "0" else deriv1(root, p)
          chain <- vapply(variables, function(k) {
            sv <- paste0(k, ".", p)
            if (!(sv %in% newvariables)) return("0")
            prod_term(dr_dx[[k]], sv)
          }, character(1))
          numerator <- sum_terms(c(dr_dp, chain))
          if (is_zero(numerator) || is_zero(denom)) return("0")
          paste0("-(", numerator, ") / (", denom, ")")
        }, character(1)), sens_pars)
      } else {
        dtau_dp <- stats::setNames(vapply(sens_pars, function(p) {
          if (is_init[[p]]) "0" else deriv1(tau_expr, p)
        }, character(1)), sens_pars)
      }

      # Emit one "replace" event per tracked (state i, parameter p). Records
      # with value equivalent to the variable itself (no-op) are filtered out.
      records <- list()
      for (ii in variables) {
        for (p in sens_pars) {
          sensvar <- paste0(ii, ".", p)
          if (!(sensvar %in% newvariables)) next

          if (ii == xone) {
            ja_terms <- vapply(variables, function(j) {
              sv <- paste0(j, ".", p)
              if (!(sv %in% newvariables)) return("0")
              prod_term(dg_dx[[j]], sv)
            }, character(1))
            salt_term <- prod_term(salt[[xone]], dtau_dp[[p]])
            value_expr <- sum_terms(c(ja_terms, dg_dp[[p]], salt_term))
          } else {
            salt_term <- prod_term(salt[[ii]], dtau_dp[[p]])
            if (is_zero(salt_term)) next
            value_expr <- sum_terms(c(sensvar, salt_term))
          }

          if (identical(value_expr, sensvar)) next  # no-op

          records[[length(records) + 1L]] <- data.frame(
            var           = sensvar,
            time          = tau_expr,
            value         = value_expr,
            root          = root,
            method        = "replace",
            use_pre_state = TRUE,
            stringsAsFactors = FALSE
          )
        }
      }

      if (length(records) == 0L) {
        return(data.frame(
          var           = character(0),
          time          = character(0),
          value         = character(0),
          root          = character(0),
          method        = character(0),
          use_pre_state = logical(0),
          stringsAsFactors = FALSE
        ))
      }
      do.call(rbind, records)
    })

    # Overwrite events
    events <- events.addon

    # Add rownames that allow to trace back records in eventframe to the list of events
    for (i in seq_along(events)) {
      if (nrow(events[[i]]) > 0L)
        rownames(events[[i]]) <- paste(i, seq_len(nrow(events[[i]])), sep = "_")
    }

    # Make sure character columns stay character but preserve the logical
    # use_pre_state flag (funC checks it with as.logical()).
    eventframe <- do.call(rbind, events)
    if (!is.null(eventframe) && nrow(eventframe) > 0L) {
      for (k in seq_along(eventframe)) {
        col <- eventframe[[k]]
        if (!is.logical(col)) eventframe[[k]] <- as.character(col)
      }
    }

    # Reduce events (by name matching with eventframe)
    events <- lapply(events, function(e) {
      if (nrow(e) == 0L) return(e)
      e[rownames(e) %in% rownames(eventframe), , drop = FALSE]
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
