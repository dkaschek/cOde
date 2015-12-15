#' Compute sensitivity equations of a function symbolically
#' 
#' @param f named vector of type character, the functions
#' @param states Character vector. Sensitivities are computed with respect to initial
#' values of these states
#' @param parameters Character vector. Sensitivities are computed with respect to initial
#' values of these parameters
#' @param inputs Character vector. Input functions or forcings. They are excluded from
#' the computation of sensitivities.
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
#' 
#' 

library(cOde)
library(dMod)

f <- c(
  A = "k1 * B * C - k2 * A",
  B = "-k1 * B * C + k2 * A",
  C  = "-k1 * B * C + k2 * A"
)

f <- c(
  A = "k1 * B  - k2 * A^2",
  B = "-k1 * B + k2 * A^2"
)


#zum Testen von Reduce Sensitivities
f <- c(
  A = "k1 * B  - k2 * A^2",
  B = "-k1 * B "
)


states = names(f)
parameters = NULL
inputs = NULL
reduce = TRUE

expand.grid.alt <- function(seq1,seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}




# secondsensitivitiesSymb <- function(f, states = names(f), parameters = NULL, inputs = NULL, reduce = FALSE) {
  variables <- names(f)
  states <- states[!states%in%inputs]
  
  if(is.null(parameters)) {
    pars <- getSymbols(f, exclude=c(variables, inputs, "time"))
  } else {
    pars <- parameters[!parameters%in%inputs]
  }
  
  Dyf <- jacobianSymb(f, variables)
  Dpf <- jacobianSymb(f, pars)
  
  df <- length(f)
  dv <- length(variables)
  ds <- length(states)
  dp <- length(pars)
  
  # generate sensitivity variable names and names with zero entries
  mygridY0 <- expand.grid.alt(variables, states)
  mygridP <- expand.grid.alt(variables, pars)
  sensParVariablesY0 <- apply(mygridY0, 1, paste, collapse = ".")
  sensParVariablesP <- apply(mygridP, 1, paste, collapse = ".")
  
  # Write sensitivity equations in matrix form
  Dy0y <- matrix(sensParVariablesY0, ncol = ds, nrow = dv)
  Dpy <- matrix(sensParVariablesP, ncol = dp, nrow = dv)
  
  
  gl <- c(as.vector(prodSymb(matrix(Dyf, ncol = dv), Dy0y)), 
          as.vector(sumSymb(prodSymb(matrix(Dyf, ncol = dv), Dpy), matrix(Dpf, nrow = dv))))
  
  newfun <- gl
  newvariables.grid <- expand.grid.alt(variables, c(states, pars))
  newvariables <- apply(newvariables.grid, 1, paste, collapse=".")
  names(newfun) <- newvariables
  
  # Reduce the sensitivities
  vanishing <- c(sensParVariablesY0, sensParVariablesP[Dpf == "0"])
  if(reduce) {
    newfun <- reduceSensitivities(newfun, vanishing)
    is.zero.sens <- names(newfun) %in% attr(newfun,"is.zero")
  } else {
    is.zero.sens <- rep(FALSE, length(newfun))
  }
  newfun <- newfun[!is.zero.sens]
  
  
  # Append initial values
  initials <- rep(0, length(newfun))
  names(initials) <- newvariables[!is.zero.sens]
  ones <- which(apply(newvariables.grid, 1, function(row) row[1] == row[2]))
  initials[newvariables[ones]] <- 1
  
  
  #Second Sensitivies of ode with respect to parameters
  # g <- newfun[(dv*ds+1):length(newfun)] #first sensitivities

  # calculate the three summands of the 2nd degree sensitivities and add them
  
  Dyg <- jacobianSymb(newfun[(sum(!is.zero.sens[1:(dv*ds)], na.rm=TRUE)+1):length(newfun)], variables) #hier muss man evtl vanishing sensitivites noch beachten!!!
  SecondSensSum1 <- as.vector(prodSymb(matrix(Dyg,ncol=dv),Dpy)) 
  
  mygridPP<-expand.grid.alt(sensParVariablesP,pars)
  sensParVariablesPP <- apply(mygridPP, 1, paste, collapse=".")
  SecondSensSum2 <- matrix(prodSymb(matrix(Dyf,ncol=dv),matrix(sensParVariablesPP, nrow = dv)), ncol=dp) 
  SecondSensSum2 <- as.vector(SecondSensSum2[!is.zero.sens[(dv*ds+1):length(is.zero.sens)],]) # to consider reduced symmetries
  
  SecondSensSum3 <- jacobianSymb(newfun[(sum(!is.zero.sens[1:(dv*ds)], na.rm=TRUE)+1):length(newfun)], pars) 
  
  SecondSens <- as.vector(sumSymb(sumSymb(SecondSensSum1,SecondSensSum2),SecondSensSum3))
  names(SecondSens) <- names(SecondSensSum3)  

  # consider symmetries of derivative indices ie d^2x_i/(dp_j dp_k) = d^2x_i/(dp_k dp_j)
  # der Vektor secondSens sieht folgendermaßen aus: (ijk), wobei erst i durchgefahren wird, dann j, dann k
  # auffassen als vektor von matrizen jk, in denen nur die diagonalelemente genommen werden
  
  secondsenssplit <-strsplit(names(SecondSens), ".", fixed= TRUE)
  secondsenssplit.1<-unlist(lapply(secondsenssplit, function(v) v[1]))
  secondsenssplit.2<-unlist(lapply(secondsenssplit, function(v) v[2]))
  secondsenssplit.3<-unlist(lapply(secondsenssplit, function(v) v[3]))
  
  permutedjk<-paste(secondsenssplit.1,secondsenssplit.3,secondsenssplit.2, sep=".") #vertausche j und k
  pairs<-unlist(lapply(permutedjk, function(v) match(v, names(SecondSens)))) #finde j<->k-Paare und eventuelle reduced Sensitivities
  
  unique.inz<-unlist(lapply((1:length(pairs)), function(s) {
    ifelse(s<=pairs[s], s, NA)
  })) # nimmt von jeder j-k-Permutation nur diejenige, bei der j<=k ist und setzt alle anderen gleich NA,
  # alle Paare, von denen eine Permutation schon wegen reduced Symmetries geflogen ist, werden auch gleich NA gesetzt 
  
  # if a combination like "yi.pj.pj" is zero but it still comes up as a symbol in other second sensitivities, set it to zero
  secondZerosSame <- strsplit(sensParVariablesPP, ".", fixed= TRUE)
  secondZerosSame <- secondZerosSame[unlist(lapply(secondZerosSame, function(v) v[2]==v[3]))]
  secondZerosSame <- unlist(lapply(secondZerosSame, function(v) paste(v[1],v[2],v[3],sep=".")))
  secondZerosSame <- secondZerosSame[!(secondZerosSame %in% names(SecondSens))]
  
  # change symbolic derivatives like "yi.pj.pk" to "yi.pk.pj" if the first index combination is being eliminated due to symmetry

  origSecondNames<-c(names(SecondSens[is.na(unique.inz)]),permutedjk[is.na(pairs)],secondZerosSame)
  
  secondZerosNames<-names(SecondSens[pairs[is.na(unique.inz)]])
  secondZerosNames[is.na(secondZerosNames)]<-0
  secondZerosNames<-c(secondZerosNames, rep_len("0", length.out = sum(is.na(pairs), length(secondZerosSame))))
  
  
  replaceSymbols(origSecondNames, secondZerosNames ,SecondSens)
  
  
  unique.inz<-unique.inz[!is.na(unique.inz)] 
  
  SecondSens<-SecondSens[unique.inz]
  
  #initials for second sensitivities (sind natürlich auch null)
  Secondyini<-rep(0,length(SecondSens))
  names(Secondyini) <- names(SecondSens)
  


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
  attr(newfun, "outputs") <- structure(rep(0, length(which(is.zero.sens))), names = newvariables[is.zero.sens])
  attr(newfun, "forcings") <- c(statesD, weightsD)
  attr(newfun, "yini") <- initials
  attr(newfun, "secondsens") <-SecondSens
  attr(newfun, "secondyini") <-Secondyini
  
  return(newfun)
  
}


# Christoffelsymbole aus den zweiten Sensitivitäten berechnen

# get parameters and variables
getParameters <- function(f, states = names(f), parameters = NULL, inputs = NULL, reduce = FALSE) {
  # to get the parameters of the sensitivity equations
  variables <- names(f)
  states <- states[!states %in% inputs]
  
  if (is.null(parameters)) {
    pars <- getSymbols(f, exclude = c(variables, inputs, "time"))
  } else {
    pars <- parameters[!parameters %in% inputs]
  }
  
  return(pars)
}

pars<-getParameters(f)
variables <- names(f)

# if you have some predicition, extract the columns of the prediction with the second derivatives
mygridP <- expand.grid.alt(variables, pars)
sensParVariablesP <- apply(mygridP, 1, paste, collapse = ".")
mygridPP<-expand.grid.alt(sensParVariablesP,pars)
sensParVariablesPP <- apply(mygridPP, 1, paste, collapse=".")

firstDeriv<-prediction[,sensParVariablesP]
secondDeriv<- prediction[,sensParVariablesPP]

# Was bringen mir die zweiten Ableitungen der Variablen? Ich bräuchte eher die zweiten Ableitungen der Residuen.
# Wie stehen diese im Zusammenhang? <- Die einzelnen Variablen sind Indizes eines einzelnen Residuums. 
# Das heißt, der Residuenvektor (des Papers) besteht aus Residuen, die selbst Vektoren sind.


