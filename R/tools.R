#' D.character
#' @export
D.character <- function(what, by) paste(deparse(D(parse(text = what), by)), collapse="")


getBlocks <- function(subsample, target) {
  
  
  
  mysub <- subset(subsample, Target == target)
  N <- dim(mysub)[1]
  Experiments <- unique(mysub$Experiment)
  Conditions <- unique(paste(mysub$min, mysub$Condition))
  
  M <- matrix(0, ncol=length(Experiments), nrow=length(Conditions))
  
  for(i in 1:N) {
    myrow <- which(Conditions == paste(mysub$min, mysub$Condition)[i])
    mycol <- which(Experiments == mysub$Experiment[i])
    M[myrow, mycol] <- 1
  }
  
  return(analyzeBlocks(M))
  
  
}

## Analyze matrix with zero and one entries for connection components
analyzeBlocks <- function(M) {
  out <- which(apply(M, 1, sum)==0)
  if(length(out) > 0) {
    M <- M[-out,]
    cat("matrix contains zero rows which have been eliminated\n")
  }
  n <- dim(M)[1]
  rcomponents <- list()
  ccomponents <- list()
  
  counter <- 0
  while(length(unlist(rcomponents))< n) {
    
    counter <- counter + 1
    
    if(length(unlist(rcomponents)) == 0) {
      w <- 1
    } else {
      mysample <- (1:n)[-unlist(rcomponents)]
      w <- min(mysample)
    }
    
    repeat {  
      
      v <- unique(rapply(as.list(w), function(i) which(M[i,] == 1)))
      wnew <- unique(rapply(as.list(v), function(j) which(M[,j] == 1)))
      if(length(wnew) == length(w)) break
      w <- wnew
    }
    rcomponents[[counter]] <- w
    ccomponents[[counter]] <- v
    
    
  }
  
  #print(M[unlist(rcomponents), unlist(ccomponents)])
  
  
  return(rcomponents)
  
}



#' reduceSensitivities
#' @export
#' @examples
#' f <- c(C = "-1*k3*C +1*k2*B", B = "-1*k2*B +1*k1*A*A", A = "-2*k1*A*A +1*buildA")
#' mysens <- sensitivitiesSymb(f)
#' red <- reduceSensitivities(mysens, attr(mysens, "is.zero.Dpf"))
reduceSensitivities <- function(sens, is.zero.Dpf) {
  
  
  sensvar <- names(sens)
    
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) v[2]))
  
  ini.zero <- which(senssplit.1 != senssplit.2)
  
  ## Neuer Ansatz (auskommentiert)
  
#   sens.jac <- jacobianSymb(sens)
#   sens.jac[sens.jac!="0"] <- "1"
#   sens.jac <- matrix(as.numeric(sens.jac), length(sens), length(sens))
#   
#   depends.on <- apply(sens.jac, 2, function(v) which(v == 1))
#   repeat {
#     new.depends.on <- lapply(depends.on, function(v) unique(unlist(depends.on[v])))
#     new.diff <- unlist(lapply(1:length(depends.on), function(i) setdiff(depends.on[[i]], new.depends.on[[i]])))      
#     if(length(new.diff) == 0) break 
#     depends.on <- new.depends.on
#   }
#   
#   lapply(depends.on, function(v) {
#     
#   })
  
  exit <- FALSE
  while(!exit) {
  
    # Check if sens depends only on itself
    is.sep <- unlist(lapply(sensvar[ini.zero], function(s) {
      condition1 <- grepl(s, sens[s])
      condition2 <- s %in% is.zero.Dpf
      condition3 <- !any(sensvar[sensvar != s]%in%getSymbols(sens[s]))
      #condition3 <- !any(sapply(sensvar[sensvar != s], function(s2) grepl(s2, sens[s])))
      condition1 & condition2 & condition3
    }))
    
    if(any(is.sep)) {    
      sens[ini.zero][is.sep] <- "0"
      sens <- replaceSymbols(sensvar[ini.zero][is.sep], "0", sens)
    } else {
      exit <- TRUE
    }
    
  }
  
  return(sens)
  
}