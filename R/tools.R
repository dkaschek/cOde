#' D.character
#' @export
D.character <- function(what, by) paste(deparse(D(parse(text = what), by)), collapse="")

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
  ini.nonzero <- which(senssplit.1 == senssplit.2) 
  
  exit <- FALSE
  sensvar.zero <- sensvar[ini.zero]
  sensvar.nonzero <- sensvar[ini.nonzero]
  while(!exit){
    find_nonzero <- unlist(lapply(sensvar.zero, function(s){
      allSyb <- sensvar[sensvar %in%getSymbols(sens[s])]
      nDpf <- setdiff(allSyb,is.zero.Dpf)
      nIni <- intersect(sensvar.nonzero,allSyb)
      nonzero <- (length(nDpf)+length(nIni) > 0)
    }))
    if(any(find_nonzero)){
      sensvar.nonzero <- c(sensvar.nonzero,sensvar.zero[find_nonzero])
      sensvar.zero <- sensvar.zero[!find_nonzero]
    }else{
      exit <- TRUE
    }
  }
  sens[sensvar.zero] <- "0"
  sens <- replaceSymbols(sensvar.zero, "0", sens)
  attr(sens,"is.zero") <- sensvar.zero
  #return(list(sens, sensvar.zero,sensvar.nonzero))
  return(sens)
}
