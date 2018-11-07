context("sensitivitiesSymb")
test_that("sensitivitiesSymb_works_for_all_possible_combinations_of_states_and_events", {
  
  #-!Start example code
  # library(dMod)

  # reactions <- NULL
  # reactions <- dMod::addReaction(reactions, "A", "B", "kcat*INPUT*A/(Km+A)", 1)
  # reactions <- dMod::addReaction(reactions, "B", "", "-k*B", 2)
  # reactions <- dMod::addReaction(reactions, "", "INPUT", "0", 3)
  # f <- unclass(dMod::as.eqnvec(reactions))
  f <- c(A = "-1*(kcat*INPUT*A/(Km+A))", B = "1*(kcat*INPUT*A/(Km+A)) -1*(-k*B)", INPUT = "1*(0)")
  
  events <- expand.grid(var = "INPUT", 
                        time = c(1, "ton"), 
                        value = c(1, "xon"), 
                        method = c("replace", "add" # , "multiply" # multiply not yet implemented
                        ), stringsAsFactors = F)
  states <- list(character(0), "DUMMY1", c("A"), c("A", "B"), c("A", "INPUT"), c("A", "B", "INPUT"), "INPUT")
  parameters <- list(character(0), "DUMMY2", c("k"), c("k", "kon"), c("k", "ton"), "kon", "ton")
  
  
  out <- lapply(seq_len(nrow(events)), function(event_row__) {
    events__ <- events[event_row__,]
    lapply(states, function(states__) {
      lapply(parameters, function(parameters__) {
        if(length(parameters__) == 0 & length(states__) == 0)
          return()
        # cat(paste0("event: ", event_row__, "\t", "states: ", paste0(states__, collapse = ", ") , "\t", "parameters: ", paste0(parameters__, collapse = ", ") , "\n"))
        sensitivitiesSymb(f, states = states__, parameters = parameters__, events = events__, reduce = T)
      })
    })
  })
  #-!End example code
  

  
  
  # myevent <- events[2,]
  # mystates <- states[[1]]
  # myparameters <- "ton"
  # # debugonce(sensitivitiesSymb)
  # sensitivitiesSymb(f, mystates, myparameters, events = myevent, reduce = T)
  
  parsed <- try(lapply(unlist(out, recursive = T), function(i) parse(text = i)))

  
  
  # Define your expectations here
  expect_false(inherits(parsed, "try-error"))

})