
IMScore <- function(matrices, scores, penalty) {
  

  
  #print(typeof(length(scores)))
  m <- length(scores)
  #print(m)
  #print(dim(matrices[[1]]))
  #n <- dim(matrices[[1]])[[2]]
  #print(typeof(n))
  #print(n)
  sum <- 0
  k <- 5 #filler value
  
  for (i in 1:length(scores)) {
    #sum <- sum + scores[[i]]$global.score(scores[[i]]$create.dag()) + ((penalty * k) * log(n))
    sum <- sum + scores[[i]]$global.score(scores[[i]]$create.dag())
  }
  
  #print(sum)
  imscore = ((-2/m) *  sum )+ ((penalty * k) * log(n))
  print(imscore)
  return(imscore)
}



call_forward <- function(
  score, 
  labels = score$getNodes(), 
  targets = score$getTargets(),
  fixedGaps = NULL, 
  adaptive = c("none", "vstructures", "triples"), 
  phase = c("forward"),
  iterate = length(phase) > 1,
  turning = NULL, 
  maxDegree = integer(0),
  verbose = FALSE, 
  ...) {
  
  #' Creates a list of options for the C++ function "causalInference";
  #' internal function
  causal.inf.options = function(
    caching = TRUE,
    phase = c("forward", "backward", "turning"),
    iterate = length(phase) > 1,
    maxDegree = integer(0),
    maxSteps = 0,
    childrenOnly = integer(0),
    fixedGaps = NULL,
    adaptive = c("none", "vstructures", "triples"),
    verbose = 0,
    p = 0) {
    # Check for deprecated calling convention and issue a warning
    if (p > 0) {
      warning(paste("Argument 'p' is deprecated in calls of ges() or gies",
                    "and will be disabled in future package versions;",
                    "please refer to the corresponding help page.", sep = " "))
    }
    
    # Error checks for supplied arguments
    # TODO extend!
    if (is.logical(adaptive)) {
      adaptive <- ifelse(adaptive, "vstructures", "none")
      warning(paste("The parameter 'adaptive' should not be provided as logical anymore;",
                    "cf. ?ges or gies", sep = " "))
    }
    phase <- match.arg(phase, several.ok = TRUE)
    stopifnot(is.logical(iterate))
    adaptive <- match.arg(adaptive)
    if (is.null(fixedGaps)) {
      adaptive <- "none"
    }
    list(caching = caching,
         phase = phase,
         iterate = iterate,
         maxDegree = maxDegree,
         maxSteps = maxSteps,
         childrenOnly = childrenOnly,
         fixedGaps = fixedGaps,
         adaptive = adaptive,
         DEBUG.LEVEL = as.integer(verbose))
  }
  
  in_edges = replicate(length(labels), integer(0), simplify = FALSE)
  alg.name = "GIES-F"
  
  new.graph <- .Call("causalInference",
                     in_edges,
                     score$pp.dat,
                     alg.name,
                     score$c.fcn,
                     causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, phase=phase ...),
                     PACKAGE = "imagestest")
  
  return(new.graph)
  #return(list("graph" = new.graph, "in_edges" = new.graph$in.edges))
}

call_backward <- function(
  score,
  forward,
  labels = score$getNodes(), 
  targets = score$getTargets(),
  fixedGaps = NULL, 
  adaptive = c("none", "vstructures", "triples"), 
  phase = c("forward"),
  iterate = length(phase) > 1,
  turning = NULL, 
  maxDegree = integer(0),
  verbose = FALSE, 
  ...) {
  
  #' Creates a list of options for the C++ function "causalInference";
  #' internal function
  causal.inf.options = function(
    caching = TRUE,
    phase = c("forward", "backward", "turning"),
    iterate = length(phase) > 1,
    maxDegree = integer(0),
    maxSteps = 0,
    childrenOnly = integer(0),
    fixedGaps = NULL,
    adaptive = c("none", "vstructures", "triples"),
    verbose = 0,
    p = 0) {
    # Check for deprecated calling convention and issue a warning
    if (p > 0) {
      warning(paste("Argument 'p' is deprecated in calls of ges() or gies",
                    "and will be disabled in future package versions;",
                    "please refer to the corresponding help page.", sep = " "))
    }
    
    # Error checks for supplied arguments
    # TODO extend!
    if (is.logical(adaptive)) {
      adaptive <- ifelse(adaptive, "vstructures", "none")
      warning(paste("The parameter 'adaptive' should not be provided as logical anymore;",
                    "cf. ?ges or gies", sep = " "))
    }
    phase <- match.arg(phase, several.ok = TRUE)
    stopifnot(is.logical(iterate))
    adaptive <- match.arg(adaptive)
    if (is.null(fixedGaps)) {
      adaptive <- "none"
    }
    list(caching = caching,
         phase = phase,
         iterate = iterate,
         maxDegree = maxDegree,
         maxSteps = maxSteps,
         childrenOnly = childrenOnly,
         fixedGaps = fixedGaps,
         adaptive = adaptive,
         DEBUG.LEVEL = as.integer(verbose))
  }
  
  in_edges = forward$in.edges
  alg.name = "GIES-B"
  
  new.graph <- .Call("causalInference",
                     in_edges,
                     score$pp.dat,
                     alg.name,
                     score$c.fcn,
                     causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, ...),
                     PACKAGE = "imagestest")
  
  return(new.graph)
}

call_turn <- function(
  score,
  backward,
  labels = score$getNodes(), 
  targets = score$getTargets(),
  fixedGaps = NULL, 
  adaptive = c("none", "vstructures", "triples"), 
  phase = c("forward"),
  iterate = length(phase) > 1,
  turning = NULL, 
  maxDegree = integer(0),
  verbose = FALSE, 
  ...) {
  
  #' Creates a list of options for the C++ function "causalInference";
  #' internal function
  causal.inf.options = function(
    caching = TRUE,
    phase = c("forward", "backward", "turning"),
    iterate = length(phase) > 1,
    maxDegree = integer(0),
    maxSteps = 0,
    childrenOnly = integer(0),
    fixedGaps = NULL,
    adaptive = c("none", "vstructures", "triples"),
    verbose = 0,
    p = 0) {
    # Check for deprecated calling convention and issue a warning
    if (p > 0) {
      warning(paste("Argument 'p' is deprecated in calls of ges() or gies",
                    "and will be disabled in future package versions;",
                    "please refer to the corresponding help page.", sep = " "))
    }
    
    # Error checks for supplied arguments
    # TODO extend!
    if (is.logical(adaptive)) {
      adaptive <- ifelse(adaptive, "vstructures", "none")
      warning(paste("The parameter 'adaptive' should not be provided as logical anymore;",
                    "cf. ?ges or gies", sep = " "))
    }
    phase <- match.arg(phase, several.ok = TRUE)
    stopifnot(is.logical(iterate))
    adaptive <- match.arg(adaptive)
    if (is.null(fixedGaps)) {
      adaptive <- "none"
    }
    list(caching = caching,
         phase = phase,
         iterate = iterate,
         maxDegree = maxDegree,
         maxSteps = maxSteps,
         childrenOnly = childrenOnly,
         fixedGaps = fixedGaps,
         adaptive = adaptive,
         DEBUG.LEVEL = as.integer(verbose))
  }
  
  in_edges = backward$in.edges
  alg.name = "GIES-T"
  
  new.graph <- .Call("causalInference",
                     in_edges,
                     score$pp.dat,
                     alg.name,
                     score$c.fcn,
                     causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, ...),
                     PACKAGE = "imagestest")
  
  return(new.graph)
}

update_globals <- function(results) {
  
  scores = list()
  
  for (i in 1:length(results)) {
    scores[[i]] <- new("GaussL0penObsScore", results[[i]])#$in.edge)s
  }
  
}


IMaGES <- function(matrices, scores, penalty = 1.5) {
  
  scores = list()
  
  for (i in 1:length(matrices)) {
    scores[[i]] <- new("GaussL0penObsScore", matrices[[i]])
    print(scores[[i]]$global.score(scores[[i]]$create.dag()))
    #scores[[i]] <- ges(scr)
    #scores[[i]] <- scr
  }
  
  #print(scores[[1]].format)
  
  imscore = IMScore(matrices, scores, penalty)
  
  in_edges = list()
  
  # for (i in 1:length(scores)) {
  #   in_edges[[i]] <- scores[[i]]$getNodes()
  #   
  # }
  
  forwards <- list()
  
  for (i in 1:length(scores)) {
    #forwards[[i]] <- ges(scores[[i]])#, phase=c("forward"))
    forwards[[i]] <- call_forward(scores[[i]])
    #res <- call_forward(scores[[i]])
    #forwards[[i]] <- res$graph
    #in_edges[[i]] <- res$in_edges
  }
  
  #print(forwards[[1]])

  #return(forwards)
  
  backwards <- list()
  
  for (i in 1:length(forwards)) {
    backwards[[i]] <- call_backward(scores[[i]], forwards[[i]])  
    
  }
  
  #return(backwards)
  
  turning <- list()
  
  for (i in 1:length(backwards)) {
    turning[[i]] <- call_turn(scores[[i]], backwards[[i]])
    
  }
  
  return(turning)
  
}


