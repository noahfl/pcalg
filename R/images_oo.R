IMaGES <- setRefClass("IMaGES",
  fields = list(matrices="list", penalty="numeric", .rawscores="list", .graphs="list", imscore = "numeric", results="list"),
  
  # validity <- function(object) {
  #   return(TRUE)
  # },
  
  
  #TODO: figure out why this isn't being called
  methods = list(
    create.graph = function(
      score, 
      labels = score$getNodes(), 
      targets = score$getTargets(),
      fixedGaps = NULL, 
      #adaptive = c("none", "vstructures", "triples"), 
      #phase = c("forward", "backward", "turning"),
      #iterate = length(phase) > 1,
      turning = NULL, 
      maxDegree = integer(0),
      verbose = FALSE, 
      ...)
    {
      # Catch calling convention of previous package versions:
      # ges(p, targets, score, fixedGaps = NULL, ...)
      # If this calling convention is used, issue a warning, but adjust the 
      # arguments
      if (is.numeric(score) && is.list(labels) && inherits(targets, "Score")) {
        score <- targets
        targets <- labels
        labels <- as.character(1:length(score$getNodes()))
        warning(paste("You are using a deprecated calling convention for gies()",
                      "which will be disabled in future versions of the package;",
                      "cf. ?gies.", sep = " "))
      }
      # If the old calling convention was used with named arguments, "p = ..."
      # would assign a numerical value to "phase" (expanding arguments...)
      # if (is.numeric(phase)) {
      #   phase <- c("forward", "backward", "turning")
      #   warning(paste("You are using a deprecated calling convention for gies()",
      #                 "which will be disabled in future versions of the package;",
      #                 "cf. ?gies.", sep = " "))
      # }
      
      # Issue warning if argument 'turning' was used
      if (!missing(turning)) {
        stopifnot(is.logical(turning))
        warning(paste0("The argument 'turning' is deprecated; please use 'phase'",
                       "instead (cf. ?ges)"))
        
      #   if (turning) {
      #     phase <- c("forward", "backward", "turning")
      #     iterate <- FALSE
      #   } else {
      #     phase <- c("forward", "backward")
      #     iterate <- FALSE
      #   }
      }
      
      # Error checks
      if (!inherits(score, "Score")) {
        stop("Argument 'score' must be an instance of a class inherited from 'Score'.")
      }
      #phase <- match.arg(phase, several.ok = TRUE)
      # TODO extend...
    
      # Catching error occurring when a user called one of the causal 
      # inference algorithms using the old calling conventions: try to
      # rearrange passed arguments, print a warning
      #
      # NOTE: old calling conventions were
      # (algorithm, p, targets, score) for caus.inf
      # (p, targets, score) for all functions allowing interventional data
      # (p, score) for GES
      if (is.numeric(score)) {
        # This happens when the old calling convention is used with all 
        # mandatory arguments unnamed
        p <- score
        if (is.list(labels) && is(targets, "Score")) {
          score <- targets
          targets <- labels
          labels <- as.character(1:p)
          warning(paste("You are using a DEPRECATED calling convention for",
                        "gies(), gds() or simy(); please refer to the documentation",
                        "of these functions to adapt to the new calling conventions."))
        } else if (is(labels, "Score")) {
          score <- labels
          labels <- as.character(1:p)
          warning(paste("You are using a DEPRECATED calling convention for",
                        "ges(); please refer to the documentation",
                        "to adapt to the new calling convention."))
        }
      } else if (is.numeric(labels) && length(labels) == 1) {
        # This happens when the old calling convention is used with only the
        # 'score' argument named
        labels <- as.character(1:labels)
        warning(paste("You are using a DEPRECATED calling convention for",
                      "gies(), ges(), gds() or simy(); please refer to the documentation",
                      "of these functions to adapt to the new calling conventions."))
      }
      
      if (!is(score, "Score")) {
        stop("'score' must be of a class inherited from the class 'Score'.")
      }
      if (!is.character(labels)) {
        stop("'labels' must be a character vector.")
      }
      if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
        stop("'targets' must be a list of integer vectors.")
      }
      
      
      #print("You made it this far")
      imgraph <- new("IMGraph", nodes = labels, targets = targets, score = score)
      return(imgraph)
      # if (essgraph$caus.inf(algorithm, ...)) {
      #   if (algorithm == "GIES") {
      #     ## GIES yields an essential graph; calculate a representative thereof
      #     list(essgraph = essgraph, repr = essgraph$repr())
      #   } else {
      #     ## GDS and SiMy yield a DAG; calculate the corresponding essential graph,
      #     ## although calculations may come from a model class where Markov equivalence
      #     ## does not hold!
      #     list(essgraph = dag2essgraph(essgraph$repr(), targets = targets),
      #          repr = essgraph$repr())
      #   }
      # } else stop("invalid 'algorithm' or \"EssGraph\" object")
    },
    
    run_phase = function(phase="forward", j) {
      
      alg.name <- ""
      
      #phase <- match.arg(phase)
      
      if (phase == "forward") {
        alg.name <- "GIES-F"
      }
      else if (phase == "backward") {
        alg.name <- "GIES-B"
      }
      else if (phase == "turning") {
        alg.name <- "GIES-T"
      }
      else {
        stop("incorrect algorithm name")
      }
      
      
      #for (i in 1:length(.graphs)) {
      if (!.graphs[[j]]$greedy.step(alg.name=alg.name, direction = phase, verbose = FALSE)) {
        #stop("something happened")
        #print("SOMETHING HAPPENED. probably that number thing")
        print(.graphs[[j]]$.nodes)
      }
      else {
        #print("good")
        print(.graphs[[j]]$.nodes)
      }
      #}
      
    },
    
    update_score = function() {
      
      imscore <<- IMScore()
      #imscore <<- IMScore()
      
      print("UPDATED")
      print(imscore)
      
       for (i in 1:length(.graphs)) {
         #print("INITIAL")
         #print(.graphs[[i]]$.score$.imscore)
         .graphs[[i]]$.score$.imscore = imscore
         
         #print("AFTER")
         #print(.graphs[[i]]$.score$.imscore)
         .graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=imscore)
       }
    },
    
    
    run = function() {
      update_score()
      phases = list("forward", "backward", "turning")
      
      for (i in 1:length(phases)) {
        for (j in 1:length(.graphs)) {
          #print(toString(j))
          run_phase(phases[[i]], j)
          update_score()

        }
      }
    print("HERE?")
    },
    
    #TODO: address snowball effect for IMScore
    IMScore = function() {
      
      #graphs <- .graphs
      #penalty <- get('penalty', envir=globalenv())
      penalty
      matrices <- get('matrices', envir=globalenv())
      
      #print(typeof(length(scores)))
      m <- length(.graphs)
      #print(m)
      #print(dim(matrices[[1]]))
      n <- dim(matrices[[1]])[[2]]
      #print("N")
      #print(n)
      #print(typeof(n))
      #print(n)
      sum <- 0
      k <- .graphs[[1]]$.score$pp.dat$lambda
      #print("VARS")
      #print(list(penalty, m, n, sum, k))
      
      for (i in 1:length(scores)) {
        #sum <- sum + scores[[i]]$global.score(scores[[i]]$create.dag()) + ((penalty * k) * log(n))
        # print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=.graphs[[i]]$.score$.imscore))
        # sum <- sum + .graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=.graphs[[i]]$.score$.imscore)
        
        print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=NULL))
        sum <- sum + .graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=NULL)
      }
      
      #print("SUM")
      #print(sum)
      #print(sum)
      imscore = ((-2/m) *  sum) + ((penalty * k) * log(n))
      #imscore <<- imscore
      #print(paste("imscore in function: ", imscore, "\n"))
      print("IMSCORE:")
      print(imscore)
      return(imscore)
    }
  )
)

IMaGES$methods(
  initialize = function(matrices = NULL, penalty = 1.5, imscore = NULL) {
    #images <-
    #print("initializing")
    #imscore = 0
    rawscores <- list()
    for (i in 1:length(matrices)) {
      #print("adding score")
      rawscores[[i]] <- new("GaussL0penObsScore", matrices[[i]])
    }
    penalty <<- penalty
    .rawscores <<- rawscores
    
    graphs <- list()
    
    for (i in 1:length(.rawscores)) {
      #print("creating graph")
      graphs[[i]] <- create.graph(rawscores[[i]])
    }
    
    .graphs <<- graphs
    
    imscore <<- IMScore()
    #imscore <<- IMScore()
    
    print("---------------")
    
    run()
    
    print("---------------")
    
    #imscore <<- IMScore()
    
    results = list()
    
    imscore <- IMScore()
    print("FINAL SCORES")
    for (i in 1:length(.graphs)) {
      print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=imscore))
      #print(list(.graphs[[i]], .graphs[[i]]$repr()))
      results[[i]] = list(.graphs[[i]], .graphs[[i]]$repr())
    }
    
    results <<- results
    
  }
)

#' Interventional essential graph for IMaGES
setRefClass("IMGraph",
            fields = list(
              .nodes = "vector",
              .in.edges = "list",
              .targets = "list",
              .score = "Score"
            ),
            
            validity = function(object) {
              ## Check in-edges
              if (!all(sapply(object$.in.edges, is.numeric))) {
                return("The vectors in 'in.edges' must contain numbers.")
              }
              if (!all(unique(unlist(object$.in.edges)) %in% 1:object$node.count())) {
                return(sprintf("Invalid edge source(s): edge sources must be in the range 1:%d.",
                               object$node.count()))
              }
              
              ## Check targets
              if (anyDuplicated(object$.targets)) {
                return("Targets are not unique.")
              }
              if (!all(unique(unlist(object$.targets)) %in% 1:object$node.count())) {
                return(sprintf("Invalid target(s): targets must be in the range 1:%d.",
                               object$node.count()))
              }
              
              return(TRUE)
            },
            
            methods = list(
              #' Constructor
              initialize = function(nodes,
                                    in.edges = replicate(length(nodes), integer(0), simplify = FALSE),
                                    targets = list(integer(0)),
                                    score = NULL) {
                ## Store nodes names
                if (missing(nodes)) {
                  stop("Argument 'nodes' must be specified.")
                }
                .nodes <<- as.character(nodes)
                
                ## Store in-edges
                stopifnot(is.list(in.edges) && length(in.edges) == length(nodes))
                # More error checking is done in validity check
                .in.edges <<- in.edges
                names(.in.edges) <<- NULL
                
                ## Store targets
                setTargets(targets)
                
                ## Store score
                setScore(score)
              },
              
              #' Yields the number of nodes
              node.count = function() {
                length(.nodes)
              },
              
              #' Yields the total number of edges in the graph
              edge.count = function() {
                sum(vapply(.in.edges, length, 1L))
              },
              
              #' Getter and setter functions for score object
              getScore = function() {
                .score
              },
              
              setScore = function(score) {
                if (!is.null(score)) {
                  .score <<- score
                }
              },
              
              #' Getter and setter functions for targets list
              getTargets = function() {
                .targets
              },
              
              setTargets = function(targets) {
                .targets <<- lapply(targets, sort)
              },
              
              #' Creates a list of options for the C++ function "causalInference";
              #' internal function
              causal.inf.options = function(
                caching = TRUE,
                phase = c("forward"),
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
                #phase <- match.arg(phase, several.ok = TRUE)
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
              },
              
              #' Performs one greedy step
              greedy.step = function(alg.name="GIES-F", direction = c("forward"), verbose = FALSE, ...) {
                stopifnot(!is.null(score <- getScore()))
                
                ## Cast direction
                #direction <- match.arg(direction)
                #alg.name <- switch(direction,
                #                   forward = "GIES-F",
                #                   backward = "GIES-B",
                #                   turning = "GIES-T")
                #alg.name <- match.arg(alg.name)
                print("calling GES")
                new.graph <- .Call("causalInference",
                                   .in.edges,
                                   score$pp.dat,
                                   alg.name,
                                   score$c.fcn,
                                   causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, phase=direction, ...),
                                   PACKAGE = "imagestest")
                #if (identical(new.graph, "interrupt"))
                #  return(FALSE)
                
                # if (new.graph$steps > 0) {
                #   .in.edges <<- new.graph$in.edges
                #   names(.in.edges) <<- .nodes
                # }
                
                print("BEFORE")
                print(.in.edges)
                .in.edges <<- new.graph$in.edges
                print("AFTER")
                print(.in.edges)
                
                names(.in.edges) <<- .nodes
                
                return(new.graph$steps == 1)
              },
              
              
              
              #' greedy.search = function(direction = c("forward", "backward", "turning")) {
              #'   stopifnot(!is.null(score <- getScore()))
              #'   
              #'   ## Cast direction
              #'   direction <- match.arg(direction)
              #'   alg.name <- switch(direction,
              #'                      forward = "GIES-F",
              #'                      backward = "GIES-B",
              #'                      turning = "GIES-T")
              #'   
              #'   new.graph <- .Call("causalInference",
              #'                      .in.edges,
              #'                      score$pp.dat,
              #'                      alg.name,
              #'                      score$c.fcn,
              #'                      causal.inf.options(caching = FALSE),
              #'                      PACKAGE = "imagestest")
              #'   if (identical(new.graph, "interrupt"))
              #'     return(FALSE)
              #'   
              #'   if (new.graph$steps > 0) {
              #'     .in.edges <<- new.graph$in.edges
              #'     names(.in.edges) <<- .nodes
              #'   }
              #'   
              #'   return(new.graph$steps)
              #' },
              #' 
              #' #' Performs a causal inference from an arbitrary start DAG
              #' #' with a specified algorithm
              #' caus.inf = function(algorithm = c("GIES", "GIES-F", "GIES-B", "GIES-T", "GIES-STEP",
              #'                                   "GDS", "SiMy"), ...) {
              #'   stopifnot(!is.null(score <- getScore()))
              #'   algorithm <- match.arg(algorithm)
              #'   
              #'   print("entering")
              #'   
              #'   new.graph <- .Call("causalInference",
              #'                      .in.edges,
              #'                      score$pp.dat,
              #'                      algorithm,
              #'                      score$c.fcn,
              #'                      causal.inf.options(...),
              #'                      PACKAGE = "imagestest")
              #'   
              #'   if (identical(new.graph, "interrupt"))
              #'     return(FALSE)
              #'   else {
              #'     .in.edges <<- new.graph$in.edges
              #'     names(.in.edges) <<- .nodes
              #'     return(TRUE)
              #'   }
              #' },
              
              #' Yields a representative (estimating parameters via MLE)
              repr = function() {
                stopifnot(!is.null(score <- getScore()))
                
                result <- score$create.dag()
                result$.in.edges <- .Call("representative", .in.edges, PACKAGE = "imagestest")
                result$.params <- score$global.fit(result)
                
                return(result)
              },
              
              #' Calculates an optimal intervention target
              #'
              #' @param   max.size    maximum target size; allowed values: 1, p (= # nodes)
              ## TODO document that function... or better: provide a documented wrapper function
              opt.target = function(max.size) {
                .Call("optimalTarget", .in.edges, max.size, PACKAGE = "imagestest")
              }
            ))


