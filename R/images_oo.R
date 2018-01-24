IMaGES <- setRefClass("IMaGES",
                      fields = list(matrices="list", penalty="numeric", .rawscores="list", .graphs="list", imscore = "numeric", results="list", scores = "list"),
                      
                      # validity <- function(object) {
                      #   return(TRUE)
                      # },
                      
                      
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
                            #print(score)
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
                        
                        run_phase = function(phase="GIES-F", j) {
                          
                          alg.name <- phase
                          
                          # alg.name <- ""
                          # 
                          # #phase <- match.arg(phase)
                          # 
                          # if (phase == "forward") {
                          #   alg.name <- "GIES-F"
                          # }
                          # else if (phase == "backward") {
                          #   alg.name <- "GIES-B"
                          # }
                          # else if (phase == "turning") {
                          #   alg.name <- "GIES-T"
                          # }
                          # else {
                          #   print(paste("alg name: ", alg.name))
                          #   stop("incorrect algorithm name")
                          # }
                          
                          print("ALGORITHM")
                          print(alg.name)
                          
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
                        
                        #only finds first max. look into this
                        find_opt = function(opt_phases) {
                          phases <- list()
                          phase_counts = list()
                          for (i in 1:length(opt_phases)) {
                            if ((length(phases) == 0) || !(is.element(opt_phases[[i]], phases))) {
                              phases[[length(phases) + 1]] <- opt_phases[[i]]
                              phase_counts[[length(phase_counts) + 1]] = 1
                            }
                            else {
                              phase_counts[[match(opt_phases[[i]], phases)]] = phase_counts[[match(opt_phases[[i]], phases)]] + 1
                            }
                            
                          }
                          max_val = max(unlist(phase_counts))
                          index = match(max_val, phase_counts)
                          
                          print(paste("max val: ", max_val, " index: ", index, " phase[[index]]: ", phases[[index]]))
                          
                          return(phases[[index]])
                        },
                        
                        update_score = function() {
                          
                          imscore <- IMScore()
                          #imscore <<- IMScore()
                          #trueIM <<- imscore
                          assign("score", imscore, envir=trueIM)
                          #unlockBinding("trueIM", environment())
                          
                          
                          print("UPDATED")
                          print(get("score", envir=trueIM))
                          
                          # for (i in 1:length(.graphs)) {
                          #   #print("INITIAL")
                          #   #print(.graphs[[i]]$.score$.imscore)
                          #   .graphs[[i]]$.score$.imscore <<- imscore
                          #   
                          #   #print("AFTER")
                          #   #print(.graphs[[i]]$.score$.imscore)
                          #   .graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=imscore)
                          # }
                        },
                        
                        mode = function(x) {
                          ux <- unique(x)
                          ux[which.max(tabulate(match(x, ux)))]
                        },
                        
                        run = function() {
                          print("start of run")
                          update_score()
                          print("score updated")
                          phases = list("forward", "backward", "turning")
                          
                          opt_phases = list()
                          
                          #for (i in 1:length(phases)) {
                          for (j in 1:length(.graphs)) {
                            #print(toString(j))
                            #print("MADE IT HERE")
                            #opt_phases[[i]] = .Call("greedyStepRFunc", .graphs[[i]$.in.edges, PACKAGE = "imagestest")
                            opt_phases[[j]] <- .Call("greedyStepRFunc",
                                                     .graphs[[j]]$.in.edges,
                                                     .graphs[[j]]$.score$pp.dat,
                                                     .graphs[[j]]$.score$c.fcn,
                                                     .graphs[[j]]$causal.inf.options(caching = FALSE, maxSteps = 1),
                                                     PACKAGE = "imagestest")
                          }
                          #}
                          
                          opt <- mode(opt_phases)
                          
                          str_opt <- ''
                          
                          if (opt == 1) {
                            str_opt <- 'GIES-F'
                          }
                          else if (opt == 2) {
                            str_opt <- 'GIES-B'
                          }
                          else if (opt == 3) {
                            str_opt <- 'GIES-T'
                          }
                          else if (opt == 0) {
                            str_opt <- 'none'
                          }                          
                          
                          
                          
                          # str_phases <- list()
                          # for (i in 1:length(opt_phases)) {
                          #   if (opt_phases[[i]] == 1) {
                          #     str_phases[[i]] <- 'GIES-F'
                          #   }
                          #   else if (opt_phases[[i]] == 2) {
                          #     str_phases[[i]] <- 'GIES-B'
                          #   }
                          #   else if (opt_phases[[i]] == 3) {
                          #     str_phases[[i]] <- 'GIES-T'
                          #   }
                          #   else if (opt_phases[[i]] == 0) {
                          #     str_phases[[i]] <- 'none'
                          #   }
                          # }
                          # 
                          # print(paste("phases: ", str_phases))
                          # opt = find_opt(str_phases)
                          
                          if (!(str_opt == "none")) {
                            for (j in 1:length(.graphs)) {
                              print(paste("opt_phase: ", opt))
                              #print(toString(j))
                              run_phase(phase=str_opt, j)
                              update_score()
                              
                            }
                          }
                          
                          
                          
                          print("HERE?")
                        },
                        
                        IMScore = function() {
                          
                          #graphs <- .graphs
                          #penalty <- get('penalty', envir=globalenv())
                          penalty <<- penalty
                          #matrices <- get('matrices', envir=globalenv())
                          
                          #print(typeof(length(scores)))
                          m <- length(.graphs)
                          #print(m)
                          #print(dim(matrices[[1]]))
                          n <- ncol(.graphs[[1]]$.score$pp.dat$data) * nrow(.graphs[[1]]$.score$pp.dat$data)
                          #print("N")
                          #print(n)
                          #print(typeof(n))
                          #print(n)
                          sum <- 0
                          k <- nrow(.graphs[[1]]$.score$pp.dat$data)
                          #print("VARS")
                          #print(list(penalty, m, n, sum, k))
                          assign("isLocalIM", FALSE, envir=trueIM)
                          
                          
                          
                          for (i in 1:length(.graphs)) {
                            #sum <- sum + scores[[i]]$global.score(scores[[i]]$create.dag()) + ((penalty * k) * log(n))
                            # print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=.graphs[[i]]$.score$.imscore))
                            # sum <- sum + .graphs[[i]]$.score$global.score(.graphs[[i]]$.score$create.dag(), .imscore=.graphs[[i]]$.score$.imscore)
                            
                            print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.current_repr))
                            sum <- sum + .graphs[[i]]$.score$global.score(.graphs[[i]]$.current_repr)
                          }
                          
                          print(paste("n: ", n, " k: ", k, "length: ", length(.graphs), "sum: ", sum))
                          
                          #print("SUM")
                          #print(sum)
                          #print(sum)
                          #imscore <- ((-2/m) *  sum) + ((penalty * k) * log(n))
                          imscore <- ((2/m) *  sum) + ((penalty * k) * log(n))
                          
                          if (!(length(.graphs) == 1)) {
                            assign("isLocalIM", TRUE, envir=trueIM) 
                          }
                          
                          #imscore <<- imscore
                          #print(paste("imscore in function: ", imscore, "\n"))
                          print("IMSCORE:")
                          print(imscore)
                          return(imscore)
                        }
                      )
)

IMaGES$methods(
  initialize = function(matrices = NULL, scores = NULL, penalty = 1.5, imscore = NULL) {
    #images <-
    #print("initializing")
    #imscore = 0
    rawscores <- list()
    
    test.env <- new.env()
    test.env$tst <- "TEST"
    #print(paste("penalty: ", penalty))
    
    if (!is.null(matrices)) {
      assign("numDatasets", length(matrices), envir=trueIM)
      for (i in 1:length(matrices)) {
        
        #print("adding score")
        rawscores[[i]] <- new("GaussL0penObsScore", matrices[[i]])
      }  
    }
    else {
      assign("numDatasets", length(scores), envir=trueIM)
      for (i in 1:length(scores)) {
        print(paste("length of scores: ", length(scores)))
        print(paste("num rows: ", nrow(scores[[i]]$pp.dat$data)))
        rawscores[[i]] <- scores[[i]]
      }
    }
    
    penalty <<- penalty
    .rawscores <<- rawscores
    
    graphs <- list()
    
    for (i in 1:length(.rawscores)) {
      #print("creating graph")
      graphs[[i]] <- create.graph(rawscores[[i]])
    }
    
    .graphs <<- graphs
    
    #imscore <<- IMScore()
    #imscore <<- IMScore()
    
    print("---------------")
    
    for (i in 1:(ncol(.graphs[[1]]$.score$pp.dat$data) * ncol(.graphs[[1]]$.score$pp.dat$data))) {
      #print("test")
      run()      
    }
    
    
    print("---------------")
    
    #imscore <<- IMScore()
    
    results <<- list()
    
    imscore <<- IMScore()
    print("FINAL SCORES")
    for (i in 1:length(.graphs)) {
      print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.current_repr))
      #print(list(.graphs[[i]], .graphs[[i]]$repr()))
      results[[i]] <<- list(.graphs[[i]], .graphs[[i]]$repr())
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
              .score = "Score",
              .current_repr = "list"
              
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
                .current.repr <<- new("GaussParDAG", nodes=.nodes, in.edges=.in.edges)
                #names(.in.edges) <<- NULL
                
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
                .current_repr$.in.edges <<- .Call("representative", new.graph$in.edges, PACKAGE = "imagestest")
                names(.in.edges) <<- .nodes
                print("AFTER")
                print(.in.edges)
                
                #names(.in.edges) <<- .nodes
                
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


### GLOBAL VAR

#trueIM <- 100
trueIM <- new.env()
assign("score", 100, env=trueIM)
assign("isLocalIM", FALSE, env=trueIM)
assign("numDatasets", 1, env=trueIM)