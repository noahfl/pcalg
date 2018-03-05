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
                          
                          print(paste("IM before: ", trueIM$score))
                          
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
                          #ux[which.max(tabulate(match(x, ux)))]
                          tab <- tabulate(match(x, ux)); ux[tab == max(tab)]
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
                          #invisible(readline(prompt="Press [enter] to continue"))
                          #}
                          
                          opt <- mode(opt_phases)
                          print(opt_phases)
                          print(paste("opt_phase(s): ", opt))
                          #invisible(readline(prompt="Press [enter] to continue"))
                          
                          str_opt <- ''
                          
                          for (i in 1:length(opt_phases)) {
                            if (opt_phases[[i]] == 1) {
                              opt_phases[[i]] <- 'GIES-F'
                            }
                            else if (opt_phases[[i]] == 2) {
                              opt_phases[[i]] <- 'GIES-B'
                            }
                            else if (opt_phases[[i]] == 3) {
                              opt_phases[[i]] <- 'GIES-T'
                            }
                            else if (opt_phases[[i]] == 0) {
                              opt_phases[[i]] <- 'none'
                            }
                          }
                          
                          if (opt == 1) {
                            str_opt <- 'GIES-F'
                          }
                          else if (opt == 2) {
                            str_opt <- 'GIES-B'
                          }
                          else if (opt == 3) {
                            str_opt <- 'GIES-T'
                            print("Detected TURN")
                            #print("Something messed up here")
                            invisible(readline(prompt="Press [enter] to continue"))
                          }
                          else if (opt == 0) {
                            str_opt <- 'none'
                            print("Detected NONE")
                            #print("Something messed up here")
                            #invisible(readline(prompt="Press [enter] to continue"))
                          }
                          
                          temp.scores <- vector()
                          if (!(str_opt == "none")) {
                            for (j in 1:length(.graphs)) {
                              print(paste("opt_phase: ", str_opt))
                              #print(toString(j))
                              
                              
                              #run_phase(phase=str_opt, j)
                              run_phase(phase=str_opt, j)
                              
                              
                              #TODO: run, get score, roll back for each
                              #then implement one with best score update for global graph
                              temp.scores[[j]] <- IMScore()
                              .graphs[[j]]$undo.step()
                            }
                            find.best.step(temp.scores)
                          }

                          #best.graph.index <- which(temp.scores == max(temp.scores))
                          #, best.graph.index)
                          # print(paste("Best index: ", best.graph.index))
                          # print(best.graph.index)
                          # 
                          # #re-enable all steps after seeing which best impacts IMScore
                          # for (j in 1:length(.graphs)) {
                          #   .graphs[[j]]$redo.step()
                          # }
                          # #update IMScore with all changes
                          # update_score()
                          # 
                          # #just use first index of best.graph.index for now
                          # #TODO: find smarter way to do this
                          # if (!update.global(.graphs[[best.graph.index[[1]]]]$.edge.change)) {
                          #   
                          # }
                              
                              #return(TRUE)
                          # else {
                          #   return(FALSE)
                          # }

  
                        },
                        
                        find.best.step = function(score.list) {
                          # for (i in 1:length(score.list)) {
                          #   .graphs[[i]]$undo.step()
                          # }
                          inf <- 0
                          for (i in 1:length(score.list)) {
                            if (is.infinite(score.list[[i]])) {
                              inf <- inf + 1
                            }
                          }
                          if (inf == length(score.list)) {
                            return()
                          }
                          print(paste("SCORES: ", score.list))
                          best.graph.index <- which(score.list == max(score.list))
                          print(paste("Best index: ", best.graph.index))
                          print(best.graph.index)
                          
                          #re-enable all steps after seeing which best impacts IMScore
                          for (j in 1:length(.graphs)) {
                            .graphs[[j]]$redo.step()
                          }
                          #update IMScore with all changes
                          update_score()
                          
                          #just use first index of best.graph.index for now
                          #TODO: find smarter way to do this
                          if (!update.global(.graphs[[best.graph.index[[1]]]]$.edge.change)) {
                            #make the score lower than the rest so it considers the next highest score
                            #better way to do this? probably
                            score.list[[best.graph.index[[1]]]] <- -Inf
                            find.best.step(score.list)
                          }
                        },
                        
                        #insert edge into global graph, where dst contains the list of edges going towards it
                        insert.global = function(src, dst) {
                          insert <- src
                          print(paste("dst: ", dst))
                          print(trueIM$global.edges)
                          #find index where edge should be inserted
                          insert.point <- which(order(c(insert, trueIM$global.edges[[dst]]))==1)
                          #insert the edge into the edgelist for that vertex
                          trueIM$global.edges[[dst]] <- append(trueIM$global.edges[[dst]], insert, insert.point - 1)
                        },
                        
                        #remove edge from global graph, where *dst* contains the list of edges going towards it
                        remove.global = function(src, dst) {
                          #remove edge by reassigning edge list to itself, where none of the values are *src*
                          print(paste("src: ", src))
                          print(paste("dst: ", dst))
                          print(paste("trueIM$global.edges[[dst]]: ", trueIM$global.edges[[dst]]))
                          #print(paste("trueIM$global.edges[[dst]][[trueIM$global.edges[[dst]]: ", trueIM$global.edges[[dst]]))
                          trueIM$global.edges[[dst]] <- trueIM$global.edges[[dst]][trueIM$global.edges[[dst]] != src]
                        },
                        
                        #perform swap on *src* and *dst* to imitate turning method
                        turn.global = function(src, dst) {
                          print(trueIM$global.edges)
                          remove.global(src, dst)
                          print(trueIM$global.edges)
                          insert.global(dst,src)
                          print(trueIM$global.edges)
                          invisible(readline(prompt="Press [enter] to continue"))
                          
                        },
                        
                        edge.exists = function(src, dst) {
                          return(src %in% trueIM$global.edges[[dst]])

                        },
                        
                        is.legal.edge = function(src, dst) {
                          print(paste("SRC: ", src, " DST: ", dst))
                          if (src > 0 && src <= ncol(.graphs[[1]]$.score$pp.dat$data)) {
                            if (dst > 0 && dst <= ncol(.graphs[[1]]$.score$pp.dat$data)) {
                              return(TRUE)
                            }
                          }
                          print("not true")
                          return(FALSE)
                        },
                        
                        #handles updating of global graph. calls insert.global, remove.global,
                        #or turn.global depending on what the edge change specifies
                        update.global = function(edge.change) {
                          src <- edge.change[[1]]
                          dst <- edge.change[[2]]
                          dir <- edge.change[[3]]
                          if (dir == 'GIES-F') {
                            #insert
                            if (is.legal.edge(src,dst) && !(edge.exists(src, dst))) {
                              print("Inserting edge")
                              insert.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              return(FALSE)
                            }
                          }
                          else if (dir == 'GIES-B') {
                            #remove
                            if (is.legal.edge(src,dst) && edge.exists(src, dst)) {
                              print("Removing edge")
                              remove.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              return(FALSE)
                            }
                          }
                          else if (dir == 'GIES-T') {
                            #turn
                            print("Made it to TURN")
                            invisible(readline(prompt="Press [enter] to continue"))
                            if (is.legal.edge(src,dst) && edge.exists(src, dst)) {
                              print("Turning edge")
                              turn.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              print("Something messed up here")
                              invisible(readline(prompt="Press [enter] to continue"))
                              return(FALSE)
                              
                            }
                          }
                          else {
                            return(TRUE)
                          }
                        },
                        
                        initialize.global = function() {
                          edges <- list()
                          
                          for (i in 1:ncol(.graphs[[1]]$.score$pp.dat$data)) {
                            edges[[i]] <- vector()
                          }
                          
                          assign("global.edges", edges, env=trueIM)
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
                            
                            print(.graphs[[i]]$.score$global.score.int(.graphs[[i]]$.in.edges))
                            #sum <- sum + .graphs[[i]]$.score$global.score(.graphs[[i]]$.current_repr)
                            sum <- sum + .graphs[[i]]$.score$global.score.int(.graphs[[i]]$.in.edges)
                          }
                          
                          print(paste("n: ", n, " k: ", k, "length: ", length(.graphs), "sum: ", sum))
                          
                          #print("SUM")
                          #print(sum)
                          #print(sum)
                          #imscore <- ((-2/m) *  sum) + ((penalty * k) * log(n))
                          imscore <- ((2/m) *  sum) + ((penalty * k) * log(n))
                          
                          # if (!(length(.graphs) == 1)) {
                          #   assign("isLocalIM", TRUE, envir=trueIM) 
                          # }
                          
                          #imscore <<- imscore
                          #print(paste("imscore in function: ", imscore, "\n"))
                          print("IMSCORE:")
                          print(imscore)
                          return(imscore)
                        },
                        
                        convert = function(from) {
                          
                          edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
                          names(edgeList) <- from$.nodes
                          print("edgeList")
                          print(edgeList)
                          result <- new("graphNEL",
                                        nodes = from$.nodes,
                                        edgeL = edgeList,
                                        edgemode = "directed")
                          return(reverseEdgeDirections(result))
                          
                          #edges<-IMtest$results$.in.edges
                        },
                        
                        #does structural equation modeling on graphNEL
                        #object (needs dataset)
                        apply.sem = function(converted, dataset) {
                          #print("made it here")
                          graph.nodes <- igraph.from.graphNEL(converted)
                          edge.list <- get.edgelist(graph.nodes)
                          #print(edge.list)
                          model <- paste(edge.list[,1], "~", edge.list[,2])
                          #print(model)
                          fit <- sem(model, data=data.frame(dataset))
                          estimate <- partable(fit)$est
                          estimate <- round(estimate,2)
                          #print(edgeNames(converted))
                          #print(estimate)
                          names(estimate) <- edgeNames(converted)
                          return(estimate)
                        },
                        
                        average.sem = function(params.list) {
                          print(length(params.list))
                          base <- params.list[[1]]
                          for (i in 1:length(params.list)) {
                            print(length(params.list[[i]]))
                            for (j in 1:length(params.list[[i]])) {
                              base[j] <- base[j] + params.list[[i]][j]
                            }
                          }
                          
                          for (i in 1:length(params.list)) {
                            base[i] <- round((base[i] / length(params.list)),2)
                          }
                          return(base)
                        }

                      )
)

IMaGES$methods(
  initialize = function(matrices = NULL, scores = NULL, penalty = 3, imscore = NULL) {
    #images <-
    #print("initializing")
    #imscore = 0
    rawscores <- list()
    
    #print(paste("penalty: ", penalty))
    
    if (!is.null(matrices)) {
      print(paste("LENGTH: ", length(matrices)))
      assign("numDatasets", length(matrices), envir=trueIM)
      for (i in 1:length(matrices)) {
        
        #print("adding score")
        rawscores[[i]] <- new("GaussL0penObsScore", matrices[[i]], lambda = penalty)
        print(rawscores[[i]]$pp.dat$lambda)
      }  
    }
    else {
      assign("numDatasets", length(scores), envir=trueIM)
      for (i in 1:length(scores)) {
        print(paste("length of scores: ", length(scores)))
        print(paste("num rows: ", nrow(scores[[i]]$pp.dat$data)))
        #scores[[i]]$pp.dat$lambda <- scores[[i]]$pp.dat$lambda * penalty
        scores[[i]]$pp.dat$lambda <- penalty
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
    
    initialize.global()
    
    #imscore <<- IMScore()
    #imscore <<- IMScore()
    
    print("---------------")
    
    
    for (i in 1:(2*ncol(.graphs[[1]]$.score$pp.dat$data) * ncol(.graphs[[1]]$.score$pp.dat$data))) {
      # if (!run()) {
      #   break
      # }
      run()
    }
    
    # for (i in 1:ncol(.graphs[[1]]$.score$pp.dat$data) * ncol(.graphs[[1]]$.score$pp.dat$data)) {
    #   temp.scores <- vector()
    # 
    #   for (j in 1:length(.graphs)) {
    # 
    #     #run_phase(phase=str_opt, j)
    #     run_phase(phase="GIES-T", j)
    #     
    #     
    #     #TODO: run, get score, roll back for each
    #     #then implement one with best score update for global graph
    #     temp.scores[[j]] <- IMScore()
    #     .graphs[[j]]$undo.step()
    #   }
    #   
    #   #best.graph.index <- which(temp.scores == max(temp.scores))
    #   find.best.step(temp.scores)
    # }
    
    single.graphs <- list()
    params.list <- list()
    converted <- convert(list(.in.edges = trueIM$global.edges, .nodes = .graphs[[1]]$.nodes))
    print(.graphs[[1]]$.in.edges)
    alt_converted <- convert(list(.in.edges = .graphs[[1]]$.in.edges, .nodes = .graphs[[1]]$.nodes))
    for (i in 1:length(.graphs)) {
      #create .in.edges structure and convert it to graphNEL object
      #converted <- convert(list(.in.edges = .graphs[[i]]$.in.edges, .nodes = .graphs[[i]]$.nodes))
      #print("Type of converted: ")
      #print(converted)
      
      params <- apply.sem(converted, .graphs[[i]]$.score$pp.dat$data)
      params.list[[i]] <- params
      single.converted <- convert(list(.in.edges = .graphs[[i]]$.in.edges, .nodes = .graphs[[i]]$.nodes))
      single.graphs[[i]] <-list(.graph = single.converted, .params = params)
    }
    
    print(single.graphs[[i]]$.params)
    
    print("---------------")
    
    #imscore <<- IMScore()
    
    #switched these for the time being
    global <- list(.graph = alt_converted, .params = average.sem(params.list))
    alt <- list(.graph = converted, .params = average.sem(params.list))
    
    results <<- list(.global = global, .single.graphs = single.graphs, .alt = alt)
    
    #results$.in.edges <- trueIM$global.edges
    #results$.nodes <- .graphs[[1]]$.nodes
    return(results)
    
    #imscore <<- IMScore()
    #print("FINAL SCORES")
    #for (i in 1:length(.graphs)) {
      #print(.graphs[[i]]$.score$global.score(.graphs[[i]]$.current_repr))
      #print(list(.graphs[[i]], .graphs[[i]]$repr()))
      #results[[i]] <<- list(.graphs[[i]], .graphs[[i]]$repr())
    #}
    
    #results <<- results
    
  }
)

#' Interventional essential graph for IMaGES
setRefClass("IMGraph",
            fields = list(
              .nodes = "vector",
              .in.edges = "list",
              .targets = "list",
              .score = "Score",
              #.current_repr = "list",
              .old.edges = "list",
              .edge.change = "list"
              
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
                new.graph <- .Call("causalInferenceEdge",
                                   #.current_repr$.in.edges,
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
                #print("BEFORE")
                #print(.in.edges)
                
                
                .old.edges <<- .in.edges
                .in.edges <<- new.graph$in.edges
                #.in.edges <<- .Call("representative", new.graph$in.edges, PACKAGE = "imagestest")
                .saved.edges <<- new.graph$in.edges
                #test <<- .Call("representative", new.graph$in.edges, PACKAGE = "imagestest")
                #print(test)
                names(.in.edges) <<- .nodes
                .edge.change <<- list(new.graph$edge.change$src, new.graph$edge.change$dst, new.graph$edge.change$alg)
                #print("AFTER")
                #print(.in.edges)
                
                #names(.in.edges) <<- .nodes
                
                #return(new.graph$steps == 1)
                return(FALSE)
              },
              
              undo.step = function() {
                .in.edges <<- .old.edges
              },
              
              redo.step = function() {
                # print("SAME?: ")
                # print(.in.edges)
                # print("---------------")
                # print(.saved.edges)
                #invisible(readline(prompt="Press [enter] to continue"))
                .in.edges <<- .saved.edges
              },
              
              
              #' Yields a representative (estimating parameters via MLE)
              repr = function() {
                stopifnot(!is.null(score <- getScore()))
                
                result <- new("GaussParDAG", nodes =.nodes)
                # result$.in.edges <- .current_repr$.in.edges
                result$.in.edges <- .in.edges
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

#finds lowest two values that, when multiplied, create square-ish grid for
#plotting graphs
find_dimensions = function(number) {
  val1 <- floor(sqrt(number))
  val2 <- ceiling(sqrt(number))
  #res1 <- val1
  #res2 <- val2
  
  while (val1*val2 < number) {
    val2 <- val2 + 1
    if (val1 * val2 < number) {
      val1 <- val1 - 1
    }
  }
  
  return(c(val1, val2))
  
}




#plot function for graphs generated by IMaGES
plotIMGraph = function(graph.list, title="Global") {
  graph <- graph.list$.graph
  params <- graph.list$.params
  finalgraph <- agopen(graph, "", attrs=list(node=list(shape="ellipse")), edgeAttrs=list(label=params))
  plot(finalgraph, main=title)
}

#plots all 
plotAll = function(im.fits) {
  single.length <- length(im.fits$results$.single.graphs)
  plot.vals <- find_dimensions(single.length + 1)
  par(mfrow=plot.vals)
  plotIMGraph(im.fits$results$.global)
  
  for (i in 1:single.length) {
    plotIMGraph(im.fits$results$.single.graphs[[i]], title=paste("Graph ", i))
  }
  
}

### GLOBAL VAR

#trueIM <- 100
trueIM <- new.env()
assign("score", 100, env=trueIM)
assign("isLocalIM", FALSE, env=trueIM)
assign("numDatasets", 1, env=trueIM)
assign("global.edges", list(), env=trueIM)