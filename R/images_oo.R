IMaGES <- setRefClass("IMaGES",
                      fields = list(matrices="list", penalty="numeric", .rawscores="list", .graphs="list", imscore = "numeric", results="list", scores = "list", num.markovs = "numeric"),
                      
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
                          
                          #print("ALGORITHM")
                          #print(alg.name)
                          
                          #for (i in 1:length(.graphs)) {
                          if (!.graphs[[j]]$greedy.step(alg.name=alg.name, direction = phase, verbose = FALSE)) {
                            #stop("something happened")
                            #print("SOMETHING HAPPENED. probably that number thing")
                            #print(.graphs[[j]]$.nodes)
                          }
                          else {
                            #print("good")
                            #print(.graphs[[j]]$.nodes)
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
                          
                          #print(paste("max val: ", max_val, " index: ", index, " phase[[index]]: ", phases[[index]]))
                          
                          return(phases[[index]])
                        },
                        
                        update_score = function() {
                          
                          #print(paste("IM before: ", trueIM$score))
                          
                          imscore <- IMScore()
                          #imscore <<- IMScore()
                          #trueIM <<- imscore
                          assign("score", imscore, envir=trueIM)
                          #unlockBinding("trueIM", environment())
                          
                          
                          #print("UPDATED")
                          #print(get("score", envir=trueIM))
                          
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
                          #print("start of run")
                          update_score()
                          #print("score updated")
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
                                                     PACKAGE = "IMaGES")
                          }
                          #invisible(readline(prompt="Press [enter] to continue"))
                          #}
                          
                          opt <- mode(opt_phases)[[1]]
                          #print(opt_phases)
                          #print(paste("opt_phase(s): ", opt))
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
                            #print("Detected TURN")
                            #print("Something messed up here")
                            #invisible(readline(prompt="Press [enter] to continue"))
                          }
                          else if (opt == 0) {
                            str_opt <- 'none'
                            #print("Detected NONE")
                            #print("Something messed up here")
                            #invisible(readline(prompt="Press [enter] to continue"))
                          }
                          
                          temp.scores <- vector()
                          if (!(str_opt == "none")) {
                            for (j in 1:length(.graphs)) {
                              #print(paste("opt_phase: ", str_opt))
                              #print(toString(j))
                              
                              
                              #run_phase(phase=str_opt, j)
                              run_phase(phase=str_opt, j)
                              
                              
                              #TODO: run, get score, roll back for each
                              #then implement one with best score update for global graph
                              temp.scores[[j]] <- IMScore()
                              .graphs[[j]]$undo.step()
                            }
                            
                            #re-enable all steps after seeing which best impacts IMScore
                            for (j in 1:length(.graphs)) {
                              .graphs[[j]]$redo.step()
                              #print(paste("graph ", j))
                              update.markovs(.graphs[[j]], temp.scores[[j]])
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
                          #print(paste("SCORES: ", score.list))
                          best.graph.index <- which(score.list == max(score.list))
                          
                          #print(paste("Best index: ", best.graph.index))
                          #print(best.graph.index)
                          

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
                          #print(paste("dst: ", dst))
                          #print(trueIM$global.edges)
                          #find index where edge should be inserted
                          insert.point <- which(order(c(insert, trueIM$global.edges[[dst]]))==1)
                          #insert the edge into the edgelist for that vertex
                          trueIM$global.edges[[dst]] <- append(trueIM$global.edges[[dst]], insert, insert.point - 1)
                        },
                        
                        #remove edge from global graph, where *dst* contains the list of edges going towards it
                        remove.global = function(src, dst) {
                          #remove edge by reassigning edge list to itself, where none of the values are *src*
                          #print(paste("src: ", src))
                          #print(paste("dst: ", dst))
                          #print(paste("trueIM$global.edges[[dst]]: ", trueIM$global.edges[[dst]]))
                          #print(paste("trueIM$global.edges[[dst]][[trueIM$global.edges[[dst]]: ", trueIM$global.edges[[dst]]))
                          trueIM$global.edges[[dst]] <- trueIM$global.edges[[dst]][trueIM$global.edges[[dst]] != src]
                        },
                        
                        #perform swap on *src* and *dst* to imitate turning method
                        turn.global = function(src, dst) {
                          #print(trueIM$global.edges)
                          remove.global(src, dst)
                          #print(trueIM$global.edges)
                          insert.global(dst,src)
                          #print(trueIM$global.edges)
                          #invisible(readline(prompt="Press [enter] to continue"))
                          
                        },
                        
                        edge.exists = function(src, dst) {
                          return(src %in% trueIM$global.edges[[dst]])

                        },
                        
                        is.legal.edge = function(src, dst) {
                          #print(paste("SRC: ", src, " DST: ", dst))
                          if (src > 0 && src <= ncol(.graphs[[1]]$.score$pp.dat$data)) {
                            if (dst > 0 && dst <= ncol(.graphs[[1]]$.score$pp.dat$data)) {
                              return(TRUE)
                            }
                          }
                          #print("not true")
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
                              #print("Inserting edge")
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
                              #print("Removing edge")
                              remove.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              return(FALSE)
                            }
                          }
                          else if (dir == 'GIES-T') {
                            #turn
                            #print("Made it to TURN")
                            #invisible(readline(prompt="Press [enter] to continue"))
                            if (is.legal.edge(src,dst) && edge.exists(src, dst)) {
                              #print("Turning edge")
                              turn.global(src, dst)
                              return(TRUE)
                            }
                            else {
                              #print("Something messed up here")
                              #invisible(readline(prompt="Press [enter] to continue"))
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
                            
                            #print(.graphs[[i]]$.score$global.score.int(.graphs[[i]]$.in.edges))
                            #sum <- sum + .graphs[[i]]$.score$global.score(.graphs[[i]]$.current_repr)
                            sum <- sum + .graphs[[i]]$.score$global.score.int(.graphs[[i]]$.in.edges)
                          }
                          
                          #print(paste("n: ", n, " k: ", k, "length: ", length(.graphs), "sum: ", sum))
                          
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
                          #print("IMSCORE:")
                          #print(imscore)
                          return(imscore)
                        },
                        
                        turn.step = function() {
                          old.graph <- .graphs[[1]]$.in.edges
                          
                          best.graph <- old.graph
                          trueIM$global.edges <- old.graph
                          
                          
                          #edge destination
                          for (i in 1:length(old.graph)) {
                            #edge source
                            for (j in 1:length(old.graph[[i]])) {
                              #modify the graphs
                              #old.graph <- trueIM$global.edges
                              
                              #it can't point to itself
                              if (i == j) {
                                next
                              }
                              
                              #skip if no edges coming in
                              if (length(old.graph[[i]]) == 0) {
                                next
                              }

                              
                              #print("BEFORE")
                              #print(trueIM$global.edges)
                              if (is.legal.edge(j, i) && edge.exists(j, i) && (edge.exists(i,j))) {
                                print("removing because reverse edge exists")
                                remove.global(i, j)
                              }
                              
                              if (is.legal.edge(j, i) && edge.exists(j, i) && !(edge.exists(i,j))) {
                                turn.global(j, i)
                              }
                              #print("AFTER")
                              #print(trueIM$global.edges)
                              for (k in 1:length(.graphs)) {
                                .graphs[[k]]$.in.edges <- trueIM$global.edges
                              }
                              
                              temp.score <- IMScore()
                              
                              if (temp.score > trueIM$score) {
                                #for adding the original graph into the 
                                #MEC if it doesn't have the highest score
                                #print("updating BEST - turn")
                                update.markovs(old.graph, trueIM$score)
                                best.graph <- .graphs[[1]]$.in.edges
                              }
                              else {
                                #print("updating MEC - turn")
                                #print(paste("src: ", j, "dst: ", i))
                                update.markovs(.graphs[[1]]$.in.edges, temp.score)
                                #revert to keep scoring baseline consistent
                                turn.global(i, j)
                              }
                              
                              
                              if (is.legal.edge(j, i) && edge.exists(j, i)) {
                                #print("Removing edge")
                                remove.global(j, i)
                              }
                              
                              for (k in 1:length(.graphs)) {
                                .graphs[[k]]$.in.edges <- trueIM$global.edges
                              }
                              
                              temp.score <- IMScore()
                              
                              if (temp.score > trueIM$score) {
                                #for adding the original graph into the 
                                #MEC if it doesn't have the highest score
                                #print("updating BEST - remove")
                                update.markovs(old.graph, trueIM$score)
                                best.graph <- .graphs[[1]]$.in.edges
                              }
                              else {
                                print("updating MEC - remove")
                                update.markovs(.graphs[[1]]$.in.edges, temp.score)
                                #revert to keep scoring baseline consistent
                                insert.global(j, i)
                              }
                            

                            }
                          }
                          
                          for (i in 1:length(.graphs)) {
                            .graphs[[i]]$.in.edges <- best.graph
                          }
                          
                        },
                        
                        #fixes edge lists so there are double-directed edges
                        #
                        #
                        fix.edges = function(graph) {
                          new.graph <- rep(list(c()), length(graph))
                          for (i in 1:length(graph)) {
                            #if there are edges going to this vertex
                            if (length(graph[[i]]) > 0) {
                              #for each edge going to vertex i
                              for (j in (1:length(graph[[i]]))) {
                                if (!((i %in% graph[[graph[[i]][[j]]]]) && (i < graph[[i]][[j]]))) {
                                  #adding edges that only point forward -- hacky fix for now
                                  new.graph[[i]] <- append(new.graph[[i]], graph[[i]][[j]])
                                }
                              }
                            }
                          }
                          new.graph
                        },
                        
                        
                        #compares graphs and returns true if they're identical 
                        #and false otherwise
                        is.identical = function(graph1, graph2) {
                          for (i in 1:length(graph1)) {
                            comp.bools <- graph1[[i]] == graph2[[i]]
                            #print(comp.bools)
                            
                            if (!(length(graph1[[i]] == 0))) {
                              if (length(graph2[[i]]) == 0) {
                                #print("graph1 not null, graph2 null")
                                #invisible(readline(prompt="Press [enter] to continue"))
                                return(FALSE)
                              }
                            }
                            if (!(length(graph2[[i]]) == 0)) {
                              if (length(graph1[[i]] == 0)) {
                                #print("graph1 null, graph2 not null")
                                #invisible(readline(prompt="Press [enter] to continue"))
                                return(FALSE)
                              }
                            }


                            # if ((typeof(graph1[[i]]) != "logical" && typeof(graph2[[i]]) != "logical" && 
                            #      graph1[[i]] != integer(0) && graph2[[i]] != integer(0) && 
                            #      !is.null(graph1[[i]]) && !is.null(graph2[[i]]))) {
                            if (length(graph1[[i]]) != 0 && length(graph2[[i]]) != 0) {
                              #print("logical true")
                              #print(graph1[[i]])
                              #print(graph2[[i]])
                              #invisible(readline(prompt="Press [enter] to continue"))
                              #print(comp.bools)
                              if (!all(comp.bools)) {
                                #print("not equal")
                                #invisible(readline(prompt="Press [enter] to continue"))
                                return(FALSE)
                                
                              }
                            }
                            # if (!all(comp.bools)) {
                            #   return(FALSE)
                            # }
                          }
                          #print("TRUE")
                          #invisible(readline(prompt="Press [enter] to continue"))
                          return(TRUE)
                        },
                        
                        update.markovs = function(graph, score) {
                          #iterate backwards to find lowest score that it's higher than
                          sum <- 0
                          
                          for (i in 1:length(trueIM$markovs)) {
                            # if it is higher, bump lowest markov and sort in the new one
                            #exclude graphs that are identical since that's not too helpful
                            #if (score > trueIM$markovs[[i]]$.score) {
                            #print("-------")
                            #print(trueIM$markovs[[i]]$.score)
                            #print(paste(i, ": ", typeof(trueIM$markovs[[i]])))
                            #print("here")
                            if (is.identical(graph$.in.edges, trueIM$markovs[[i]]$.graph)) {
                              #print("made it here")
                              #sum <- sum + 1
                              return()
                            }
                            #if (sum == length(trueIM$markovs)) {
                            #  return()
                            #}
                          }
                          
                          
                          for (i in 1:length(trueIM$markovs)) {
                            # if it is higher, bump lowest markov and sort in the new one
                            #exclude graphs that are identical since that's not too helpful
                            #if (score > trueIM$markovs[[i]]$.score) {
                            #print(trueIM$markovs)
                            #print(paste("score: ", score))
                            #print(paste("MScore: ", trueIM$markovs[[i]]$.score))
                            res <- (score > trueIM$markovs[[i]]$.score)
                            #print(res)
                            if (score > trueIM$markovs[[i]]$.score && !is.null(score) && !is.null(trueIM$markovs[[i]]$.score)) {
                              # for (k in length(trueIM$markovs) - 1:i + 1) {
                              #   trueIM$markovs[[k + 1]] <- trueIM$markovs[[k]]
                              # }
                              #print("made it into edit")
                              converted <- convert(list(.in.edges = graph$.in.edges, .nodes = graph$.nodes))
                              markov <- list(.graph=graph$.in.edges, .score=score, .data=graph$.score$pp.dat$data)#,
                                                          #.params=apply.sem(converted, graph$.score$pp.dat$data))
                              
                              #print(trueIM$markovs[[i]], )
                              num.markovs <<- num.markovs
                              
                              for (k in num.markovs - 1:i + 1) {
                                if (i == num.markovs) {
                                  trueIM$markovs[[i]] <- markov
                                  break
                                }
                                #print(paste("i: ", i))
                                #print(paste(k, k - 1))
                                trueIM$markovs[[k]] <- trueIM$markovs[[k - 1]]
                              }
                              trueIM$markovs[[i]] <- markov
                              #trueIM$markovs <-  append(trueIM$markovs, markov, i - 1)
                              #print(trueIM$markovs)
                              # for (k in num.markovs + 1:length(trueIM$markovs)) {
                              #   trueIM$markovs[[i]] <- NULL
                              # }
                              
                              #trueIM$markovs <- trueIM$markovs[1:trueIM$num.markovs]
                              #trueIM$markovs <- trueIM$markovs[1:5]
                              # markovs <- list()
                              # for (j in 1:trueIM$num.markovs) {
                              #   markovs[[i]] <- trueIM$markovs[[i]]
                              # }
                              # trueIM$markovs <- markovs
                              # print(trueIM$markovs)
                              break
                            }
                          }
                          
                        },
                        
                        convert = function(from) {
                          
                          edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
                          names(edgeList) <- from$.nodes
                          #print("edgeList")
                          #print(edgeList)
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
                          #print(length(params.list))
                          base <- params.list[[1]]
                          if (length(params.list) == 1) {
                            return(base)
                          }
                          for (i in 2:length(params.list)) {
                            #print(length(params.list[[i]]))
                            for (j in 1:length(params.list[[i]])) {
                              base[j] <- base[j] + params.list[[i]][j]
                            }
                          }
                          
                          for (i in 1:length(params.list)) {
                            #print(paste("base[", i,"] before: ", base[i]))
                            base[i] <- round((base[i] / length(params.list)),2)
                            #print(paste("base[", i,"] after: ", base[i]))
                          }
                          return(base)
                        }

                      )
)

IMaGES$methods(
  initialize = function(matrices = NULL, scores = NULL, penalty = 3, imscore = NULL, num.markovs=5) {
    #images <-
    #print("initializing")
    #imscore = 0
    rawscores <- list()
    
    #print(paste("penalty: ", penalty))
    
    if (!is.null(matrices)) {
      #print(paste("LENGTH: ", length(matrices)))
      assign("numDatasets", length(matrices), envir=trueIM)
      for (i in 1:length(matrices)) {
        
        #print("adding score")
        rawscores[[i]] <- new("GaussL0penObsScore", matrices[[i]], lambda = penalty)
      }  
    }
    else {
      assign("numDatasets", length(scores), envir=trueIM)
      for (i in 1:length(scores)) {
        scores[[i]]$pp.dat$lambda <- penalty
        rawscores[[i]] <- scores[[i]]
      }
    }
    
    #trueIM$num.markovs <- num.markovs
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
    
    print("Running...")
    
    #create list of size num.markovs
    #initialize to lowest signed int value because you can't compare
    #ints and NULL
    trueIM$markovs <- rep(list(list(.graph=NULL, .score=-2147483648)), num.markovs)
    #print(trueIM$markovs)
    
    for (i in 1:(ncol(.graphs[[1]]$.score$pp.dat$data) * ncol(.graphs[[1]]$.score$pp.dat$data))) {
      #run IMaGES
      run()
    }
    
    #####turn.step()
    .graphs[[1]]$.in.edges <- fix.edges(.graphs[[1]]$.in.edges)
    
    if (length(.graphs) > 1) {
      for (i in 2:length(.graphs)) {
        .graphs[[i]]$.in.edges <- .graphs[[1]]$.in.edges
      }
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
    converted <- convert(list(.in.edges = graphs[[1]]$.in.edges, .nodes = .graphs[[1]]$.nodes))
    #print(.graphs[[1]]$.in.edges)
    #alt_converted <- convert(list(.in.edges = .graphs[[1]]$.in.edges, .nodes = .graphs[[1]]$.nodes))
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
    
    #print(single.graphs[[i]]$.params)
    
    print("Done.")
    
    #imscore <<- IMScore()
    
    #switched these for the time being
    global <- list(.graph = converted, .params = average.sem(params.list))
    #alt <- list(.graph = converted, .params = average.sem(params.list))
    
    markovs <- list()
    
    #print(trueIM$markovs)
    
    for (i in 1:num.markovs) {
      #print(trueIM$markovs)
      
      attempt <- tryCatch(
        {
          converted.markov <- convert(list(.in.edges = fix.edges(trueIM$markovs[[i]]$.graph), .nodes = .graphs[[1]]$.nodes))
          markovs[[i]] <- list(.graph=converted.markov, .params = apply.sem(converted.markov, trueIM$markovs[[i]]$.data))
        },
        error = function(e) {
          markovs[[i]] <- NULL
        }
      )

    }
    
    results <<- list(.global = global, .single.graphs = single.graphs, .markovs = markovs)
    
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
                #print("calling GES")
                new.graph <- .Call("causalInferenceEdge",
                                   #.current_repr$.in.edges,
                                   .in.edges,
                                   score$pp.dat,
                                   alg.name,
                                   score$c.fcn,
                                   causal.inf.options(caching = FALSE, maxSteps = 1, verbose = verbose, phase=direction, ...),
                                   PACKAGE = "IMaGES")
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
                .Call("optimalTarget", .in.edges, max.size, PACKAGE = "IMaGES")
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

#plots all 
plotMarkovs = function(im.fits) {
  single.length <- length(im.fits$results$.markovs)
  print(single.length)
  plot.vals <- find_dimensions(single.length + 1)
  par(mfrow=plot.vals)
  plotIMGraph(im.fits$results$.global)
  
  for (i in 1:single.length) {
    plotIMGraph(im.fits$results$.markovs[[i]], title=paste("MEC Graph ", i))
  }
  
}

### GLOBAL VAR

#trueIM <- 100
options(warn = 0)
trueIM <- new.env()
assign("num.markovs", 5, env=trueIM)
assign("best", list(), env=trueIM)
assign("score", 100, env=trueIM)
assign("isLocalIM", FALSE, env=trueIM)
assign("numDatasets", 1, env=trueIM)
assign("global.edges", list(), env=trueIM)
assign("markovs", list(), env=trueIM)