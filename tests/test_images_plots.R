library(IMaGES)
library(sfsmisc)
library(graph)

get_gmg <- function() {
  set.seed(40)
  p <- 8
  n <- 5000
  ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- randomDAG(p, prob = 0.3, V = vars)
  gmG  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
  gmG8 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)
  return(gmG8)
}

## Define the score (BIC)
#score <- new("GaussL0penObsScore", gmG8$x)#, lambda=2)

randomDAG <- function (n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
{
  stopifnot(n >= 2, is.numeric(prob), length(prob) == 1,
            0 <= prob, prob <= 1,
            is.numeric(lB), is.numeric(uB), lB <= uB)
  edL <- vector("list", n)
  nmbEdges <- 0L
  for (i in seq_len(n - 2)) {
    listSize <- rbinom(1, n - i, prob)
    nmbEdges <- nmbEdges + listSize
    edgeList <- sample(seq(i + 1, n), size = listSize)
    weightList <- runif(length(edgeList), min = lB, max = uB)
    edL[[i]] <- list(edges = edgeList, weights = weightList)
  }
  ## i=n-1 separately
  ## (because of sample(7,1) is actually sample(1:7,1) and not 7)
  listSize <- rbinom(1, 1, prob)
  if (listSize > 0) {
    nmbEdges <- nmbEdges + 1
    edgeList <- n
    weightList <- runif(1, min = lB, max = uB)
  } else {
    edgeList <- integer(0)
    weightList <- numeric(0)
  }
  edL[[n-1]] <- list(edges = edgeList, weights = weightList)
  if (nmbEdges > 0) {
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    names(edL) <- V
    new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  }
  else
    new("graphNEL", nodes = V, edgemode = "directed")
}

rmvDAG <-
  function(n, dag,
           errDist = c("normal", "cauchy", "t4", "mix", "mixt3", "mixN100"),
           mix = 0.1, errMat = NULL, back.compatible = FALSE,
           use.node.names = !back.compatible)
  {
    ## Purpose: Generate data according to a given DAG (with weights) and
    ## given node distribution (rows: number of samples; cols: node values in
    ## topological order)
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## - n     : Number of samples
    ## - dag   : Graph object containing the DAG and weights
    ## - errDist: "normal" or "mix" for pure standard normal node distribution
    ##           or mixing with standard cauchy
    ## - mix   : Percentage of points sampled from standard cauchy
    ## ----------------------------------------------------------------------
    ## Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler
    
    ## check input &  initialize variables
    stopifnot(is(dag, "graph"),
              (p <- length(nodes(dag))) >= 2)
    
    ##  as(.,"matrix") now {for some versions of 'graph' pkg} is 0/1
    ## weightMatrix <- t(as(dag,"matrix"))
    weightMatrix <- if(back.compatible) wgtMatrix.0(dag) else wgtMatrix(dag)
    
    ## check if top. sorted
    nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
    if (nrow(nonZeros) > 0) {
      if (any(nonZeros[,1] - nonZeros[,2] < 0) || any(diag(weightMatrix) != 0))
        stop("Input DAG must be topologically ordered!")
    }
    
    errDist <- match.arg(errDist)
    if(grepl("^mix", errDist))
      eMat <- function(outs) { # (n,p)
        X <- c(rnorm(n*p - length(outs)), outs)
        matrix(sample(X), nrow = n)
      }
    if(is.null(errMat)) {
      ## generate errors e_i
      errMat <-
        switch(errDist,
               "normal" = matrix(rnorm  (n*p),  nrow = n),
               "cauchy" = matrix(rcauchy(n*p),  nrow = n),
               "t4" =     matrix(rt(n*p, df = 4), nrow = n),
               "mix"    = eMat(rcauchy(round(mix*n*p))),
               "mixt3"  = eMat(     rt(round(mix*n*p), df = 3)),
               "mixN100"= eMat(  rnorm(round(mix*n*p), sd = 10)))
    }
    else { ## check & use 'errMat' argument:
      stopifnot(!is.null(dim.eM <- dim(errMat)),
                dim.eM == c(n,p), is.numeric(errMat))
    }
    if(use.node.names)
      colnames(errMat) <- nodes(dag) # == colnames(weightMatrix)
    
    ## compute X matrix X_i
    if (sum(weightMatrix) > 0) {
      X <- errMat
      for (j in 2:p) { ## uses X[*, 1:(j-1)] -- "recursively" !
        ij <- 1:(j-1)
        X[,j] <- X[,j] + X[, ij, drop = FALSE] %*% weightMatrix[j, ij]
      }
      X
    }
    else
      errMat
  }

wgtMatrix.0 <- function(g, transpose = TRUE)
{
  ## Purpose: work around "graph" package's  as(g, "matrix") bug
  ## ----------------------------------------------------------------------
  ## ACHTUNG: mat_[i,j]==1 iff j->i,
  ## whereas with as(g,"matrix") mat_[i,j]==1 iff i->j
  ## ----------------------------------------------------------------------
  ## Arguments: g: an object inheriting from (S4) class "graph"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006
  
  ## MM: another buglet for the case of  "no edges":
  if(numEdges(g) == 0) {
    p <- length(nd <- nodes(g))
    return( matrix(0, p,p, dimnames = list(nd, nd)) )
  }
  ## Usual case, when there are edges:
  if(!("weight" %in% names(edgeDataDefaults(g))))
    edgeDataDefaults(g, "weight") <- 1L
  w <- unlist(edgeData(g, attr = "weight"))
  ## we need the *transposed* matrix typically:
  tm <- if(transpose) t(as(g, "matrix")) else as(g, "matrix")
  ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
  if(any(w != 1)) ## fix it
    tm[tm != 0] <- w
  ## tm_[i,j]==1 iff i->j
  tm
}

wgtMatrix <- function(g, transpose = TRUE) {
  res <- as(g, "matrix") # from 'graph' package, now reliable (we hope)
  if (transpose) ## default!
    t(res) else res
}


create_im_dags <- function(num_sets, noise) {
  gmG8 <- get_gmg()
  start_seed <- 3.1415926
  data_list <- list()
  #data_list[[1]] <- gmG8
  
  for (i in 1:num_sets) {
    set.seed(start_seed)
    p <- 8
    n <- 5000
    ## true DAG:
    vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
    gGtrue <- gmG8$g
    #s2  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
    set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(40000,0,noise),5000,8),                       g = gGtrue)
    print(dim(set8$x))
    data_list[[i]] <- set8
    start_seed = start_seed + 2000
  }
  return(data_list)
}


create_scores <- function(datasets) {
  scores <- list()
  for (i in 1:length(datasets)) {
    scores[[i]] <- new("GaussL0penObsScore", datasets[[i]]$x)
  }
  return(scores)
}

run_im <- function(datasets) {
  print(length(datasets))
  results <- new("IMaGES", scores = datasets, penalty=3)
  return(results)
}

plot_graph <- function(fit) {
  gmG8 <- get_gmg()
  if (require(Rgraphviz)) {
    par(mfrow=c(1,2))
    plot(fit, main = "Estimated CPDAG")
    plot(gmG8$g, main = "True DAG")
  } else { ## alternative:
    str(fit, max=2)
    as(as(fit$essgraph,"graphNEL"),"Matrix")
  }
}

#' Auxiliary function reading an edge list (as used in the constructors
#' of DAGs) out of an adjacency matrix or a graphNEL object
#' @param from adjacency matrix, graphNEL object, or object inherited
#'  from ParDAG
#' @return list of in-edges; length of list = number of vertices,
#' entries for i-th vertex = indices sources of in-edges
inEdgeList <- function(from)
{
  if (is.matrix(from)) {
    p <- nrow(from)
    stopifnot(p == ncol(from))
    lapply(1:p, function(i) which(from[, i] != 0))
  } else if (class(from) == "graphNEL") {
    nodeNames <- graph::nodes(from)
    edgeList <- lapply(graph::inEdges(from), function(v) match(v, nodeNames))
    names(edgeList) <- NULL
    edgeList
  } else if (length(grep(".*ParDAG", class(from)) == 1)) {
    from$.in.edges
  }else {
    stop(sprintf("Input of class '%s' is not supported.", class(from)))
  }
}

find_error <- function(graph) {
  
  gmG8 <- get_gmg()  
  true <- inEdgeList(gmG8$g)
  print(true)
  print(graph)
  graph <- inEdgeList(graph)
  #print(graph)
  
  positives <- 0
  true_positives <- 0
  false_positives <- 0
  false_negatives <- 0
  
  #true_list <- list()
  
  for (i in 1:length(graph)) {
    #print(length(graph))
    #if (length(graph$.in.edges[[i]]) > 0) {
    for (j in 1:length(graph[[i]])) {
      #print(paste("graph[[i]]: ", graph[[i]], "\n"))
      # if !(true$.in.edges[[i]][[j]] %in% graph$.in.edges[[i]]) {
      #   
      # }
      #print(paste("length of graph edgelist: ", length(graph$.in.edges[[i]])))
      #print(paste("length of true edgelist: ", length(true[[i]])))
      #print("HERE")
      
      #print(paste("len graph[[i]]: ", length(graph[[i]]), "len true[[i]]: ", length(true[[i]])))
      if (((length(graph[[i]]) > 0) && (length(true[[i]]) > 0)) && (graph[[i]][[j]] %in% true[[i]])) {
        #if (!(is.null(graph[[i]])) && !(is.null(true[[i]])) && (graph[[i]][[j]] %in% true[[i]])) {
        #if ((length(true[[i]]) > 0) && (graph$.in.edges[[i]][[j]] %in% true[[i]])) {
        #print("positive")
        positives <- positives + 1
        
      }
      
      else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(graph[[i]][[j]] %in% true[[i]])) {
        #else if ((length(true[[i]]) > 0) && !(graph$.in.edges[[i]][[j]] %in% true[[i]])) {
        #print("false")
        false_positives <- false_positives + 1
      }
    }
    #}
    #print("here")
    if (length(true[[i]]) > 0) {
      #print("HERE")
      for (k in 1:length(true[[i]])) {
        if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && (true[[i]][[k]] %in% graph[[i]])) {
          #print("true pos")
          true_positives <- true_positives + 1
        }
        else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(true[[i]][[k]] %in% graph[[i]])) {
          #print("false neg")
          false_negatives <- false_negatives + 1
        }
      }
    }
  }
  
  #print(paste("Positives: ", positives, "True positives: ", true_positives, "False positives: ", false_positives, "False negatives: ", false_negatives))
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f_measure <- 2 * ((precision * recall) / (precision + recall))
  return(f_measure)
}


#plot_error <- function(results, fname, num_sets) {
#  inv_measures <- list()
#  
#  #print(paste("length: ", length(results)))
#  for (i in 1:length(results)) {
#    f_list <- list()
#    #print(paste("results: ", length(results[[i]])))
#    for (k in 1:length(results[[i]]$results)) {
#
#      #f_list[[k]] <- find_error(results[[i]]$results[[k]][[2]]$.in.edges)
#      f_list[[k]] <- 1 - find_error(results[[i]]$results[[k]][[2]]$.in.edges)
#
#    }
#    #print(f_list)
#    inv_measures[[i]] <- mean(unlist(f_list))
#  }
#
#  print(inv_measures)
#  
#  png(filename=fname)
#  plot_measures <- unlist(inv_measures)
#  
#  plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
#  axis(1, at=1:num_sets)
#  dev.off()
#}
plot_error <- function(results, fname, num_sets, noise) {
  inv_measures <- list()
  
  #print(paste("length: ", length(results)))
  for (i in 1:length(results)) {
    f_list <- list()
    #print(paste("results: ", length(results[[i]])))
    for (k in 1:length(results[[i]]$results)) {
      
      #f_list[[k]] <- find_error(results[[i]]$results[[k]][[2]]$.in.edges)
      f_list[[k]] <- 1 - find_error(results[[i]]$results$.global$.graph)
      
    }
    #print(f_list)
    inv_measures[[i]] <- mean(unlist(f_list))
  }
  
  print(inv_measures)
  
  png(filename=fname, width=800, height=400)
  plot_measures <- unlist(inv_measures)
  
  #dev.new(width=10, height=5)
  #plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
  plot(plot_measures, main=paste("IMaGES Error for noise value", noise), type="o", col='blue', ylim=c(0,0.5), xlab='', ylab='')
  at <- seq(from=0, to=num_sets, by=num_sets/20)
  title(xlab="Number of datasets")
  title(ylab="Error")
  axis(side = 1, at = at)
  #axis(1, at=1:num_sets)
  dev.off()
}


plot_driver <- function(num_sets, noise, fname) {
  
  #gmG8 <- get_gmg()
  
  result_sets <- list()
  
  for (k in 1:num_sets) {
    
    im_run_dags <- create_im_dags(k, noise)
    
    im_run_scores <- create_scores(im_run_dags)
    #print(scores == im_run_scores)
    im_fits <- run_im(im_run_scores)
    result_sets[[k]] <- im_fits
    
    
  }
  
  saveRDS(result_sets, paste(noise, "_", num_sets, "results.rds", sep=""))
  
  plot_error(result_sets, fname, num_sets, noise)
  
  print(length(result_sets))
  
  #for (k in 1:length(result_sets)) {
  #  for (i in 1:length(result_sets[[k]])) {
  #    plot_graph(result_sets[[k]]$results[[i]][[2]])
  #  }
  #}
  
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("Please supply noise value", call.=FALSE)
} else if (length(args) > 1) {
  stop("Too many arguments!")
}

noise <- as.numeric(args[[1]])
filename <- paste("../plot_noise_", noise, ".png", sep="")
plot_driver(10, noise, filename)

