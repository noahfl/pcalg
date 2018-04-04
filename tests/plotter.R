library(sfsmisc)
library(graph)
library(igraph)
library(Rgraphviz)


## Define the score (BIC)
#score <- new("GaussL0penObsScore", gmG8$x)#, lambda=2)




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


find_error <- function(graph, seed, prob) {
  gmG8 <- get_gmg(seed, prob)
  target <- as_adjacency_matrix(igraph.from.graphNEL(graph), names=FALSE)
  true <- as_adjacency_matrix(igraph.from.graphNEL(gmG8$g), names=FALSE)
  
  retrieved <- sum(target)
  precision <- sum(target & true) / retrieved
  recall <- sum(target & true) / sum(true)
  f_measure <- 2 * (precision * recall) / (precision + recall)
  return(f_measure)

}


plot_error <- function(error_sets, fname, num_sets, prob) {

  #print(error_sets)

  averaged <- list()

  for (i in 1:num_sets) {
    current_sum <- 0
    for (k in 1:length(error_sets)) {
      current_sum <- current_sum + error_sets[[k]][[i]]
    }
    averaged[[i]] <- current_sum / length(error_sets)
  }
  
  png(filename=paste(fname,".png", sep=""), width=800, height=400)
  plot_measures <- unlist(averaged)
  print(plot_measures) 
  #dev.new(width=10, height=5)
  #plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
  plot(plot_measures, main=paste("IMaGES Averaged Error", sep=""), type="o", col='blue', ylim=c(0,1), xlab='Number of datasets', ylab='Error')
  at <- seq(from=0, to=num_sets, by=num_sets/20)
  #title(xlab="Number of datasets")
  #title(ylab="Error")
  axis(side = 1, at = at)
  #axis(1, at=1:num_sets)
  abline(lsfit(1:num_sets,plot_measures),lwd=4,col="red")
  dev.off()
}



plotter_driver <- function(num_sets, prob, fname) {
 
  sapply(list.files(pattern="[.]rds$", path="R/", full.names=TRUE), source)

  filenames <- list.files(".", pattern=paste("poster_", num_sets, "_errors", sep=""), full.names=TRUE)

  error_sets <- list()
  
  for (k in 1:length(filenames)) {
    
    
    error_sets[[k]] <- readRDS(filenames[[k]])
    
    
  }
  
  
  plot_error(error_sets, fname, num_sets, prob)
  
  #print(length(result_sets))
  
  #for (k in 1:length(result_sets)) {
  #  for (i in 1:length(result_sets[[k]])) {
  #    plot_graph(result_sets[[k]]$results[[i]][[2]])
  #  }
  #}
  
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("Please supply arguments: num_runs", call.=FALSE)
} else if (length(args) > 2) {
  stop("Too many arguments!")
}

num_runs <- as.numeric(args[[1]])
prob <- as.numeric(args[[2]])

filename <- paste("../poster_plot_", num_runs, "_", prob, sep="")
plotter_driver(num_runs, prob, filename)


