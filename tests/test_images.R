## Load predefined data
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

make_data <- function(prob) {
  
  set.seed(40)
  p <- 8
  n <- 5000
  ## true DAG:
  vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
  gGtrue <- randomDAG(p, prob = prob, V = vars)
  dataset <- list(x = rmvDAG(n, gGtrue), g = gGtrue)
  #dataset<-dataset$x
}


#create the DAGs generated by gmG8 plus added noise
create_im_dags <- function(num_sets) {
  gmG8 <- get_gmg()
  #initial seed for generation of dataset
  start_seed <- 3.1415926
  data_list <- list()
  
  for (i in 1:num_sets) {
    set.seed(start_seed)
    p <- 8
    n <- 5000
    ## true DAG:
    vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
    gGtrue <- gmG8$g
    #inject noise into DAGs using rnorm
    set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(40000,0,0.1),5000,8), g = gGtrue)
    data_list[[i]] <- set8
    #increment start seed
    start_seed = start_seed + 2000
  }
  return(data_list)
}

#create and return list of score objects to be passed into IMaGES
create_scores <- function(datasets) {
  scores <- list()
  for (i in 1:length(datasets)) {
    scores[[i]] <- new("GaussL0penObsScore", datasets[[i]]$x)
  }
  return(scores)
}

#function used for evaluating original, GES-style runs on individual datasets
run_originals <- function(datasets) {
  fits = list()
  for (i in 1:length(datasets)) {
    images <- new("IMaGES", scores = list(datasets[[i]]), penalty=3)
    fits[[i]] <- images
  }
  return(fits)
}


## Plot the estimated essential graph and the true DAG
plot_graph <- function(fit) {
  if (require(Rgraphviz)) {
    par(mfrow=c(1,2))
    plot(fit, main = "Estimated CPDAG")
    plot(gmG8$g, main = "True DAG")
  } else { ## alternative:
    str(fit, max=2)
    as(as(fit$essgraph,"graphNEL"),"Matrix")
  }
}


#calculates precision, recall, and f-measure of graphs by finding:
#    true positives (edges in both graphs)
#    false positives (edges in the generated graph but not in the true DAG)
#    false negatives (edges that should have been in the generated graph but weren't)
#and using these for the calculations
find_error <- function(graph) {
  
  true <- inEdgeList(gmG8$g)
  
  positives <- 0
  true_positives <- 0
  false_positives <- 0
  false_negatives <- 0
  
  for (i in 1:length(graph)) {
      for (j in 1:length(graph[[i]])) {

        if (((length(graph[[i]]) > 0) && (length(true[[i]]) > 0)) && (graph[[i]][[j]] %in% true[[i]])) {
          positives <- positives + 1
        }
        else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(graph[[i]][[j]] %in% true[[i]])) {
          false_positives <- false_positives + 1
        }
      }

    if (length(true[[i]]) > 0) {
      for (k in 1:length(true[[i]])) {
        if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && (true[[i]][[k]] %in% graph[[i]])) {
          true_positives <- true_positives + 1
        }
        else if ((length(graph[[i]]) > 0) && (length(true[[i]]) > 0) && !(true[[i]][[k]] %in% graph[[i]])) {
          false_negatives <- false_negatives + 1
        }
      }
    }
  }
  
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  
  f_measure <- 2 * ((precision * recall) / (precision + recall))
  return(f_measure)
}


#######
#create accuracy measure
#precision, recall, f-measure
#0 or 1 for each edge/direction

#function used for plotting the error across multiple iterations of IMaGES and showing
#how it decreases due to the "IMaGES effect"
plot_error <- function(results) {
  
  inv_measures <- list()
  
  for (i in 1:length(results)) {
    #calculate f-measures for each graph in the set
    f_list <- list()
    for (k in 1:length(results[[i]]$results)) {
      #use 1 minus error for the sake of presentation
      f_list[[k]] <- 1 - find_error(results[[i]]$results[[k]][[2]]$.in.edges)

    }
    #print(f_list)
    #use mean error (although should all be the same) for graph
    inv_measures[[i]] <- mean(unlist(f_list))
  }

  #cannot be in list form for plotting
  plot_measures <- unlist(inv_measures)
  
  plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
  axis(1, at=1:5)
}

#driver for individual GES-like runs
driver <- function() {
  #change to how many graphs you want
  num_sets <- 3
  
  gmG8 <- get_gmg()
  
  #generate DAGS
  dags <- create_im_dags(num_sets)
  #create score objects
  scores <- create_scores(dags)
  #find GES-like fits using IMaGES
  orig_fits <- run_originals(scores)
  
  #print(orig_fits[[1]][[1]][[2]])
  
  #plot in.edges for each graph
  for (i in 1:length(orig_fits)) {
    plotIMGraph(orig_fits[[i]]$results$.global)
  }
  
  #now do the same thing for IMaGES
  
  #create DAGS
  im_run_dags <- create_im_dags(num_sets)
  
  #create score objects
  im_run_scores <- create_scores(im_run_dags)
  #run IMaGES
  im_fits <- new("IMaGES", scores = im_run_scores, penalty=3)
  
  #plot results
  par(mfrow=c(1,2))
  plotIMGraph(im_fits$results$.global)
}

driver_prob <- function() {
  #change to how many graphs you want
  num_sets <- 3
  
  #now do the same thing for IMaGES
  
  #create DAGS
  #im_run_dags <- create_im_dags(num_sets)
  
  dataset1 <- make_data(0.3)
  dataset2 <- make_data(0.3)
  dataset3 <- make_data(0.3)
  
  #create score objects
  #im_run_scores <- create_scores(list(dataset1,dataset2,dataset3))
  im_run_scores <- create_scores(list(dataset1))
  #run IMaGES
  im_fits <- new("IMaGES", scores = im_run_scores, penalty=3)
  
  
  plotIMGraph(im_fits$results$.global)
  
  #plot results
  # for (i in 1:length(im_fits$results)) {
  #   plot_graph(im_fits$results[[i]][[2]])
  # }
}

#driver for calculation of errors across runs of increasing size
plot_driver <- function() {
  #change to number of sets to iterate up to
  num_sets <- 10
  
  #generate gmG8 data
  gmG8 <- get_gmg()

  #stores fits for each set size
  result_sets <- list()
  
  for (k in 1:num_sets) {
    #create DAGs
    im_run_dags <- create_im_dags(k)
    #create score objects
    im_run_scores <- create_scores(im_run_dags)
    #run IMaGES
    im_fits <- new("IMaGES", scores = im_run_scores, penalty=3)
    #append results to result_sets
    result_sets[[k]] <- im_fits
    
    
  }
  #calculates errors for each of the result sets
  plot_error(result_sets)
  
  #plots individual sets (might creash computer as it's a lot of plots)
  for (k in 1:length(result_sets)) {
    for (i in 1:length(result_sets[[k]])) {
      plot_graph(result_sets[[k]]$results[[i]][[2]])
    }
  }
  
  
}
#driver for running IMaGES on autism data. works but the plot still isn't showing up properly
#it might be due to the fact that the labels aren't included?
autism_driver <- function() {
  #get file locations
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames 
  filenames <- list.files("test/steve", pattern="autism*", full.names=TRUE)
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
  #for (i in 1:2) {
    #this might be causing the plotting issue but i'm not yet sure of a workaround
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=TRUE))
    #test1 <- as.matrix(read.table(filenames[[i]]))#, skip=1))
    #test2 <- as.matrix(read.table(filenames[[i]]), skip=1)
    
  }
  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=5)
  
  plotIMGraph(results$results$.global)
  plotIMGraph(results$results$.alt)
  
  plotAll(results)
  
  #plot resulting DAGs
  # for (i in 1:length(results)) {
  #   par(mfrow=c(1,2))
  #   plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  # }

}

powerball_driver <- function() {
  #get file locations
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames 
  filenames <- list.files("test/powerball", pattern="pb*", full.names=TRUE)
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
    #for (i in 1:2) {
    #this might be causing the plotting issue but i'm not yet sure of a workaround
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=FALSE))
    #test1 <- as.matrix(read.table(filenames[[i]]))#, skip=1))
    #test2 <- as.matrix(read.table(filenames[[i]]), skip=1)
    
  }
  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=1)
  
  plotIMGraph(results$results$.global)
  
  
  #plot resulting DAGs
  # for (i in 1:length(results)) {
  #   par(mfrow=c(1,2))
  #   plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  # }
  
}

test_dataset <- function() {
  #get file locations
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  #get filenames c
  filenames <- list.files("test/d9", pattern="dataset*", full.names=TRUE)
  matrices = list()
  
  #import data
  for (i in 1:length(filenames)) {
    #for (i in 1:2) {
    #this might be causing the plotting issue but i'm not yet sure of a workaround
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], header=FALSE))
    #test1 <- as.matrix(read.table(filenames[[i]]))#, skip=1))
    #test2 <- as.matrix(read.table(filenames[[i]]), skip=1)
    
  }
  
  #run IMaGES on data
  results = new("IMaGES", matrices = matrices, penalty=3)
  
  
  
  plotIMGraph(results$results$.global)
  plotAll(results)
  
  #plot resulting DAGs
  # for (i in 1:length(results)) {
  #   par(mfrow=c(1,2))
  #   plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  # }
  
}

convert <- function(from) {
  edgeList <- lapply(from$.in.edges, function(v) from$.nodes[v])
  names(edgeList) <- from$.nodes
  result <- new("graphNEL",
                nodes = from$.nodes,
                edgeL = edgeList,
                edgemode = "directed")
  return(reverseEdgeDirections(result))

}

