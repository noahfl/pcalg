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

## Define the score (BIC)
#score <- new("GaussL0penObsScore", gmG8$x)#, lambda=2)

create_im_dags <- function(num_sets) {
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
    set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(40000,0,0.1),5000,8),                       g = gGtrue)
    print(dim(set8$x))
    data_list[[i]] <- set8
    start_seed = start_seed + 2000
  }
  return(data_list)
}

# create_orig_dags <- function(num_sets) {
#   data(gmG)
#   start_seed <- 50
#   data_list <- list()
#   data_list[[1]] <- gmG8
#   
#   for (i in 2:num_sets) {
#     set.seed(start_seed)
#     p <- 8
#     n <- 5000
#     ## true DAG:
#     vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
#     #gGtrue <- gmG8$g
#     gGtrue <- randomDAG(p, prob = 0.3, V = vars)
#     #s2  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
#     set8 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)
#     
#     data_list[[i]] <- set8
#     start_seed = start_seed + 10
#   }
#   return(data_list)
# }

create_scores <- function(datasets) {
  scores <- list()
  for (i in 1:length(datasets)) {
    scores[[i]] <- new("GaussL0penObsScore", datasets[[i]]$x)
  }
  return(scores)
}

# #set.seed(40)
# set.seed(50)
# p <- 8
# n <- 5000
# ## true DAG:
# vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
# gGtrue <- gmG8$g
# s2  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
# s28 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)

#score_alt <- new("GaussL0penObsScore", s28$x)#, lambda=2)

#print(log(nrow(gmG8$x)))

## Estimate the essential graph

run_originals <- function(datasets) {
  fits = list()
  for (i in 1:length(datasets)) {
    images <- new("IMaGES", scores = list(datasets[[i]]), penalty=0)
    fits[[i]] <- images$results
  }
  return(fits)
}

run_im <- function(datasets) {
  print(length(datasets))
  results <- new("IMaGES", scores = datasets, penalty=0)
  return(results)
}

# images <- new("IMaGES", scores = list(score, score_alt), penalty=0)
# 
# images.fit <- images$results[[1]]
# images2.fit <- images$results[[2]]
# print(images.fit)



#print(str(images.fit))

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



find_error <- function(graph) {
  
  true <- inEdgeList(gmG8$g)
  #print(true)
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


#######
#create accuracy measure
#precision, recall, f-measure
#0 or 1 for each edge/direction


plot_error <- function(results) {
  inv_measures <- list()
  
  #print(paste("length: ", length(results)))
  for (i in 1:length(results)) {
    f_list <- list()
    #print(paste("results: ", length(results[[i]])))
    for (k in 1:length(results[[i]]$results)) {

      #f_list[[k]] <- find_error(results[[i]]$results[[k]][[2]]$.in.edges)
      f_list[[k]] <- 1 - find_error(results[[i]]$results[[k]][[2]]$.in.edges)

    }
    #print(f_list)
    inv_measures[[i]] <- mean(unlist(f_list))
  }

  print(inv_measures)
  
  plot_measures <- unlist(inv_measures)
  
  plot(plot_measures, type="o", col="blue", main="Error", ylim=c(0,0.5))
  axis(1, at=1:5)
}

driver <- function() {
  num_sets <- 2
  
  gmG8 <- get_gmg()
  
  dags <- create_im_dags(num_sets)
  scores <- create_scores(dags)
  orig_fits <- run_originals(scores)
  
  print(orig_fits[[1]][[1]][[2]])
  
  for (i in 1:length(orig_fits)) {
    plot_graph(orig_fits[[i]][[1]][[2]])
  }
  
  im_run_dags <- create_im_dags(num_sets)
  
  im_run_scores <- create_scores(im_run_dags)
  #print(scores == im_run_scores)
  im_fits <- run_im(im_run_scores)
  
  print(im_fits$results[[1]][[2]])
  
  par(mfrow=c(1,6))
  par(mar=c(1,1,1,1))
  for (i in 1:length(im_fits$results)) {
    plot_graph(im_fits$results[[i]][[2]])
  }
}

plot_driver <- function() {
  num_sets <- 50
  
  gmG8 <- get_gmg()

  result_sets <- list()
  
  for (k in 1:num_sets) {
    
    im_run_dags <- create_im_dags(k)
    
    im_run_scores <- create_scores(im_run_dags)
    #print(scores == im_run_scores)
    im_fits <- run_im(im_run_scores)
    result_sets[[k]] <- im_fits
    
    
  }

  plot_error(result_sets)
  
  print(length(result_sets))
  
  for (k in 1:length(result_sets)) {
    for (i in 1:length(result_sets[[k]])) {
      plot_graph(result_sets[[k]]$results[[i]][[2]])
    }
  }
  
  
}

autism_driver <- function() {
  sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);
  
  
  filenames <- list.files("test/steve", pattern="autism*", full.names=TRUE)
  matrices = list()
  
  #for (i in 1:3) {
  for (i in 1:length(filenames)) {
    matrices[[i]] <- as.matrix(read.table(filenames[[i]], skip=1))
  }
  
  results = new("IMaGES", matrices = matrices, penalty=0)
  
  for (i in 1:length(results)) {
    par(mfrow=c(1,2))
    plot(results$results[[i]][[2]], main = "Estimated CPDAG")
  }

}

# if (require(Rgraphviz)) {
#   par(mfrow=c(1,2))
#   plot(images2.fit[[2]], main = "Estimated CPDAG")
#   plot(s28$g, main = "True DAG")
# } else { ## alternative:
#   str(images2.fit, max=2)
#   as(as(images2.fit$essgraph,"graphNEL"),"Matrix")
# }

# data(gmG)
# 
# score2 <- new("GaussL0penObsScore", gmG8$x, lambda=2)
# 
# ges.fit <- ges(score2)
# 
# ## Plot the estimated essential graph and the true DAG
# if (require(Rgraphviz)) {
#   par(mfrow=c(1,2))
#   plot(ges.fit[[2]], main = "Estimated CPDAG")
#   plot(gmG8$g, main = "True DAG")
# } else { ## alternative:
#   str(ges.fit, max=2)
#   as(as(ges.fit$essgraph,"graphNEL"),"Matrix")
# }

#driver()
