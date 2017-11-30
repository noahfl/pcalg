## Load predefined data
data(gmG)

## Define the score (BIC)
#score <- new("GaussL0penObsScore", gmG8$x)#, lambda=2)

create_im_dags <- function(num_sets) {
  data(gmG)
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
    set8 <- list(x = rmvDAG(n, gGtrue)+ matrix(rnorm(40000,0,1),5000,8),                       g = gGtrue)
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
    fits[[i]] <- images$results[[1]]
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
    plot(fit[[2]], main = "Estimated CPDAG")
    plot(gmG8$g, main = "True DAG")
  } else { ## alternative:
    str(fit, max=2)
    as(as(fit$essgraph,"graphNEL"),"Matrix")
  }
}

#######
#create accuracy measure
#precision, recall, f-measure
#0 or 1 for each edge/direction


driver <- function() {
  num_sets <- 3
  
  dags <- create_im_dags(num_sets)
  scores <- create_scores(dags)
  orig_fits <- run_originals(scores)
  
  for (i in 1:length(orig_fits)) {
    plot_graph(orig_fits[[i]])
  }
  
  im_run_dags <- create_im_dags(num_sets)
  
  im_run_scores <- create_scores(im_run_dags)
  #print(scores == im_run_scores)
  im_fits <- run_im(im_run_scores)
  
  
  par(mfrow=c(1,6))
  par(mar=c(1,1,1,1))
  for (i in 1:length(im_fits$results)) {
    plot_graph(im_fits$results[[i]])
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

score2 <- new("GaussL0penObsScore", gmG8$x, lambda=2)

ges.fit <- ges(score2)

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(ges.fit[[2]], main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
} else { ## alternative:
  str(ges.fit, max=2)
  as(as(ges.fit$essgraph,"graphNEL"),"Matrix")
}

