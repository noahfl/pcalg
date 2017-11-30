## Load predefined data
data(gmG)

## Define the score (BIC)
score <- new("GaussL0penObsScore", gmG8$x)#, lambda=2)


#set.seed(40)
set.seed(50)
p <- 8
n <- 5000
## true DAG:
vars <- c("Author", "Bar", "Ctrl", "Goal", paste0("V",5:8))
gGtrue <- randomDAG(p, prob = 0.3, V = vars)
gmA  <- list(x = rmvDAG(n, gGtrue, back.compatible=TRUE), g = gGtrue)
gmA8 <- list(x = rmvDAG(n, gGtrue),                       g = gGtrue)

score_alt <- new("GaussL0penObsScore", gmA8$x)#, lambda=2)

#print(log(nrow(gmG8$x)))

## Estimate the essential graph
images <- new("IMaGES", scores = list(score, score_alt), penalty=0)

images.fit <- images$results[[1]]
images2.fit <- images$results[[2]]
print(images.fit)



print(str(images.fit))

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(images.fit[[2]], main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
} else { ## alternative:
  str(images.fit, max=2)
  as(as(images.fit$essgraph,"graphNEL"),"Matrix")
}

if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(images2.fit[[2]], main = "Estimated CPDAG")
  plot(gmA8$g, main = "True DAG")
} else { ## alternative:
  str(images2.fit, max=2)
  as(as(images2.fit$essgraph,"graphNEL"),"Matrix")
}

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

