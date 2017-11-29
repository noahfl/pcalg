## Load predefined data
data(gmG)

## Define the score (BIC)
score <- new("GaussL0penObsScore", gmG8$x, lambda=2)

print(log(nrow(gmG8$x)))

## Estimate the essential graph
images <- new("IMaGES", scores = list(score), penalty=0)

images.fit <- images$results[[1]]

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

