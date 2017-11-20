## Load predefined data
data(gmG)

## Define the score (BIC)
score <- new("GaussL0penObsScore", gmG8$x, lambda=2)

print(log(nrow(gmG8$x)))

## Estimate the essential graph
images <- new("IMaGES", scores = list(score), penalty=0)

ges.fit <- images$results[[1]]



print(str(ges.fit))

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(ges.fit[[2]], main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
} else { ## alternative:
  str(ges.fit, max=2)
  as(as(ges.fit$essgraph,"graphNEL"),"Matrix")
}

score <- new("GaussL0penObsScore", gmG8$x, lambda=2)

ges2.fit <- ges(score)

## Plot the estimated essential graph and the true DAG
if (require(Rgraphviz)) {
  par(mfrow=c(1,2))
  plot(ges2.fit[[2]], main = "Estimated CPDAG")
  plot(gmG8$g, main = "True DAG")
} else { ## alternative:
  str(ges2.fit, max=2)
  as(as(ges2.fit$essgraph,"graphNEL"),"Matrix")
}

