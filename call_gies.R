#source("https://bioconductor.org/biocLite.R")
#biocLite()

#loadNamespace("NAMESPACE")
#library(pcalg)
#library(Rcpp)
#library(methods)
#print("1")
#library(graph)
#print("2")
#library(utils)
#print("3")



sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);

#print(lst);

filenames <- list.files("test/steve", pattern="autism*", full.names=TRUE)
matrices = list()

for (i in 1:3) {
  #for (i in 1:length(filenames)) {
        matrices[[i]] <- as.matrix(read.table(filenames[[i]], skip=1))
}

tst <- new("GaussL0penObsScore", matrices[[1]])

#print(tst$global.score(tst$create.dag()))

scores = list()

for (i in 1:length(matrices)) {
  scr <- new("GaussL0penObsScore", matrices[[i]])
  #scores[[i]] <- ges(scr)
  scores[[i]] <- scr
}

#tst <- ges(scores[[1]])

images = new("IMaGES", matrices = matrices, penalty=1.5)

print(images$.graphs[[1]])

#res <- IMaGES(matrices,scores)

#print(scores[[1]]$global.score(scores[[1]]$create.dag()))

# IMScore <- function(matrices, scores, penalty) {
#   #print(typeof(length(scores)))
#   m <- length(scores)
#   #print(m)
#   #print(dim(matrices[[1]]))
#   #n <- dim(matrices[[1]])[[2]]
#   #print(typeof(n))
#   #print(n)
#   sum <- 0
#   #k <- 5
#   
#   for (i in 1:length(scores)) {
#     #sum <- sum + scores[[i]]$global.score(scores[[i]]$create.dag()) + ((penalty * k) * log(n))
#     sum <- sum + scores[[i]]$global.score(scores[[i]]$create.dag())
#   }
#   
#   print(sum)
#   imscore = (-2/m) *  sum
#   print(imscore)
#   return(imscore)
# }



#ges(tst)

