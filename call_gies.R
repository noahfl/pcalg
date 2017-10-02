#source("https://bioconductor.org/biocLite.R")
#biocLite()

#loadNamespace("NAMESPACE")
#library(pcalg)
library(Rcpp)
library(methods)
print("1")
library(graph)
print("2")
library(utils)
print("3")



sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);

#print(lst);

filenames <- list.files("test/5-10", pattern="*.txt", full.names=TRUE)
matrices = list()

for (i in 1:length(filenames)) {
        matrices[[i]] <- as.matrix(read.table(filenames[[i]], skip=1))
}

tst <- new("GaussL0penObsScore", matrices[[1]])

ges(tst)

