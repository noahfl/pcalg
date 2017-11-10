

sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source);


filenames <- list.files("test/steve", pattern="autism*", full.names=TRUE)
matrices = list()

#for (i in 1:3) {
for (i in 1:length(filenames)) {
  matrices[[i]] <- as.matrix(read.table(filenames[[i]], skip=1))
}
images = new("IMaGES", matrices = matrices, penalty=1.5)

