## MRes Project 2013
## Compare phylogenies to references
## 5_screen, x_referencephylogenies
## 26/11/2013

## Libraries
source(file.path('functions','EcoDataTools.R'))

## Dirs
input.dirs <- c("5_screen", "x_referencephylogenies")

## Print stage
cat("\n\nComparing phylogenies to references\n")

## Produced phylogenies
cat('\nIdentifying studies with phylogenies ...')
phylo.studies.files <- list.files(path = input.dirs[1], pattern = "^.*\\.tre$")
phylo.studies <- unique(sub("_phylo\\.tre$", "", phylo.studies.files))
cat(paste0('\nDone. Found [', length(phylo.studies), '] studies with phylogenies.'))

## Reference phylogenies
cat('\nIdentifying reference phylogenies ...')
phylo.refs.files <- list.files(path = input.dirs[2], pattern = "^.*\\.tre$")
phylo.refs <- unique(sub("\\.tre$", "", phylo.refs.files))
cat(paste0('\nDone. Found [', length(phylo.refs), '] studies with phylogenies.'))

## Compare
for (i in 1:length(phylo.studies)) {
  cat(paste0("\nWorking on ", phylo.studies[i], "... "))
  phylos <- read.tree(file.path(input.dirs[1], phylo.studies.files[i]))
  refi <- match(phylo.studies[i], phylo.refs)
  ref <- read.tree(file.path(input.dirs[2], phylo.refs.files[refi]))
  score.dist <- topo.dist <- tree.sizes <- rep(NA, length(phylos))
  for (j in 1:length(phylos)) {
    temp.ref <- drop.tip(ref, ref$tip.label[!ref$tip.label %in% phylos[[i]]$tip.label])
    topo.dist[j] <- dist.topo(phylos[[i]], temp.ref, method="PH85")
    score.dist[j] <- dist.topo(phylos[[i]], temp.ref, method="score")
    tree.sizes[j] <- length(phylos[[i]]$tip.label)
  }
  cat(paste0("\nMean tree size: ", mean(tree.sizes)))
  cat(paste0("\nMean topo. dist: ", mean(topo.dist)))
  cat(paste0("\nMean score dist: ", mean(score.dist)))
}