## MRes Project 2013
## Compare phylogenies to references
## 5_screen, x_referencephylogenies
## 26/11/2013

## Libraries
source(file.path('functions','EcoDataTools.R'))

## Dirs
input.dirs <- c("4_phylogenies", "x_referencephylogenies")

## Print stage
cat("\n\nComparing phylogenies to references\n")

## Produced phylogenies
cat('\nIdentifying studies with phylogenies ...')
phylo.studies.files <- list.files(path = input.dirs[1], pattern = "^.*\\.tre$")
phylo.studies <- unique(sub("_gene_.*\\.tre$", "", phylo.studies.files))
cat(paste0('\nDone. Found [', length(phylo.studies), '] studies with phylogenies.'))

## Reference phylogenies
cat('\nIdentifying reference phylogenies ...')
phylo.refs.files <- list.files(path = input.dirs[2], pattern = "^.*\\.tre$")
phylo.refs <- unique(sub("\\.tre$", "", phylo.refs.files))
cat(paste0('\nDone. Found [', length(phylo.refs), '] studies with phylogenies.'))

## Compare
pdf("consensus_phylogeneies.pdf", w = 14)
for (i in 1:length(phylo.studies)) {
  cat(paste0("\nWorking on ", phylo.studies[i], "... "))
  phylos <- list()
  for (each in phylo.studies.files) {
    if (grepl(paste0("^", phylo.studies[i]), each)) {
      phylos <- c(phylos, list(read.tree(file.path(input.dirs[1], each))))
    }
  }
  class(phylos) <- 'multiPhylo'
  #phylo.consensus <- consensus(phylos, p = 0.5)
  #plot(phylo.consensus, main = paste0(phylo.studies[i], " Consensus p = 0.5"))
  refi <- match(phylo.studies[i], phylo.refs)
  ref <- read.tree(file.path(input.dirs[2], phylo.refs.files[refi]))
  if (class(ref) == "multiPhylo") {
    ref <- ref[[1]]
  }
  ref <- drop.tip(ref, ref$tip.label[!ref$tip.label %in% phylo.consensus$tip.label])
  plot(ref, main = paste0(phylo.studies[i], " Reference"))
  score.dist <- topo.dist <- tree.sizes <- rep(NA, length(phylos))
  for (j in 1:length(phylos)) {
    if (length(phylos[[j]]$tip.label) > 10) {
      temp.ref <- drop.tip(ref, ref$tip.label[!ref$tip.label %in% phylos[[j]]$tip.label])
      #plot(temp.ref)
      #plot(phylos[[j]])
      topo.dist[j] <- dist.topo(phylos[[j]], temp.ref, method="PH85")
      score.dist[j] <- dist.topo(phylos[[j]], temp.ref, method="score")
      tree.sizes[j] <- length(phylos[[j]]$tip.label)
    }
  }
  cat(paste0("\nMean tree size: ", mean(tree.sizes, na.rm = TRUE)))
  cat(paste0("\nMean topo. dist: ", mean(topo.dist, na.rm = TRUE)))
  cat(paste0("\nMean score dist: ", mean(score.dist, na.rm = TRUE)))
}
dev.off()