## MRes Project 2013
## Stage 5: Screen phylogenies
## In: 2_taxids, 5_phylogenies | Out: 5_screen
## 28/08/2013

## Print stage
cat("\n\nThis is stage 5: screen\n")

## Parameters
min.rtt.multiplier <- 100 # how much larger than max rtt dist is to the median
pblack <- 0.9 # the proportion of trees dropped due to long branches across a study that
              # cause the responsible ids and genes to be black listed when the process is
              # re-run.

## Libraries
require("ape")

## Dirs
input.dirs <- c("1_taxids", "4_phylogenies")
output.dir <- "5_screen"
if(!file.exists(output.dir))
{
  dir.create(output.dir)
}

## Phylogenies
cat('\nIdentifying studies with phylogenies ...')
phylo.files <- list.files(path = input.dirs[2], pattern = "^.*\\.tre$")
phylo.studies <- unique(sub("_gene_.*_[0-9]{1,3}\\.tre$", "", phylo.files))
cat(paste0('\nDone. Found [', length(phylo.studies), '] studies with phylogenies.'))

## Names
cat('\nImporting qnames and taxids and writing to taxa.names.obj...')
taxa.names.obj <- list()
for (i in 1:length(phylo.studies)) {
  qname.file <- paste0(phylo.studies[i], '_qnames.txt')
  taxids.file <- paste0(phylo.studies[i], '_taxids.txt')
  qnames <- read.delim(file.path(input.dirs[1], qname.file), header = FALSE,
                       stringsAsFactors = FALSE)[,1]
  taxids <- read.delim(file.path(input.dirs[1], taxids.file), header = FALSE,
                       stringsAsFactors = FALSE)[,1]
  names(taxids) <- qnames
  taxa.names.obj <- c(taxa.names.obj, list(taxids))
}
cat('\nDone.')

cat('\nMerging names that have been resolved to genus level ...')
resolved.names <- unresolved.names <- vector()
for (i in 1:length(taxa.names.obj)) {
  taxids <- taxa.names.obj[[i]]
  if (any(duplicated(taxids))) {
    utaxids <- unique(taxids)
    for (each in utaxids) {
      pull <- taxids %in% each
      if (sum(pull) > 1) {
        temp.unresolved.names <- names(taxids)[pull]
        temp.resolved.name <- strsplit(temp.unresolved.names[1], "\\s")[[1]]
        temp.resolved.name <- paste0(temp.resolved.name[1], "_merged")
        temp.resolved.names <- rep(temp.resolved.name, sum(pull))
        names(taxids)[pull] <- temp.resolved.names
        resolved.names <- c(resolved.names, temp.resolved.names)
        unresolved.names <- c(unresolved.names, temp.unresolved.names)
      }
    }
    taxa.names.obj[[i]] <- taxids[!duplicated(taxids)]
  }
}
name.modifiers <- list(unresolved.names, resolved.names)
save(name.modifiers, file = file.path(output.dir, "name_modifiers.RData"))
cat('\nDone.')

## Screen the phylogenies
cat("\n\nLooping through study phylogenies ...")
cat(paste0("\n... dropping all phylogenies with root to tip lengths greater than [",
       min.rtt.multiplier, "] times the median root to tip length ..."))
counter <- 0
nphylos.all <- 0
black.studies <- black.ids <- black.genes <- vector()
for (i in 1:length(phylo.studies)) {
  # empty vectors and numbers for study
  cat(paste0("\nWorking on [", phylo.studies[i], "] ..."))
  temp.black.ids <- temp.black.genes <- vector()
  temp.phylo.files <- list.files(input.dirs[2],
                                 pattern = paste0('^', phylo.studies[i]))
  temp.phylos <- rep(NA, length(temp.phylo.files))
  cat(paste0("\n... Found [", length(temp.phylo.files), "] phylogenies ..."))
  taxon <- names(taxa.names.obj[[i]])
  taxids <- taxa.names.obj[[i]]
  no.outgroup <- incorrect.names <- poor.phylo <- 0
  index1 <- regexpr("_gene_", temp.phylo.files) + 6
  index2 <- regexpr("._[0-9]*\\.tre$", temp.phylo.files)
  genes.used <- substring(temp.phylo.files, index1, index2)
  cat("\n...for genes: [", unique(genes.used), "] ...")
  ntips <- rep(NA, length(temp.phylo.files))
  # loop through phylos
  for (j in 1:length(temp.phylo.files)) {
    temp.phylo <- read.tree(file.path(input.dirs[2],temp.phylo.files[j]))
    ntips[j] <- length(temp.phylo$tip.label)
    if ("outgroup" %in% temp.phylo$tip.label){
      temp.phylo <- unroot(temp.phylo)
      temp.phylo <- root(temp.phylo, outgroup = "outgroup")
      temp.phylo <- drop.tip(temp.phylo, "outgroup")
      rtt.dists <- sort(diag(vcv.phylo(temp.phylo)))
      if (max(rtt.dists) < median(rtt.dists)*min.rtt.multiplier) {
        tip.labels <- as.numeric(sub("tx", "", temp.phylo$tip.label))
        tip.labels <- taxon[match(tip.labels, taxids)]
        if (length(tip.labels) == length(temp.phylo$tip.label)) {
          temp.phylo$tip.label <- tip.labels
          temp.phylos[j] <- list(temp.phylo)
        } else {
          incorrect.names <- incorrect.names + 1
        }
      } else {
        # add black ids
        temp.black.id <- names(rtt.dists)[rtt.dists > median(rtt.dists)*min.rtt.multiplier]
        temp.black.id <- as.numeric(sub("^tx", "", temp.black.id))
        temp.black.gene <- rep(genes.used[j], length(temp.black.id))
        temp.black.ids <- c(temp.black.ids, temp.black.id)
        temp.black.genes <- c(temp.black.genes, temp.black.gene)
        poor.phylo <- poor.phylo + 1
      }
    } else {
      no.outgroup <- no.outgroup + 1
    }
  }
  # Black listing:
  #  check if the proportion of trees dropped is greater than pblack for each gene
  #  check if it's worth keeping the gene and dropping certain ids or dropping the gene
  if (poor.phylo > 1) {
    phylos.by.gene <- table(genes.used)
    dphylos.by.gene <- tapply(as.numeric(is.na(temp.phylos)), genes.used, sum)
    ntips.by.gene <- tapply(ntips, genes.used, mean)
    bids.table <- table(temp.black.genes,temp.black.ids)
    nbids.by.gene <- rowSums(bids.table > 1)
    prop.dropped <- dphylos.by.gene/phylos.by.gene
    for (j in 1:length(prop.dropped)) {
      if (prop.dropped[j] > pblack){
        black.gene <- names(prop.dropped)[j]
        pull <- rownames(bids.table) == black.gene
        if (nbids.by.gene[pull]/ntips.by.gene[j] < 0.5) {
          # if the prop of black ids is less than 0.5, remove the black ids
          black.temp.ids.gene <- colnames(bids.table)[bids.table[pull,] > 1]
          black.ids <- c(black.ids, black.temp.ids.gene)
          black.genes <- c(black.genes, rep(black.gene, length(black.temp.ids.gene)))
          black.studies <- c(black.studies, rep(phylo.studies[i],
                                                length(black.temp.ids.gene)))
          cat("\n[", length(black.temp.ids.gene), "] ids black listed for gene [",
              black.gene, "] ...", sep = "")
          
        } else {
          # else remove the gene
          black.ids <- c(black.ids, "all")
          black.genes <- c(black.genes, black.gene)
          black.studies <- c(black.studies, phylo.studies[i])
        }
      }
    }
  }
  # write out trees
  temp.phylos <- temp.phylos[!is.na(temp.phylos)]
  nphylos <- length(temp.phylos)
  ndropped <- length(temp.phylo.files) - nphylos
  cat(paste0("\n[", ndropped, "] dropped ([", no.outgroup, "] no outgroup, [",
             incorrect.names, "] incorrect names, [", poor.phylo, "] poor phylogenies)."))
  if (nphylos > 1) {
    class(temp.phylos) <- 'multiPhylo'
  } else if (length(temp.phylos) == 1) {
    temp.phylos <- temp.phylos[[1]]
  } else {
    cat("\nAll phylogenies dropped.\nDropping study.")
    next
  }
  cat(paste0("\nKept [", length(temp.phylos),"] phylogenies."))
  write.tree(phy = temp.phylos, file = file.path(output.dir, paste0(
    phylo.studies[i], "_phylo.tre")), append = TRUE)
  counter <- counter + 1
  nphylos.all <- nphylos.all + nphylos
}

## Write out black list
cat("\n\nWriting out black list ...")
if (is.logical(black.genes)){
  cat("\nNo ids written to black list...")
} else {
  black.list <- data.frame(study = black.studies, genes = black.genes, taxids = black.ids)
  black.list <- black.list[!duplicated(black.list), ]
  write.csv(x = black.list, file = file.path("blacklist.csv"), row.names = FALSE)
  cat("\nA black list was created for [", nrow(black.list),"] genes and ids!", sep = "")
  cat("\nTo produce better phylogenies delete step folders and re-run steps alignment to screen with black list")
  cat("\n  in 0_names folder!! ")
}
cat("\n\nStage finished. [", nphylos.all, "] phylogenies across [", counter, "] studies.",
    sep = "")
