### Anolis data processing
lib <- c("ape", "phytools", "phyloch","bite", "devtools")
sapply(lib, library, character.only = T)

tree <- read.nexus("data/data-raw/Anolis_tree.nex")
tree$tip.label <- sprintf("A._%s", tree$tip.label)
data <- read.table("data/data-raw/Anolis.txt", header = T, sep = "\t", stringsAsFactors = F)
data[,1] <- gsub("A.", "A._", data[,1])

Anolis_tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% unique(data[,1])])

Anolis_map <- matrix(0, nrow = nrow(Anolis_tree$edge), ncol = 2, dimnames = list(apply(Anolis_tree$edge, 1, paste0, collapse = ","), c("Hispaniola", "Cuba")))
for(i in 1:nrow(Anolis_tree$edge)){
  island <- unique(data[data[,1] %in% Anolis_tree$tip.label[descendants(Anolis_tree, Anolis_tree$edge[i,2])],2])
  Anolis_map[i,island == c("cybotes", "sagrei")] <- 1
}
Anolis_map <- Anolis_map * Anolis_tree$edge.length

Anolis_traits <- data[,c(1,7,6)]

use_data(Anolis_tree, Anolis_map, Anolis_traits, overwrite = T)

