
##########################################################################
## Impact of life cycle variation on feeding system musculature in Caudata
## Morgane Taillades
## Phylogeny work
##########################################################################

### Libraries 
library(readxl)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(gridExtra)
library(ape)
library(geiger)

# Load the dataset
DataBrut <- read_excel("Mean_by_species.xlsx", col_names = TRUE)

# Load the phylogenetic tree 
phy <- read.tree("BigTree.tre")

# Keep only the species present in the dataset
phy_filtered <- drop.tip(phy,tip = setdiff(phy$tip.label, unique(DataBrut$Espece)))
phy <- phy_filtered
plot(phy)

edgelabels()
nodelabels()

# This function adds a new species (tip) to a phylogenetic tree,

tip_bind <- function(tree, sp_names = NULL, interactive = TRUE, position = 0,
                     where = "root", tipPos = NULL, rootPos = NULL) {
    require(ape)
  require(phytools)

  if (!is.null(rootPos)) {
    interactive <- FALSE
    where <- "root" }
  
  tip <- list(
    edge = matrix(c(2, 1), 1, 2),
    tip.label = sp_names,
    edge.length = 0.0,
    Nnode = 1)
  class(tip) <- "phylo"
 
   plot(tree)
  tree_temp <- bind.tree(tree, tip,
                         interactive = interactive,
                         position = position,
                         where = where)
  
  node_heights <- nodeHeights(tree_temp)
  tip_index <- which(tree_temp$edge[, 2] == which(tree_temp$tip.label == sp_names))
  
  diff_length <- max(node_heights) - node_heights[tip_index, 2]
  
  if (!is.null(tipPos)) {
    diff_length <- diff_length - tipPos
  }
  
  tree_temp$edge.length[tip_index] <-
    tree_temp$edge.length[tip_index] + diff_length

  if (!is.null(rootPos)) {
    tree_temp <- root(tree_temp,
                      outgroup = sp_names,
                      resolve.root = TRUE)
    
    index_root <- which(tree_temp$edge[, 1] == Ntip(tree_temp) + 1)
    tree_temp$edge.length[index_root] <-
      tree_temp$edge.length[index_root] + rootPos
  }
  
  tree_temp <- untangle(tree_temp)
  
  return(tree_temp)
}

# Add metamorphic Ambystoma mexicanum
phy_modifie <- tip_bind(phy,sp_names = "Ambystoma_mexicanum_metamorph", interactive = FALSE,
  position = 0.5, where = 6)

# Add metamorphic Ambystoma andersoni
phy_modifie <- tip_bind(phy_modifie,sp_names = "Ambystoma_andersoni_metamorph", interactive = FALSE,
  position = 0.5, where = 5)

plot(phy_modifie)

#Export the new phylogeny
write.nexus(phy_modifie, file = "PhySalamanders.nex")
