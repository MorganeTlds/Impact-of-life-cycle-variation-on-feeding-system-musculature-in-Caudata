
##########################################################################
## Impact of life cycle variation on feeding system musculature in Caudata
## Morgane Taillades
## PCA plot
##########################################################################

### Libraries 
library(ggplot2)
library(ggrepel)

# Generic function to produce a phylomorphospace plot
plot_phylomorphospace <- function(scores, df, title_text) {
  
  # Data preparation
  components <- data.frame(
    PC1 = scores[, 1],
    PC2 = scores[, 2],
    Adult_habitat = df$Adult_habitat,
    T_eco = df$T_eco,
    T_morpho = df$T_morpho
  )
  row.names(components)<-row.names(scores)
  
  # Convex hull polygons for each Adult_habitat group
  splitData <- split(components[, c("PC1", "PC2", "Adult_habitat")],
                     components$Adult_habitat)
  
  polyData <- do.call(rbind, lapply(splitData, function(x) x[chull(x$PC1, x$PC2), ]))
  
  # Plot
  ggplot(components, aes(x = PC1, y = PC2)) +
    geom_polygon(
      data = polyData,
      aes(fill = Adult_habitat, group = Adult_habitat),
      alpha = 0.3
    ) +
    geom_point(aes(color = T_eco, shape = T_morpho), size = 3) +
    geom_text_repel(
      aes(label = row.names(components)),
      size = 4.5, color = "black",
      max.overlaps = Inf, force = 2, force_pull = 0
    ) +
    scale_fill_manual(values = c("A" = "blue", "SA" = "cyan", "T" = "red")) +
    theme_classic() +
    labs(
      title = title_text,
      x = paste0("PC1"),
      y = paste0("PC2")
    )
}

# --- Run the two plots ---


plotPCAVol <- plot_phylomorphospace(
  scores = as.data.frame(ACP_Vol$scores),
  df = df,
  title_text = "Phylomorphospace of PCA muscle volume")

plotPCApcsa <- plot_phylomorphospace(
  scores = as.data.frame(ACP_pcsa$scores),
  df = df,
  title_text = "Phylomorphospace of PCA muscle PCSA")


# Display

plot(plotPCAVol)
plot(plotPCApcsa)




