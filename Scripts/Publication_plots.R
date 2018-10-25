library("grid")
library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

load("~/Koop_Domaszewska/Data/trainingMDS.RDa")
expr_raw = trainingMDS$E
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "")
source("~/Koop_Domaszewska/Scripts/Variance_selection.R")

dim(expr_raw)
expr_raw = expr_raw[,sample(size = 300, x = 1:ncol(expr_raw))]
expr_raw[1:5,1:5]

### Prep

meta_info = trainingMDS$targets
rownames(meta_info) = meta_info$donor
meta_data = meta_info[ colnames(expr_raw), ]
table(is.na(meta_data))

##
aka3 = list(
  Study   = c(
    Anderson = "BLACK",
    Berry = "Orange",
    Berry_2010 = "Yellow",
    Bloom = "purple",
    Kaforou_2013 = "Darkgreen"),
  study   = c(
    GSE19491 = "BLACK",
    GSE28623 = "Orange",
    GSE34608 = "Yellow",
    GSE39941 = "purple",
    GSE42834 = "Darkgreen",
    GSE47673 = "Brown",
    GSE73408 = "gray"
    ),
  activeTB = c(
      no = "green",
      TB = "red"
  ),
  IFN_I = c(
    low = "green",
    medium = "yellow",
    high = "red"
  )
)
###

load("~/Koop_Domaszewska/Data/myIFN_I_set.RDa")
sad_genes = as.character(unlist(myIFN_I_set$MODULES2GENES[1]))
sad_genes = sad_genes[ sad_genes != ""]
expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
expr[1:5,1:5]
dim(expr)
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

## Figure 1

# Plot 1

pheatmap::pheatmap(
    cor_mat,
    annotation_col = meta_data[c("activeTB","IFN_I","study")],
    annotation_colors = aka3,
    show_rownames = F,
    show_colnames = F,
    treeheight_col = 0,
    treeheight_row = 0,
    legend = F,
    fontsize_col = 7,
    clustering_method = "ward.D2"
)

## Figure 2

# Plot 1

#meta_data$Deco_type = factor( meta_data$Deco_type, levels = c("Alpha","Beta","Gamma","Delta", "Ductal", "Acinar") )

ggbiplot::ggbiplot(
  pcr,
  choices = c(1,2),
  obs.scale = .75,
  groups = as.character(meta_data[colnames(expr),"activeTB"]),
  #shape = as.character(meta_data$Grading),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F
)
p

library("umap")

### umap

sad_genes = as.character(unlist(myIFN_I_set$MODULES2GENES[5]))
sad_genes = sad_genes[ sad_genes != ""]
expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]

umap_plot = umap::umap(t(expr))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")
dist_mat = dist((vis_data))
p = ggplot2::qplot(
    x = vis_data$x,
    y = vis_data$y,
    color = meta_data[colnames(expr),"activeTB"],
    geom=c("point", "smooth"),
    xlab = "Umap dim 1",
    ylab = "Umap dim 2"
) + guides(color=guide_legend(title="Active TB"), position = "top")
p

