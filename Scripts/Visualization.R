library("Rtsne") # Load package

#meta_info = read.table("~/Koop_Domaszewska/Misc/Meta_Information.tsv", sep ="\t", header =  T, stringsAsFactors = F)
#meta_info = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW/Annotation.tsv", sep ="\t", header =  T, stringsAsFactors = F)
#meta_data = meta_info[ match(colnames(bam_counts), meta_info$Sample, nomatch = 0), ]
#rownames(meta_data) = meta_data$Sample

load("~/Koop_Domaszewska/Data/myIFN_I_set.RDa")

myIFN_I_set$MODULES
sad_genes = myIFN_I_set["LI.M75_I"]$GENES$ID # CXL10
sad_genes = myIFN_I_set["DC.M5.12_I"]$GENES$ID # ACTA2
sad_genes = myIFN_I_set["DC.M1.2_I"]$GENES$ID # OTOF
sad_genes = myIFN_I_set["DC.M3.4_I"]$GENES$ID
sad_genes = myIFN_I_set["ifn_I"]$GENES$ID

# illumina_humanht_12_v4
#Gene_mapping = getBM( attributes=c("hgnc_symbol","illumina_humanht_12_v4"),filters="illumina_humanht_12_v4",values = rownames(bam_counts), mart = ensembl)

#sad_genes
expr = expr_raw[rownames(expr_raw) %in% sad_genes,]

### VIZ
tsne_out <- Rtsne( t(expr), perplexity = 5) # Run TSNE
#plot( tsne_out$Y, col = meta_data$characteristics_ch1.1)
plot( tsne_out$Y, col = meta_data$characteristics_ch1.5)

cor_mat = cor(expr)
#pdf( "~/Koop_Domaszewska/Results/QA/Heatmap_raw_300.pdf" , onefile=FALSE)
pheatmap::pheatmap(
  cor_mat,
  #annotation_col = meta_data[c("characteristics_ch1","characteristics_ch1.1","characteristics_ch1.2")],
  annotation_col = meta_data[c("characteristics_ch1.3","characteristics_ch1.5")],
  show_rownames = F,
  show_colnames = F#,
  #annotation_colors = aka3,
  #color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
)
#dev.off()

pcr = prcomp(t(cor_mat))
#pdf( "~/Koop_Domaszewska/Results/QA/PCA_raw_300.pdf" , onefile=FALSE)
ggbiplot::ggbiplot(
  pcr,
  obs.scale = 1,
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = meta_data$characteristics_ch1.5,
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)
#dev.off()

###

cor_mat = cor(data_matrix)
pheatmap::pheatmap(
  cor_mat,
  show_colnames = F,
  show_rownames = F
)