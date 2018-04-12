library("Rtsne") # Load package

load("~/Koop_Domaszewska/Data/MDS.RDa")
meta_info = read.table("~/Koop_Domaszewska/Misc/Meta_Information.tsv", sep ="\t", header =  T, stringsAsFactors = F)

expr_raw = read.table("~/Koop_Domaszewska/Data/Illumina_merged.HGNC.non_normalized.tsv",sep ="\t", stringsAsFactors = F, header = T)

match = match( colnames(expr_raw), rf_meta$targets$donor, nomatch = 0 )

expr_raw = expr_raw[,match != 0]
meta_data = meta_info[ match( colnames(expr_raw), meta_info$Sample,  nomatch = 0),]
rownames(meta_data) = meta_data$Sample
meta_data = cbind(meta_data, rf_meta$targets[match(rownames(meta_data), rf_meta$targets$donor, nomatch = 0),])

load("~/Koop_Domaszewska/Data/myIFN_I_set.RDa")

sad_genes = myIFN_I_set["LI.M75_I"]$GENES$ID # CXL10
sad_genes = myIFN_I_set["DC.M5.12_I"]$GENES$ID # ACTA2
sad_genes = myIFN_I_set["DC.M1.2_I"]$GENES$ID # OTOF
sad_genes = myIFN_I_set["DC.M3.4_I"]$GENES$ID
sad_genes = myIFN_I_set["ifn_I"]$GENES$ID

# illumina_humanht_12_v4
#Gene_mapping = getBM( attributes=c("hgnc_symbol","illumina_humanht_12_v4"),filters="illumina_humanht_12_v4",values = rownames(bam_counts), mart = ensembl)

#sad_genes

expr = expr_raw[rownames(expr_raw) %in% sad_genes,]
#expr = expr_raw
row_var = apply(expr, FUN = var, MARGIN = 1)
col_var = apply(expr, FUN = var, MARGIN = 2)
expr = expr[row_var > 0, col_var > 0]
cor_mat = cor( expr[,sample(1:ncol(expr), 300)])

pheatmap::pheatmap(
  cor_mat,
  show_rownames = F,
  show_colnames = F,
  annotation_col = meta_data[c("Study","activeTB","HIV","rf_nonTB","rf","NK","IFN")],
  #annotation_col = meta_data[c("Study","LTB","healthy","OD","status","IFN_I","IFN_II","IFN_I_TB")],
  treeheight_row = 0,
  treeheight_col = 0
)

### VIZ
#tsne_out <- Rtsne( t(expr), perplexity = 5) # Run TSNE
#plot( tsne_out$Y, col = meta_data$characteristics_ch1.1)
#plot( tsne_out$Y, col = meta_data$characteristics_ch1.5)

#pdf( "~/Koop_Domaszewska/Results/QA/Heatmap_raw_300.pdf" , onefile=FALSE)
#dev.off()

pcr = prcomp(t(cor_mat))
#pdf( "~/Koop_Domaszewska/Results/QA/PCA_raw_300.pdf" , onefile=FALSE)
ggbiplot::ggbiplot(
  pcr,
  obs.scale = 1,
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = meta_data$Study[match( colnames(cor_mat), meta_data$donor, nomatch = 0) ],
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