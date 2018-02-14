library("stringr")
library("limma")

load("~/Koop_Domaszewska/MDSraw.RDa")

# init

raw_data = MDSraw$E
cutoff_20 = as.integer(quantile(rowSums(raw_data), seq(0,1,by=.1))[3])
raw_data_flt = raw_data[ rowSums(raw_data) >= cutoff_20,]
dim(raw_data_flt)

neg_val = which(as.double(apply(raw_data, MARGIN = 2, FUN = min)) < 0)
MDSraw$targets$study[ match( colnames(raw_data[,neg_val]), rownames(MDSraw$targets)) ]
table(MDSraw$targets$study[ match( colnames(raw_data[,neg_val]), rownames(MDSraw$targets)) ])

# create stichprobe

meta_match = match(colnames(raw_data_flt), rownames(MDSraw$targets), nomatch = 0)
r = sample(1:ncol(raw_data_flt),300,replace = F)
meta_match[-r] = 0

# create subsampled datasets

meta_data = as.data.frame( MDSraw$targets[ meta_match ,] )
raw_data_flt = raw_data_flt[ , meta_match != 0]
cor_mat = cor(raw_data_flt)
dim(meta_data)
dim(raw_data_flt)
dim(cor_mat)

cbind(rownames(meta_data), colnames(cor_mat)) # <- check if ordered correctly

colsums = colSums(raw_data_flt)


neg_val = which(as.double(apply(raw_data, MARGIN = 2, FUN = min)) < 0)
#agg_exp = aggregate( raw_data_flt, by = list(meta_data$study), FUN = rowSums)

#abs_plot = ggplot( data = , aes( Library, log2( Count ) ) )
#abs_plot = abs_plot + geom_boxplot(aes(fill = Library))

pdf( "~/Koop_Domaszewska/Results/QA/Heatmap_raw_300.pdf" , onefile=FALSE)
  pheatmap::pheatmap(
    cor_mat,
    annotation_col = meta_data[c("study")],
    show_rownames = F,
    show_colnames = F#,
    #annotation_colors = aka3,
    #color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
  )
dev.off()

pcr = prcomp(t(cor_mat))
pdf( "~/Koop_Domaszewska/Results/QA/PCA_raw_300.pdf" , onefile=FALSE)
  ggbiplot::ggbiplot(
    pcr,
    obs.scale = 1,
    var.scale = 1, 
    labels.size = 4,
    alpha = 1,
    groups = meta_data$study,
    ellipse = TRUE, 
    circle = TRUE,
    var.axes = F
  )
dev.off()

sort(table(MDSraw$targets$study), decreasing = T)
cumsum(sort(table(MDSraw$targets$study), decreasing = F)) / sum(sort(table(MDSraw$targets$study), decreasing = F)) * 100
# losing 47% of samples

passing_studies = c("Anderson","kaforou")
meta_data_filt = meta_data[meta_data$study %in% passing_studies,]
raw_filt = raw_data[ , colnames( raw_data) %in% rownames(meta_data_filt) ]
dim(raw_filt)

cor_mat = cor(raw_filt)
pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data_filt[c("study")],
  show_rownames = F,
  show_colnames = F#,
  #annotation_colors = aka3,
  #color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
)

ggbiplot::ggbiplot(
  pcr,
  obs.scale = 1,
  var.scale = 1, 
  labels.size = 4,
  alpha = 1,
  groups = MDSraw$targets$study[r],
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)