library("stringr")
library("beadarray")

# load hard

i_files = list.files( "~/Koop_Domaszewska/Data/GSE37250_RAW/", pattern = ".txt", full.names = T)
data = readIllumina(
  dir = "~/Koop_Domaszewska/Data/GSE37250_RAW/",
  #useImages = F,
  illuminaAnnotation = "Humanv4"
)

# init

load("~/Koop_Domaszewska/Data/MDSraw.RDa")
raw_data = MDSraw$E
raw_data[1:5,1:5]

# sample wise scaling

scale_mat = apply(raw_data, MARGIN = 2, FUN = function(vec){return(scale(vec))})
rownames(scale_mat) = rownames(raw_data)
scale_mat[1:5,1:5]

###

cutoff_20 = as.integer(quantile(rowSums(raw_data), seq(0,1,by=.1))[3])
raw_data_flt = raw_data[ rowSums(raw_data) >= cutoff_20,]
dim(raw_data_flt)
min(as.double(apply(raw_data, MARGIN = 2, FUN = min)))

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
pdf( "~/Koop_Domaszewska/Results/QA/Hist_ColSums_raw_300.pdf" , onefile=FALSE)
  hist(colsums)
dev.off()

neg_val = which(as.double(apply(raw_data, MARGIN = 2, FUN = min)) < 0)
agg_exp = aggregate( raw_data_flt, by = list(meta_data$study), FUN = rowSums)
gather_matrix <- reshape2:::melt.matrix(as.matrix(raw_data_flt))
gather_matrix[1:5,]





