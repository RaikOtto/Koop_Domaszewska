library("stringr")

meta_file = read.table("~/Koop_Domaszewska/Misc/Meta_Information.tsv",sep ="\t", header = T)

###

bam_counts_gse_37250 = read.table("~/Koop_Domaszewska/Data/GSE37250_RAW/GSE37250_non-normalized.txt",sep ="\t", header = T,  fill = T)
rownames(bam_counts_gse_37250) = bam_counts_gse_37250[,1]
bam_counts_gse_37250 = bam_counts_gse_37250[, str_detect(colnames(bam_counts_gse_37250),pattern = "AVG_Signal")]
colnames(bam_counts_gse_37250) = sapply(colnames(bam_counts_gse_37250), FUN = function(vec){return(
  head(as.character(unlist(str_split(vec, pattern = "_"))),1)
)})

na_detect = apply(bam_counts_gse_37250, MARGIN = 1, FUN=function(vec){return(sum(is.na(vec))>0)})
dim(bam_counts_gse_37250)
bam_counts_gse_37250 = bam_counts_gse_37250[!na_detect,]
colnames(bam_counts_gse_37250) = read.table("~/Koop_Domaszewska/Data/GSE37250_RAW/Mapping_table.tsv",sep ="\t", header = F, stringsAsFactors = F)[,1]
bam_counts_gse_37250[1:5,1:5]

bam_counts_gse_39939 = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW/GSE39939_non-normalized_data.txt",sep ="\t", header = T, fill = T)
rownames(bam_counts_gse_39939) = bam_counts_gse_39939[,1]
bam_counts_gse_39939 = bam_counts_gse_39939[, str_detect(colnames(bam_counts_gse_39939),pattern = "AVG_Signal")]
colnames(bam_counts_gse_39939) = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW//Mapping_Table39.tsv",sep ="\t", header = F, stringsAsFactors = F)[,1]
bam_counts_gse_39939[1:5,1:5]

bam_counts_gse_39940 = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW/GSE39940_GSE39938_Sample_Probe_Profile_with_no_BG_subtracted_and_not_normalised_GEO_upload_CT_MLW_07_08_12.txt",sep ="\t", header = T, row.names = 1, fill = T)
bam_counts_gse_39940 = bam_counts_gse_39940[,! str_detect(colnames(bam_counts_gse_39940),pattern = "\\.Pval")]
colnames(bam_counts_gse_39940) = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW//Mapping_Table40.tsv",sep ="\t", header = F, stringsAsFactors = F)[,1]
bam_counts_gse_39940[1:5,1:5]

bam_counts_gse_19491 = read.table("~/Koop_Domaszewska/Data/GSE19491_RAW/Merged.txt",sep ="\t", header = T, row.names = 1, fill = T)
bam_counts_gse_19491 = bam_counts_gse_19491[,! str_detect(colnames(bam_counts_gse_19491),pattern = "\\.Pval")]
colnames(bam_counts_gse_19491) = read.table("~/Koop_Domaszewska/Data/GSE19491_RAW//Mapping_table.tsv",sep ="\t", header = F, stringsAsFactors = F)[,1]
bam_counts_gse_19491[1:5,1:5]

bam_counts_gse_42834 = read.table("~/Koop_Domaszewska/Data/GSE42834_RAW/Merged.txt",sep ="\t", header = T, row.names = 1, fill = T)
bam_counts_gse_42834[1:5,1:5]

merged_probes = Reduce( intersect, list(
  rownames(bam_counts_gse_37250),
  rownames(bam_counts_gse_39939),
  rownames(bam_counts_gse_39940),
  rownames(bam_counts_gse_19491),
  rownames(bam_counts_gse_42834)
) )

bam_counts_gse_37250 = bam_counts_gse_37250[ rownames(bam_counts_gse_37250) %in% merged_probes, ]
bam_counts_gse_39939 = bam_counts_gse_39939[ rownames(bam_counts_gse_39939) %in% merged_probes, ]
bam_counts_gse_39940 = bam_counts_gse_39940[ rownames(bam_counts_gse_39940) %in% merged_probes, ]
bam_counts_gse_19491 = bam_counts_gse_19491[ rownames(bam_counts_gse_19491) %in% merged_probes, ]
bam_counts_gse_42834 = bam_counts_gse_42834[ rownames(bam_counts_gse_42834) %in% merged_probes, ]

new_mat = merged_probes
new_mat = cbind(new_mat,bam_counts_gse_37250[match(rownames(bam_counts_gse_37250), merged_probes, nomatch = 0),])
new_mat = cbind(new_mat,bam_counts_gse_39939[match(rownames(bam_counts_gse_39939), merged_probes, nomatch = 0),])
new_mat = cbind(new_mat,bam_counts_gse_39940[match(rownames(bam_counts_gse_39940), merged_probes, nomatch = 0),])
new_mat = cbind(new_mat,bam_counts_gse_19491[match(rownames(bam_counts_gse_19491), merged_probes, nomatch = 0),])
new_mat = cbind(new_mat,bam_counts_gse_42834[match(rownames(bam_counts_gse_42834), merged_probes, nomatch = 0),])
new_mat = new_mat[,-1]
new_mat = matrix( as.double(as.character(unlist(new_mat))), ncol = ncol(new_mat), nrow = nrow(new_mat) )
new_mat[1:5,1:5]

write.table(new_mat, "~/Koop_Domaszewska/Data/Illumina_merged.naked.tsv",sep="\t",quote =F, row.names = F)
