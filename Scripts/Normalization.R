
#eset = threestep( raw_data, background.method = "RMA.2", normalize.method="quantile", summary.method="median.polish")
#eset = rma( raw_data, target='extended', normalize = T)
#eset = rma( raw_data, target='core', normalize = T)

bam_counts = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW/GSE39941_series_matrix.txt",sep ="\t", header = T, row.names = 1)
rn = rownames(bam_counts)
bam_counts = apply( bam_counts, FUN = as.double, MARGIN = 2)
rownames(bam_counts) = rn
bam_counts[1:5,1:5]

dim(bam_counts)

row_var = apply( bam_counts, FUN = var, MARGIN = 1 )
col_var = apply( bam_counts, FUN = var, MARGIN = 2 )
table(row_var <= 1)
table(col_var <= 1)
