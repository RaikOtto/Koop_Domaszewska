library("GEOquery")

gse_37250 <- getGEO("GSE37250",GSEMatrix=T) # Kaforou_2013 Illumina HumanHT-12 V4.0 expression beadchip
gse_39941 <- getGEO("GSE39941",GSEMatrix=T) # Anderson     Illumina HumanHT-12 V4.0 expression beadchip
gse_42834 <- getGEO("GSE42834",GSEMatrix=T) # Bloom        Illumina HumanHT-12 V4.0 expression beadchip
gse_19491 <- getGEO("GSE19491",GSEMatrix=T) # Kaforou_2011 Illumina HumanHT-12 V3.0 expression beadchip
gse_34608 <- getGEO("GSE34608",GSEMatrix=T) # Maerzdorf_2012 Agilent-014850 Whole Human Genome Microarray 4x44K G4112F
gse_28623 <- getGEO("GSE28623",GSEMatrix=T) # Maerzdorf_2011 Agilent-014850 Whole Human Genome Microarray 4x44K G4112F
gse_73408 <- getGEO("GSE73408",GSEMatrix=T) # Walter [HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array [transcript (gene) version]

#mat_37250 = exprs(gse_37250$)
mat_gse_19491 = exprs(gse_19491$GSE19491_series_matrix.txt.gz)

t_mat = mat_gse_19491[1:5,1:5] # test_mat
apply(t_mat, FUN= function(vec){return(2**vec)}, MARGIN = 2)

hist(mat_gse_19491)

pheatmap::pheatmap(cor(mat_gse_19491))

p_data_37250 = pData(gse_37250[[1]])
p_data_37250$data_processing
p_data_39941 = pData(gse_39941[[1]])
p_data_39941$data_processing
colnames(pData(gse_42834[[1]]))
colnames(pData(gse_19491[[1]]))
colnames(pData(gse_34608[[1]]))
colnames(pData(gse_28623[[1]]))
colnames(pData(gse_73408[[1]]))

p_data_19491 = pData(gse_19491[[1]])
p_data_19491$data_processing
