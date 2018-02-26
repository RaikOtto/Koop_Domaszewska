library("GEOquery")

gse_37250 <- getGEO("GSE37250",GSEMatrix=TRUE) # Kaforou_2013
gse_73408 <- getGEO("GSE73408",GSEMatrix=TRUE) # Walter
gse_39941 <- getGEO("GSE39941",GSEMatrix=TRUE) # Anderson
gse_42834 <- getGEO("GSE42834",GSEMatrix=TRUE) # Bloom
gse_34608 <- getGEO("GSE34608",GSEMatrix=TRUE) # Maerzdorf_2012
gse_28623 <- getGEO("GSE28623",GSEMatrix=TRUE) # Maerzdorf_2011
gse_19491 <- getGEO("GSE19491",GSEMatrix=TRUE) # Kaforou_2011

exprs(gse$GSE73408_series_matrix.txt.gz)

head(pData(phenoData(gse[[1]])))
