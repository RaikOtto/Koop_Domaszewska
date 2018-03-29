library("stringr")
library("GEOquery")

gse_37250 <- getGEO("GSE37250",GSEMatrix=T) # Kaforou_2013 Illumina HumanHT-12 V4.0 expression beadchip

gse_39939 <- getGEO("GSE39939",GSEMatrix = T) # Anderson     Illumina HumanHT-12 V4.0 expression beadchip
gse_39940 <- getGEO("GSE39940",GSEMatrix = T) # Anderson     Illumina HumanHT-12 V4.0 expression beadchip
gse_39941 <- getGEO("GSE39941",GSEMatrix = T) # Anderson     Illumina HumanHT-12 V4.0 expression beadchip

gse_42834 <- getGEO("GSE42834",GSEMatrix=T) # Bloom        Illumina HumanHT-12 V4.0 expression beadchip

gse_19491 <- getGEO("GSE19491",GSEMatrix=T) # Berry Illumina HumanHT-12 V3.0 expression beadchip

gse_34608 <- getGEO("GSE34608",GSEMatrix=T) # Maerzdorf_2012 Agilent-014850 Whole Human Genome Microarray 4x44K G4112F
gse_28623 <- getGEO("GSE28623",GSEMatrix=T) # Maerzdorf_2011 Agilent-014850 Whole Human Genome Microarray 4x44K G4112F
gse_73408 <- getGEO("GSE73408",GSEMatrix=T) # Walter [HuGene-1_1-st] Affymetrix Human Gene 1.1 ST Array [transcript (gene) version]


###

bam_counts = read.table("~/Koop_Domaszewska/Data/GSE39941_RAW/GSE39939_non-normalized_data.txt",sep ="\t", header = T, row.names = 1)
colnames(bam_counts) = str_replace(colnames(bam_counts), pattern = "^X", "")
bam_counts = bam_counts[,seq(1,ncol(bam_counts), by = 2)]
dim(bam_counts)
bam_counts[1:5,1:5]

meta_data = meta_data[match(
  str_replace( colnames(bam_counts), pattern = "\\.AVG_Signal",""), 
  meta_data$`barcode:ch1`
),]

#filePaths = getGEOSuppFiles("gse_39941")

meta_data = meta_data[match(
  str_replace( colnames(bam_counts), pattern = "\\.AVG_Signal",""), 
  meta_data$`barcode:ch1`
),]

#filePaths = getGEOSuppFiles("GSE37250")

row_var = apply( bam_counts, FUN = var, MARGIN = 1)
table(row_var <= 1)
bam_counts = bam_counts[row_var >= 1,]
dim(bam_counts)

###

bam_counts = exprs(gse_39941$GSE39941_series_matrix.txt.gz)
meta_data = as.data.frame( pData(gse_39941[[1]]) )

bam_counts = exprs(gse_37250$GSE37250_series_matrix.txt.gz)
meta_data = as.data.frame( pData(gse_37250[[1]]) )
  
bam_counts = exprs(gse_42834$GSE42834_series_matrix.txt.gz)
meta_data = as.data.frame( pData(gse_42834[[1]]) )

bam_counts = exprs(gse_19491$GSE19491_series_matrix.txt.gz)
meta_data = as.data.frame( pData(gse_19491[[1]]) )

colnames(bam_counts) = str_replace(colnames(bam_counts), pattern = "^X", "")
