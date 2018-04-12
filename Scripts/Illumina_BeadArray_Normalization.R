library(knitr)
library("beadarray")
#library(illuminaHumanv3.db)

expr_raw = read.table("~/Koop_Domaszewska/Data/Illumina_merged.HGNC.non_normalized.tsv",sep ="\t", stringsAsFactors = F, header = T)
expr_raw_new = matrix(as.double(as.character(unlist(expr_raw))), ncol=ncol(expr_raw), nrow = nrow(expr_raw))
colnames(expr_raw_new) = colnames(expr_raw)
rownames(expr_raw_new) = rownames(expr_raw)

expr_raw = expr_raw_new

eset_m = new(
  Class ="ExpressionSet",
  exprs = expr_raw
)

eset_n = beadarray::normaliseIllumina(eset_m)
colnames(eset_n) = colnames(expr_raw)
rownames(eset_n) = rownames(expr_raw)
exprs(eset_n)[1:5,1:5]
expr_raw = exprs(eset_n)
expr_raw = scale(expr_raw)

design <- model.matrix(~0 + meta_data$Study)

DGE = edgeR::DGEList(expr_raw)
DGE = edgeR::calcNormFactors(DGE,method =c("TMM"))
expr_raw = DGE$counts
v = limma::voom(DGE,design)

expr_raw = v$E
expr_raw = scale(expr_raw)
expr_raw[1:5,1:5]

### sva combat

expr_raw = sva::ComBat(expr_raw, batch = meta_data$Study)
boxplot(expr_raw[,sample(1:ncol(expr), 100)])
