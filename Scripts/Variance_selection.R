#gene_mapping = read.table("~/Koop_Domaszewska/Misc/Gene_mapping_Illumina.tsv", sep ="\t", header = T, stringsAsFactors = F)
### hgnc

####

hgnc_list = rf_meta$genes$hgnc_symbol[  rf_meta$genes$ensembl_gene_id == rownames(expr_raw)]
hgnc_list_uni = unique(hgnc_list)

expr_raw = expr_raw[ ! is.na(hgnc_list),]
hgnc_list = hgnc_list[!is.na(hgnc_list)]
expr_raw = expr_raw[ hgnc_list != "",]
hgnc_list = hgnc_list[ hgnc_list != ""]
hgnc_list_uni = as.character( unique(hgnc_list) )

max_list = as.integer( sapply( hgnc_list_uni, FUN = function(gene){
  
  var_match = which( hgnc_list == gene )
  if (length(var_match) > 1){
    row_var_max = which.max( apply( as.matrix(expr_raw[var_match,]), FUN = var, MARGIN = 1)  );
  } else{
    row_var_max = 1
  }
  return( var_match[row_var_max] )
}))

expr_raw = expr_raw[max_list,]
rownames(expr_raw) = hgnc_list[max_list]
dim(expr_raw)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X", "" )
expr_raw[1:5,1:5]

#write.table(expr_raw,"~/Koop_Domaszewska/Data/Illumina_merged.HGNC.tsv",sep="\t",quote = F)
