library("stringr")

raw_files = list.files("~/Koop_Domaszewska/Data/GSE42834_RAW//", pattern = ".txt", full.names = T)
raw_files_short = list.files("~/Koop_Domaszewska/Data/GSE42834_RAW//", pattern = ".txt", full.names = F)
#annot_file = read.table( "~/Koop_Domaszewska/Data/GSE19491_RAW/GPL6947_HumanHT-12_V3_0_R1_11283641_A.bgx", sep ="\t", header = T, stringsAsFactors = F , skip = 8, fill = T)
#annot_file[1,]

bam_counts <<- read.table(raw_files[1], sep ="\t", header = T)[,1]
#48803
for ( i in 1:length(raw_files)){
    print(i)
    #input = read.table(raw_files[i], sep ="\t", header = T, colClasses = c("character","numeric","NULL"))
    input = read.table(raw_files[i], sep ="\t", header = T)[,-1]
    
    na_detect = apply(input, MARGIN = 2, FUN=function(vec){return(sum(is.na(vec))>0)})
    input = input[,!na_detect]
    input = input[,! str_detect(colnames(input),pattern = "Detection")]
    
    bam_counts <<- cbind(bam_counts, input)
    print(ncol(input))
}

mapping_table = read.table("~/Koop_Domaszewska/Data/GSE42834_RAW/mapping_table.tsv",sep ="\t", header = F, stringsAsFactors = F)
colnames(bam_counts) = c("Probe_ID", mapping_table[,1])

dim(bam_counts)
colnames(bam_counts) = sapply(colnames(bam_counts), FUN = function(vec){return(
    head(as.character(unlist(str_split(vec, pattern = "_"))),1)
)})
colnames(bam_counts)[1] = "Probe_ID"
bam_counts[1:5,1:5]

write.table(bam_counts,"~/Koop_Domaszewska/Data/GSE42834_RAW//Merged.txt", sep ="\t", quote = F , row.names = F)
