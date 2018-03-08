library("stringr")

raw_files = list.files("~/Koop_Domaszewska/Data/GSE19491_RAW/", pattern = ".txt", full.names = T)
raw_files_short = list.files("~/Koop_Domaszewska/Data/GSE19491_RAW/", pattern = ".txt", full.names = F)

bam_counts <<- matrix(as.integer(),ncol = 1, nrow = 48803)
#48803
for ( i in 1:length(raw_files)){
    print(i)
    input = read.table(raw_files[i], sep ="\t", header = T, row.names = 1, colClasses = c("character","double","NULL"))
    colnames(input) = raw_files_short[i]
    bam_counts <<- cbind(bam_counts, input)
}

bam_counts = bam_counts[,-1]
dim(bam_counts)
bam_counts[1:5,1:5]

t = read.table("/home/ottoraik/Koop_Domaszewska/Data/GSE19491_RAW//GSM484599_LTB_LON_test009.txt", sep ="\t", header = T, row.names = 1)
