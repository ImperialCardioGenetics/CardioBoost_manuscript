
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
input_file=args[1]
output_file=args[2]

#input file is output from annovar annotation: "XX_variants.vcf.hg19_multianno.txt"
file=input_file
NAMES <- read.table(file, nrow = 1, stringsAsFactors = FALSE, sep = "\t")
DATA <- read.table(file, skip = 1, stringsAsFactors = FALSE, sep = "\t")
#NAMES[,which(NAMES=="fathmm-MKL_coding_score")]="fathmm_MKL_coding_score"
#NAMES[,which(NAMES=="GERP++_RS")]="GERP"
NAMES<-NAMES[,1:length(NAMES)-1]
DATA <- DATA[, 1:length(NAMES)]
names(DATA) <- NAMES


write.table(DATA,output_file,quote=FALSE,sep="\t",col.names = T,row.names = F)
