library('tidyverse')
library("gplots")
library('pheatmap')

#################### Heat map (3-subtypes) #################### 
snv_summary_w_length
columns_of_interest <- c("Sample_ID","NRAS","KRAS","HRAS","BRAF","NF1")
snv_test <- as.data.frame(snv_summary_count[,columns_of_interest])
rownames(snv_test) <- snv_test$Sample_ID
snv_test %>% 
  as.matrix(.) %>% 
  pheatmap(., colorRampPalette(c('white','red'))(100), labels_row = rep('',dim(snv_test)[1]))

snv_test %>% 
  as.matrix(.) %>% 
  heatmap(.)


#################### Gene enrichment gene-expression #################### 
# Data filtration
count_data_raw <- read.csv("./data/_raw/TCGA-SKCM.htseq_counts.tsv.gz",
                           sep="\t",
                           fill = 0)


count_data <- count_data_raw[,2:473]
rownames(count_data) <- count_data_raw$Ensembl_ID
colnames(count_data) <- str_replace_all(colnames(count_data), "\\.", "-")

survival_death <- as.data.frame(survival_filter) 
survival_death <- survival_death[,-1]
rownames(survival_death) = survival_filter$sample


intersect_patients <- intersect(rownames(survival_death),colnames(count_data))

survival_death <- survival_death[intersect_patients,]$OS
count_data <- count_data[,intersect_patients]


mel_expr_cpm <- apply(count_data, 
                      2, 
                      function(x) x/sum(x)*1000000)

# Finding enriched genes

# Make an empty list to store your signatures in
gene_signatures <- list()
# Make a function to do a MWW U test
row_mww <- function(x) {
  subtype <- x[survival_death == i]
  rest <- x[!survival_death == i]
  res <- wilcox.test(subtype, rest)
  return(res$p.value)
}
# Make a function to calculate log2 fc
row_fc <- function(x) {
  subtype <- x[survival_death == i]
  rest <- x[!survival_death == i]
  res <- log2(median(subtype)/median(rest))
  return(res)
}

for (i in c(1)) {
  pvals <- apply(mel_expr_cpm, 1, FUN = row_mww)
  padj <- p.adjust(pvals, method = "BY")
  sig_probes <- rownames(mel_expr_cpm)[padj<0.05]
  subtype_sig <- mel_expr_cpm[rownames(mel_expr_cpm) %in% sig_probes,]
  log2fc <- apply(subtype_sig, 1, FUN = row_fc)
  
  print(i)
  print(length(gene_signatures))
}

gene_signatures
