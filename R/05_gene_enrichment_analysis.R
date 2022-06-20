library('tidyverse')
library("gplots")
library('pheatmap')
library("DESeq2")

#################### Gene enrichment gene-expression #################### 
# Data filtration
survival_death <- survival_filter$OS

# Finding enriched genes

# Make an empty list to store your signatures in
gene_signatures <- list()
# Make a function to do a MWW U test
row_mww <- function(x) {
  subtype <- x[survival_death == 1]
  rest <- x[!survival_death == 0]
  res <- wilcox.test(subtype, rest)
  return(res$p.value)
}
# Make a function to calculate log2 fc
row_fc <- function(x) {
  subtype <- x[survival_death == 1]
  rest <- x[!survival_death == 0]
  res <- log2(median(subtype)/median(rest))
  return(res)
}
dim(mel_cpm)


pvals <- apply(mel_cpm, 1, FUN = row_mww)


padj <- p.adjust(pvals, method = "BY")

sig_probes <- rownames(mel_cpm)[padj<0.05]
subtype_sig <- mel_cpm[rownames(mel_cpm) %in% sig_probes,]
log2fc <- apply(subtype_sig, 1, FUN = row_fc)


subtype_sig

survival_filter[,'OS'] <- factor(survival_filter[,'OS'])

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(mel_counts)),
                              colData = survival_filter[,-2],
                              design = ~ OS)


dds_analysis <- DESeq(dds)
dds_analysis

resultsNames(dds_analysis)


sorted_padj <- results(dds_analysis,
                       name = "OS_1_vs_0"
) %>%
  as.data.frame() %>%
  mutate(Significance = case_when(
    padj <= 0.05 ~ "Significant",
    padj > 0.05 ~ "Not significant"
  )) %>%
  drop_na(
    padj,
    log2FoldChange
  ) %>%
  arrange(., padj)

significant_genes <- sorted_padj[sorted_padj$Significance == "Significant",]

significant_genes <- rownames(significant_genes[abs(significant_genes$log2FoldChange) > 0.5,])
significant_genes


mel_cpm_sig <- mel_cpm[rownames(mel_cpm) %in% significant_genes,]

pca <- prcomp(t(mel_cpm_sig))

pca %>% broom::tidy(matrix = "eigenvalues")


pca %>% broom::augment((survival_filter)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = factor(OS))) + # Skal farves på survival
  geom_point() + 
  stat_ellipse()

