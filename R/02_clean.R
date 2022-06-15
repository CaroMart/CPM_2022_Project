library(tidyverse)
library(Hmisc)
### SNVs
snvs_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.mutect2_snv.tsv.gz")

snvs <- snvs_raw %>% 
  filter(filter == "PASS",
         !grepl("synonymous_variant", effect))

snvs[rev(order(snvs$dna_vaf)),]

gene_lengths = read_delim(file = "data/gene_lengths.txt", delim = " ", col_names = c("gene", "length"))
gene_lengths

max(gene_lengths$length)

snv_summary_w_length = snvs %>% 
  group_by(Sample_ID, gene) %>% 
  summarise(count = n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = Sample_ID,
              values_from = count) %>% 
  mutate_at(-1, replace_na, replace = 0) %>%
  left_join(gene_lengths, by="gene") %>% 
  mutate_at(vars(matches("length")), 
            replace_na, 
            replace = max(gene_lengths$length)) %>% 
  mutate_at(
    vars(-matches(c("gene", "length"))),
    .funs = ~ . / length
  ) %>% 
  pivot_longer(-gene) %>% 
  pivot_wider(names_from = gene,
              values_from = value)

snv_summary = snv_summary_w_length %>% 
  filter(name != "length")
snv_summary

table(snv_summary[,"BRAF"])

snv_summary[dim(snv_summary)[1],]
snv_summary[-dim(snv_summary)[1],-1]

snv_summary_count = snvs %>% 
  group_by(Sample_ID, gene) %>% 
  summarise(count = n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = gene,
              values_from = count) %>% 
  mutate_at(-1, replace_na, replace = 0) 

sum_gene_freq_no_correction <- colSums(snv_summary_count[,-1])

ggplot(mapping = aes(x = sum_gene_freq_no_correction, y =as.numeric(snv_summary_w_length[dim(snv_summary_w_length)[1],-1]))) +
  geom_point() +
  
length(sum_gene_freq_no_correction)
length(snv_summary_w_length[dim(snv_summary_w_length)[1],-1])
as.numeric(snv_summary_w_length[dim(snv_summary_w_length)[1],-1])


names(avg_gene_freq[order(avg_gene_freq, decreasing=TRUE)][1:25]) %>% 
  cat(sep="\n")

snv_summary

?rowMedians

avg_gene_freq = colMeans(snv_summary[-dim(snv_summary)[1],-1])

avg_gene_freq["DEFB115"]

names(avg_gene_freq[order(avg_gene_freq, decreasing=TRUE)][1:25]) %>% 
  cat(sep="\n")

ggplot(mapping = aes(x =names(gene_counts[order(gene_counts, decreasing=TRUE)][1:25]), y = gene_counts[order(gene_counts, decreasing=TRUE)][1:25]))+
    geom_col()

### CNVs
cnvs_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.masked_cnv.tsv.gz")

cnvs <- cnvs_raw %>% 
  filter(abs(value) > 0.3) %>% 
  mutate(cn_type = case_when(
    value > 0 ~ "gain",
    TRUE ~ "loss"
  ),
  chrom_cn_type = paste(Chrom, cn_type, sep="_"))

cnv_summary <- cnvs %>% 
  group_by(sample, chrom_cn_type) %>% 
  summarise(count = n()) %>% 
  pivot_wider(names_from = chrom_cn_type,
              values_from = count) %>% 
  mutate_at(vars(matches("^[0-9]")), replace_na, replace = 0)

cnv_summary

### Counts
counts <- read_tsv(file = "data/_raw/TCGA-SKCM.htseq_counts.tsv.gz")
counts

### Phenotype
pheno <- read_tsv(file = "data/_raw/TCGA-SKCM.GDC_phenotype.tsv.gz")

### Methylation
# methylation_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.methylation450.tsv.gz")


### Survival
survival <- read_tsv(file = "data/_raw/TCGA-SKCM.survival.tsv")
