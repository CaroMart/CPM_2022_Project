library(tidyverse)

### SNVs
snvs_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.mutect2_snv.tsv.gz")

snvs <- snvs_raw %>% 
  filter(filter == "PASS",
         !grepl("synonymous_variant", effect))

snvs[rev(order(snvs$dna_vaf)),]

snv_summary = snvs %>% 
  group_by(Sample_ID, gene) %>% 
  summarise(count = 1) %>% 
  ungroup() %>% 
  pivot_wider(names_from = gene,
            values_from = count) %>%
  mutate_at(-1, replace_na, replace = 0)

snv_summary[-1,]

gene_counts = colSums(snv_summary[,-1])

names(gene_counts[order(gene_counts, decreasing=TRUE)][1:25]) %>% 
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
pheno

### Methylation
#methylation_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.methylation450.tsv.gz")


### Survival
survival <- read_tsv(file = "data/_raw/TCGA-SKCM.survival.tsv")
