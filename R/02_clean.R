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
View(pheno)
pheno$malignant_neoplasm_mitotic_count_rate
features_of_interest <- c("age_at_initial_pathologic_diagnosis","breslow_depth_value","age_at_initial_pathologic_diagnosis","days_to_submitted_specimen_dx","malignant_neoplasm_mitotic_count_rate") # submitter_id.samples
hist.data.frame(pheno[,features_of_interest])
pheno$new_tumor_dx_prior_submitted_specimen_dx

ggplot(data = pheno, aes(x=factor(history_of_neoadjuvant_treatment))) + geom_bar()
ggplot(data = pheno, aes(x=factor(melanoma_clark_level_value))) + geom_bar()
ggplot(data = pheno, aes(x=factor(melanoma_origin_skin_anatomic_site))) + geom_bar()
ggplot(data = pheno, aes(x=factor(melanoma_ulceration_indicator))) + geom_bar()
ggplot(data = pheno, aes(x=factor(new_tumor_dx_prior_submitted_specimen_dx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(other_dx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(pathologic_M))) + geom_bar()
ggplot(data = pheno, aes(x=factor(pathologic_N))) + geom_bar()
ggplot(data = pheno, aes(x=factor(pathologic_T))) + geom_bar()
ggplot(data = pheno, aes(x=factor(person_neoplasm_cancer_status))) + geom_bar()
ggplot(data = pheno, aes(x=factor(postoperative_rx_tx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(primary_melanoma_at_diagnosis_count))) + geom_bar()
ggplot(data = pheno, aes(x=factor(primary_neoplasm_melanoma_dx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(primary_tumor_multiple_present_ind))) + geom_bar()
ggplot(data = pheno, aes(x=factor(prior_radiation_therapy))) + geom_bar()
ggplot(data = pheno, aes(x=factor(prior_systemic_therapy))) + geom_bar()
ggplot(data = pheno, aes(x=factor(prior_systemic_therapy_type))) + geom_bar()
ggplot(data = pheno, aes(x=factor(radiation_therapy))) + geom_bar()
ggplot(data = pheno, aes(x=factor(radiation_therapy_to_primary))) + geom_bar()
ggplot(data = pheno, aes(x=factor(submitted_tumor_location))) + geom_bar() 
ggplot(data = pheno, aes(x=factor(subsequent_primary_melanoma_during_followup))) + geom_bar() 
ggplot(data = pheno, aes(x=factor(subsequent_primary_melanoma_during_followup))) + geom_bar() 

pheno_factor <- apply(pheno,1, as.factor)
class(pheno_factor)
phen


View(head(pheno))

colnames(pheno_factor)


names(pheno_factor)

for (feature in names(pheno_factor)){
  plot <- ggplot(data = pheno_factor, aes_string(x=feature)) + 
    geom_bar()
  ggsave(filename = paste("./tmp/",feature,"_bar",".png",sep=""),plot)
}

for (feature in colnames(pheno)){
  plot <- ggplot(data = pheno, aes_string(x=factor(pheno[,feature]))) + 
    geom_histogram()
  ggsave(filename = paste("./tmp/",feature,"_hist",".png",sep=""),plot)
}





write.csv(pheno,file = "./tmp/pheno.csv")

### Methylation
methylation_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.methylation450.tsv.gz")


### Survival
survival <- read_tsv(file = "data/_raw/TCGA-SKCM.survival.tsv")
