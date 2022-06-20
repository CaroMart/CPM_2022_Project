library(tidyverse)
library(Hmisc)
library("MCPcounter")

### Phenotype
pheno <- read_tsv(file = "data/_raw/TCGA-SKCM.GDC_phenotype.tsv.gz")

pheno_filter_all <- pheno %>% 
  filter(!(submitted_tumor_location %in% c('', "Primary Tumor"))) %>% 
  mutate(sample = submitter_id.samples) %>%
  dplyr::select(sample,
         age_at_initial_pathologic_diagnosis,
         breslow_depth_value,
         days_to_collection.samples,
         days_to_submitted_specimen_dx,
         gender.demographic,
         height.exposures,
         history_of_neoadjuvant_treatment,
         icd_10_code.diagnoses,
         initial_weight.samples,
         malignant_neoplasm_mitotic_count_rate,
         melanoma_clark_level_value,
         melanoma_ulceration_indicator,
         pathologic_M,
         pathologic_N,
         pathologic_T,
         primary_diagnosis.diagnoses,
         prior_malignancy.diagnoses,
         prior_systemic_therapy_type, #NA = None
         race.demographic,
         ethnicity.demographic,
         site_of_resection_or_biopsy.diagnoses,
         submitted_tumor_location,
         system_version,
         weight) %>%
  replace_na(list(prior_systemic_therapy_type = "None")) %>%
  mutate_at(vars(gender.demographic,
                 history_of_neoadjuvant_treatment,
                 icd_10_code.diagnoses,
                 melanoma_clark_level_value,
                 melanoma_ulceration_indicator,
                 pathologic_M,
                 pathologic_N,
                 pathologic_T,
                 primary_diagnosis.diagnoses,
                 prior_malignancy.diagnoses,
                 prior_systemic_therapy_type,
                 race.demographic,
                 ethnicity.demographic,
                 site_of_resection_or_biopsy.diagnoses,
                 submitted_tumor_location,
                 system_version), factor)

pheno_filter_categorical <- pheno_filter_all %>%
  dplyr::select(sample,
         gender.demographic,
         history_of_neoadjuvant_treatment,
         icd_10_code.diagnoses,
         melanoma_clark_level_value,
         melanoma_ulceration_indicator,
         pathologic_M,
         pathologic_N,
         pathologic_T,
         primary_diagnosis.diagnoses,
         prior_malignancy.diagnoses,
         prior_systemic_therapy_type,
         race.demographic,
         ethnicity.demographic,
         site_of_resection_or_biopsy.diagnoses,
         submitted_tumor_location,
         system_version) %>% 
  column_to_rownames('sample')

pheno_filter_numeric <- pheno_filter_all %>%
  dplyr::select(-gender.demographic,
         -history_of_neoadjuvant_treatment,
         -icd_10_code.diagnoses,
         -melanoma_clark_level_value,
         -melanoma_ulceration_indicator,
         -pathologic_M,
         -pathologic_N,
         -pathologic_T,
         -primary_diagnosis.diagnoses,
         -prior_malignancy.diagnoses,
         -prior_systemic_therapy_type,
         -race.demographic,
         -ethnicity.demographic,
         -site_of_resection_or_biopsy.diagnoses,
         -submitted_tumor_location,
         -system_version) %>% 
  column_to_rownames('sample')


pheno
pheno <- pheno %>% 
  mutate(prior_systemic_therapy_type = replace_na(prior_systemic_therapy_type, 
                                                  "None")) %>% 
  column_to_rownames('submitter_id.samples')

survival <- read_tsv("./data/_raw/TCGA-SKCM.survival.tsv")
survival_filter <- survival %>% 
  filter(sample %in% pull(pheno_filter_all,sample)) %>% 
  column_to_rownames('sample')


### SNVs
snvs_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.mutect2_snv.tsv.gz")

snvs <- snvs_raw %>% 
  filter(filter == "PASS",
         !grepl("synonymous_variant", effect))


gene_lengths = read_delim(file = "data/gene_lengths.txt", delim = " ", col_names = c("gene", "length"))


snv_summary_binary <- snvs %>% 
  group_by(Sample_ID, gene) %>% 
  summarise(count = 1) %>% 
  ungroup() %>%
  pivot_wider(names_from = Sample_ID,
              values_from = count) %>% 
  mutate_at(-1,replace_na, replace = 0) %>% 
  column_to_rownames('gene')

snv_summary_binary <- snv_summary_binary[rowSums(snv_summary_binary) > 4,]
dim(snv_summary_binary)

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
  mutate_at(vars(matches("^[0-9]")), replace_na, replace = 0) %>% 
  column_to_rownames('sample') %>% 
  mutate_at(-1,replace_na, replace = 0)



### Counts
counts <- read_tsv(file = "data/_raw/TCGA-SKCM.htseq_counts.tsv.gz")
counts <- counts %>% 
  column_to_rownames('Ensembl_ID')



#################### Estimating purity in samples #################### 
BiocManager::install("org.Hs.eg.db")
library('org.Hs.eg.db')
library('AnnotationDbi')
counts_for_estimate <- round(2^counts - 1,0)
counts_for_estimate
rownames(counts_for_estimate) <- gsub("\\..*","",rownames(counts))
ensemble_to_genename <- mapIds(org.Hs.eg.db, keys = rownames(counts_for_estimate), keytype = "ENSEMBL", column = "SYMBOL")
ensemble_to_genename <- as.data.frame(ensemble_to_genename)

rownames(counts_for_estimate) <- ensemble_to_genename$ensemble_to_genename
rownames(counts_for_estimate) <- make.names(ensemble_to_genename$ensemble_to_genename, unique = TRUE)

estimate_output <- LARSEstimateFunction(counts_for_estimate)
hist(purity_scores$TumorPurity)
summary(purity_scores$ESTIMATEScore[purity_scores$TumorPurity > 0 ])
estimate_output


ggplot(estimate_output, aes(x=ESTIMATEScore, y=TumorPurity)) + geom_point()



MCP_output <- MCPcounter.estimate(counts_for_estimate,featuresType = "HUGO_symbols")
MCP_output <- as.data.frame(MCP_output)



mel_cpm <- apply(counts, 
                      2, 
                      function(x) x/sum(x)*1000000)

mad_values <- apply(mel_cpm, 1, mad)
mad_mel_cpm <- mel_cpm[order(mad_values,decreasing = TRUE)[1:5000], ]
# mad_mel_cpm <- mel_cpm[mad_values > 1.696, ]
mad_mel_cpm <- as.data.frame(mad_mel_cpm)
dim(mad_mel_cpm)
mad_mel_cpm

mel_cpm <- as.data.frame(mel_cpm)

################### Creating the intersects ###################

intersect_patients <- Reduce(intersect,list(rownames(pheno),rownames(survival_filter),colnames(snv_summary_binary),rownames(cnv_summary),colnames(mel_cpm),rownames(pheno_filter_categorical),rownames(pheno_filter_numeric)))
length(intersect_patients)

pheno <- pheno[intersect_patients,]
survival_filter <- survival_filter[intersect_patients,]
snv_summary_binary <- snv_summary_binary[,intersect_patients]
cnv_summary <- cnv_summary[intersect_patients,]
mel_counts <- counts[,intersect_patients]
mel_cpm <- mel_cpm[,intersect_patients]
mad_mel_cpm <- mad_mel_cpm[,intersect_patients]
pheno_filter_categorical <- as.data.frame(pheno_filter_categorical)[intersect_patients,]
pheno_filter_numeric <- pheno_filter_numeric[intersect_patients,]
pheno_filter_numeric

####### PLOTS #######
## CNV ##
cnv_summary %>% 
  mutate(response = survival_filter$OS) %>% 
  group_by(response) %>% 
  summarise(count = n())
# summarise(.,sum_ = sum(count)/n()) %>% 


cnv_summary %>% 
  mutate(response = survival_filter$OS) %>% 
  rownames_to_column('sample') %>% 
  pivot_longer(!c(response,sample), names_to = "chromosome", values_to = "count") %>%
  group_by(chromosome,response) %>% 
  summarise(.,sum_ = sum(count)/n()) %>% 
  ggplot(.,aes(x=factor(chromosome),y=sum_,group = factor(response),fill = factor(response))) +
  geom_col(position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cnv_summary
pca <- prcomp(cnv_summary)
pca %>% broom::tidy(matrix = "eigenvalues")


pca %>% broom::augment((survival_filter)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2, 
    color = factor(OS))) + # Skal farves på survival
  geom_point() + 
  stat_ellipse()


dist_mi <- as.matrix(dist(pca$x))
heatmap(dist_mi)



##### SNV #####  
rownames(snv_summary_binary)
dim(snv_summary_binary)
snv_sums_1 <- as.data.frame(rowSums(snv_summary_binary[,survival_filter$OS == 1]/sum(survival_filter$OS == 1)))
snv_sums_0 <- as.data.frame(rowSums(snv_summary_binary[,survival_filter$OS == 0]/sum(survival_filter$OS == 0)))


difference_in_freq <- snv_sums_1 - snv_sums_0
colnames(difference_in_freq) <- c("diff_in_freq")
summary(difference_in_freq)
difference_in_freq %>% 
  ggplot(.,aes(x=diff_in_freq)) + geom_histogram() + 
  theme_classic()


upregulated_in_death <- difference_in_freq %>% 
  rownames_to_column('gene') %>% 
  filter((diff_in_freq > 0.008858))

upregulated_in_death <- difference_in_freq %>% 
  rownames_to_column('gene') %>% 
  filter((diff_in_freq < -0.014698))


differential_snvs <- difference_in_freq %>% 
  rownames_to_column('gene') %>% 
  filter((diff_in_freq < -0.014698) | (diff_in_freq > 0.008858))

differential_snvs <- differential_snvs$gene
differential_snvs

write_csv(upregulated_in_death,"./tmp/up_reg_death.csv")

summary(difference_in_freq)

survival_filter$OS

# PCA SNV
SNV_outlier <- "TCGA-FW-A3R5-06A"
pca <- prcomp(t(snv_summary_binary))
pca %>% broom::tidy(matrix = "eigenvalues")


pca %>% broom::augment((survival_filter)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2, 
    color = factor(OS))) + # Skal farves på survival
  geom_point() + 
  stat_ellipse()


dist_mi <- as.matrix(dist(pca$x))
heatmap(dist_mi)


# library(qgraph)
# dist_mi <- 1/dist_mi # one over, as qgraph takes similarity matrices as input
# jpeg('tmp/dist_map.jpg', width=1000, height=1000, unit='px')
# qgraph(dist_mi, layout='spring', vsize=3)
# dev.off()
# 
# pca$x[168,1]
# 
# 
# 
# rownames(as.data.frame(pca$x))[168]
# dim(dist(pca$x[,c(1,2)]))

##### Count data #####  

pca <- prcomp(t(mel_cpm))

pca %>% broom::tidy(matrix = "eigenvalues")


pca %>% broom::augment((survival_filter)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = factor(OS))) + # Skal farves på survival
  geom_point() + 
  stat_ellipse()


####### OLD ####### 
# 
# snv_summary_w_length = snvs %>% 
#   group_by(Sample_ID, gene) %>% 
#   summarise(count = n()) %>% 
#   ungroup() %>%
#   pivot_wider(names_from = Sample_ID,
#               values_from = count) %>% 
#   mutate_at(-1, replace_na, replace = 0) %>%
#   left_join(gene_lengths, by="gene") %>% 
#   mutate_at(vars(matches("length")), 
#             replace_na, 
#             replace = max(gene_lengths$length)) %>% 
#   mutate_at(
#     vars(-matches(c("gene", "length"))),
#     .funs = ~ . / length
#   ) %>% 
#   pivot_longer(-gene) %>% 
#   pivot_wider(names_from = gene,
#               values_from = value)
# 
# snv_summary = snv_summary_w_length %>% 
#   filter(name != "length")
# snv_summary
# 
# ## snv for FactoMineR
# snv_summary_filter <- snv_summary %>%
#   filter(name %in% pull(survival_filter, sample)) %>%
#   mutate(sample = name) %>%
#   select(-name)
# 
# table(snv_summary[,"BRAF"])
# 
# snv_summary[dim(snv_summary)[1],]
# snv_summary[-dim(snv_summary)[1],-1]
# 
# snv_summary_count = snvs %>% 
#   group_by(Sample_ID, gene) %>% 
#   summarise(count = n()) %>% 
#   ungroup() %>%
#   pivot_wider(names_from = gene,
#               values_from = count) %>% 
#   mutate_at(-1, replace_na, replace = 0) 
# 
# sum_gene_freq_no_correction <- colSums(snv_summary_count[,-1])
# 
# ggplot(mapping = aes(x = sum_gene_freq_no_correction, y =as.numeric(snv_summary_w_length[dim(snv_summary_w_length)[1],-1]))) +
#   geom_point() +
#   
# length(sum_gene_freq_no_correction)
# length(snv_summary_w_length[dim(snv_summary_w_length)[1],-1])
# as.numeric(snv_summary_w_length[dim(snv_summary_w_length)[1],-1])
# 
# 
# names(avg_gene_freq[order(avg_gene_freq, decreasing=TRUE)][1:25]) %>% 
#   cat(sep="\n")
# 
# snv_summary
# 
# ?rowMedians
# 
# avg_gene_freq = colMeans(snv_summary[-dim(snv_summary)[1],-1])
# 
# avg_gene_freq["DEFB115"]
# 
# names(avg_gene_freq[order(avg_gene_freq, decreasing=TRUE)][1:25]) %>% 
#   cat(sep="\n")
# 
# ggplot(mapping = aes(x =names(gene_counts[order(gene_counts, decreasing=TRUE)][1:25]), y = gene_counts[order(gene_counts, decreasing=TRUE)][1:25]))+
#     geom_col()
# 
# ### CNVs
# cnvs_raw <- read_tsv(file = "data/_raw/TCGA-SKCM.masked_cnv.tsv.gz")
# 
# cnvs <- cnvs_raw %>% 
#   filter(abs(value) > 0.3) %>% 
#   mutate(cn_type = case_when(
#     value > 0 ~ "gain",
#     TRUE ~ "loss"
#   ),
#   chrom_cn_type = paste(Chrom, cn_type, sep="_"))
# 
# cnv_summary <- cnvs %>% 
#   group_by(sample, chrom_cn_type) %>% 
#   summarise(count = n()) %>% 
#   pivot_wider(names_from = chrom_cn_type,
#               values_from = count) %>% 
#   mutate_at(vars(matches("^[0-9]")), replace_na, replace = 0)
# 
# cnv_summary
# 
# ## for FactoMineR
# cnv_summary_filter <- cnv_summary %>%
#   filter(sample %in% pull(survival_filter, sample))
# 
# ### Counts
# counts <- read_tsv(file = "data/_raw/TCGA-SKCM.htseq_counts.tsv.gz")
# counts
# 
# 
# ### Survival
# survival <- read_tsv(file = "data/_raw/TCGA-SKCM.survival.tsv")
# 
#   
# 
# 
# 
