library("tidyverse")
SNV_data <- read.csv("./data/_raw/TCGA-SKCM.mutect2_snv.tsv.gz",sep="\t")
SNV_data_tibble <- as_tibble(SNV_data)

View(head(SNV_data_tibble))

SNV_data_tibble %>%
  filter(effect == "missense_variant") %>%
  group_by(chrom) %>% 
  summarise(n = n(),) %>% 
  arrange(desc(n)) %>% 
  ggplot(.,aes(x=chrom,y=n)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

SNV_data_tibble

SNV_data_tibble %>% 
  filter(effect == "missense_variant") %>%
  distinct(gene) %>% 
  summarise(n=n())


pheno_data <- read.csv("./data/_raw/TCGA-SKCM.GDC_phenotype.tsv.gz",sep="\t")
pheno_data_tibble <- as_tibble(pheno_data)

pheno_data_tibble %>% 
  group_by(distant_metastasis_anatomic_site) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  ggplot(.,aes(x=distant_metastasis_anatomic_site,y=log10(n))) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



colnames(pheno_data_tibble)
pheno_data_tibble %>% 
  group_by(melanoma_clark_level_value) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  ggplot(.,aes(x=melanoma_clark_level_value,y=n)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pheno_data_tibble %>% 
  group_by(submitted_tumor_location) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

View(pheno_data_tibble)



count_data <- read.csv("./data/_raw/TCGA-SKCM.htseq_counts.tsv.gz",sep="\t",fill = 0)
count_data

count_data <- count_data[,2:473]
count_data <- sapply(count_data, as.numeric)
View(head(count_data))

gbm_expr_cpm <- apply(count_data, 2, function(x) x/sum(x)*1000000)

pheno_data[,submitted_tumor_location]
pheno_data <- pheno_data[pheno_data$submitted_tumor_location != '',]
pheno_data

colnames(pheno_data)


pheno_data$submitter_id.samples 

pheno_data

gbm_expr_cpm_filtered <- gbm_expr_cpm[, str_replace_all(colnames(gbm_expr_cpm), "\\.", "-") %in% pheno_data$submitter_id.samples]
pheno_data_filtered <- pheno_data[pheno_data$submitter_id.samples %in% str_replace_all(colnames(gbm_expr_cpm_filtered), "\\.", "-"),]

count_data_tibble <- as_tibble(count_data)

gbm_expr_cpm_filtered

colnames(pheno)

pheno_data_filtered$submitter_id.samples == str_replace_all(colnames(gbm_expr_cpm_filtered), "\\.", "-")
order(pheno_data_filtered$submitter_id.samples)
str_replace_all(colnames(gbm_expr_cpm_filtered), "\\.", "-")

pheno_sort <- pheno_data_filtered[order(pheno_data_filtered$submitter_id.samples),]
pheno_sort$submitter_id.samples
expr_sort <- gbm_expr_cpm_filtered[, order(colnames(gbm_expr_cpm_filtered))]
expr_sort

dim(gbm_expr_cpm_filtered)
dim(pheno_data_filtered)

pca <- prcomp(t(expr_sort))

pca %>% broom::tidy(matrix = "eigenvalues")

pca %>% broom::augment((pheno_sort)) %>% 
  ggplot(mapping = aes(
  x = .fittedPC1,
  y = .fittedPC2,
  color = submitted_tumor_location)) +
  geom_point() +
  stat_ellipse()
