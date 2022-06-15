#-------------------------------------------------------------------------------
# Load libraries
library("tidyverse")

#-------------------------------------------------------------------------------
# Load SNV data
SNV_data <- read.csv("./data/_raw/TCGA-SKCM.mutect2_snv.tsv.gz",
                     sep="\t")
SNV_data_tibble <- as_tibble(SNV_data)

#-------------------------------------------------------------------------------
# Plotting number of missense variants pr. chromosome
SNV_data_tibble %>%
  filter(effect == "missense_variant") %>%
  group_by(chrom) %>% 
  summarise(n = n(),) %>% 
  arrange(desc(n)) %>% 
  ggplot(.,aes(x = chrom,
               y = n)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1))

# Unique genes with missense variants
SNV_data_tibble %>% 
  filter(effect == "missense_variant") %>%
  distinct(gene) %>% 
  summarise(n=n())

#-------------------------------------------------------------------------------
# Load pheno data
pheno_data <- read.csv("./data/_raw/TCGA-SKCM.GDC_phenotype.tsv.gz",
                       sep="\t")
pheno_data_tibble <- as_tibble(pheno_data)

#-------------------------------------------------------------------------------

# Plotting distant metastasis sites
pheno_data_tibble %>% 
  group_by(distant_metastasis_anatomic_site) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  ggplot(.,aes(x=distant_metastasis_anatomic_site,
               y=log10(n))) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))



colnames(pheno_data_tibble)

# Plotting frequency of cancer stages (?)
pheno_data_tibble %>% 
  group_by(melanoma_clark_level_value) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n)) %>% 
  ggplot(.,aes(x=melanoma_clark_level_value,
               y=n)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust=1))

# frequency of submitted tumor location
pheno_data_tibble %>% 
  group_by(submitted_tumor_location) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

# Plotting Age distribution
age_distribution_plot <- pheno_data_tibble %>% 
  ggplot(mapping = aes(
    x = age_at_initial_pathologic_diagnosis,
    fill = gender.demographic)) +
  geom_histogram(bins = 30) +
  facet_wrap(~gender.demographic) +
  labs(title = "Distribution of Age at diagnosis stratified on sex",
       x = "Age at initial diagnosis") +
  theme_classic()
  
age_distribution_plot

# Plotting race distribution
race_distribution_plot <- pheno_data_tibble %>% 
  ggplot(mapping = aes(
    x = race.demographic)) +
  geom_bar() +
  #facet_wrap(~gender.demographic) +
  labs(title = "Distribution of race",
       x = "Race" ) +
  theme_classic()

race_distribution_plot

# Plotting cancer stage distribution
stage_distribution_plot <- pheno_data_tibble %>% 
  ggplot(mapping = aes(
    x = tumor_stage.diagnoses)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Distribution of tumor stages",
       x = "Tumor stage")

stage_distribution_plot

#-------------------------------------------------------------------------------
# Load count data
count_data_raw <- read.csv("./data/_raw/TCGA-SKCM.htseq_counts.tsv.gz",
                       sep="\t",
                       fill = 0)


count_data <- count_data_raw[,2:473]

count_data <- sapply(count_data, 
                     as.numeric)
View(head(count_data))

#-------------------------------------------------------------------------------

mel_expr_cpm <- apply(count_data, 
                      2, 
                      function(x) x/sum(x)*1000000)

#pheno_data[,submitted_tumor_location]
dim(mel_expr_cpm)

pheno_data <- pheno_data[!(pheno_data$submitted_tumor_location %in% c('', "Primary Tumor")),]
pheno_data

colnames(pheno_data)


pheno_data$submitter_id.samples 

pheno_data

mel_expr_cpm_filtered <- mel_expr_cpm[, str_replace_all(colnames(mel_expr_cpm), "\\.", "-") %in% pheno_data$submitter_id.samples]
pheno_data_filtered <- pheno_data[pheno_data$submitter_id.samples %in% str_replace_all(colnames(mel_expr_cpm_filtered), "\\.", "-"),]

count_data_tibble <- as_tibble(count_data)

mel_expr_cpm_filtered

colnames(pheno)


#-------------------------------------------------------------------------------
# PCA
pheno_data_filtered$submitter_id.samples == str_replace_all(colnames(mel_expr_cpm_filtered), "\\.", "-")
order(pheno_data_filtered$submitter_id.samples)
str_replace_all(colnames(mel_expr_cpm_filtered), "\\.", "-")

pheno_sort <- pheno_data_filtered[order(pheno_data_filtered$submitter_id.samples),]
pheno_sort$submitter_id.samples
expr_sort <- mel_expr_cpm_filtered[, order(colnames(mel_expr_cpm_filtered))]

expr_sort

mad_values = apply(expr_sort, 1, mad)
gene_names_mad = count_data_raw[order(mad_values)[1:5000],1]
expr_mad = expr_sort[order(mad_values)[1:5000], ]
expr_mad


colnames(expr_mad)

dim(expr_mad)

dim(mel_expr_cpm_filtered)
dim(pheno_data_filtered)

pca <- prcomp(t(expr_mad))

pca %>% broom::tidy(matrix = "eigenvalues")
pheno_sort$pathologic_T

pheno_sort[c("breslow_depth_value", "submitted_tumor_location")]

pca %>% broom::augment((pheno_sort)) %>% 
  ggplot(mapping = aes(
  x = .fittedPC1,
  y = .fittedPC2,
  color = pathologic_N)) +
  geom_point() + 
  stat_ellipse()

pca %>% broom::augment((pheno_sort)) %>% 
  ggplot(mapping = aes(
    y = .fittedPC1,
    x = submitted_tumor_location)) +
  geom_boxplot()

wilcox.test(data=(pca %>% broom::augment((pheno_sort))),
              .fittedPC1~(submitted_tumor_location=="Primary Tumor"))


survival_filter = survival[survival$sample %in% pheno_sort$submitter_id.samples,]

new_survival = data.frame(sample = pheno_sort$submitter_id.samples, OS = rep(NA, length(pheno_sort$submitter_id.samples)))

new_survival[new_survival$sample %in% survival_filter$sample,2] = survival_filter$OS



table(new_survival$OS)
survival_filter
pca %>% broom::augment((new_survival)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color =  as.factor(OS))) +
  geom_point(alpha = 0.8)

survival
