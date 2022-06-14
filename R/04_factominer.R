library("FactoMineR")
library("missMDA")
library("tidyverse")


# ncp is eqivalent to number of PC's in a PCA
# num.group.sup can be used to indicate which groups should NOT be included in MFA

data("wine")
dim(wine)
res <- MFA(wine, group=c(2,5,3,10,9,2), type=c("n",rep("s",5)),
           ncp=5, name.group=c("orig","olf","vis","olfag","gust","ens"),
           num.group.sup=c(1,6))
summary(res)
barplot(res$eig[,1],main="Eigenvalues",names.arg=1:nrow(res$eig))

res$global.pca

data(orange)
res.impute <- imputeMFA(orange, group=c(5,3), type=rep("s",2),ncp=2) 
res.mfa <- MFA(res.impute$completeObs,group=c(5,3),type=rep("s",2)) 
summary(res.mfa)
barplot(res.mfa$eig[,1],main="Eigenvalues",names.arg=1:nrow(res.mfa$eig))

#### Testing factoMineR on small dataframe ####

## Creating small dataframe
## Only keeping samples for which we have survival data 
snv_summary
# snv_summary_df <- as.data.frame(snv_summary)
# rownames(snv_summary_df) <- snv_summary_df$name
# 
# snv_summary_df <- snv_summary_df[,2:1000]
# snv_summary_df_sample <- snv_summary_df
# 
# dim(snv_summary_df_sample)
# 
# 
# res <- MFA(snv_summary_df_sample, group=c(499,), type=c("s"),
#            ncp=5, name.group=c("snv"))
snv_summary

snvs_factominer_test <- snv_summary %>%
  mutate(sample = name) %>%
  select(-name) %>%
  inner_join(survival_filter[1], by = "sample")

cnvs_factominer_test <- cnv_summary %>%
  inner_join(survival_filter[1], by = "sample")

mel_expr_cpm_t <- t(mel_expr_cpm)
mel_expr_cpm_t <- as.data.frame(mel_expr_cpm_t)
mel_expr_cpm_t <- str_replace_all(rownames(mel_expr_cpm_t), "\\.", "-")

rownames(mel_expr_cpm_t) <- (str_replace_all(rownames(mel_expr_cpm_t), "\\.", "-"))
mel_expr_cpm_t


expr_mad_t <- t(expr_mad)
colnames(expr_mad_t) <- gene_names_mad
View(head(expr_mad_t))
expr_mad_t <- as.data.frame(expr_mad_t)
rownames(expr_mad_t) <- str_replace_all(rownames(expr_mad_t), "\\.", "-")
View(head(expr_mad_t))


expr_mad_tibble <- as.tibble(expr_mad_t) %>% 
  mutate(sample = rownames(expr_mad_t))
expr_mad_tibble

pheno_factominer_test <- pheno_data_filtered %>% 
  mutate(sample = submitter_id.samples) %>%
  select(sample, 
         ethnicity.demographic,
         gender.demographic,
         race.demographic,
         year_of_birth.demographic,
         primary_diagnosis.diagnoses,
         prior_malignancy.diagnoses,
         prior_treatment.diagnoses,
         tumor_stage.diagnoses,
         name.tissue_source_site) %>%
  inner_join(survival_filter[1], by = "sample")

data_factominer_test <- cnvs_factominer_test %>% 
  inner_join(snvs_factominer_test, by = "sample") %>%
  inner_join(pheno_factominer_test, by = "sample") %>%
  inner_join(expr_mad_tibble, by = "sample") %>%
  column_to_rownames(., var = "sample") %>%
  drop_na() %>%
  mutate_at(vars(ethnicity.demographic,
                 gender.demographic,
                 race.demographic,
                 year_of_birth.demographic,
                 primary_diagnosis.diagnoses,
                 prior_malignancy.diagnoses,
                 prior_treatment.diagnoses,
                 tumor_stage.diagnoses,
                 name.tissue_source_site), factor)

data_factominer_test

dim(data_factominer_test)

dim(expr_mad_tibble)

dim(data_factominer_test)


### running imputation to fill out e.g. BMI 
#res.impute <- imputeMFA(data_factominer_test,  group=c(46,19210,9), type=c(rep("s",2), "n"),ncp=5)

res <- MFA(data_factominer_test, group=c(46,19210,9,5000), type=c(rep("s",2), "n","s"),
    ncp=5, name.group=c("cnv","snv","pheno","gene_expr"))

summary(res)
PC_matrix <- res$ind
PC_matrix_df <- as.data.frame(PC_matrix$coord)


pheno_data_filtered_subset <- pheno_data_filtered[pheno_data_filtered$submitter_id.samples %in% rownames(PC_matrix_df),]
rownames(pheno_data_filtered_subset) <- pheno_data_filtered_subset$submitter_id.samples
pheno_data_filtered_subset <- pheno_data_filtered_subset[rownames(PC_matrix_df),]
pheno_data_filtered_subset
PC_matrix_df$feature_pheno <- pheno_data_filtered_subset$age_at_initial_pathologic_diagnosis

survival_df <- as.data.frame(survival)
rownames(survival_df) <- survival_df$sample
survival_df <- survival_df[,-1]
survival_subset <- survival_df[rownames(PC_matrix_df),]
PC_matrix_df$feature_survival <- survival_subset$OS


ggplot(PC_matrix_df,aes(x=Dim.1, y=Dim.2, color = factor(feature_survival))) + 
  geom_point() + 
  stat_ellipse()

