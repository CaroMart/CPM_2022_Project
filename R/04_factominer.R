library("FactoMineR")
library("missMDA")
library("tidyverse")



#### Testing factoMineR on small dataframe ####

snv_summary_binary_t <- as.data.frame(t(snv_summary_binary))

survival_filter_tibble <- survival_filter %>% 
  rownames_to_column('sample')
  


snvs_factominer_test <- snv_summary_binary_t[,differential_snvs] %>%
  rownames_to_column('sample') %>% 
  as_tibble(.)


cnvs_factominer_test <- cnv_summary %>%
  rownames_to_column('sample')


pheno_filter_categorical <- pheno_filter_categorical %>%
  rownames_to_column('sample') %>% 
  as_tibble(.)

pheno_filter_numeric <- pheno_filter_numeric %>%
  rownames_to_column('sample')

estimate_output_test <- estimate_output %>%
  rownames_to_column('sample') %>% 
  dplyr::select(c('sample','TumorPurity'))



MCP_output <- t(MCP_output)
MCP_output <- as.data.frame(MCP_output)
MCP_output <- MCP_output %>% 
  rownames_to_column('sample') %>% 
  as_tibble(.)
  

mel_cpm
mel_expr_cpm_t <- as.data.frame(t(mel_cpm))
mel_expr_cpm_t <- mel_expr_cpm_t[,significant_genes]
dim(mel_expr_cpm_t)
mel_expr_cpm_t
expr_mad_tibble <- mel_expr_cpm_t %>% 
  rownames_to_column('sample')


pheno_factominer_test <- pheno %>% 
  rownames_to_column('sample') 


## Running MFA for real 

## The data frames - put in clean or augment

samples_in_all_data <- as.data.frame(intersect_patients)
samples_in_all_data
colnames(samples_in_all_data) <- c('sample')
samples_in_all_data <- as_tibble(samples_in_all_data)
samples_in_all_data

impute_data_numeric <-  cnvs_factominer_test %>% 
  inner_join(snvs_factominer_test, by = "sample") %>%
  inner_join(pheno_filter_numeric, by = "sample") %>%
  inner_join(expr_mad_tibble, by = "sample") %>%
  inner_join(estimate_output_test, by = "sample") %>%
  column_to_rownames(., var = "sample")

dim(snvs_factominer_test)
dim(estimate_output_test)
dim(impute_data_numeric)

impute_data_categorical <- pheno_filter_categorical %>%
  column_to_rownames(., var = "sample")  

survival_of_351 <- samples_in_all_data %>% 
  left_join(survival_filter_tibble, by = "sample")
dim(survival_of_351)
## Imputing NA values for continuous values 
# should continuous variables be s or c? (that is, are they scaled or not). 
# impute seems to only work when they are c
dim(expr_mad_tibble)


dim(impute_data_numeric)



# res.imputeMFA <- imputeMFA(impute_data_numeric,group=c(46,13317,8,136), type=c("c","c","c","c"),ncp=2)
res.imputeMFA <- imputeMFA(impute_data_numeric,group=c(46,6517,8,136,1), type=c("c","c","c","c","c"),ncp=2)
res.imputeMCA <- imputeMCA(impute_data_categorical, ncp=2)
dim(res.imputeMFA$completeObs)
dim(res.imputeMCA$completeObs)



data_factominer <- samples_in_all_data %>%
  bind_cols(res.imputeMCA$completeObs) %>%
  bind_cols(res.imputeMFA$completeObs) 

dim(res.imputeMFA$completeObs)

data_factominer <- data_factominer %>% 
  as_tibble() %>% 
  column_to_rownames('sample')

dim(data_factominer)

## Running MFA

# res <- MFA(data_factominer, group=c(16,46,13317,8,60488), type=c("n","s","s","s","s"),
#            ncp=5, name.group=c("pheno_cat","cnv","snv","pheno_num","gene_expr"),graph = FALSE)
res <- MFA(data_factominer, group=c(16,46,6517,8,136,1), type=c("n","s","s","s","s","s"),
           ncp=5, name.group=c("pheno_cat","cnv","snv","pheno_num","gene_expr","purity"))

saveRDS(res, file = "/tmp/res_object_sig_snv_sig_gene_purity.rds")
## Running HMFA

summary(res)
## Plotting result
PC_matrix <- res$ind
PC_matrix_df <- as.data.frame(PC_matrix$coord)

ggplot(PC_matrix_df,aes(x=Dim.1, y=Dim.2, color = factor(survival_of_351$OS))) + 
  geom_point() + 
  stat_ellipse()
