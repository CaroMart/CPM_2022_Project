library('tidyverse')
library("MCPcounter")
library("tidyverse")

source("R_new/estimate_function.R")

####### TPM data ####### 
estimate_output <- myEstimateFunction(tpm_ensg)
MCP_output <- MCPcounter.estimate(tpm_ensg,featuresType = "HUGO_symbols")
MCP_output

## Fitting linear models to estimate output
round(estimate_output$TumorPurity, )

tpm_ensg_t <- t(tpm_ensg)
tpm_ensg_t <- as.data.frame(tpm_ensg_t)

purity_scores <- estimate_output[rownames(tpm_ensg_t),]$TumorPurity
tpm_ensg_t$purity <- purity_scores

tpm_ensg_t

tpm_ensg_t_long <- tpm_ensg_t %>% 
  as_tibble() %>% 
  mutate(sample_nr = rownames(tpm_ensg_t)) %>% 
  pivot_longer(!c("sample_nr","purity"),names_to = "genes",values_to = "counts") 

tpm_ensg_t_grouped_data <- tpm_ensg_t_long %>% 
  group_by(genes) %>% 
  nest()  %>% 
  ungroup()

tpm_ensg_t_grouped_data <- tpm_ensg_t_grouped_data %>%
  mutate(mdl = map(data, ~glm(counts ~ purity,
                              data = .x)))  

#tpm_ensg_t_grouped_data$data
tpm_ensg_t_grouped_data_calc <- tpm_ensg_t_grouped_data %>% 
  mutate(max_purity_counts = map(mdl, function(x){
    tibble(estimated_purity = sum(coef(x)*c(1,1)) + residuals(x),
           sample_nr = rownames(tpm_ensg_t))
  }))

#tpm_ensg_t_grouped_data_calc
tpm_ensg_t_grouped_data_calc_wider <- tpm_ensg_t_grouped_data_calc %>% 
  unnest(max_purity_counts) %>% 
  select(-data, -mdl) %>% 
  pivot_wider(names_from = genes, values_from = estimated_purity)

tpm_ensg_t_grouped_data_calc_wider

g_t_grouped_data_calc[1,4][[1]]

grouped_data_sample_tidy <- tpm_ensg_t_grouped_data %>% 
  mutate(tidy = map(mdl, broom::tidy)) %>% 
  unnest(tidy)

ggplot(data = estimate_output, aes(x=ESTIMATEScore)) + geom_histogram()


####### LOF data ####### 
LOF_mut_t <- t(LOF_mut)

LOF_sums_0 <- as.data.frame(colSums(LOF_mut_t[response$response == 0,]/sum(response$response == 0)))
names(LOF_sums_0)
LOF_sums_0 %>% 
  rownames_to_column(.,"gene") %>% 
  rename(., sum_count = 'colSums(LOF_mut_t[response$response == 0, ]/sum(response$response == 0))') %>% 
  top_n(., 10, sum_count) %>% 
  arrange(.,desc(sum_count)) %>% 
  ggplot(.,aes(x=fct_reorder(gene,sum_count),y=sum_count)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


LOF_sums_1 <- as.data.frame(colSums(LOF_mut_t[response$response == 1,]/sum(response$response == 1)))
names(LOF_sums_1)
LOF_sums_1 %>% 
  rownames_to_column(.,"gene") %>% 
  rename(., sum_count = 'colSums(LOF_mut_t[response$response == 1, ]/sum(response$response == 1))') %>% 
  top_n(., 10, sum_count) %>% 
  arrange(.,desc(sum_count)) %>% 
  ggplot(.,aes(x=fct_reorder(gene,sum_count),y=sum_count)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




lof_1_clust <- hclust(dist(LOF_mut_t[response$response == 1,]))

plot(lof_1_clust)
