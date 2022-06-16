library('tidyverse')
library("MCPcounter")

####### TPM data ####### 
estimate_output <- myEstimateFunction(tpm_ensg)
MCP_output <- MCPcounter.estimate(tpm_ensg,featuresType = "HUGO_symbols")
MCP_output
## Fitting linear models to estimate output
estimate_output
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
  mutate(mdl = map(data, ~glm(purity ~ counts,
                              data = .x)))  

grouped_data_sample_tidy <- tpm_ensg_t_grouped_data %>% 
  mutate(tidy = map(mdl, broom::tidy)) %>% 
  unnest(tidy)

ggplot(data = estimate_output, aes(x=ESTIMATEScore)) + geom_histogram()


####### LOF data ####### 


