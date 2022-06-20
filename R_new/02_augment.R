# DATA AUGMENTATION

#-------------------------------------------------------------------------------
# Load Libraries
library("MCPcounter")
library("tidyverse")

source("R_new/estimate_function.R")
source("R_new/LARS_estimate_function.R")

#-------------------------------------------------------------------------------
####### TPM data ####### 
estimate_output_lars <- LARSEstimateFunction(tpm_ensg)
ggplot(estimate_output_lars,aes(x=ESTIMATEScore, y=TumorPurity)) + geom_point()
summary(estimate_output_lars$ESTIMATEScore)

estimate_output <- myEstimateFunction(tpm_ensg)

hist(estimate_output$TumorPurity)

MCP_output <- MCPcounter.estimate(tpm_ensg,featuresType = "HUGO_symbols")
MCP_output


summary(estimate_output$ESTIMATEScore)

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



mad_values <- apply(tpm_ensg, 1, mad)
# mad_tpm_ensg <- tpm_ensg[order(mad_values)[1:5000], ]
mad_tpm_ensg <- tpm_ensg[mad_values > 1.696, ]


mad_values_df %>% 
  rownames_to_column(.)  %>% 
  filter(mad_values > 0) %>% 
  ggplot(.,aes(x=log10(mad_values))) + geom_density() + geom_vline(xintercept = 1.696, linetype="dotted", 
                                                                   color = "blue")

summary(mad_values_df)
clustering_mad_tpm <- hclust(dist(t(mad_tpm_ensg)))
plot(clustering_mad_tpm)

#-------------------------------------------------------------------------------
####### LOF data ####### 
LOF_mut_t <- t(LOF_mut)
LOF_mut_t_filtered <- LOF_mut_t[,colSums(LOF_mut_t) > 1]
dim(LOF_mut_t_filtered)
clustering_LOF_filtered <- hclust(dist(LOF_mut_t_filtered))
plot(clustering_LOF_filtered)
heatmap(LOF_mut_t_filtered)

LOF_sums_0 <- as.data.frame(colSums(LOF_mut_t_filtered[response$response == 0,]/sum(response$response == 0)))
names(LOF_sums_0)
LOF_sums_0 %>% 
  rownames_to_column(.,"gene") %>% 
  rename(., sum_count = 'colSums(LOF_mut_t_filtered[response$response == 0, ]/sum(response$response == 0))') %>% 
  top_n(., 10, sum_count) %>% 
  arrange(.,desc(sum_count)) %>% 
  ggplot(.,aes(x=fct_reorder(gene,sum_count),y=sum_count)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


LOF_sums_1 <- as.data.frame(colSums(LOF_mut_t_filtered[response$response == 1,]/sum(response$response == 1)))
names(LOF_sums_1)
LOF_sums_1 %>% 
  rownames_to_column(.,"gene") %>% 
  rename(., sum_count = 'colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1))') %>% 
  top_n(., 10, sum_count) %>% 
  arrange(.,desc(sum_count)) %>% 
  ggplot(.,aes(x=fct_reorder(gene,sum_count),y=sum_count)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


as.data.frame(LOF_sums_1 - LOF_sums_0) %>% 
  rownames_to_column("gene") %>% 
  rename(diff_freq = 'colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1))') %>% 
  ggplot(.,aes(x=diff_freq)) + geom_density() + geom_vline(xintercept = -0.12, linetype="dotted", 
                                                           color = "blue", size=1.5)

more_common_LOF_res <- as.data.frame(LOF_sums_1 - LOF_sums_0) %>% 
  rownames_to_column("gene") %>% 
  rename(diff_freq = 'colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1))') %>% 
  filter(diff_freq > 0.15)

dim(more_common_LOF_res)
write.table(more_common_LOF_res$gene,"./tmp/more_common_LOF_res.csv")

more_common_LOF_nonres <- as.data.frame(LOF_sums_1 - LOF_sums_0) %>% 
  rownames_to_column("gene") %>% 
  rename(diff_freq = 'colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1))') %>% 
  filter((diff_freq < -0.12))

more_common_LOF_nonres

write.table(common_LOF_res_non_res$gene,"./tmp/common_LOF_res_non_res.csv")


common_LOF_res_non_res <- as.data.frame(LOF_sums_1 - LOF_sums_0) %>% 
  rownames_to_column("gene") %>% 
  rename(diff_freq = 'colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1))') %>% 
  mutate(freq1 = colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1)),
         freq0 = colSums(LOF_mut_t_filtered[response$response == 0, ]/sum(response$response == 0))) %>% 
  filter((diff_freq > -0.12) & 
           (diff_freq < 0.12) & 
           (freq1 > 0.0) & 
           (freq0 > 0.0))

LOF_data <- as.data.frame(LOF_sums_1 - LOF_sums_0) %>% 
  rownames_to_column("gene") %>% 
  rename(diff_freq = 'colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1))') %>% 
  mutate(freq1 = colSums(LOF_mut_t_filtered[response$response == 1, ]/sum(response$response == 1)),
         freq0 = colSums(LOF_mut_t_filtered[response$response == 0, ]/sum(response$response == 0)))


LOF_data


common_LOF_res_non_res[common_LOF_res_non_res$diff_freq > 0.11,]

write.table(common_LOF_res_non_res$gene,"./tmp/common_LOF_res_non_res.csv")


row_mww_wilcox <- function(x) {
  subtype <- x[response$response == 1]
  rest <- x[!response$response == 1]
  res <- wilcox.test(subtype, rest)
  return(res$p.value)
}

wilcox_teston_lof <- apply(LOF_mut,1,row_mww_wilcox)
dim(wilcox_teston_lof)
response$response

View(head(LOF_mut_t))


dim(LOF_mut_t_filtered)
heatmap(LOF_mut_t_filtered)

#-------------------------------------------------------------------------------
####### Response data ####### 

response_char <- response %>% 
  mutate(response = case_when(response == 1 ~ "Response",
                              response == 0 ~ "No Response"),
         response_RECIST = case_when(RECIST == "CR" ~ "Response",
                                     RECIST == "PR" ~ "Response",
                                     RECIST == "SD" ~ "No Response",
                                     RECIST == "PD" ~ "No Response"))
