library("FactoMineR")
library("missMDA")
library("tidyverse")
library("dplyr")

## For augment maybe

# Making mean-absolute-deviation (mad) cleansing on tpm_ensg data 
mad_values <- apply(tpm_ensg, 1, mad)
mad_tpm_ensg <- tpm_ensg[order(mad_values)[1:5000], ]

# Subsetting to only include genes that are also in mad_tpm_ensg

purity_corrected_tpm <- tpm_ensg_t_grouped_data_calc_wider %>%
  mutate(patient = sample_nr) %>%
  select(-sample_nr) %>%
  select(patient, row.names(mad_tpm_ensg))

## Making data frame for MFA ##

# First, transpose all data frames so there is one patient/sample pr. row
# and so "patient" is a column in the data frame (for joining)

MCP_transposed <- as.data.frame(t(MCP_output))
MCP_transposed <- tibble::rownames_to_column(MCP_transposed, "patient")

mad_tpm_ensg_transposed <- as.data.frame(t(mad_tpm_ensg))
mad_tpm_ensg_transposed <- tibble::rownames_to_column(mad_tpm_ensg_transposed, "patient")

LOF_transposed <- as.data.frame(t(LOF_mut))
LOF_transposed <- tibble::rownames_to_column(LOF_transposed, "patient")

estimate_output <- tibble::rownames_to_column(estimate_output, "patient")

# Now, merge them all together
# MCPcount udregner kompositionen af celler i kræften baseret på gene expression
# ESTIMATE udregner bl.a. purity, som er andelen af stromal og immune cells
# LOF angiver om et givent gen har mistet funktionen (loss of function) i patienten 
# tpm_ensg transcript pr million 
 
 
data_factominer <- MCP_transposed %>%
  left_join(mad_tpm_ensg_transposed, by = "patient") %>%
  left_join(LOF_transposed, by = "patient") %>%
  left_join(purity_corrected_tpm, by = "patient") %>% 
  left_join(response[,c("patient", "mut_load")], by = "patient") %>%
  cbind(., purity_scores) %>%
  column_to_rownames(., var = "patient")

## MFA ##
dim(data_factominer)
res_MFA <- MFA(data_factominer, group=c(10,5000,3468,5000,1,1), type=c(rep("c", 6)),
           ncp=20, name.group=c("MCP_count","mad_tpm_ensg","LOF","purity_corrected_tpm","TMB", "purity_score"))

## Plotting result
PC_matrix <- res_MFA$ind
PC_matrix_df <- as.data.frame(PC_matrix$coord)

ggplot(PC_matrix_df,aes(x=Dim.1, y=Dim.2, color = response$response)) + 
  geom_point() + 
  stat_ellipse()

fosmp = c()
fosmp = c(fosmp, rownames(PC_matrix_df[PC_matrix_df$Dim.1>(3),]))
fosmp = c(fosmp, rownames(PC_matrix_df[PC_matrix_df$Dim.2>2.5,]))

fosmp

# Running again, having removed three outliers
# OBS look at plot before to determine new outliers!
data_factominer_removed_fosmp = data_factominer[!(rownames(data_factominer) %in% fosmp),]

res_MFA_fosmp <- MFA(data_factominer_removed_fosmp, group=c(10,5000,3468,5000,1,1), type=c(rep("c", 6)),
               ncp=20, name.group=c("MCP_count","mad_tpm_ensg","LOF","purity_corrected_tpm","TMB", "purity_score"))


## Plotting result
PC_matrix_fosmp <- res_MFA_fosmp$ind
PC_matrix_fosmp_df <- as.data.frame(PC_matrix_fosmp$coord)

ggplot(PC_matrix_fosmp_df,aes(x=Dim.1, y=Dim.2, color = factor(response$response[!(rownames(data_factominer) %in% fosmp)]))) + 
  geom_point() + 
  stat_ellipse()

heatmap(PC_matrix_fosmp$coord)


res_MFA_fosmp$ind$coord
res_MFA_fosmp

response$response

## Hierarchical ##
#hierar <- list(c(10,35808,3468,4), c(3,1))
#res_HMFA <- HMFA(data_factominer, H = hierar, type=c(rep("s",4)), graph = TRUE)

PC_matrix_fosmp
heatmap(PC_matrix_fosmp$coord)
clusters = hclust(d=dist(PC_matrix_fosmp$coord))
clusters$labels <- response$response
plot(clusters)
summary(res_MFA_fosmp)
