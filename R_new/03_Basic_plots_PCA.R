# DATA VISUALISATION

#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)

#-------------------------------------------------------------------------------
# Basic Plots

# Response
response_dist <- response_char %>% 
  ggplot(mapping = aes(x = response,
                       fill = response)) +
  geom_bar(color = "black") +
  theme_minimal() +
  labs(title = "Distribution of response vs. no response")

ggsave(
  filename = "03_response_dist.png",
  plot = response_dist,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# Purity
TumorPurity_dens <- estimate_output %>% 
  ggplot(mapping = aes(x = TumorPurity,
                       fill = "blue")) +
  geom_density(bw = 0.03) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Tumor Purity") +
  ylab("Density") + 
  labs(title = "Tumor Purity")

TumorPurity_dens

ggsave(
  filename = "03_TumorPurity.png",
  plot = TumorPurity_dens,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# Mut load (LIGEGYLDIGT)
Mut_load_dens <- response %>% 
  ggplot(mapping = aes(x = mut_load,group = factor(response) ,fill=factor(response))) +
  geom_histogram() + 
  theme_minimal() +
  labs(title = "Mutational load pr. patient")

ggsave(
  filename = "03_MutLoad_dens.png",
  plot = Mut_load_dens,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# Purity for non-responders and responders
library("ggpubr")


# bleu: #53BCC2
# rot: #E87E72

dim(estimate_output)

estimate_response <- estimate_output %>% 
  rownames_to_column() %>% 
  inner_join(response_char, by = c("rowname" = "patient")) %>% 
  mutate(response = case_when(response == "Response" ~ "Response (n=10)",
                              response == "No Response" ~ "No Response (n=12)"))

p <- ggboxplot(estimate_response, x = "response", y = "TumorPurity",
               color = "response", palette = c("#E87E72", "#53BCC2"),
               add = "jitter")
#  Add p-value
p + stat_compare_means()

p <- ggboxplot(estimate_response, x = "response", y = "mut_load",
               color = "response", palette = c("#E87E72", "#53BCC2"),
               add = "jitter", xlab = "Response to ACT", ylab = "Tumor mutational burden")
#  Add p-value
ggpar(p, legend.title = "Response", legend = "none") + stat_compare_means()

Purity_respond_dens <- estimate_output %>% 
  ggplot(mapping = aes(x = TumorPurity,
                       fill = "blue")) +
  geom_histogram(bins = 15, color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~response_char$response) +
  labs(title = "Tumor purity stratified on responders vs. non-responders")

ggsave(
  filename = "03_Purity_respond.png",
  plot = Purity_respond_dens,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

grouped_data_sample_tidy_beta <- tpm_ensg_t_grouped_data %>% 
  mutate(tidy = map(mdl, function(x) {broom::tidy(x)[2,]})) %>% 
  unnest(tidy)

grouped_data_sample_tidy_beta %>% top_n(10, -p.value)
grouped_data_sample_tidy_beta %>% filter(estimate > 0) %>% top_n(10, -p.value)

one_gene_glm_model = tpm_ensg_t_grouped_data_calc %>% filter(genes == "MLANA") %>% pluck(3) %>% .[[1]]
one_gene_data = tpm_ensg_t_grouped_data_calc %>% filter(genes == "MLANA") %>% pluck(2) %>% .[[1]]

one_gene_data %>% ggplot(mapping = aes(x = purity, y = counts, fill = "MLANA")) +
  geom_point() +
  geom_smooth(method = "glm", col = "red", se = FALSE) +
  geom_hline(yintercept = sum(one_gene_glm_model$coef), linetype="dashed") +
  theme_minimal() + 
  xlab("Purity") +
  ylab("Transcripts per million of MLANA") + 
  labs(fill = "Gene GLM") +
  theme(legend.position = "bottom")

one_gene_glm_model = tpm_ensg_t_grouped_data_calc %>% filter(genes == "PLAGL1") %>% pluck(3) %>% .[[1]]
one_gene_data = tpm_ensg_t_grouped_data_calc %>% filter(genes == "PLAGL1") %>% pluck(2) %>% .[[1]]

one_gene_data %>% ggplot(mapping = aes(x = purity, y = counts, fill = "PLAGL1")) +
  geom_point() +
  geom_smooth(method = "glm", col = "red", se = FALSE) +
  geom_hline(yintercept = sum(one_gene_glm_model$coef), linetype="dashed") +
  theme_minimal() + 
  xlab("Purity") +
  ylab("Transcripts per million of PLAGL1") + 
  labs(fill = "Gene GLM") +
  theme(legend.position = "bottom")

one_gene_glm_model = tpm_ensg_t_grouped_data_calc %>% filter(genes == "CARD11") %>% pluck(3) %>% .[[1]]
one_gene_data = tpm_ensg_t_grouped_data_calc %>% filter(genes == "CARD11") %>% pluck(2) %>% .[[1]]

one_gene_data %>% ggplot(mapping = aes(x = purity, y = counts, fill = "CARD11")) +
  geom_point() +
  geom_smooth(method = "glm", col = "red", se = FALSE) +
  geom_hline(yintercept = sum(one_gene_glm_model$coef), linetype="dashed") +
  theme_minimal() + 
  xlab("Purity") +
  ylab("Transcripts per million of CARD11") + 
  labs(fill = "Gene GLM") +
  theme(legend.position = "bottom")

# Mut_load for non-responders and responders
Mut_load_respond_dens <- response_char %>% 
  ggplot(mapping = aes(x = mut_load,
                       fill = "blue")) +
  geom_boxplot(color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~response) +
  labs(title = "Tumor Mutational Burden",
       x = "Tumor Mutational Burden (TMB)")

# Same men som boxplot, skal laves som specielt boxplot
Mut_load_respond_dens <- response_char %>% 
  ggplot(mapping = aes(x = response,
                       y = mut_load)) +
  geom_boxplot(color = "black", fill = "#e87e72") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Tumor Mutational Burden",
       y = "Tumor Mutational Burden (TMB)")

Mut_load_respond_dens

ggsave(
  filename = "03_Mut_load_respond_dens.png",
  plot = Mut_load_respond_dens,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

#-------------------------------------------------------------------------------
# PCA

# tpm
pca_tpm <- prcomp(t(tpm_ensg))

pca_tpm %>% broom::tidy(matrix = "eigenvalues")


pca_tpm_plot <- pca_tpm %>% 
  broom::augment((response_char)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = response)) + 
  geom_point() + 
  theme_minimal() +
  labs(title = "PCA for tpm",
       x = "PC1 [42%]",
       y = "PC2 [29.4%]") +
  stat_ellipse()

pca_tpm_plot
ggsave(
  filename = "03_PCA_tpm.png",
  plot = pca_tpm_plot,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# Outlier analysis
pca_tpm$rotation %>%
  as.data.frame() %>% 
  ggplot(aes(x = log10(PC1))) +
  geom_histogram()

ct <- hclust(dist(pca_tpm$x))

plot(ct)

tpm_PC1 <- as.data.frame(pca_tpm$rotation[,2])

cdotpro <- apply(tpm_ensg,2, function(x) {x*tpm_PC1})

cdotpro <- as.data.frame(cdotpro)

colnames(cdotpro) <- colnames(tpm_ensg)


rownames(cdotpro$`pca_tpm$rotation[, 1]`)[order(abs(cdotpro$`pca_tpm$rotation[, 1]`), decreasing = TRUE)]

mean_expr_no <- cdotpro %>% 
  rownames_to_column("gene") %>% 
  select(-c("MM909_14","MM909_22","MM909_26")) %>%
  pivot_longer(!gene, names_to = "patients", values_to = "score") %>% 
  group_by(gene) %>% 
  summarise(mean_ = mean(score))

cdotpro %>% 
  rownames_to_column("gene") %>% 
  select(c("MM909_14","MM909_22","MM909_26")) %>%
  mutate(across(-matches("gene"),~. -mean_expr_no$mean_)) %>% 
  mutate(gene = mean_expr_no$gene) %>% 
  pivot_longer(!gene, names_to = "outlier", values_to = "score") %>% 
  ggplot(aes(x = log10(score), color = outlier)) + geom_density()

o <- c("MM909_14","MM909_22","MM909_26")


colorforfun <- c("lightpink" , "lightpink","lightpink","lightpink","black","lightpink","lightpink","lightpink","black", "lightpink","black","lightpink","lightpink","lightpink","lightpink","lightpink","lightpink","lightpink","lightpink","lightpink","lightpink","lightpink")

cdotpro %>% 
  rownames_to_column("gene") %>%
  pivot_longer(!gene, names_to = "patients", values_to = "score") %>%
  mutate(abs_val = abs(score)) %>% 
  mutate(outlier = case_when(patients %in% o ~ 1,
                   !(patients %in% o) ~ 0)) %>% 
  filter(abs_val>1) %>% 
  ggplot(aes(x = log10(abs_val), 
             color = patients)) +
  geom_density() +
  theme_minimal() +
  scale_color_manual(values = colorforfun)


# estimate output
pca_estimate <- prcomp(estimate_output)

pca_estimate %>% broom::tidy(matrix = "eigenvalues")

pca_estimate_plot <- pca_estimate %>% 
  broom::augment((response_char)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = response)) + 
  geom_point() + 
  theme_minimal() +
  labs(title = "PCA for estimate",
       x = "PC1 [97.8%]",
       y = "PC2 [2.8%]") +
  stat_ellipse()

ggsave(
  filename = "03_PCA_estimate.png",
  plot = pca_estimate_plot,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# LOF

pca_LOF <- prcomp(t(LOF_mut))

pca_LOF %>% broom::tidy(matrix = "eigenvalues")

PCA_LOF_plot <- pca_LOF %>% 
  broom::augment((response_char)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = response)) + 
  geom_point() + 
  theme_minimal() +
  labs(title = "PCA for LOF",
       x = "PC1 [19.9%]",
       y = "PC2 [16.5%]") +
  stat_ellipse()

ggsave(
  filename = "03_PCA_LOF.png",
  plot = PCA_LOF_plot,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# MCP counter

MCP_output %>% t %>% 
  as.data.frame %>% 
  rownames_to_column %>% 
  left_join((response_char %>% select(patient, response)), by = c("rowname" = "patient")) %>% 
  pivot_longer(-c(rowname, response)) %>% 
  mutate(value = log(value))%>% 
  ggplot(mapping = aes(x = name, y = value, color=response, group=rowname)) +
    geom_line() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("log(MCP-counter)") +
    xlab("Cell type") +
    labs(color="Patient")

#Removing _026, as it is an outlier
MCP_output_pca <- MCP_output %>% 
  subset(., select = -c(MM909_26))

response_char_MCP <- response_char %>% 
  filter(patient != "MM909_26")

pca_MCP <- prcomp(t(MCP_output_pca))

pca_MCP %>% broom::tidy(matrix = "eigenvalues")

pca_MCP_plot <- pca_MCP %>% 
  broom::augment((response_char_MCP)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = response)) + 
  geom_point() + 
  theme_minimal() +
  labs(title = "PCA for MCP",
       x = "PC1 [%]",
       y = "PC2 [%]") +
  stat_ellipse()

ggsave(
  filename = "03_PCA_MCP.png",
  plot = pca_MCP_plot,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)
#-------------------------------------------------------------------------------