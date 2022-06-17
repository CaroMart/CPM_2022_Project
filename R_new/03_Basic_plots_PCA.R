# DATA VISUALISATION

#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)

#-------------------------------------------------------------------------------
# Basic Plots

# Response
response_char <- response %>% 
  mutate(response = case_when(response == 1 ~ "Response",
                              response == 0 ~ "No Response"),
         response_RECIST = case_when(RECIST == "CR" ~ "Response",
                                     RECIST == "PR" ~ "Response",
                                     RECIST == "SD" ~ "No Response",
                                     RECIST == "PD" ~ "No Response"))


response_dist <- response_char %>% 
  ggplot(mapping = aes(x = response,
                       fill = response)) +
  geom_bar() +
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
  ggplot(mapping = aes(x = TumorPurity)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Tumor Purity")

ggsave(
  filename = "03_TumorPurity_dens.png",
  plot = TumorPurity_dens,
  path = "results/",
  scale = 1,
  width = 5,
  height = 5,
  units = "in",
  dpi = 300
)

# Mut load
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
Purity_respond_dens <- estimate_output %>% 
  ggplot(mapping = aes(x = TumorPurity)) +
  geom_density() +
  theme_minimal() +
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

# Mut_load for non-responders and responders
Mut_load_respond_dens <- response_char %>% 
  ggplot(mapping = aes(x = mut_load)) +
  geom_density() +
  theme_minimal() +
  facet_wrap(~response) +
  labs(title = "Mut_load")

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

# tpm (Kører på RECIST)
pca_tpm <- prcomp(t(tpm_ensg))

pca_tpm %>% broom::tidy(matrix = "eigenvalues")

pca_tpm_plot <- pca_tpm %>% 
  broom::augment((response_char)) %>% 
  ggplot(mapping = aes(
    x = .fittedPC1,
    y = .fittedPC2,
    color = response_RECIST)) + 
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

