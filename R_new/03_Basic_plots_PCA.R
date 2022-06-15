# DATA VISUALISATION

#-------------------------------------------------------------------------------
# Load libraries
library(tidyverse)

#-------------------------------------------------------------------------------
# Basic Plots

# Response
response_char <- response %>% 
  mutate(response = case_when(response == 1 ~ "Response",
                              response == 0 ~ "No Response"))
response_char %>% 
  ggplot(mapping = aes(x = response,
                       fill = response)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Distribution of response vs. no response")

# Purity
estimate_output %>% 
  ggplot(mapping = aes(x = TumorPurity)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Tumor Purity")

# Mut load
response %>% 
  ggplot(mapping = aes(x = mut_load)) +
  geom_density() + 
  theme_minimal() +
  labs(title = "Mutational load pr. patient")

#-------------------------------------------------------------------------------
# PCA

# tpm
pca_tpm <- prcomp(t(tpm_ensg))

pca_tpm %>% broom::tidy(matrix = "eigenvalues")

pca_tpm %>% 
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

# estimate output
pca_estimate <- prcomp(estimate_output)

pca_estimate %>% broom::tidy(matrix = "eigenvalues")

pca_estimate %>% 
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

# LOF

pca_LOF <- prcomp(t(LOF_mut))

pca_LOF %>% broom::tidy(matrix = "eigenvalues")

pca_LOF %>% 
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
