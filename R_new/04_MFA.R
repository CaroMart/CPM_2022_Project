library("FactoMineR")
library("missMDA")
library("tidyverse")
library("dplyr")

## For augment maybe

## Making data frame for MFA ##

# First, transpose all data frames so there is one patient/sample pr. row
# and so "patient" is a column in the data frame (for joining)

MCP_transposed <- as.data.frame(t(MCP_output))
MCP_transposed <- tibble::rownames_to_column(MCP_transposed, "patient")

tpm_ensg_transposed <- as.data.frame(t(tpm_ensg))
tpm_ensg_transposed <- tibble::rownames_to_column(tpm_ensg_transposed, "patient")

LOF_transposed <- as.data.frame(t(LOF_mut))
LOF_transposed <- tibble::rownames_to_column(LOF_transposed, "patient")

estimate_output <- tibble::rownames_to_column(estimate_output, "patient")

# Now, merge them all together

data_factominer <- MCP_transposed %>%
  left_join(tpm_ensg_transposed, by = "patient") %>%
  left_join(LOF_transposed, by = "patient") %>%
  left_join(estimate_output, by = "patient") %>%
  column_to_rownames(., var = "patient")

## MFA ##
data_factominer
View(head(data_factominer))
dim(data_factominer)
res_MFA <- MFA(data_factominer, group=c(10,35808,3468,4), type=c(rep("c", 4)),
           ncp=20, name.group=c("MCP_count","tpm_ensg","LOF","estimate"))

## Plotting result
PC_matrix <- res_MFA$ind
PC_matrix_df <- as.data.frame(PC_matrix$coord)

ggplot(PC_matrix_df,aes(x=Dim.1, y=Dim.2, color = response$response)) + 
  geom_point() + 
  stat_ellipse()

fosmp = c()
fosmp = c(fosmp, rownames(PC_matrix_df[PC_matrix_df$Dim.1<(-2),]))
fosmp = c(fosmp, rownames(PC_matrix_df[PC_matrix_df$Dim.2>2,]))

fosmp

data_factominer_removed_fosmp = data_factominer[!(rownames(data_factominer) %in% fosmp),]

res_MFA <- MFA(data_factominer_removed_fosmp, group=c(10,35808,3468,4), type=c(rep("c", 4)),
               ncp=20, name.group=c("MCP_count","tpm_ensg","LOF","estimate"))

## Plotting result
PC_matrix <- res_MFA$ind
PC_matrix_df <- as.data.frame(PC_matrix$coord)

ggplot(PC_matrix_df,aes(x=Dim.1, y=Dim.2, color = factor(response$response[!(rownames(data_factominer) %in% fosmp)]))) + 
  geom_point() + 
  stat_ellipse()

heatmap(PC_matrix$coord)


res_MFA$ind$coord
res_MFA

response$response

## Hierarchical ##
#hierar <- list(c(10,35808,3468,4), c(3,1))
#res_HMFA <- HMFA(data_factominer, H = hierar, type=c(rep("s",4)), graph = TRUE)

PC_matrix
heatmap(PC_matrix$coord)
clusters = hclust(dist(PC_matrix$coord))
clusters$labels <- response$RECIST
plot(clusters)
summary(res_MFA)
