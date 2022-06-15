pheno_factor <- pheno %>% 
  mutate_all(factor)



pheno$malignant_neoplasm_mitotic_count_rate
features_of_interest <- c("age_at_initial_pathologic_diagnosis","breslow_depth_value","age_at_initial_pathologic_diagnosis","days_to_submitted_specimen_dx","malignant_neoplasm_mitotic_count_rate") # submitter_id.samples
hist.data.frame(pheno[,features_of_interest])
pheno$new_tumor_dx_prior_submitted_specimen_dx

ggplot(data = pheno, aes(x=factor(history_of_neoadjuvant_treatment))) + geom_bar()
ggplot(data = pheno, aes(x=factor(melanoma_clark_level_value))) + geom_bar()
ggplot(data = pheno, aes(x=factor(melanoma_origin_skin_anatomic_site))) + geom_bar()
ggplot(data = pheno, aes(x=factor(melanoma_ulceration_indicator))) + geom_bar()
ggplot(data = pheno, aes(x=factor(new_tumor_dx_prior_submitted_specimen_dx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(other_dx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(pathologic_M))) + geom_bar()
ggplot(data = pheno, aes(x=factor(pathologic_N))) + geom_bar()
ggplot(data = pheno, aes(x=factor(pathologic_T))) + geom_bar()
ggplot(data = pheno, aes(x=factor(person_neoplasm_cancer_status))) + geom_bar()
ggplot(data = pheno, aes(x=factor(postoperative_rx_tx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(primary_melanoma_at_diagnosis_count))) + geom_bar()
ggplot(data = pheno, aes(x=factor(primary_neoplasm_melanoma_dx))) + geom_bar()
ggplot(data = pheno, aes(x=factor(primary_tumor_multiple_present_ind))) + geom_bar()
ggplot(data = pheno, aes(x=factor(prior_radiation_therapy))) + geom_bar()
ggplot(data = pheno, aes(x=factor(prior_systemic_therapy))) + geom_bar()
ggplot(data = pheno, aes(x=factor(prior_systemic_therapy_type))) + geom_bar()
ggplot(data = pheno, aes(x=factor(radiation_therapy))) + geom_bar()
ggplot(data = pheno, aes(x=factor(radiation_therapy_to_primary))) + geom_bar()
ggplot(data = pheno, aes(x=factor(submitted_tumor_location))) + geom_bar() 
ggplot(data = pheno, aes(x=factor(subsequent_primary_melanoma_during_followup))) + geom_bar() 
ggplot(data = pheno, aes(x=factor(subsequent_primary_melanoma_during_followup))) + geom_bar() 

pheno_factor <- apply(pheno,2, as.factor)
class(pheno_factor)
View(pheno_factor)
colnames(pheno_factor) <- colnames(pheno)
pheno_factor <- as.data.frame(pheno_factor)
class(pheno_factor$batch_number)
summary(pheno_factor)

class(pheno)


factor(pheno)

View(head(pheno))

colnames(pheno_factor)


names(pheno_factor)

for (feature in names(pheno_factor)){
  plot <- ggplot(data = pheno_factor, aes_string(x=feature)) + 
    geom_bar()
  ggsave(filename = paste("./tmp/",feature,"_bar",".png",sep=""),plot)
}

for (feature in colnames(pheno)){
  plot <- ggplot(data = pheno, aes_string(x=factor(pheno[,feature]))) + 
    geom_histogram()
  ggsave(filename = paste("./tmp/",feature,"_hist",".png",sep=""),plot)
}
