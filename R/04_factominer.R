library("FactoMineR")
library("missMDA")


# ncp is eqivalent to number of PC's in a PCA
# num.group.sup can be used to indicate which groups should NOT be included in MFA

data("wine")
res <- MFA(wine, group=c(2,5,3,10,9,2), type=c("n",rep("s",5)),
           ncp=5, name.group=c("orig","olf","vis","olfag","gust","ens"),
           num.group.sup=c(1,6))
summary(res)
barplot(res$eig[,1],main="Eigenvalues",names.arg=1:nrow(res$eig))


data(orange)
res.impute <- imputeMFA(orange, group=c(5,3), type=rep("s",2),ncp=2) 
res.mfa <- MFA(res.impute$completeObs,group=c(5,3),type=rep("s",2)) 
summary(res.mfa)
barplot(res.mfa$eig[,1],main="Eigenvalues",names.arg=1:nrow(res.mfa$eig))

#### Testing factoMineR on small dataframe ####

## Creating small dataframe
## Only keeping samples for which we have survival data 

snvs_factominer_test <- snv_summary %>% 
  mutate(sample = Sample_ID) %>%
  select(-Sample_ID) %>%
  inner_join(survival_filter[1], by = "sample")

cnvs_factominer_test <- cnv_summary %>% 
  inner_join(survival_filter[1], by = "sample")

pheno_factominer_test <- pheno_data_filtered %>% 
  mutate(sample = submitter_id.samples) %>%
  select(sample, 
         ethnicity.demographic,
         gender.demographic,
         race.demographic,
         year_of_birth.demographic,
         age_at_diagnosis.diagnoses,
         primary_diagnosis.diagnoses,
         prior_malignancy.diagnoses,
         prior_treatment.diagnoses,
         tumor_stage.diagnoses,
         name.tissue_source_site) %>%
  inner_join(survival_filter[1], by = "sample")

data_factominer_test <- cnvs_factominer_test %>% 
  inner_join(snvs_factominer_test, by = "sample") %>%
  inner_join(pheno_factominer_test, by = "sample") %>%
  column_to_rownames(., var = "sample")


### running imputation to fill out e.g. BMI 
#res.impute <- imputeMFA(data_factominer_test, group=c(19210,46,11), type=c(rep("s",19210),rep("s",46),rep("n",11)),ncp=5)

res <- MFA(data_factominer_test, group=c(46,19210,10), type=c(rep("s",46), rep("s",19210), rep("n",4), "s", rep("n",5)),
    ncp=5, name.group=c("cnv","snv","pheno"))
