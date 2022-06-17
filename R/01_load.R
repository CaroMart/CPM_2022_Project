# Set timeout for downloads
options(timeout = 10000000)

# Download data ----------------------------------------------------------------
## SNVs
download.file(
  url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.mutect2_snv.tsv.gz",
  destfile = "data/_raw/TCGA-SKCM.mutect2_snv.tsv.gz"
)

## Methylation
# download.file(
#   url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.methylation450.tsv.gz",
#   destfile = "data/_raw/TCGA-SKCM.methylation450.tsv.gz"
# )

## Counts
download.file(
  url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.htseq_counts.tsv.gz",
  destfile = "data/_raw/TCGA-SKCM.htseq_counts.tsv.gz"
)

## CNVs
download.file(
  url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.masked_cnv.tsv.gz",
  destfile = "data/_raw/TCGA-SKCM.masked_cnv.tsv.gz"
)

## Phenotype
download.file(
  url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.GDC_phenotype.tsv.gz",
  destfile = "data/_raw/TCGA-SKCM.GDC_phenotype.tsv.gz"
)

## Survival
download.file(
  url = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-SKCM.survival.tsv",
  destfile = "data/_raw/TCGA-SKCM.survival.tsv"
)
