# Set timeout for downloads
options(timeout = 100000)

# Download data ----------------------------------------------------------------
download.file(
  url = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.GBM.sampleMap%2FHiSeqV2.gz",
  destfile = "data/_raw/TCGA.GBM.sampleMap-HiSeqV2.gz"
)