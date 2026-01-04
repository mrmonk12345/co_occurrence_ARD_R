
## ---- Install required packages ----

# Install CRAN packages
cran_packages <- c(
  "tidyverse",
  "data.table",
  "Matrix",
  "igraph"
)

install.packages(setdiff(cran_packages, installed.packages()[, "Package"]))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "phyloseq",
  "Biostrings",
  "decontam"
)

BiocManager::install(setdiff(bioc_packages, installed.packages()[, "Package"]))

# Install SpiecEasi (CRAN + GitHub fallback)
if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
  tryCatch({
    install.packages("SpiecEasi")
  }, error = function(e) {
    message("Installing SpiecEasi from GitHub...")
    BiocManager::install("zdk123/SpiecEasi")
  })
}

## ---- Load all libraries ----

library(tidyverse)
library(data.table)
library(phyloseq)
library(Biostrings)
library(decontam)
library(SpiecEasi)
library(Matrix)
library(igraph)
