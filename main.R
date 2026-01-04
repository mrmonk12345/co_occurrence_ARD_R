# import libraries

# Core data structures and I/O
library(tidyverse)
library(data.table)
library(phyloseq)     # integrated microbiome objects
library(Biostrings)   # sequences when needed

## Decontamination
#library(decontam)     # frequency- and prevalence-based contaminant detection

# Network inference
library(SpiecEasi)    # SPIEC-EASI main methods (mb/glasso) + pulsar stability
library(Matrix)
library(igraph)


# Replace with your file paths / objects
asv_tab_path   <- "data/raw/project_otu.csv"      # rows = samples, cols = ASVs
tax_tab_path   <- "data/raw/taxon_tab.csv"        # ASV, Kingdom,..., Genus, Species
metadata_path  <- "data/raw/16s_sample_metadata_file.tsv" # includes: SampleID, DNA_conc, is_neg (TRUE/FALSE)

# Read
asv_df <- read.csv(asv_tab_path, header = TRUE, row.names = 1, check.names = FALSE)
taxon_df <- read.csv(tax_tab_path, header = TRUE, row.names = 1, check.names = FALSE)

meta <- read.delim(
  metadata_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1
)


asv_matrix <- asv_df %>% as.matrix()
tax_matrix <- taxon_df %>% as.matrix()


# Build a phyloseq object
OTU  <- otu_table(asv_matrix, taxa_are_rows = TRUE)  # samples are rows
TAX  <- tax_table(tax_matrix)
SAM  <- sample_data(meta)
ps0  <- phyloseq(OTU, TAX, SAM)
ps0  <- prune_samples(sample_sums(ps0) > 0, ps0)

# Decontaminate
ps1<-ps0 #not doing this yet
ps2<-ps1

# Parameters to tune

# Set your threshold: keep taxa present in ≥5% of samples
min_preval <- 0.05

# Extract the OTU/ASV abundance table
otu <- otu_table(ps2)

# Compute prevalence per taxon (fraction of samples with non-zero counts)
# Respect 'taxa_are_rows(ps2)' to choose the right margin:
prev <- apply(otu,
              MARGIN = ifelse(taxa_are_rows(ps2), 1, 2),
              FUN    = function(x) mean(x > 0))

# Make sure the prevalence vector is named by taxa IDs (not sample IDs)
names(prev) <- taxa_names(ps2)

# Option A (safest): pass a CHARACTER vector of taxa IDs to keep
keep_ids <- names(prev)[prev >= min_preval]
ps2 <- prune_taxa(keep_ids, ps2)

# Option B: if you prefer a LOGICAL vector, align its order to taxa_names(ps2)
# logic_keep <- (prev >= min_preval)[ taxa_names(ps2) ]
# ps2 <- prune_taxa(logic_keep, ps2)

# Sanity check
stopifnot(ntaxa(ps2) > 0, nsamples(ps2) > 0)
cat(sprintf("After prevalence filtering: %d samples, %d taxa\n", nsamples(ps2), ntaxa(ps2)))


ntaxa(ps0)
length(prev)
head(taxa_names(ps2))
head(names(prev))



