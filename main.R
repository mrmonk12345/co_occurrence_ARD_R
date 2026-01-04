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
ps_all  <- phyloseq(OTU, TAX, SAM)
ps0 <- subset_samples(ps_all, treatment == "BjSa")
ps0  <- prune_samples(sample_sums(ps0) > 0, ps0)

# Decontaminate
ps1<-ps0 #not doing this yet


# Parameters to tune

min_depth  <- 0    # minimum library size per sample (example)
min_total  <- 0       # ASV must have at least X total counts across cohort


ps2 <- ps1 %>%
  prune_samples(sample_sums(.) >= min_depth, .) %>%
  prune_taxa(taxa_sums(.) >= min_total, .)


min_preval <- 0.25 # Set your threshold: keep taxa present in ≥5% of samples

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


# SPIEC-EASI

set.seed(1)

# SPIEC-EASI expects taxa as columns and samples as rows; phyloseq-friendly wrapper handles that
# Parameters to tune:
nlambda          <- 30           # length of regularization path
lambda_min_ratio <- 1e-2         # smallest lambda as a fraction of max
rep_num          <- 50           # stability selection resamples (increase for larger data)
sel_crit         <- "stars"      # stability criterion; alternatives: "bstars"
mb_alpha         <- 0.05         # target instability threshold for STARS (default 0.05)

se.mb <- spiec.easi(ps2,
                    method          = "mb",
                    sel.criterion   = sel_crit,
                    pulsar.params   = list(rep.num = rep_num, thresh = mb_alpha),
                    nlambda         = nlambda,
                    lambda.min.ratio= lambda_min_ratio,
                    verbose         = TRUE)

# Extract adjacency and edge weights
adj  <- getRefit(se.mb)                       # 0/1 adjacency (symmetric)
beta <- symBeta(getOptBeta(se.mb), mode="maxabs")  # weighted edges (partial associations)

# Build igraph
g_asv <- adj2igraph(adj, vertex.attr = list(name = taxa_names(ps2)))
E(g_asv)$weight <- beta[upper.tri(beta)][which(adj[upper.tri(adj)] == 1)]
cat(sprintf("ASV network: %d nodes, %d edges\n", gorder(g_asv), gsize(g_asv)))

