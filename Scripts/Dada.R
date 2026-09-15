# -------------------------------
# Load packages
# -------------------------------
library(dada2)
library(phyloseq)
library(tidyverse)

# -------------------------------
# Define paths
# -------------------------------
path <- "~/RNA16S/fastq_files/"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# -------------------------------
# Filter and trim
# -------------------------------
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

sample.names <- sub("_S[0-9]+_.*", "", basename(fnFs))

all.equal(
  sub("_R1_001.fastq$", "", basename(fnFs)),
  sub("_R2_001.fastq$", "", basename(fnRs))
)
stopifnot(length(fnFs) == 194)
stopifnot(all.equal(
  sub("_R1_001.fastq$", "", basename(fnFs)),
  sub("_R2_001.fastq$", "", basename(fnRs))
))

plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen=c(240,160), # adjust based on your quality plots
  maxN=0,
  maxEE=c(2,2),
  truncQ=2,
  rm.phix=TRUE,
  compress=TRUE,
  multithread=TRUE
)

head(out)

# Check for failed samples
failed <- rownames(out)[out[,2] == 0]
if(length(failed) > 0) {
  message("Samples with 0 reads after filtering: ", paste(failed, collapse=", "))
}
keep <- out[,2] > 0
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]
sample.names <- sample.names[keep]
# -------------------------------
# Learn error rates
# -------------------------------
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Optional: plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# -------------------------------
# Denoise sequences
# -------------------------------
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
# -------------------------------
# Merge paired reads
# -------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])
# -------------------------------
# Make ASV table
# -------------------------------
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
table(nchar(getSequences(seqtab)))
#Remove the Sequences that are much longer or shorter than expected 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
dim(seqtab2)
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab.nochim0 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)  # fraction of reads retained

# Check

# -------------------------------
# Assign taxonomy (genus+species)
# -------------------------------
taxa <- assignTaxonomy(
  seqtab.nochim,
  "silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  multithread = TRUE
)

# Check taxonomy
head(taxa)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa)
)

summary(sample_sums(ps))

save.image(file = "dada_out_default.RData")


saveRDS(ps, file="phyloseq_default.rds")

# -------------------------------
# Create Phyloseq object
# -------------------------------

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa)
)
# Optional: save Phyloseq object
saveRDS(ps, file="~/RNA16S/phyloseq_otu_240.rds")

# -------------------------------
# Ready for downstream analysis
# -------------------------------
ps  # inspect object

load("dada_out_default.RData")
library(DECIPHER); packageVersion("DECIPHER")

#http://DECIPHER.codes/Downloads.html
load("silva/SILVA_SSU_r138.2_v2.RData")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs

ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

ps_DECIPHER <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxid)
)
# Optional: save Phyloseq object
saveRDS(ps_DECIPHER, file="~/RNA16S/phyloseq_DECIPHER.rds")
