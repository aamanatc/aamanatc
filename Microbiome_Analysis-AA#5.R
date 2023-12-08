library(dada2)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(microbiome)
library(vegan)
library(picante)
library(DT)
library(ggpubr)
# Obtain the data:
path <- "~/Downloads/MiSeq_SOP" # CHANGE ~/Downloads download to the directory including fastq file 
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plot F and R
plotQualityProfile(fnFs[1:2])
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# 1- Quality Control
#1-Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#evaluating the error rate:
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plot the error rate:
plotErrors(errF, nominalQ=TRUE)
#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#Inspecting the returned dada-class object:
dadaFs[[1]]
# merge the reads:
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
# make the sequence table:
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# removing chimers:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #0.9640374, here chimeras make up about 21% of the merged sequence variants, but when we account for the abundances of those variants we see they account for only about 4% of the merged sequence reads.
#Track reads through the pipeline
  
  getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write_csv(as.data.frame(track), "Track_raw_reads.csv")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2- Assign Taxonomy:
# the taxanomy take time to run 
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
# print the taxonomic assignments:
  taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#Import into phyloseq:
theme_set(theme_bw())
#construct a simple sample data.frame from the information encoded in the filenames. Usually this step would instead involve reading the sample data in from a file.
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$State <- "Controls"
samdf$State[samdf$Day>100] <- "Patients"
rownames(samdf) <- samples.out
samdf$Season <- samdf$When
#We now construct a phyloseq object directly from the dada2 outputs.

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
# clean the taxa name and use ASVs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

######################################################3
# Taxonomy analysis:
# add a tree to the phyloseq:
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
#merge the tree with the phyloseq
physeq = merge_phyloseq(ps, random_tree)
#plot the tree:
plot_tree(physeq, color="State", ladderize="left") + coord_polar(theta="y")
plot_tree(physeq, color="State", label.tips="Genus")
#Relative abundance analysis:
top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~State, scales="free_x")
######################################################################
# Diversity Analysis
# 1- alpha-diversity:
# Does the diversity analysis differ if we use different normalization methods?
# Alpha diversity analysis using rarification method
set.seed(129) # keep result reproductive
ps.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=1000, replace=F)  # two samples F3D142 and F3D143 were removed due to low read

# alpha diversity
rich = estimate_richness(ps.rarefied, measures = "Observed")
meta <- meta(ps.rarefied)
meta$Observed <- rich$Observed
kruskal.test(Observed ~ State, data=meta) # p-value = 0.03034 controls had sig more richness and observed taxa than patients
# plot alpha diversity observed index
p.observed <- boxplot_alpha(ps.rarefied, 
                           index = "Observed",
                           x_var = "State")
p.observed <- p.shannon + theme_minimal() + 
  labs(x="State", y="Observed diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.observed
# alpha without rarification
set.seed(1610)
rich = estimate_richness(physeq, measures = "Observed") # using the original reads
meta <- meta(physeq)
meta$Observed <- rich$Observed
kruskal.test(Observed ~ State, data=meta) # p-value = 0.1205 # without rarification there was no differences between the two groups. 
# plot alpha diversity observed index
p.observed <- boxplot_alpha(physeq, 
                            index = "Observed",
                            x_var = "State")
p.observed <- p.shannon + theme_minimal() + 
  labs(x="State", y="Observed diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
p.observed

# 2- Beta diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="State", title="Bray NMDS")
p1 = plot_ordination(ps.prop, ord.nmds.bray, type="taxa", color="Phylum", title="taxa")# display by taxa
print(p1)

# 3- Phylogenetics analysis
set.seed(7171)
#Using the rarifired phyloseq object : ps.rarefied
ps0.rar <- ps.rarefied
ps0.rar.asvtab <- as.data.frame(ps0.rar@otu_table)
ps0.rar.tree <- ps0.rar@phy_tree
# make sure if the tree is rooted or not 
ps0.rar@phy_tree
# it is a rooted tree
df.pd <- pd(ps0.rar.asvtab, ps0.rar.tree,include.root=T) # t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
datatable(df.pd)
# We will add the results of PD to this file and then plot.
newmeta <- meta(ps0.rar)
newmeta$Phylogenetic_Diversity <- df.pd$PD
pd.plot <- ggboxplot(newmeta,
                     x = "State",
                     y = "Phylogenetic_Diversity",
                     fill = "State",
                     palette = "jco",
                     ylab = "Phylogenetic Diversity",
                     xlab = "State",
                     legend = "right"
                     
                     
)
pd.plot <- pd.plot + rotate_x_text()
pd.plot
# test if there is phylogeny difference between taxa from Controls and patients group
kruskal.test(Phylogenetic_Diversity ~ State, data=newmeta) # sig  p-value = 0.02224 Control taxa had a higher phylogenetics alpha diversity





