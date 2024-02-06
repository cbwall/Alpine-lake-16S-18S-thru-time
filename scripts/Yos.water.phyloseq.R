########### ########### investigate the final 16S and 18S library ########### ########### 

# load packages
if (!require("pacman")) install.packages("pacman") # for rapid install if not in library

devtools::install_github("benjjneb/dada2", ref="v1.20") # update to most recent dada2
devtools::install_github("zdk123/SpiecEasi")
remotes::install_github("microbiome/microbiome")


# use pacman to load CRAN packages missing
pacman::p_load('knitr', 'microbiome', 'phyloseq', 'tidyr', 'tidyverse', 'knitr', 'magrittr', 'effects', 'devtools',
               'stringi', 'dplyr', "ggplot2", "gridExtra", "dada2", "phyloseq", "vegan", "cowplot", 'doBy', 'ecodist',
               'glue', 'geosphere', 'data.table', 'patchwork', 'car', 'ggcorrplot', 'FactoMineR', 'devtools', 'reshape',
               'lattice',  'plyr', 'magrittr', 'factoextra', 'multcompView', 'decontam', 'factoextra', 'car', 'mia', 
               'ade4', 'fossil', 'picante', 'reshape', 'readr', 'corrr', 'scipplot', 'Hmisc', 
               "decontam","BiocManager", 'ggpubr', 'ggmap', "ggordiplots", "fossil", "SpiecEasi", "igraph", "huge", "MASS")




######### 16S
# ps.prune is the final phyloseq object from pipeline, can load it in here...
### or load in the raw data and re-assemble

ps.prune.16S.og<-readRDS("output/16S_Yos_data/ps.prune.16S.RDS")

#######
# idiot check: make sure to remove all mitochondria and chloroplast taxonomic IDs
ps.prune.16S.og <- subset_taxa(ps.prune.16S.og, Family!= "Mitochondria" | 
                                is.na(Family) & Class!="Chloroplast" | is.na(Class))

##### Read counts and rarefaction curves
# Make a data frame with a column for the read counts of each sample
sample_sum_df<-as.data.frame(sample_sums(ps.prune.16S.og))
colnames(sample_sum_df)<-"read.sum"
sample_sum_df$sampleNames<-rownames(sample_sum_df)

# bring out the metadata
metaD = microbiome::meta(ps.prune.16S.og)

# are row names the same? if so, add them in
identical(rownames(sample_sum_df), rownames(metaD))
metaD$read.sum<-sample_sum_df$read.sum

# bring in the env data
MD.ENV.data<-read.csv("data/metadata/metadata_eDNA_16S.csv")

sub.ID<-MD.ENV.data$submission_sample_ID
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames.p2 <- sapply(strsplit(sub.ID, "_"), `[`, 2) # extract sample names
sampleNames.p3 <- sapply(strsplit(sub.ID, "_"), `[`, 3) # extract the run # sample
sampleNames<-paste(sampleNames.p2,sampleNames.p3) # compile
sampleNames<-gsub(" ", "_", sampleNames) # remove space and add an underscore

row.names(MD.ENV.data)<-sampleNames

MD.ENV.df<- MD.ENV.data %>%
  dplyr::select(-sample_number, -sample_year_extract_ID, -year, -sample_type, -gene, -pH.method,
                -alkalinity.drops, -salinity, -tow.depth)

# if need to merge by a column, say if sequences not all in a single run or separated for some reason...
run.metaD.merge <- merge(metaD, MD.ENV.df, by="submission_sample_ID", all.x=TRUE)
rownames(run.metaD.merge)<-run.metaD.merge$sampleNames

#remove unwanted columns (year is a duplicate and the 'site.x' is riddled with poor names)
run.metaD.merge<- run.metaD.merge %>%
  dplyr::select(-site.x)

# make an above and below treeline factor
run.metaD.merge$treeline<-ifelse(run.metaD.merge$elevation.m > 2900, "above", "below")

# bring it back in
sample_data(ps.prune.16S.og)<-run.metaD.merge

saveRDS(ps.prune.16S.og, "output/final phyloseq/ps.prune.16S.og.16S.RDS") # save phyloseq


########## EXPLORE PLOTS
#### inspect plots
pdf(file="figures/16S/reads.year.pdf", height=4, width=7)
ggplot(metaD, aes(x=sample_type, y=read.sum, color=year)) + geom_boxplot() + theme_bw()
dev.off() 

pdf(file="figures/16S/read.site.pdf", height=4, width=12)
ggplot(metaD, aes(x=site, y=read.sum, color=site)) + geom_boxplot() + theme_bw()
dev.off() 

#############

# Histogram of sample read counts
hist.depth <- ggplot(metaD, aes(read.sum)) + 
  geom_histogram(color="gold4", fill= "gold2", binwidth = 3000) + 
  xlab("Read counts") +
  ggtitle("Sequencing Depth") + theme_classic()

hist.depth
dev.copy(pdf, "figures/16S/hist.depth.pdf", height=7, width=9)
dev.off() 

# a different view
read_depth_violin <- ggplot(metaD, aes(x=1, y=read.sum)) +
  geom_violin() + ggtitle("Distribution of sample sequencing depth") + scale_y_log10() +
  xlab("Samples") + ylab("Read Depth") + geom_boxplot(width=0.1)

read_depth_violin
dev.copy(pdf, "figures/16S/read_depth_violin.pdf", height=4, width=5)
dev.off()

####################
# Rarefaction curve
count_table_filt <- as.data.frame(otu_table(ps.prune.16S.og))
rarecurve((count_table_filt), step=50, cex=0.5, xlim=c(0,10000), ylab = "ASVs", label=FALSE)
#abline(v = 500, lty = "dotted", col="red", lwd=2)

dev.copy(pdf, "figures/16S/rare.raw.pdf", height=6, width=8)
dev.off() 

################
###################################

###### subset the PS objects for needs...
# final and NOT rarefied data: 
ps.prune.16S.og

# remove samples with < 5000 reads
ps.prune.16S <- prune_samples(sample_sums(ps.prune.16S.og) > 5000, ps.prune.16S.og)
ps.prune.16S # = 88 samples

# rarified to 6800 reads
PS.rar = rarefy_even_depth(ps.prune.16S, rngseed=111, sample.size=0.9*min(sample_sums(ps.prune.16S)), replace=F) 
table(sample_data(PS.rar)$sample_type) 
# 88 samples, 2440 taxa
# 449 OTUs were removed

######################
# rarify, and hellinger
PS.rar.hell = transform_sample_counts(PS.rar, function(x) x^0.5)

###### ###### ###### 
# Compare full data vs. rarefied data beta diversity
####### full data
ord.full <- ordinate(ps.prune.16S.og, method = 'PCoA', distance = 'bray')
plot_ordination(ps.prune.16S.og, ord.full, color = "year")

# what samples have low reads
plot_ordination(ps.prune.16S.og, ord.full, shape = "year") +
  geom_point(aes(color = read.sum > 6800))

Rare.culling<-plot_ordination(ps.prune.16S.og, ord.full, color = "year") +
  geom_point(aes(shape = read.sum > 6800, size= read.sum > 6800)) +
  scale_size_manual(values= c(6,2)) +
  scale_shape_manual(values=c(1, 16))+
  ggtitle("Samples culled by rarefaction") + theme_classic()

# plot the PCoA for non-rarified
Non.rar.NMDS<-plot_ordination(
  physeq = ps.prune.16S.og,                                                   
  ordination = ord.full) +                                                
  geom_point(aes(color = year), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=year)) +
  ggtitle("ALL taxa-all data") +
  theme_classic() 

###### if rarefied
ord.rar <- ordinate(PS.rar, method = 'PCoA', distance = 'bray')
plot_ordination(PS.rar, ord.rar, color = "sample_type")

Rar.NMDS<- plot_ordination(
  physeq = PS.rar,                                                   
  ordination = ord.rar) +                                                
  geom_point(aes(color = year), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=year)) +
  ggtitle("ALL taxa- 6800reads") +
  theme_classic() 

#### rarefied, but hellinger transformed
ord.hell <- ordinate(PS.rar.hell, method = 'PCoA', distance = 'bray')

Hell.NMDS<- plot_ordination(
  physeq=PS.rar.hell, 
  ordination = ord.hell) +
  geom_point(aes(color = year), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=year)) +
  ggtitle("ALL taxa - hellinger-6800") +
  theme_classic() 


NMDS.tests<-plot_grid(Rare.culling + guides(color="none"), 
                      Non.rar.NMDS +  guides(color="none"), 
                      Rar.NMDS +  guides(color="none"), 
                      Hell.NMDS, ncol=4, rel_widths=c(4,4,4,5.5))
NMDS.tests
dev.copy(pdf, "figures/16S/NMDS.tests.pdf", height=6, width=20)
dev.off()
#######


Hell.NMDS.treeline<- plot_ordination(
  physeq=PS.rar.hell, 
  ordination = ord.hell) +
  geom_point(aes(color = treeline), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=treeline)) +
  ggtitle("ALL taxa - hellinger-6800") +
  theme_classic() 

