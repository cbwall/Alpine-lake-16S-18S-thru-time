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



######### ######### ######### 
######### load in phloseq objects for 16 and 18S
# ps.prune is the final phyloseq object from pipeline, can load it in here...

ps.prune.16S.og<-readRDS("output/16S_Yos_data/ps.prune.16S.RDS")
ps.prune.18S.og<-readRDS("output/18S_Yos_data/ps.prune.18S.RDS")

#######
# idiot check: make sure to remove all mitochondria and chloroplast taxonomic IDs
ps.prune.16S.og <- subset_taxa(ps.prune.16S.og, Family!= "Mitochondria" | 
                                is.na(Family) & Class!="Chloroplast" | is.na(Class))

ps.prune.18S.og <- subset_taxa(ps.prune.18S.og, Family!= "Mitochondria" | is.na(Family) 
                               & Class!="Chloroplast" | is.na(Class) &
                                 Kingdom!="Bacteria" | is.na(Kingdom))

table(tax_table(ps.prune.16S.og)[, "Kingdom"], exclude = NULL)
table(tax_table(ps.prune.18S.og)[, "Kingdom"], exclude = NULL)


########## ##### #####  Read counts and rarefaction curves, make df with a column for the read counts of each sample
##### sample read sums for 16S
sample_sum_df.16S<-as.data.frame(sample_sums(ps.prune.16S.og))
colnames(sample_sum_df.16S)<-"read.sum"
sample_sum_df.16S$sampleNames<-rownames(sample_sum_df.16S)

# bring out the metadata
metaD.16S = microbiome::meta(ps.prune.16S.og)

# are row names the same? if so, add them in
identical(rownames(sample_sum_df.16S), rownames(metaD.16S))
metaD.16S$read.sum<-sample_sum_df.16S$read.sum


##### sample read sums for 18S
sample_sum_df.18S<-as.data.frame(sample_sums(ps.prune.18S.og))
colnames(sample_sum_df.18S)<-"read.sum"
sample_sum_df.18S$sampleNames<-rownames(sample_sum_df.18S)

# bring out the metadata
metaD.18S = microbiome::meta(ps.prune.18S.og)

# are row names the same? if so, add them in
identical(rownames(sample_sum_df.18S), rownames(metaD.18S))
metaD.18S$read.sum<-sample_sum_df.18S$read.sum


# bring in the env data
MD.ENV.data.16S<-read.csv("data/metadata/metadata_eDNA_16S.csv")
MD.ENV.data.18S<-read.csv("data/metadata/metadata_eDNA_18S.csv")


########## ##### #####  16S metadata smithing
sub.ID.16S<-MD.ENV.data.16S$submission_sample_ID
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames.p2.16S <- sapply(strsplit(sub.ID.16S, "_"), `[`, 2) # extract sample names
sampleNames.p3.16S <- sapply(strsplit(sub.ID.16S, "_"), `[`, 3) # extract the run # sample
sampleNames.16S<-paste(sampleNames.p2.16S,sampleNames.p3.16S) # compile
sampleNames<-gsub(" ", "_", sampleNames.16S) # remove space and add an underscore

row.names(MD.ENV.data.16S)<-sampleNames

MD.ENV.df.16S<- MD.ENV.data.16S %>%
  dplyr::select(-sample_number, -sample_year_extract_ID, -year, -sample_type, -gene, -pH.method,
                -alkalinity.drops, -salinity, -tow.depth)

# if need to merge by a column, say if sequences not all in a single run or separated for some reason...
run.metaD.merge.16S <- merge(metaD.16S, MD.ENV.df.16S, by="submission_sample_ID", all.x=TRUE)
rownames(run.metaD.merge.16S)<-run.metaD.merge.16S$sampleNames

#remove unwanted columns (year is a duplicate and the 'site.x' is riddled with poor names)
run.metaD.merge.16S<- run.metaD.merge.16S %>%
  dplyr::select(-site.x)

# adjust name
run.metaD.merge.16S<-run.metaD.merge.16S %>% 
  dplyr::rename(
    site = site.y)

# make an above and below treeline factor
run.metaD.merge.16S$treeline<-ifelse(run.metaD.merge.16S$elevation.m > 2900, "above", "below")

# adjust factors
# format metadata (have to redo this if loading back in)
make.fac<-c("gene", "site", "fish", "treeline")
run.metaD.merge.16S[make.fac]<-lapply(run.metaD.merge.16S[make.fac], factor) # make all these factors

# reorder lake by relative elevation
run.metaD.merge.16S$site<-reorder(run.metaD.merge.16S$site, run.metaD.merge.16S$elevation.m)

# bring it back in, write export to final
sample_data(ps.prune.16S.og)<-run.metaD.merge.16S
saveRDS(ps.prune.16S.og, "output/final phyloseq/ps.prune.16S.og.RDS") # save phyloseq


########## ##### #####  18S metadata smithing
sub.ID.18S<-MD.ENV.data.18S$submission_sample_ID
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sampleNames.p2.18S <- sapply(strsplit(sub.ID.18S, "_"), `[`, 2) # extract sample names
sampleNames.p3.18S <- sapply(strsplit(sub.ID.18S, "_"), `[`, 3) # extract the run # sample
sampleNames.18S<-paste(sampleNames.p2.18S,sampleNames.p3.18S) # compile
sampleNames<-gsub(" ", "_", sampleNames.18S) # remove space and add an underscore

row.names(MD.ENV.data.18S)<-sampleNames

MD.ENV.df.18S<- MD.ENV.data.18S %>%
  dplyr::select(-sample_number, -sample_year_extract_ID, -year, -sample_type, -gene, -pH.method,
                -alkalinity.drops, -salinity, -tow.depth)

# if need to merge by a column, say if sequences not all in a single run or separated for some reason...
run.metaD.merge.18S <- merge(metaD.18S, MD.ENV.df.18S, by="submission_sample_ID", all.x=TRUE)
rownames(run.metaD.merge.18S)<-run.metaD.merge.18S$sampleNames

#remove unwanted columns ('site.x' is riddled with poor names)
run.metaD.merge.18S<- run.metaD.merge.18S %>%
  dplyr::select(-site.x)

# adjust name
run.metaD.merge.18S<-run.metaD.merge.18S %>% 
  dplyr::rename(
    site = site.y)

# make an above and below treeline factor
run.metaD.merge.18S$treeline<-ifelse(run.metaD.merge.18S$elevation.m > 2900, "above", "below")

# adjust factors
# format metadata (have to redo this if loading back in)
make.fac<-c("gene", "site", "fish", "treeline")
run.metaD.merge.18S[make.fac]<-lapply(run.metaD.merge.18S[make.fac], factor) # make all these factors

# reorder lake by relative elevation
run.metaD.merge.18S$site<-reorder(run.metaD.merge.18S$site, run.metaD.merge.18S$elevation.m)

# bring it back in, write export to final
sample_data(ps.prune.18S.og)<-run.metaD.merge.18S
saveRDS(ps.prune.18S.og, "output/final phyloseq/ps.prune.18S.og.RDS") # save phyloseq

########## ########## ########## ########## 

########## EXPLORE PLOTS

##### 16S
#### inspect plots
read.yr.16S<-ggplot(run.metaD.merge.16S, aes(x=sample_type, y=read.sum, color=year)) + geom_boxplot() +
  theme_bw() + ylim(0,60000)+
  ggtitle("16S: Reads by year")

read.site.16S<-ggplot(run.metaD.merge.16S, aes(x=site, y=read.sum, color=treeline)) + geom_boxplot() + 
  theme_bw()+ ylim(0,60000)+ theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  ggtitle("16S: Reads by site")

#############
# Histogram of sample read counts 16S
hist.depth.16S <- ggplot(run.metaD.merge.16S, aes(read.sum)) + 
  geom_histogram(color="gold4", fill= "gold2", binwidth = 3000) +
  xlab("Read counts") + ylim(0,16)+
  ggtitle("16S: Sequencing Depth") + theme_classic()

# a different view
read_depth_violin.16S <- ggplot(run.metaD.merge.16S, aes(x=1, y=read.sum)) +
  geom_violin(color="gold4", fill= "gold2", alpha=0.2) + scale_y_log10(limits = c(10e2,10e4)) +
  xlab("Samples") + ylab("Read Depth") + geom_boxplot(width=0.1)+ theme_classic() +
  ggtitle("16S: Distribution of sample sequencing depth")


####################
# Rarefaction curve 16S
count_table_filt.16S <- as.data.frame(otu_table(ps.prune.16S.og))

pdf(file="figures/16S/rarefac.16S.pdf", height=6, width=8)
rarecurve(count_table_filt.16S, col="gold4", step=50, cex=0.5, xlim=c(0,10000), ylab = "16S ASVs", label=FALSE)
#abline(v = 500, lty = "dotted", col="red", lwd=2)
dev.off()

################
## 18S 
#### inspect plots 18S
read.yr.18S<-ggplot(run.metaD.merge.18S, aes(x=sample_type, y=read.sum, color=year)) + geom_boxplot() + theme_bw() + 
  ylim(0,60000) +
  ggtitle("18S: Reads by year")

read.site.18S<-ggplot(run.metaD.merge.18S, aes(x=site, y=read.sum, color=treeline)) + geom_boxplot() + theme_bw() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylim(0,60000) +
  ggtitle("18S: Reads by site")

 #############

# Histogram of sample read counts 18S
hist.depth.18S <- ggplot(run.metaD.merge.18S, aes(read.sum)) + 
  geom_histogram(color="coral4", fill= "coral", binwidth = 3000) + 
  xlab("Read counts") + ylim(0,16)+
  ggtitle("18S: Sequencing Depth") + theme_classic()

# a different view
read_depth_violin.18S <- ggplot(run.metaD.merge.18S, aes(x=1, y=read.sum)) +
  geom_violin(color="coral", fill= "coral3", alpha=0.2) + scale_y_log10(limits = c(10e2,10e4)) +
  xlab("Samples") + ylab("Read Depth") + geom_boxplot(width=0.1) + theme_classic() +
  ggtitle("18S: Distribution of sample sequencing depth")


####################
# Rarefaction curve
count_table_filt.18S <- as.data.frame(otu_table(ps.prune.18S.og))

pdf(file="figures/18S/rarefac.18S.pdf", height=6, width=8)
rarecurve(count_table_filt.18S, col="coral4", step=50, cex=0.5, xlim=c(0,10000), ylab = "18S ASVs", label=FALSE)
#abline(v = 500, lty = "dotted", col="red", lwd=2)
dev.off()


########## combine and export plots ################
inspect.raw.16S.18S<- plot_grid(read.yr.16S, read.site.16S, read_depth_violin.16S, hist.depth.16S,
                                read.yr.18S, read.site.18S, read_depth_violin.18S, hist.depth.18S,
                                ncol=4, rel_widths = c(5,12,5,5))
inspect.raw.16S.18S
dev.copy(pdf, "figures/inspect.raw.16S.18S.pdf", height=8, width=20)
dev.off()




################### 16S ################### ################### ################### 

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
dev.copy(pdf, "figures/16S/16S.NMDS.tests.pdf", height=6, width=20)
dev.off()
#######

# test effects by treeline
Hell.NMDS.treeline<- plot_ordination(
  physeq=PS.rar.hell, 
  ordination = ord.hell) +
  geom_point(aes(color = treeline), size = 2) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=treeline)) +
  ggtitle("ALL taxa - hellinger-6800") +
  theme_classic() 

Hell.NMDS.treeline
dev.copy(pdf, "figures/16S/16S.NMDS.tree.rar.hell.pdf", height=6, width=8)
dev.off()
