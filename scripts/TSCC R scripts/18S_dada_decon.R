############ 16S part 4
# decontam

# remember to run through the 4.1.2 version of R instead of default

#install.packages("ggplot2")
library("ggplot2")

#install.packages("decontam")
library("decontam")
packageVersion("decontam")

library("BiocManager")

#BiocManager::install("phyloseq")
library("phyloseq")
packageVersion("phyloseq")

#BiocManager::install("Biostrings", force=TRUE)
library("Biostrings")
packageVersion("Biostrings")

### the below help with formatting, editing, exporting phyloseq
#install.packages("remotes")
#remotes::install_github("adrientaudiere/MiscMetabar")
#remotes::install_github("peterolah001/BiMiCo")
library("MiscMetabar")
library("BiMiCo")


start_time <- Sys.time() # track timing

setwd('/projects/ps-shurinlab/users/cbwall/Yos_water_18S')

## required load in
###################################################################################################
# read in formatted metaData
run.metaD<-read.csv("output/run.metaD.edit.csv")
run.metaD<-run.metaD[-1] # remove the "X" column

# format metadata (have to redo this if loading back in)
make.fac<-c("submission_sample_ID", "sample_year_extract_ID", "year", "location", "site", "sample_type", "sample_number", "sample_control")
run.metaD[make.fac]<-lapply(run.metaD[make.fac], factor) # make all these factors

# Sample Names
sampleNames<-run.metaD$sampleNames

taxa <- readRDS("output/taxa.rds")
seqtab.nochim <- readRDS("output/seqtab.nochim.rds")

###################################################################################################
### end required load in 

# load the three files generated in tables.R
counts <- read.table(file = 'output/ASV_counts.tsv', sep = '\t', header = TRUE, row.names = 1)
tax_table <- read.table(file = 'output/ASV_taxonomy.tsv', sep = '\t', header = TRUE, row.names = 1)
asv_fasta <- readRDS("output/ASV_fasta.rds")

# sanity check to make sure the tables loaded right
# counts colnames should be a list of the samples, tax_table colnames should be kingdom phylum class etc

print("Counts columns")
colnames(counts)

# need this to remove "X" if reading back in
colnames(counts) <- gsub(x = colnames(counts), pattern = "X", replacement = "")
colnames(counts)

# now the "X..." is removed
print("Tax columns")
colnames(tax_table)


## sample data
# metadata is run.metaD
all(rownames(seqtab.nochim) %in% run.metaD$sampleNames)

rownames(run.metaD) <- run.metaD$sampleNames

ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
              sample_data(run.metaD), 
              tax_table(taxa))


# make a string of DNA names and add to phyloseq
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

###########################
# Show available ranks in the dataset
rank_names(ps)
# 78 bacteria, 8694 Eukaryota, 81 NA

table(tax_table(ps)[, "Kingdom"], exclude = NULL)
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# remove NAs in taxonomic table
ps <- subset_taxa(ps, !is.na(Kingdom) & !Kingdom %in% c("", "Bacteria"))
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Chloroplast"))

# just in case! 
ps<- subset_taxa(ps, Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class)) 

# check and see that the Bacteria are gone
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
# only 5953 eukaryota now

#191 samples

#######
#remove Colombia samples
ps<- subset_samples(ps, !(location %in% "colombia")) # removes colombia "real" samples

# controls from Colombia to remove 
remove.Col.cont<- c("colombia_pcr_NEG_2", "colombia_pcr_POS_16S", "colombia_pcr_NEG_1")

ps.Yos<- subset_samples(ps, !(sample_year_extract_ID %in% remove.Col.cont)) 


#######
# revise run.metaD to remove those samples too
run.metaD<-run.metaD[(!run.metaD$location=="colombia"),]
run.metaD<-run.metaD[!(run.metaD$sample_year_extract_ID=="colombia_pcr_NEG_3" |
                         run.metaD$sample_year_extract_ID=="colombia_pcr_POS_16S" |
                         run.metaD$sample_year_extract_ID=="colombia_pcr_NEG_4"),]
Yos.metaD.18S<-run.metaD
write.csv(Yos.metaD.18S, "output/Yos.metaD.18S.csv")

# re-examine table, NAs gone, and all Colombia gone
table(tax_table(ps.Yos)[, "Phylum"], exclude = NULL)
sample_data(ps.Yos)
######

# now PS object has all Sierra + controls

###########################
## ID contaminants
sample_data(ps.Yos)$is.neg <- sample_data(ps.Yos)$sample_control == "neg.controls" 
contamdf.prev <- isContaminant(ps.Yos, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant) # which are contaminants? 
head(which(contamdf.prev$contaminant))
# 5947 not, 6 contams.

###########################
### prune controls and low reads
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps.Yos)
ps.noncontam

# 5947 in 108 samples

#make sure those negative and pos controls are out! 
ps.noncontam.controls.out<- subset_samples(ps.noncontam, 
                                           !(sample_control %in% "neg.controls"))
ps.noncontam.controls.out<- subset_samples(ps.noncontam.controls.out, 
                                           !(sample_control %in% "pos.controls"))


###########################
# prune those not in at least 1 sample
ps.prune <- prune_taxa(taxa_sums(ps.noncontam.controls.out) > 1, ps.noncontam.controls.out) 


# remove samples with < 100 reads
ps.prune <- prune_samples(sample_sums(ps.prune) > 100, ps.prune)
ps.prune
# 4514 taxa in 96 samples (with colombia and controls removed)


###########################
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps.prune),
               MARGIN = ifelse(taxa_are_rows(ps.prune), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps.prune),
                    tax_table(ps.prune))

prev.ASVs<-plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
colnames(prev.ASVs)<- c("Phylum", "mean.prevalence", "sum.prevalence")
write.csv(prev.ASVs, "output/prev.ASVs.csv")


rich<-estimate_richness(ps.prune, split = TRUE, measures = NULL)
write.csv(rich, "output/richness.table.csv")

########### let's inspect
df.ps.prune <- as.data.frame(sample_data(ps.prune))
df.ps.prune$LibrarySize <- sample_sums(ps.prune) # this is the # of reads
df.ps.prune <- df.ps.prune[order(df.ps.prune$LibrarySize),]
df.ps.prune$Index <- seq(nrow(df.ps.prune))

# library size / number of reads
df.ps.prune$LibrarySize
#6600 - 49000 size

########### 
saveRDS(ps.prune, "output/ps.prune.18S.RDS") # save phyloseq
ps.prune.18S<-readRDS("output/ps.prune.18S.RDS") # bring it back

# save ps as 4 csvs
# One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv) 
write_phyloseq(ps.prune.18S, path = "output")


##### END!!!!
