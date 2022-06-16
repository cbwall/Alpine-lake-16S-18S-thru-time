# edit metadata
md<-read.csv("data/metadata/lake_eDNA.csv")
md$year <- sapply(strsplit(md$sample_year_extract_ID, "_"), `[`, 1) # extract sample names to YEARS
md$sample.ID <- sapply(strsplit(md$sample_year_extract_ID, "_"), `[`, 2) # extract sample names to YEARS

# these are extraction controls
md$year[md$year=='sierra']<- NA
md$year[md$year=='colombia']<- NA

# make a column for location
md$location <- md$sample_type
md$sample_type[md$sample_type=='blank' | md$sample_type=='pcr_pos' | md$sample_type=='pcr_neg']<- "NA"

# change sample type to reflect eDNA or controls
md$sample_type[md$sample_type=='colombia' | md$sample_type=='sierra']<- "eDNA"


run.metaD<- md %>% 
  dplyr::select(submission_sample_ID, gene, year, location, site, 
                sample_type, sample_number, sample_year_extract_ID)

# export clean metadata
write.csv(run.metaD, "data/metadata/run.metaD.csv")
write.csv(run.metaD[c(1:10),], "data/metadata/16Stest.run.metaD.csv")

# filter and trim}

miseq_path<-"data/Yos_test" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)

## Filter and Trim
### remove low quality reads, trim to consistent length

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames.p2 <- sapply(strsplit(fnFs, "_"), `[`, 2) # extract sample names
sampleNames.p3 <- sapply(strsplit(fnFs, "_"), `[`, 3) # extract the run # sample
sampleNames<-paste(sampleNames.p2,sampleNames.p3) # compile
sampleNames<-gsub(" ", "_", sampleNames) # remove space and add an underscore

# add this SampleNames to the metadata file
run.metaD<-run.metaD[c(1:10),]
run.metaD$sampleNames<-sampleNames

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]

# quality score plot for forward reads
plotQualityProfile(fnFs[c(1,10)])

# quality score plot for reverse reads
plotQualityProfile(fnRs[c(2,8)])

# for test data could use 270/250

# We define the filenames for the filtered fastq.gz files:

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# Trimming and filtering is performed on paired reads jointly, i.e. both reads must pass the filter for the pair to pass.

# trim primer for forward (52) and reverse (54)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,250), #trimLeft=c(52,54)
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

### Plot error rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

### Derep
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

dadaFs <-dada(derepFs, err=errF, multithread=2)
dadaRs <-dada(derepRs, err=errF, multithread=2)

# inspect data
dadaFs[[1]]

# get sequences
head(getSequences(dadaFs[[2]]))

# merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtab <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtab)))

# remove chimera
seqtab.nochim <-removeBimeraDenovo(seqtab, method="consensus", multithread=2, verbose=TRUE)

# Set multithread=TRUE to use all cores
sum(seqtab.nochim)/sum(seqtab) # 99% of samples kept

getN <-function(x)sum(getUniques(x))
track <-cbind(out,sapply(dadaFs, getN), 
              sapply(dadaRs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(track) <-c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
head(track)

# assign taxonomy
fastaRef <- "data/rdp_train_set_16.fa.gz"
taxTab <- assignTaxonomy(seqtab.nochim, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

## sample data
# metadata is run.metaD
all(rownames(seqtab.nochim) %in% run.metaD$sampleNames)

rownames(run.metaD) <- run.metaD$sampleNames

ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
              sample_data(run.metaD), 
              tax_table(taxTab))
              
ps


## let's inspect
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps) # this is the # of reads
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))


########### plot the inspected figures
# figure formatting conditions
Fig.formatting<-(theme_classic()) +
  theme(text=element_text(size=10),
        axis.line=element_blank(),
        legend.text.align = 0,
        legend.text=element_text(size=10),
        #legend.title = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        aspect.ratio=1, 
        axis.ticks.length=unit(0.25, "cm"),
        axis.text.y=element_text(
          margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"), colour="black", size=10), 
        axis.text.x=element_text(
          margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"), colour="black", size=8)) +
  theme(legend.key.size = unit(0.4, "cm")) +
  theme(aspect.ratio=1.3) +
  theme(panel.spacing=unit(c(0, 0, 0, 0), "cm"))

#####
plot.inspect.reads.type<-ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type)) + 
  geom_point()+
  scale_color_manual(name="Sample Type", values = c("darkgoldenrod1", "cornflowerblue")) +
  xlab("Sample Index") + ylab("Library Size") +
  ylim(0, 20000) +
  Fig.formatting


# Show available ranks in the dataset
rank_names(ps)

table(tax_table(ps)[, "Phylum"], exclude = NULL)

#remove NAs
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# re-examine table, NAs gone
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


#####

# inspect # of reads
sort(rowSums(otu_table(ps1))) #reads: 134 in PCR and 65,000 in highest
rich<-estimate_richness(ps1, split = TRUE, measures = NULL)
richness.test<-cbind(run.metaD, rich$Observed)

# plot it
plot_richness(ps1, x="species", measures=c("Observed", "Shannon"))
