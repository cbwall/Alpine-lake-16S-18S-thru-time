######### First in a series of scripts for R
# before you run make sure you close conda 
# in terminal: conda info --env 
# close open conda with 'conda deactivate'


#### load packages

#### load dada2 package
# install an older version of devtools form CRAN
dev<- "https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.4.3.tar.gz"
install.packages(dev, repos=NULL, type="source")
library("devtools")

library("dada2")
library("dplyr")
# (.packages()) to check if packages are loaded

start_time <- Sys.time() # track timing


#########

setwd('/projects/ps-shurinlab/users/cbwall/Yos_water_18S')

# edit metadata
md<-read.csv("data/lake_18S_eDNA.csv")
md$year <- sapply(strsplit(md$sample_year_extract_ID, "_"), `[`, 1) # extract sample names to YEARS
md$sample.ID <- sapply(strsplit(md$sample_year_extract_ID, "_"), `[`, 2) # extract sample names to YEARS

# these are extraction controls
md$year[md$year=='sierra']<- NA
md$year[md$year=='colombia']<- NA

# make a column for location
md$location <- md$sample_type
md$location[md$location=='blank' | 
              md$location=='pcr_pos' | 
              md$location=='pcr_neg' | 
              md$location=='ext.blank'] <- "NA"

# change sample type to reflect eDNA or controls
md$sample_type[md$sample_type=='colombia' | md$sample_type=='sierra']<- "eDNA"

# made a sample or control factor
md$sample_control<- md$sample_type
md$sample_control[md$sample_control=='eDNA']<- "samples"

# ID the - controls
md$sample_control[md$sample_control=='blank' |
                    md$sample_control=='pcr_neg' | 
                    md$sample_control=='ext.blank'] <- "neg.controls"

# ID the + controls
md$sample_control[md$sample_control=='pcr_pos'] <- "pos.controls"

#### add this SampleNames to the metadata file
S1<-sapply(strsplit(md$submission_sample_ID, "_"), `[`, 2)
S2<-sapply(strsplit(md$submission_sample_ID, "_"), `[`, 3)
S.name<-sampleNames<-paste(S1,S2)
sampleNames<-gsub(" ", "_", S.name) # remove space and add an underscore

md$sampleNames<-(sampleNames)

# rearrange
run.metaD<- md %>% 
  dplyr::select(submission_sample_ID, sample_year_extract_ID, gene, sampleNames,
                year, location, site, 
                sample_type, sample_number, sample_control)

make.fac<-c("year", "location", "site")
run.metaD[make.fac] <- lapply(run.metaD[make.fac], factor) # make all these factors


# load in the cut-adapt samples in the "trimmed" folder
# run this in terminal 'gunzip data/trimmed_sequences/*'
# now unzipped... proceed

# path to folder containing demultiplexed library sequencing files
# make sure unzipped
miseq_path<-"data/trimmed" # CHANGE to the directory containing the fastq files after unzipping.
list.files(miseq_path)


################################# Filter and Trim
### remove low quality reads, trim to consistent length

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_trimmed.fastq.gz"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_trimmed.fastq.gz"))

##### ##### ##### ##### 
##### had issue here, the way #s reading in not in order with metadata sheet...
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames.p2 <- sapply(strsplit(fnFs, "_"), `[`, 2) # extract sample names
sampleNames.p3 <- sapply(strsplit(fnFs, "_"), `[`, 3) # extract the run # sample
sampleNames<-paste(sampleNames.p2,sampleNames.p3) # compile
sampleNames<-gsub(" ", "_", sampleNames) # remove space and add an underscore

#### add this SampleNames to the metadata file
run.metaD$sampleNames<-sampleNames

write.csv(run.metaD, "output/run.metaD.edit.csv")

################################ Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
#fnFs[1:3]

########################################
# create pdf of quality profiles for forward samples
# the aggregate = TRUE gives quality plot for ALL reads

pdf("output/qual_profiles_F.pdf")
plotQualityProfile(fnFs, aggregate = TRUE)
dev.off()

# create pdf of quality profiles for reverse samples
pdf("output/qual_profiles_R.pdf")
plotQualityProfile(fnRs, aggregate = TRUE)
dev.off()

# final output is two pdfs of quality score profiles - forward and reverse
#################################

Sys.time() - start_time

