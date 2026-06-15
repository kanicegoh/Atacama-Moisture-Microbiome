##################################
##### 1. Preparation of Data #####
##################################
library(readxl)

setwd("D:/atacama_project")

#1.1: Renaming files 
D_meta <- read_excel("SummaryOfPapers.xlsx", sheet = "D_rename")
for (i in 1:nrow(D_meta)){
  DRR <- D_meta$DRR[i]
  used_names <- D_meta$used_names[i]
  
  file.rename(
    from = file.path("raw_reads", paste0(DRR, "_1.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_1.fastq"))
  )
  
  file.rename(
    from = file.path("raw_reads", paste0(DRR, "_2.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_2.fastq"))
  )
}

S_meta <- read_excel("SummaryOfPapers.xlsx", sheet = "S_rename")
for (i in 1:nrow(S_meta)){
  ERR <- S_meta$ERR[i]
  used_names <- S_meta$used_names[i]
  
  file.rename(
    from = file.path("raw_reads", paste0(ERR, "_1.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_1.fastq"))
  )
  
  file.rename(
    from = file.path("raw_reads", paste0(ERR, "_2.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_2.fastq"))
  )
}

F_meta <- read_excel("SummaryOfPapers.xlsx", sheet = "F_rename")
for (i in 1:nrow(F_meta)){
  SRR <- F_meta$SRR[i]
  used_names <- F_meta$used_names[i]
  
  file.rename(
    from = file.path("raw_reads", paste0(SRR, "_1.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_1.fastq"))
  )
  
  file.rename(
    from = file.path("raw_reads", paste0(SRR, "_2.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_2.fastq"))
  )
}

N_meta <- read_excel("SummaryOfPapers.xlsx", sheet = "N_rename")
for (i in 1:nrow(N_meta)){
  ERR <- N_meta$ERR[i]
  used_names <- N_meta$used_names[i]
  
  file.rename(
    from = file.path("raw_reads", paste0(ERR, ".fastq")), 
    to = file.path("raw_reads", paste0(used_names, ".fastq"))
  )
  
  file.rename(
    from = file.path("raw_reads", paste0(ERR, ".fastq.gz")), 
    to = file.path("raw_reads", paste0(used_names, ".fastq"))
  )
}

W_meta <- read_excel("SummaryOfPapers.xlsx", sheet = "W_rename")
for (i in 1:nrow(W_meta)){
  ERR <- W_meta$ERR[i]
  used_names <- W_meta$used_names[i]
  
  file.rename(
    from = file.path("raw_reads", paste0(ERR, "_1.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_1.fastq"))
  )
  
  file.rename(
    from = file.path("raw_reads", paste0(ERR, "_2.fastq")), 
    to = file.path("raw_reads", paste0(used_names, "_2.fastq"))
  )
}

#1.2: Set-up File Paths - to raw seq
datasets.all <- c("D", "F", "S", "W", "N")
datasets.pair <- c("D", "F", "S", "W")
datasets.single <- c("N")
raw.pair.fp <- setNames(vector("list", length(datasets.pair)), datasets.pair)
raw.single.fp <- setNames(vector("list", length(datasets.single)), datasets.single)

for (ds in datasets.pair){
  num <- match(ds, datasets.pair)
  raw.pair.fp[[ds]] <- paste0("1.", num, "_", ds, "_raw_reads")
}

for (ds in datasets.single){
  num <- match(ds, datasets.single) + length(datasets.pair)
  raw.single.fp[[ds]] <- paste0("1.", num, "_", ds, "_raw_reads")
}

#1.3: Sorted Sample Names
Fs.pair.fp <- setNames(vector("list", length(datasets.pair)), datasets.pair)
Rs.pair.fp <- setNames(vector("list", length(datasets.pair)), datasets.pair)
single.fp <- setNames(vector("list", length(datasets.single)), datasets.single)

sample.names <- list()

for(ds in datasets.pair){
  Fs.pair.fp[[ds]] <- sort(list.files(raw.pair.fp[[ds]], pattern="_1.fastq", full.names = TRUE))
  Rs.pair.fp[[ds]] <- sort(list.files(raw.pair.fp[[ds]], pattern="_2.fastq", full.names = TRUE))

  sample.names[[ds]] <- sub("_1.fastq", "", basename(Fs.pair.fp[[ds]]))
}

for(ds in datasets.single){
  single.fp[[ds]] <- sort(list.files(raw.single.fp[[ds]], pattern=".fastq", full.names = TRUE))
  
  sample.names[[ds]] <- sub(".fastq", "", basename(single.fp[[ds]]))
}

#1.4: set up path to all the other folders
processed.fp <- list()
preprocess.fp <- list()
filtN.fp <- list()
trimmed.fp <- list()
filter.fp <- list()
table.fp <- list()

for(ds in datasets.all){
  
  num <- match(ds, datasets.all)
  
  # Main processed folder
  processed.fp[[ds]] <- file.path("D:/atacama_project", paste0("2.", num, "_", ds, "_processed"))
  
  # Subfolders
  preprocess.fp[[ds]] <- file.path(processed.fp[[ds]], "01_preprocess")
  filtN.fp[[ds]] <- file.path(preprocess.fp[[ds]], "filtN")
  trimmed.fp[[ds]] <- file.path(preprocess.fp[[ds]], "trimmed.cutadapt")
  filter.fp[[ds]] <- file.path(processed.fp[[ds]], "02_filter")
  table.fp[[ds]] <- file.path(processed.fp[[ds]], "03_tabletax")
  
  # Create folders
  dir.create(preprocess.fp[[ds]], recursive = TRUE, showWarnings = FALSE)
  dir.create(filtN.fp[[ds]], recursive = TRUE, showWarnings = FALSE)
  dir.create(trimmed.fp[[ds]], recursive = TRUE, showWarnings = FALSE)
  dir.create(filter.fp[[ds]], recursive = TRUE, showWarnings = FALSE)
  dir.create(table.fp[[ds]], recursive = TRUE, showWarnings = FALSE)
}

#1.5: Paths for Output 
Fs.filtN <- list()
Rs.filtN <- list()
single.filtN <- list()

for(ds in datasets.pair){
  Fs.filtN[[ds]] <- file.path(preprocess.fp[[ds]], "filtN", basename(Fs.pair.fp[[ds]]))
  Rs.filtN[[ds]] <- file.path(preprocess.fp[[ds]], "filtN", basename(Rs.pair.fp[[ds]]))
}

for(ds in datasets.single){
  single.filtN[[ds]] <- file.path(preprocess.fp[[ds]], "filtN", basename(single.fp[[ds]]))
}


###################################
##### 2. QC & Primer Trimming #####
###################################
library(dada2)
library(Biostrings)
library(ShortRead)

#2.1: Remove all N
pre_out.pair <- list()
pre_out.single <- list()

for (ds in datasets.pair){
  pre_out.pair[[ds]] <- filterAndTrim(Fs.pair.fp[[ds]], Fs.filtN[[ds]], 
                                 Rs.pair.fp[[ds]], Rs.filtN[[ds]],
                                 maxN = 0, 
                                 multithread=6, 
                                 verbose=TRUE)
}

for (ds in datasets.single){
  pre_out.single[[ds]] <- filterAndTrim(single.fp[[ds]], single.filtN[[ds]], 
                                 maxN = 0, 
                                 multithread=6,
                                 verbose=TRUE)
}

saveRDS(pre_out.pair, "R_objects/pre_out.pair.rds")
saveRDS(pre_out.pair, "R_objects/pre_out.single.rds")

#2.2 Remove Primers from Reads
#Run cutadapt
cutadapt <- "C:/Users/kanic/AppData/Local/Programs/Python/PYTHON~1/Scripts/cutadapt.exe"

#Primers Used:
FWD <- setNames(vector("list", length(datasets.all)), datasets.all)
REV <- setNames(vector("list", length(datasets.all)), datasets.all)

FWD[["D"]] <- "GTGYCAGCMGCCGCGGTAA"
REV[["D"]] <- "GGACTACNVGGGTWTCTAAT"
FWD[["F"]] <- "CCTACGGGNGGCWGCAG"
REV[["F"]] <- "GACTACHVGGGTATCTAATCC"
FWD[["S"]] <- "GTGCCAGCMGCCGCGGTAA"
REV[["S"]] <- "GGACTACHVGGGTWTCTAAT"
FWD[["W"]] <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
REV[["W"]] <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
FWD[["N"]] <- "GTGCCAGCMGCCGCGGTAA"
REV[["N"]] <- "GGACTACHVGGGTWTCTAAT"

#create all possible combinations of primers 
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), 
               Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) 
}

FWD.orients <- list()
REV.orients <- list()
for (ds in datasets.all){
  FWD.orients[[ds]] <- allOrients(FWD[[ds]])
  REV.orients[[ds]] <- allOrients(REV[[ds]])
}

FWD2 <- list()
REV2 <- list()
for (ds in datasets.all){
  FWD2[[ds]] <- FWD.orients[[ds]][["RevComp"]]
  REV2[[ds]] <- REV.orients[[ds]][["RevComp"]]
}

FWD.orients.2 <- list()
REV.orients.2 <- list()
for (ds in datasets.all){
  FWD.orients.2[[ds]] <- allOrients(FWD2[[ds]])
  REV.orients.2[[ds]] <- allOrients(REV2[[ds]])
}

#Count the number of times each primer appear
primerHits <- function(primer, fn) { 
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE) 
  return(sum(nhits > 0)) 
} 

primer.pair.results <- list()
primer.pair.results.2 <- list()
primer.single.results <- list()
primer.single.results.2 <- list()
for(ds in datasets.pair){
  primer.pair.results[[ds]] <- 
    rbind(FWD.ForwardReads = sapply(FWD.orients[[ds]], primerHits, fn = Fs.filtN[[ds]]),
          FWD.ReverseReads = sapply(FWD.orients[[ds]], primerHits, fn = Rs.filtN[[ds]]),
          REV.ForwardReads = sapply(REV.orients[[ds]], primerHits, fn = Fs.filtN[[ds]]),
          REV.ReverseReads = sapply(REV.orients[[ds]], primerHits, fn = Rs.filtN[[ds]]))
  primer.pair.results.2[[ds]] <- 
    rbind(FWD.ForwardReads = sapply(FWD.orients.2[[ds]], primerHits, fn = Fs.filtN[[ds]]),
          FWD.ReverseReads = sapply(FWD.orients.2[[ds]], primerHits, fn = Rs.filtN[[ds]]),
          REV.ForwardReads = sapply(REV.orients.2[[ds]], primerHits, fn = Fs.filtN[[ds]]),
          REV.ReverseReads = sapply(REV.orients.2[[ds]], primerHits, fn = Rs.filtN[[ds]]))
}
for(ds in datasets.single){
  primer.single.results[[ds]] <- 
    rbind(FWD.Reads = sapply(FWD.orients[[ds]], primerHits, fn = single.filtN[[ds]]),
          REV.Reads = sapply(REV.orients[[ds]], primerHits, fn = single.filtN[[ds]]))
  primer.single.results.2[[ds]] <- 
    rbind(FWD.Reads = sapply(FWD.orients.2[[ds]], primerHits, fn = single.filtN[[ds]]),
          REV.Reads = sapply(REV.orients.2[[ds]], primerHits, fn = single.filtN[[ds]]))
}

print(primer.pair.results)
print(primer.pair.results.2)
print(primer.single.results)
print(primer.single.results.2)

saveRDS(primer.pair.results, "PrimerHits/primer.pair.results.rds")
saveRDS(primer.pair.results.2, "PrimerHits/primer.pair.results.2.rds")
saveRDS(primer.single.results, "PrimerHits/primer.single.results.rds")
saveRDS(primer.single.results.2, "PrimerHits/primer.single.results.2.rds")

#create flags to remove all primers 
#for paired sequences
FWD.RC <- list()
REV.RC <- list()
FWD.RC.2 <- list()
REV.RC.2 <- list()

R1.flags <- list()
R2.flags <- list()
R1.flags.2 <- list()
R2.flags.2 <- list()

for(ds in datasets.pair){
  FWD.RC[[ds]] <- dada2:::rc(FWD[[ds]])
  REV.RC[[ds]] <- dada2:::rc(REV[[ds]])
  FWD.RC.2[[ds]] <- dada2:::rc(FWD2[[ds]])
  REV.RC.2[[ds]] <- dada2:::rc(REV2[[ds]])
  
  R1.flags[[ds]] <- paste("-g", FWD[[ds]], "-a", REV.RC[[ds]])
  R2.flags[[ds]] <- paste("-G", REV[[ds]], "-A", FWD.RC[[ds]])
  R1.flags.2[[ds]] <- paste("-g", FWD2[[ds]], "-a", REV.RC.2[[ds]])
  R2.flags.2[[ds]] <- paste("-G", REV2[[ds]], "-A", FWD.RC.2[[ds]])
}

#for single sequences
FWD.RC.single <- list()
REV.RC.single <- list()
FWD.RC.single.2 <- list()
REV.RC.single.2 <- list()

single.flags <- list()
single.flags.2 <- list()

for(ds in datasets.single){
  FWD.RC.single[[ds]] <- dada2:::rc(FWD[[ds]])
  REV.RC.single[[ds]] <- dada2:::rc(REV[[ds]])
  FWD.RC.single.2[[ds]] <- dada2:::rc(FWD2[[ds]])
  REV.RC.single.2[[ds]] <- dada2:::rc(REV2[[ds]])
  
  single.flags[[ds]] <- paste("-g", FWD[[ds]], "-a", REV.RC.single[[ds]])
  single.flags.2[[ds]] <- paste("-g", FWD2[[ds]], "-a", REV.RC.single.2[[ds]])
}

#prepare file path for cutadapt outputs
Fs.cut <- list()
Rs.cut <- list()
for(ds in datasets.pair){
  Fs.cut[[ds]] <- file.path(trimmed.fp[[ds]], basename(Fs.pair.fp[[ds]]))
  Rs.cut[[ds]] <- file.path(trimmed.fp[[ds]], basename(Rs.pair.fp[[ds]]))
}

single.cut <- list()
for (ds in datasets.single){
  single.cut[[ds]] <- file.path(trimmed.fp[[ds]], basename(single.fp[[ds]]))
}

#run cutadapt - paired
cutadapt.out <- list()
for (ds in datasets.pair){
  cutadapt.out[[ds]] <- vector("list", length(Fs.pair.fp[[ds]]))
  
  for (i in seq_along(Fs.pair.fp[[ds]])) { 
    cutadapt.out[[ds]][[i]] <- system2(cutadapt, args = c(R1.flags[[ds]], R2.flags[[ds]], R1.flags.2[[ds]], R2.flags.2[[ds]], #flags to use 
                               "-n", 5, #allow up to 5 primer trimming operations per read 
                               "-o", Fs.cut[[ds]][i], #output file
                               "-p", Rs.cut[[ds]][i], #output file
                               "--minimum-length", 150, #discard reads shorter than 150 after trimming
                               "-e", 0, #allow 0 mismatches when matching primers (perfect match)
                               "--report", "minimal", #short concise output, one line summary
                               "-j", 0, #detect the number of available cores 
                               Fs.filtN[[ds]][i], Rs.filtN[[ds]][i] #input files to look at
                               )
            )
  } 
}

#run cutadapt - single
cutadapt.out.single <- list()

for(ds in datasets.single){
  cutadapt.out.single[[ds]] <- vector("list", length(single.fp[[ds]]))
  
  for(i in seq_along(single.fp[[ds]])){
    cutadapt.out.single[[ds]][[i]] <- system2(cutadapt, args = c(single.flags[[ds]], single.flags.2[[ds]],
        "-n", 5,
        "-o", single.cut[[ds]][i],
        "--minimum-length", 150,
        "-e", 0,
        "--report", "minimal",
        "-j", 0,
        single.filtN[[ds]][i]))
    }
  }

#check if all primers removed
primer.pair.results.after <- list()
primer.pair.results.2.after <- list()
primer.single.results.after <- list()
primer.single.results.2.after <- list()
for(ds in datasets.pair){
  primer.pair.results.after[[ds]] <- 
    rbind(FWD.ForwardReads = sapply(FWD.orients[[ds]], primerHits, fn = Fs.cut[[ds]]),
          FWD.ReverseReads = sapply(FWD.orients[[ds]], primerHits, fn = Rs.cut[[ds]]),
          REV.ForwardReads = sapply(REV.orients[[ds]], primerHits, fn = Fs.cut[[ds]]),
          REV.ReverseReads = sapply(REV.orients[[ds]], primerHits, fn = Rs.cut[[ds]]))
  primer.pair.results.2.after[[ds]] <- 
    rbind(FWD.ForwardReads = sapply(FWD.orients.2[[ds]], primerHits, fn = Fs.cut[[ds]]),
          FWD.ReverseReads = sapply(FWD.orients.2[[ds]], primerHits, fn = Rs.cut[[ds]]),
          REV.ForwardReads = sapply(REV.orients.2[[ds]], primerHits, fn = Fs.cut[[ds]]),
          REV.ReverseReads = sapply(REV.orients.2[[ds]], primerHits, fn = Rs.cut[[ds]]))
}
for(ds in datasets.single){
  primer.single.results.after[[ds]] <- 
    rbind(FWD.Reads = sapply(FWD.orients[[ds]], primerHits, fn = single.cut[[ds]]),
          REV.Reads = sapply(REV.orients[[ds]], primerHits, fn = single.cut[[ds]]))
  primer.single.results.2.after[[ds]] <- 
    rbind(FWD.Reads = sapply(FWD.orients.2[[ds]], primerHits, fn = single.cut[[ds]]),
          REV.Reads = sapply(REV.orients.2[[ds]], primerHits, fn = single.cut[[ds]]))
}

print(primer.pair.results.after)
print(primer.pair.results.2.after)
print(primer.single.results.after)
print(primer.single.results.2.after)
saveRDS(primer.pair.results, "PrimerHits/primer.pair.results.after.rds")
saveRDS(primer.pair.results.2, "PrimerHits/primer.pair.results.2.after.rds")
saveRDS(primer.single.results, "PrimerHits/primer.single.results.after.rds")
saveRDS(primer.single.results.2, "PrimerHits/primer.single.results.2.after.rds")
## Primer hits for "D" dataset only 4, and original sequence length only 150bp -> likely already trimmed
## cutadapt causes all reads to be discarded -> step was skipped for this "D" dataset 

###########################################
##### 3. Visualizing Quality Profiles #####
###########################################
library(ggplot2)
library(tidyverse)
library(dplyr)

#3.1: create path and new files 
#paired
subF.fp <- list()
subR.fp <- list()

for (ds in datasets.pair){
  subF.fp[[ds]] <- file.path(filter.fp[[ds]], "preprocessed_F")
  subR.fp[[ds]] <- file.path(filter.fp[[ds]], "preprocessed_R")
  
  dir.create(subF.fp[[ds]], recursive=TRUE, showWarnings=FALSE)
  dir.create(subR.fp[[ds]], recursive=TRUE, showWarnings=FALSE)
}

#single
subSingle.fp <- list()

for (ds in datasets.single){
  subSingle.fp[[ds]] <- file.path(filter.fp[[ds]], "preprocessed_Single")
  
  dir.create(subSingle.fp[[ds]], recursive=TRUE, showWarnings=FALSE)
}

#3.2: copy trimmed cutadapt stuff
Fs.Q <- list()
Rs.Q <- list()
for (ds in datasets.pair){
  Fs.Q[[ds]] <- file.path(subF.fp[[ds]], basename(Fs.pair.fp[[ds]]))
  Rs.Q[[ds]] <- file.path(subR.fp[[ds]], basename(Rs.pair.fp[[ds]]))
  
  file.copy(from = Fs.cut[[ds]], to = Fs.Q[[ds]])
  file.copy(from = Rs.cut[[ds]], to = Rs.Q[[ds]])
}

single.Q <- list()
for (ds in datasets.single){
  single.Q[[ds]] <- file.path(subSingle.fp[[ds]], basename(single.fp[[ds]]))
  
  file.copy(from = single.cut[[ds]], to = single.Q[[ds]])
}

#Check 
fastqFs <- list()
fastqRs <- list()

for (ds in datasets.pair){
  fastqFs[[ds]] <- sort(list.files(subF.fp[[ds]], pattern="_1.fastq",full.names = TRUE))  
  fastqRs[[ds]] <- sort(list.files(subR.fp[[ds]], pattern="_2.fastq",full.names = TRUE))  
  if(length(fastqFs[[ds]]) != length(fastqRs[[ds]])) 
  stop("Forward and reverse files do not match.") 
}

fastqSingle <- list()

for (ds in datasets.single){
  fastqSingle[[ds]] <- sort(list.files(subSingle.fp[[ds]], pattern=".fastq", full.names = TRUE))
}

#3.3: Plot Quality Plots of Random Samples
#pick 3 from each dataset, only plot more if quality plots look problematic
rand_samples <- list()
for (ds in datasets.pair){
  rand_samples[[ds]] <- sample(size = 3, 1:length(fastqFs[[ds]]))
}
for (ds in datasets.single){
  rand_samples[[ds]] <- sample(size = 3, 1:length(fastqSingle[[ds]]))
}

for (ds in datasets.pair){
  fwd.qual.plot <- plotQualityProfile(fastqFs[[ds]][rand_samples[[ds]]]) + labs(x = "Sequence Position")  
  rev.qual.plot <- plotQualityProfile(fastqRs[[ds]][rand_samples[[ds]]]) + labs(x = "Sequence Position") 
  
  num <- match(ds, datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  pdf(file=file.path(loc, "03_tabletax", "fwd.qual.plot.pdf"), width=15, height=8.5)
  print(fwd.qual.plot)
  dev.off()  
  
  pdf(file=file.path(loc, "03_tabletax", "rev.qual.plot.pdf"), width=15, height=8.5)
  print(rev.qual.plot)  
  dev.off()  
  
  saveRDS(fwd.qual.plot, file.path(loc, "03_tabletax", "fwd.qual.plot.rds"))
  saveRDS(rev.qual.plot, file.path(loc, "03_tabletax", "rev.qual.plot.rds"))
}

for (ds in datasets.single){
  qual.plot <- plotQualityProfile(fastqSingle[[ds]][rand_samples[[ds]]]) + labs(x = "Sequence Position")  
  
  num <- match(ds, datasets.single) + length(datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  pdf(file=file.path(loc, "03_tabletax", "qual.plot.pdf"), width=15, height=8.5)
  print(qual.plot)  
  dev.off() 
  
  saveRDS(qual.plot, file.path(loc, "03_tabletax", "qual.plot.rds"))
}


##################################
##### 4. Trim & Filter Reads #####
##################################
#4.1: Create Output Location Folders
filtpathF <- list()
filtpathR <- list()
filtpathF.file <- list()
filtpathR.file <- list()
for (ds in datasets.pair){
  filtpathF[[ds]] <- file.path(filter.fp[[ds]], "filtered_F")
  filtpathR[[ds]] <- file.path(filter.fp[[ds]], "filtered_R")
  
  dir.create(filtpathF[[ds]]) 
  dir.create(filtpathR[[ds]])
  
  filtpathF.file[[ds]] <- file.path(filtpathF[[ds]], basename(Fs.pair.fp[[ds]])) 
  filtpathR.file[[ds]] <- file.path(filtpathR[[ds]], basename(Rs.pair.fp[[ds]]))
}


filtpath <- list()
filtpath.file <- list()
for (ds in datasets.single){
  filtpath[[ds]] <- file.path(filter.fp[[ds]], "filtered_single")
  
  dir.create(filtpath[[ds]])
  
  filtpath.file[[ds]] <- file.path(filtpath[[ds]], basename(single.fp[[ds]])) 
}

#4.2: Filter and Trim Low Quality Tails and Reads
#Used truncation length (tried to ensure that the quality score higher than 20)
#D: F - 150bp, R - 130bp (reads are too short)
#F: F - 275bp, R - 225bp
#S: F/R - 200bp
#W: F - 270bp, R - 220bp, MaxEE = c(10, 10) (way too little being retained if MaxEE is any less)
#N: F/R - 150bp (only 150bp sequenced)

#load the R environment before continuing 
filt_out_summary <- list()

for (ds in datasets.pair){
  filt_out <- filterAndTrim( 
    fastqFs[[ds]], filtpathF.file[[ds]], fastqRs[[ds]], filtpathR.file[[ds]],
    maxEE=c(2,2), #change as required
    truncQ=2, 
    maxN=0, 
    truncLen=c(200,200), #change as required
    rm.phix=TRUE,
    compress=FALSE, 
    verbose=TRUE,
    multithread=24) 
  
  #Get the Summary Before / After filtering 
  filt_out_summary[[ds]] <- filt_out %>% 
    data.frame() %>%  
    mutate(Samples = rownames(.), percent_kept = 100*(reads.out/reads.in)) %>%  
    select(Samples, everything()) 
}

for (ds in datasets.single){
  filt_out <- filterAndTrim(
    fastqSingle[[ds]], filtpath.file[[ds]],  
    maxEE=2, 
    truncQ=2, 
    maxN=0, 
    truncLen=150, #changed from 200bp to 150bp since reads are only 150bp long 
    rm.phix=TRUE, 
    compress=FALSE, 
    verbose=TRUE,
    multithread=24) 
  
  #Get the Summary Before / After filtering 
  filt_out_summary[[ds]] <- filt_out %>% 
    data.frame() %>%  
    mutate(Samples = rownames(.), percent_kept = 100*(reads.out/reads.in)) %>%  
    select(Samples, everything()) 
}

#4.3: Plot Quality Profiles of the Filtered fastq Files 
filtFs <- list()
filtRs <- list()
filtSingle <-list()

for (ds in datasets.pair){
  filtFs[[ds]] <- sort(list.files(filtpathF[[ds]], pattern="_1.fastq",full.names = TRUE))  
  filtRs[[ds]] <- sort(list.files(filtpathR[[ds]], pattern="_2.fastq",full.names = TRUE))
  
  rand_samples_filt <- sample(size = 3, 1:length(filtFs[[ds]]))   
  fwd.fil.plot <- plotQualityProfile(filtFs[[ds]][rand_samples_filt]) + labs(x = "Sequence Position") 
  rev.fil.plot <- plotQualityProfile(filtRs[[ds]][rand_samples_filt]) + labs(x = "Sequence Position") 
  
  num <- match(ds, datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  pdf(file=file.path(loc, "03_tabletax", "fwd.fil.plot.pdf"), width=15, height=8.5)
  print(fwd.fil.plot)
  dev.off()  
  
  pdf(file=file.path(loc, "03_tabletax", "rev.fil.plot.pdf"), width=15, height=8.5)
  print(rev.fil.plot)  
  dev.off()  
  
  saveRDS(fwd.fil.plot, file.path(loc, "03_tabletax", "fwd.fil.plot.rds"))
  saveRDS(rev.fil.plot, file.path(loc, "03_tabletax", "rev.fil.plot.rds"))
}

for (ds in datasets.single){
  filtSingle[[ds]] <- sort(list.files(filtpath[[ds]], pattern=".fastq",full.names = TRUE))  
  
  rand_samples_filt <- sample(size = 3, 1:length(filtSingle[[ds]]))   
  fil.plot <- plotQualityProfile(filtSingle[[ds]][rand_samples_filt]) + labs(x = "Sequence Position") 
  
  num <- match(ds, datasets.single) + length(datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  pdf(file=file.path(loc, "03_tabletax", "fil.plot.pdf"), width=15, height=8.5)
  print(fil.plot)
  dev.off()  
  
  saveRDS(fil.plot, file.path(loc, "03_tabletax", "fil.plot.rds"))
}

saveRDS(filt_out_summary, "R_objects/filt_out_summary.rds")
saveRDS(filtFs, "R_objects/filtFs.rds")
saveRDS(filtRs, "R_objects/filtRs.rds")
saveRDS(filtSingle, "R_objects/filtSingle.rds")

###################################
##### 5. Learning Error Rates #####
###################################
#Main Function 
loessErrfun_mod <- function(trans) { 
  qq <- as.numeric(colnames(trans)) 
  est <- matrix(0, nrow=0, ncol=length(qq)) 
  for(nti in c("A","C","G","T")) { 
    for(ntj in c("A","C","G","T")) { 
      if(nti != ntj) { 
        errs <- trans[paste0(nti,"2",ntj),] 
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),]) 
        rlogp <- log10((errs+1)/tot)
        rlogp[is.infinite(rlogp)] <- NA 
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp) 
        
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot)) 
        
        pred <- predict(mod.lo, qq) 
        maxrli <- max(which(!is.na(pred))) 
        minrli <- min(which(!is.na(pred))) 
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]] 
        pred[seq_along(pred)<minrli] <- pred[[minrli]] 
        est <- rbind(est, 10^pred) 
      } # if(nti != ntj) 
    } # for(ntj in c("A","C","G","T")) 
  } # for(nti in c("A","C","G","T")) 
  
  # HACKY 
  MAX_ERROR_RATE <- 0.25 
  MIN_ERROR_RATE <- 1e-7 
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE 
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE 
  
  #removed this portion as reads are short and 
  # estorig <- est 
  # est <- est %>% 
  #   data.frame() %>% 
  #   mutate_all(funs(case_when(. < X40 ~ X40, 
  #                             . >= X40 ~ .))) %>% as.matrix() 
  # rownames(est) <- rownames(estorig) 
  # colnames(est) <- colnames(estorig) 
  
  # Expand the err matrix with the self-transition probs 
  err <- rbind(1-colSums(est[1:3,]), est[1:3,], 
               est[4,], 1-colSums(est[4:6,]), est[5:6,], 
               est[7:8,], 1-colSums(est[7:9,]), est[9,], 
               est[10:12,], 1-colSums(est[10:12,])) 
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", 
                          c("A","C","G","T")) 
  colnames(err) <- colnames(trans) 
  
  return(err) 
} 

#5.1: Learn and Plot the Errors
errF <- list()
errR <- list()
errSingle <- list()

for (ds in datasets.pair){
  errF[[ds]] <- learnErrors( 
    filtFs[[ds]], 
    multithread = TRUE, 
    nbases = 1e10, 
    errorEstimationFunction = loessErrfun_mod, 
    verbose = TRUE 
  ) 
  errR[[ds]] <- learnErrors( 
    filtRs[[ds]], 
    multithread = TRUE, 
    nbases = 1e10, 
    errorEstimationFunction = loessErrfun_mod, 
    verbose = TRUE) 
  
  errF.plot <- plotErrors(errF[[ds]], nominalQ=TRUE) 
  errR.plot <- plotErrors(errR[[ds]], nominalQ=TRUE)
  
  num <- match(ds, datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  pdf(file=file.path(loc, "03_tabletax", "errF.plot.pdf"), width=15, height=8.5)
  print(errF.plot)
  dev.off()  
  
  pdf(file=file.path(loc, "03_tabletax", "errR.plot.pdf"), width=15, height=8.5)
  print(errR.plot)  
  dev.off()  
  
  saveRDS(errF.plot, file.path(loc, "03_tabletax", "errF.plot.rds"))
  saveRDS(errR.plot, file.path(loc, "03_tabletax", "errR.plot.rds"))
  
  saveRDS(errF[[ds]], file.path(loc, "03_tabletax", "errF.rds"))
  saveRDS(errR[[ds]], file.path(loc, "03_tabletax", "errR.rds"))
}

for (ds in datasets.single){
  errSingle[[ds]] <- learnErrors( 
    filtSingle[[ds]], 
    multithread = TRUE, 
    nbases = 1e10, 
    errorEstimationFunction = loessErrfun_mod, 
    verbose = TRUE 
  ) 
  
  errSingle.plot <- plotErrors(errSingle[[ds]], nominalQ=TRUE) 
  
  num <- match(ds, datasets.single) + length(datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  pdf(file=file.path(loc, "03_tabletax", "errSingle.plot.pdf"), width=15, height=8.5)
  print(errSingle.plot)
  dev.off()  

  saveRDS(errSingle.plot, file.path(loc, "03_tabletax", "errSingle.plot.rds"))
  
  saveRDS(errSingle[[ds]], file.path(loc, "03_tabletax", "errSingle.rds"))
}

saveRDS(errF, "R_objects/errF.rds")
saveRDS(errR, "R_objects/errR.rds")
saveRDS(errSingle, "R_objects/errSingle.rds")

############################
##### 6. Dereplication #####
############################
derepF <- list()
derepR <- list()
derep <- list()
names(derepF) <- list()
names(derepR) <- list()
names(derep) <- list()

for (ds in datasets.pair){
  #Keep only unique sequences and keep count of each sequence 
  derepF[[ds]] <- derepFastq(filtFs[[ds]], verbose = TRUE) 
  derepR[[ds]] <- derepFastq(filtRs[[ds]], verbose = TRUE) 
  
  #label each dereplicated sample object with correct sample name 
  names(derepF[[ds]]) <- sample.names[[ds]] 
  names(derepR[[ds]]) <- sample.names[[ds]]
}

for (ds in datasets.single){
  #Keep only unique sequences and keep count of each sequence 
  derep[[ds]] <- derepFastq(filtSingle[[ds]], verbose = TRUE) 
  
  #label each dereplicated sample object with correct sample name 
  names(derep[[ds]]) <- sample.names[[ds]] 
}

saveRDS(derepF, "R_objects/derepF.rds")
saveRDS(derepR, "R_objects/derepR.rds")
saveRDS(derep, "R_objects/derep_single.rds")


########################################
##### 7. Sequence Inference (ASVs) #####
########################################
dadaF <- list()
dadaR <- list()
dadaSingle <- list()

for (ds in datasets.pair){
  dadaF[[ds]] <- dada(derepF[[ds]], err = errF[[ds]], multithread = TRUE, pool = "pseudo") 
  dadaR[[ds]] <- dada(derepR[[ds]], err = errR[[ds]], multithread = TRUE, pool = "pseudo")
  
  dadaF[[ds]][[1]]
  dadaR[[ds]][[1]]
  
  num <- match(ds, datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  saveRDS(dadaF[[ds]], file.path(loc, "03_tabletax", "dadaF.rds"))
  saveRDS(dadaR[[ds]], file.path(loc, "03_tabletax", "dadaR.rds"))
}

for (ds in datasets.single){
  dadaSingle[[ds]] <- dada(derep_Single[[ds]], err = errSingle[[ds]], multithread = TRUE, pool = "pseudo") 
  
  dadaSingle[[ds]][[1]]
  
  num <- match(ds, datasets.single) + length(datasets.pair)
  loc <- paste0("2.", num, "_", ds, "_processed")
  saveRDS(dadaSingle[[ds]], file.path(loc, "03_tabletax", "dadaSingle.rds"))
}

saveRDS(dadaF, "R_objects/dadaF.rds")
saveRDS(dadaR, "R_objects/dadaR.rds")
saveRDS(dadaSingle, "R_objects/dadaSingle.rds")

#################################
##### 8. Merge Paired Reads #####
#################################
#8.1: Merge the Paired Reads + Construct ASV Table 
mergers <- list()
seqtab <- list()
dim.seqtab <- list()
seqlen <- list()

for (ds in datasets.pair){
  mergers[[ds]] <- mergePairs(dadaF[[ds]], derepF[[ds]], dadaR[[ds]], derepR[[ds]], verbose = TRUE)
  
  #construct ASV Tables 
  seqtab[[ds]] <- makeSequenceTable(mergers[[ds]])
  
  #Get the dimensions 
  dim.seqtab[[ds]] <- dim(seqtab[[ds]]) 
  
  #Distribution of Sequence Lengths
  seqlen[[ds]] <- table(nchar(getSequences(seqtab[[ds]])))
}

#8.2: Construct ASV Table for Single End Reads
for (ds in datasets.single){
  #construct ASV Tables 
  seqtab[[ds]] <- makeSequenceTable(dadaSingle[[ds]])
  
  #Get the dimensions 
  dim.seqtab[[ds]] <- dim(seqtab[[ds]]) 
  
  #Distribution of Sequence Lengths
  seqlen[[ds]] <- table(nchar(getSequences(seqtab[[ds]])))
}

saveRDS(mergers, "R_objects/mergers.rds")
saveRDS(seqtab, "R_objects/seqtab.rds")
saveRDS(dim.seqtab, "R_objects/dim.seqtab.rds")
saveRDS(seqlen, "R_objects/seqlen.rds")


##############################
##### 9. Remove Chimeras #####
##############################
seqtab.nochim <- list()
total.reads <- list()
total.nonchim.reads <- list()
nonchim <- list()
dim.seq.nonchim <- list()

for (ds in datasets.all){

  print(ds)
  
  #9.1: Remove Chimeras
  seqtab.nochim[[ds]] <- removeBimeraDenovo(seqtab[[ds]], method="pooled", multithread=FALSE, verbose = TRUE) 
  
  num <- match(ds, datasets.all)
  loc <- paste0("2.", num, "_", ds, "_processed")
  saveRDS(seqtab.nochim[[ds]], file.path(loc, "03_tabletax", "seqtab.nochim.rds"))
  
  #9.2: Get Dimensions and Percentages 
  total.reads[[ds]] <- sum(seqtab[[ds]])
  total.nonchim.reads[[ds]] <- sum(seqtab.nochim[[ds]])
  saveRDS(total.reads[[ds]], file.path(loc, "03_tabletax", "total.reads.rds"))
  saveRDS(total.nonchim.reads[[ds]], file.path(loc, "03_tabletax", "total.nonchim.reads.rds"))
  
  nonchim[[ds]] <- 100*sum(seqtab.nochim[[ds]])/sum(seqtab[[ds]]) 
  saveRDS(nonchim[[ds]], file.path(loc, "03_tabletax", "nonchim.rds"))
  
  dim.seq.nonchim[[ds]] <- dim(seqtab.nochim[[ds]]) 
  saveRDS(dim.seq.nonchim[[ds]], file.path(loc, "03_tabletax", "dim.seq.nonchim.rds"))
}

saveRDS(seqtab.nochim, "R_Objects/seqtab.nochim.rds")
saveRDS(total.reads, "R_Objects/total.reads.rds")
saveRDS(total.nonchim.reads, "R_Objects/total.nonchim.reads.rds")
saveRDS(nonchim, "R_Objects/nonchim.rds")
saveRDS(seqtab.nochim, "R_Objects/seqtab.nochim.rds")

#9.3: Tracking the Number of Reads 
getN <- function(x) sum(getUniques(x))
track <- list()

#Get filt_out
filt_out <- list()
for(ds in names(filt_out_summary)){
  filt_out[[ds]] <- filt_out_summary[[ds]][ , c("reads.in","reads.out")]
  rownames(filt_out[[ds]]) <- filt_out_summary[[ds]]$Samples
}

for (ds in datasets.pair){
  track[[ds]] <- cbind(pre_out.pair[[ds]], filt_out[[ds]], sapply(dadaF[[ds]], getN), 
                       sapply(dadaR[[ds]], getN), sapply(mergers[[ds]], getN), 
                       rowSums(seqtab.nochim[[ds]]))
  colnames(track[[ds]]) <- c("pre.DADA2.input", "pre.DADA2.filtered", 
                     "DADA2.input", "DADA2.filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim")
  rownames(track[[ds]]) <- sample.names[[ds]]
}

for (ds in datasets.single){
  track[[ds]] <- cbind(pre_out.single[[ds]], filt_out[[ds]], sapply(dadaSingle[[ds]], getN), 
                       rowSums(seqtab.nochim[[ds]]))
  colnames(track[[ds]]) <- c("pre.DADA2.input", "pre.DADA2.filtered", 
                       "DADA2.input", "DADA2.filtered", "denoisedSingle",
                       "nonchim")
  rownames(track[[ds]]) <- sample.names[[ds]]
}

saveRDS(track, "R_objects/track.rds")
head(track)
## Very poor reverse reads for W -> results in 0 ASVs in the end ##


##################################
##### 10. Assign Taxonomy ########
##################################
library(dada2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(phyloseq)
library(readxl)
library(fantaxtic)
library(microbiome)

taxa.nb <- list()
taxa.nb.1 <- list()
taxa.nb.print <- list()

for (ds in datasets.all){
  taxa.nb[[ds]] <- assignTaxonomy(seqtab.nochim[[ds]], 
                                  "taxonomy_reference/silva_nr99_v138.2_toSpecies_trainset.fa.gz", 
                                  tryRC = TRUE,
                                  multithread=TRUE, 
                                  outputBootstraps = TRUE, 
                                  verbose = TRUE) 
  
  taxa.nb.1[[ds]] <- addSpecies(taxa.nb[[ds]][[1]],
                                "taxonomy_reference/silva_v138.2_assignSpecies.fa.gz")
  
  num <- match(ds, datasets.all)
  loc <- paste0("2.", num, "_", ds, "_processed")
  saveRDS(taxa.nb.1[[ds]], file.path(loc, "03_tabletax", "taxa.nb.1.rds"))
  saveRDS(taxa.nb[[ds]][[2]], file.path(loc, "03_tabletax", "taxa.nb.2.rds"))
  
  #Checking the taxonomy
  taxa.nb.print[[ds]] <- taxa.nb.1[[ds]] 
  rownames(taxa.nb.print[[ds]]) <- NULL 
  head(taxa.nb.print[[ds]]) 
  saveRDS(taxa.nb.print[[ds]], file.path(loc, "03_tabletax", "taxa.nb.print.rds"))
}

saveRDS(taxa.nb, "R_Objects/tax.nb.rds")
saveRDS(taxa.nb.1, "R_Objects/tax.nb.1.rds")


#########################################
##### 11. Table of Tracked Reads ########
#########################################
track.total <- list()
track.p1 <- list()
track.p2 <- list()
track.p1.avg <- list()
track.p2.avg <- list()
track.plot <- list()

for (ds in datasets.pair){
  #11.1: Get column totals
  track.total[[ds]] <- data.frame(t(colSums(track[[ds]])))  
  
  #11.2: Get percentage retained after each step (wrt the intial)
  track.p1[[ds]] <- as.data.frame(track[[ds]]) %>% data.frame() %>% 
    dplyr::mutate (pre.DADA2.NFilt.pct = pre.DADA2.filtered/pre.DADA2.input*100, 
                   Cutadapt.Filt.pct = DADA2.input/pre.DADA2.input*100, 
                   DADA2.TrimNFilt.pct = DADA2.filtered/pre.DADA2.input*100, 
                   Denoised.F.pct = denoisedF/pre.DADA2.input*100, 
                   Denoised.R.pct = denoisedR/pre.DADA2.input*100, 
                   merged.pct = merged/pre.DADA2.input*100, 
                   nonchim.pct = nonchim/pre.DADA2.input*100) %>% 
    dplyr::select(ends_with(".pct")) 
  
  #11.3: Get percentage lost after each step (wrt to the initial)
  track.p2[[ds]] <- as.data.frame(track[[ds]]) %>% data.frame() %>% 
    dplyr::mutate (pre.DADA2.NFilt.pct = (pre.DADA2.input-pre.DADA2.filtered)/pre.DADA2.input*100, 
                   Cutadapt.Filt.pct = (pre.DADA2.filtered-DADA2.input)/pre.DADA2.filtered*100, 
                   DADA2.TrimNFilt.pct = (DADA2.input-DADA2.filtered)/DADA2.input*100, 
                   Denoised.F.pct = (DADA2.filtered-denoisedF)/DADA2.filtered*100, 
                   Denoised.R.pct = (DADA2.filtered-denoisedR)/DADA2.filtered*100, 
                   merged.pct = (((denoisedF-merged)/denoisedF*100)+((denoisedR-merged)/denoisedR*100))/2, 
                   nonchim.pct = (merged-nonchim)/merged*100)  %>% 
    dplyr::select(ends_with(".pct"))
  
  #11.4: get the average percentage retained and lost at each step 
  track.p1.avg[[ds]] <- track.p1[[ds]]  %>% summarize_at(vars(ends_with(".pct")), list(avg = mean)) 
  track.p2.avg[[ds]] <- track.p2[[ds]]  %>% summarize_at(vars(ends_with(".pct")), list(avg = mean))
  
  #11.5: create the plot for track 
  track.plot[[ds]]  <- as.data.frame(track[[ds]]) %>% data.frame() %>% #put the track data into a data frame 
    mutate(Sample = rownames(.)) %>% #add new column with the sample names 
    gather(key = "Step", value = "Reads", -Sample) %>% #reshape data from wide -> long format (required for ggplot boxplot)
    mutate(Step = factor(Step, 
                         levels = c("pre.DADA2.input", "pre.DADA2.filtered", "DADA2.input", 
                                    "DADA2.filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>% #reorder the steps 
    ggplot(aes(x = Step, y = Reads)) + #start plotting ggplot, define x and y axis
    geom_boxplot (varwidth =TRUE ,fill = "cadetblue",alpha =0.3,outlier.colour="red") + #creates a boxplot for each step 
    stat_summary(fun.y = mean, geom = "point", group = 1, color = "blue", size = 3, alpha = 0.7,shape = 18) + #add the mean read count as blue diamonds
    geom_label(data = t(track.p1.avg[[ds]][1:6]) %>% #takes first 6 percentage columns and transposes them vertically
                 data.frame() %>%
                 rename(Percent = 1) %>% #rename col to Percent
                 mutate(Step = c("pre.DADA2.filtered", "DADA2.input", "DADA2.filtered",
                                 "denoisedF", "denoisedR", "merged"),
                        Percent = sprintf("%.2f %%", Percent)), #format the labels (eg. 74.52 %)
               aes(label = Percent),
               y = 1.1 * max(track[[ds]][,3])) +
    geom_label(data = track.p1.avg[[ds]][7] %>% data.frame() %>%
                 rename(total = 1),
               aes(label = paste("Total\nRemaining:\n", #add final total remaining label
                                 round(track.p1.avg[[ds]][1,7], 2), "%")),
               y = mean(track[[ds]][,8]), x = 8.7) +
    geom_label(aes(label = paste("Total reads (%)")), y = 1.1 * #add a total reads label
                 max(track[[ds]][,3]), x = 8.7) +
    geom_label(aes(label = paste("Loss in reads (%)")), y = 1.05 * #add a loss in reads label
                 max(track[[ds]][,3]), x = 8.7) +
    geom_label(data = t(track.p2.avg[[ds]][1:6]) %>% data.frame() %>%
                 rename(Percent = 1) %>%
                 mutate(Step = c("pre.DADA2.filtered", "DADA2.input", "DADA2.filtered",
                                 "denoisedF", "denoisedR", "merged"),
                        Percent = sprintf("%.2f %%", Percent)),
               aes(label = Percent), y = 1.05 * max(track[[ds]][,3])) +
    geom_label(data = track.p2.avg[[ds]][7] %>% data.frame() %>%
                 rename(total = 1),
               aes(label = paste("Lost\nReads:\n", #add percentage lost at each stage
                                 round(track.p2.avg[[ds]][1,7], 2), "%")),
               y = (mean(track[[ds]][,8])), x = 9.3) +
    expand_limits(y = 1.1 * max(track[[ds]][,3]), x = 9.7) + #add extra empty space to fit labels nicely 
    theme_classic() #publication style theme 
  
  #11.6: Save plot 
  num <- match(ds, datasets.all)
  loc <- paste0("2.", num, "_", ds, "_processed")
  pdf(file=file.path(loc, "03_tabletax", "track.plot.pdf"), width=15, height=8.5)
  print(track.plot[[ds]])
  dev.off()
}

for (ds in datasets.single){
  #11.1: Get column totals
  track.total[[ds]] <- data.frame(t(colSums(track[[ds]])))  
  
  #11.2: Get percentage retained after each step (wrt the intial)
  track.p1[[ds]] <- as.data.frame(track[[ds]]) %>% data.frame() %>% 
    dplyr::mutate (pre.DADA2.NFilt.pct = pre.DADA2.filtered/pre.DADA2.input*100, 
                   Cutadapt.Filt.pct = DADA2.input/pre.DADA2.input*100, 
                   DADA2.TrimNFilt.pct = DADA2.filtered/pre.DADA2.input*100, 
                   Denoised.Single.pct = denoisedSingle/pre.DADA2.input*100,
                   nonchim.pct = nonchim/pre.DADA2.input*100) %>% 
    dplyr::select(ends_with(".pct")) 
  
  #11.3: Get percentage lost after each step (wrt to the previous step)
  track.p2[[ds]] <- as.data.frame(track[[ds]]) %>% data.frame() %>% 
    dplyr::mutate (pre.DADA2.NFilt.pct = (pre.DADA2.input-pre.DADA2.filtered)/pre.DADA2.input*100, 
                   Cutadapt.Filt.pct = (pre.DADA2.filtered-DADA2.input)/pre.DADA2.filtered*100, 
                   DADA2.TrimNFilt.pct = (DADA2.input-DADA2.filtered)/DADA2.input*100, 
                   Denoised.Single.pct = (DADA2.filtered-denoisedSingle)/DADA2.filtered*100, 
                   nonchim.pct = (denoisedSingle-nonchim)/denoisedSingle*100)  %>% 
    dplyr::select(ends_with(".pct"))
  
  #11.4: get the average percentage retained and lost at each step 
  track.p1.avg[[ds]] <- track.p1[[ds]]  %>% summarize_at(vars(ends_with(".pct")), list(avg = mean)) 
  track.p2.avg[[ds]] <- track.p2[[ds]]  %>% summarize_at(vars(ends_with(".pct")), list(avg = mean))
  
  #11.5: create the plot for track 
  track.plot[[ds]]  <- as.data.frame(track[[ds]]) %>% data.frame() %>% 
    mutate(Sample = rownames(.)) %>% 
    gather(key = "Step", value = "Reads", -Sample) %>% 
    mutate(Step = factor(Step, 
                         levels = c("pre.DADA2.input", "pre.DADA2.filtered", "DADA2.input", 
                                    "DADA2.filtered", "denoisedSingle", "nonchim"))) %>% 
    ggplot(aes(x = Step, y = Reads)) + 
    geom_boxplot (varwidth =TRUE ,fill = "cadetblue",alpha =0.3,outlier.colour="red") + 
    stat_summary(fun.y = mean, geom = "point", group = 1, color = "blue", size = 3, alpha = 0.7,shape = 18) + 
    geom_label(data = t(track.p1.avg[[ds]][1:4]) %>% 
                 data.frame() %>% 
                 rename(Percent = 1) %>% 
                 mutate(Step = c("pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", 
                                 "denoisedSingle"), 
                        Percent = sprintf("%.2f %%", Percent)), 
               aes(label = Percent), 
               y = 1.1 * max(track[[ds]][,3])) + 
    geom_label(data = track.p1.avg[[ds]][5] %>% data.frame() %>% 
                 rename(total = 1), 
               aes(label = paste("Total\nRemaining:\n", 
                                 round(track.p1.avg[[ds]][1,5], 2), "%")), 
               y = mean(track[[ds]][,6]), x = 6.7) + 
    geom_label(aes(label = paste("Total reads (%)")), y = 1.1 * 
                 max(track[[ds]][,3]), x = 6.7) + 
    geom_label(aes(label = paste("Loss in reads (%)")), y = 1.05 * 
                 max(track[[ds]][,3]), x = 6.7) + 
    geom_label(data = t(track.p2.avg[[ds]][1:4]) %>% data.frame() %>% 
                 rename(Percent = 1) %>% 
                 mutate(Step = c("pre.DADA2.filtered", "DADA2.input", "DADA2.filtered", 
                                 "denoisedSingle"), 
                        Percent = sprintf("%.2f %%", Percent)), 
               aes(label = Percent), y = 1.05 * max(track[[ds]][,3])) + 
    geom_label(data = track.p2.avg[[ds]][5] %>% data.frame() %>% 
                 rename(total = 1),  
               aes(label = paste("Lost\nReads:\n",  
                                 round(track.p2.avg[[ds]][1,5], 2), "%")), 
               y = (mean(track[[ds]][,6])), x = 7.3) + 
    expand_limits(y = 1.1 * max(track[[ds]][,3]), x = 7.7) + 
    theme_classic() 
  
  #11.6: Save plot 
  num <- match(ds, datasets.all)
  loc <- paste0("2.", num, "_", ds, "_processed")
  pdf(file=file.path(loc, "03_tabletax", "track.plot.pdf"), width=15, height=8.5)
  print(track.plot[[ds]])
  dev.off()
}


#########################################
##### 12. Create Phyloseq Object ########
#########################################
#12.1: Import Metadata
metadata_xls <- read_xlsx("AtacamaMetadata.xlsx")
metadata <- data.frame(metadata_xls, row.names = 1) #put into data frame 
head(metadata)

#12.2: Create Phyloseq Object
ps <- list()
meta <- list()
dna <- list()

for (ds in datasets.all){
  meta[[ds]] <- metadata[metadata$Dataset == ds, ]
  ps[[ds]] <- phyloseq(otu_table(seqtab.nochim[[ds]], taxa_are_rows=FALSE), #add ASV table 
                 sample_data(meta[[ds]]), #add metadata information (helps with classification later on)
                 tax_table(taxa.nb.1[[ds]])) #add taxonomy table (analysis related to taxonomy)
  
  #Cleaning the data
  dna[[ds]] <- Biostrings::DNAStringSet(taxa_names(ps[[ds]])) #convert taxa names (long string of ASV actual seq) into DNA seq obj
  names(dna[[ds]]) <- taxa_names(ps[[ds]]) #assign names to each DNA sequence 
  ps[[ds]] <- merge_phyloseq(ps[[ds]], dna[[ds]]) #add DNA Sequence into phyloseq object 
  taxa_names(ps[[ds]]) <- paste0("ASV", seq(ntaxa(ps[[ds]]))) #rename all taxa/ASVs in ps object to simple labels
}

#save ps object
saveRDS(ps, "R_Objects/ps.rds")
saveRDS(dna, "R_Objects/.rds")

#############################################
##### 13. Pre-processing of Phyloseq ########
#############################################
ps.edit <- list()
ASV.taxa.new <- list()
ASV.taxa.new.unassigned <- list()
Unassigned.percent <- list()
Unwanted.taxa <- list()
Mito.chloro.percent <- list()
ps.edit.melt <- list()
ps.edit.unassigned <- list()
ps.edit.chloroplast <- list()
ps.edit.mitochondria <- list()
ps.edit.unassigned.reads.percent <- list()
ps.edit.chloroplast.reads.percent <- list()
ps.edit.mitochondria.reads.percent <- list()

for (ds in datasets.all){
  
  # Check for empty ASVs
  print(ds)
  print(any(taxa_sums(ps[[ds]]) == 0))
  
  # Clean taxonomy names
  ps.edit[[ds]] <- name_na_taxa(
    ps[[ds]],
    na_label = "Unassigned <tax> (<rank>)"
  )
  
  # Extract taxonomy table
  ASV.taxa.new[[ds]] <- as.matrix(tax_table(ps.edit[[ds]]))
  
  # Replace remaining NAs
  ASV.taxa.new[[ds]][is.na(ASV.taxa.new[[ds]])] <- "Unassigned"
  
  # Put taxonomy back into phyloseq object
  tax_table(ps.edit[[ds]]) <- tax_table(ASV.taxa.new[[ds]])
  
  ####################################
  # Unassigned taxa
  ####################################
  
  ASV.taxa.new.unassigned[[ds]] <-
    as.data.frame(ASV.taxa.new[[ds]]) %>%
    filter(if_any(everything(),
                  ~ .x == "Unassigned"))
  
  Unassigned.percent[[ds]] <-
    nrow(ASV.taxa.new.unassigned[[ds]]) /
    nrow(ASV.taxa.new[[ds]]) * 100
  
  ####################################
  # Chloroplast + Mitochondria
  ####################################
  
  Unwanted.taxa[[ds]] <-
    as.data.frame(ASV.taxa.new[[ds]]) %>%
    filter_all(any_vars(. %in%
                          c("Chloroplast",
                            "Mitochondria")))
  
  Mito.chloro.percent[[ds]] <-
    nrow(Unwanted.taxa[[ds]]) /
    nrow(ASV.taxa.new[[ds]]) * 100
  
  ####################################
  # Read counts
  ####################################
  
  ps.edit.melt[[ds]] <- psmelt(ps.edit[[ds]])
  
  ps.edit.unassigned[[ds]] <-
    sum(ps.edit.melt[[ds]]$Abundance[
      ps.edit.melt[[ds]]$Kingdom == "Unassigned"
    ])
  
  ps.edit.chloroplast[[ds]] <-
    sum(ps.edit.melt[[ds]]$Abundance[
      ps.edit.melt[[ds]]$Order == "Chloroplast"
    ])
  
  ps.edit.mitochondria[[ds]] <-
    sum(ps.edit.melt[[ds]]$Abundance[
      ps.edit.melt[[ds]]$Family == "Mitochondria"
    ])
  
  ####################################
  # Percentages of reads
  ####################################
  
  total.reads <- sum(ps.edit.melt[[ds]]$Abundance)
  
  ps.edit.unassigned.reads.percent[[ds]] <-
    ps.edit.unassigned[[ds]] /
    total.reads * 100
  
  ps.edit.chloroplast.reads.percent[[ds]] <-
    ps.edit.chloroplast[[ds]] /
    total.reads * 100
  
  ps.edit.mitochondria.reads.percent[[ds]] <-
    ps.edit.mitochondria[[ds]] /
    total.reads * 100
}


##############################################
##### 13.1. Phyloseq w/o Contaminants ########
##############################################
meta <- list()

ps.no.unassign <- list()
ps.clean <- list()

check.nomito.nochloro.no.unassign <- list()

new.ps.clean <- list()

for(ds in datasets.all){
  
  meta[[ds]] <- metadata[metadata$Dataset == ds, ]
  
  # Remove unassigned taxa
  ps.no.unassign[[ds]] <-
    subset_taxa(
      ps.edit[[ds]],
      Kingdom != "Unassigned"
    )
  
  # Remove chloroplast and mitochondria
  ps.clean[[ds]] <-
    subset_taxa(
      ps.no.unassign[[ds]],
      (Order != "Chloroplast" &
         Family != "Mitochondria") |
        is.na(Family)
    )
  
  ###################################
  # Sanity check 1
  ###################################
  
  cat("\nDataset:", ds, "\n")
  
  cat("Reads before filtering:",
      sum(sample_sums(ps.edit[[ds]])),
      "\n")
  
  cat("Reads after filtering:",
      sum(sample_sums(ps.clean[[ds]])),
      "\n")
  
  ###################################
  # Sanity check 2
  ###################################
  
  check.nomito.nochloro.no.unassign[[ds]] <-
    as.data.frame(tax_table(ps.clean[[ds]])) %>%
    filter_all(
      any_vars(
        . %in% c(
          "chloroplast",
          "Chloroplast",
          "Mitochondria",
          "Unassigned"
        )
      )
    )
  
  print(
    nrow(
      check.nomito.nochloro.no.unassign[[ds]]
    )
  )
  
  ###################################
  # Save cleaned phyloseq object
  ###################################
  
  num <- match(ds, datasets.all)
  loc <- paste0("2.", num, "_", ds, "_processed")
  
  saveRDS(
    ps.clean[[ds]],
    file.path(
      loc,
      "03_tabletax",
      "ps.clean.rds"
    )
  )
  
  ###################################
  # Rebuild phyloseq object
  ###################################
  
  ASV.wide <-
    as.data.frame(otu_table(ps.clean[[ds]])) %>%
    rownames_to_column(var = "Samples") %>%
    as_tibble()
  
  ASV.wide <-
    ASV.wide %>%
    column_to_rownames(var = "Samples")
  
  ASV.wide <- as.matrix(ASV.wide)
  
  ASV.wide.new <-
    otu_table(
      ASV.wide,
      taxa_are_rows = FALSE
    )
  
  samples <-
    sample_data(
      meta[[ds]]
    )
  
  ASV.tax.wide <-
    tax_table(
      ps.clean[[ds]]
    )
  
  ref <- ps.clean[[ds]]@refseq
  
  new.ps.clean[[ds]] <-
    phyloseq(
      ASV.wide.new,
      ASV.tax.wide,
      samples,
      ref
    )
  
  saveRDS(
    new.ps.clean[[ds]],
    file.path(
      loc,
      "03_tabletax",
      "new.ps.clean.rds"
    )
  )
  
}

saveRDS(ps.clean, "R_Objects/ps.clean.rds")
saveRDS(new.ps.clean, "R_Objects/new.ps.clean.rds")


################################
##### 13.2. Rarefaction ########
################################
# Inspect the number of reads per sample 
for(ds in datasets.all){
  cat("\nDataset:", ds, "\n")
  print(sample_sums(new.ps.clean[[ds]]))
  cat("Min: ", min(sample_sums(new.ps.clean[[ds]])), "\n")
}

# List of the rarefaction depth 
# Min from each dataset was chosen
rare.depth <- list()
rare.depth[["D"]] <- 1130 # Removed 2 anomalies with only 197 and 212 reads
rare.depth[["F"]] <- 54432
rare.depth[["S"]] <- 23307
rare.depth[["W"]] <- 1968
rare.depth[["N"]] <- 15708

# Rarefy 
min.seq.depth <- list()
max.seq.depth <- list()
mean.seq.depth <- list()
median.seq.depth <- list()

rarefied.min <- list()

for (ds in datasets.all){
  min.seq.depth[[ds]] <- min(sample_sums(ps.clean[[ds]])) 
  max.seq.depth[[ds]] <- max(sample_sums(ps.clean[[ds]])) 
  mean.seq.depth[[ds]] <- mean(sample_sums(ps.clean[[ds]])) 
  median.seq.depth[[ds]] <- median(sample_sums(ps.clean[[ds]]))

  #goal: make every sample have the same number of reads
  rarefied.min[[ds]] <- rarefy_even_depth(new.ps.clean[[ds]], #input, cleaned ps object
                                          rngseed=123, #set random seed (to allow for reproducible results)
                                          sample.size=rare.depth[[ds]], #how many reads each sample should end up with (in this case the min)
                                          replace=F) #sample without replacement 
}

saveRDS(rarefied.min, "R_Objects/rarefied.min.rds")


#Check 
for (ds in datasets.all){
  summarize_phyloseq(rarefied.min[[ds]])
  sum(sample_sums(rarefied.min[[ds]])) 
}

for (ds in datasets.all){
  cat(
    "\nDataset:", ds,
    "\nSamples before:", nsamples(new.ps.clean[[ds]]),
    "\nSamples after:", nsamples(rarefied.min[[ds]]),
    "\nTotal Reads:", sum(sample_sums(rarefied.min[[ds]])),
    "\n"
  )
}

##############################################
##### 13.3. Get Final Working Dataset ########
##############################################
rarefied.min.prop <- list()
rarefied.min.int <- list()
total <- list()
subtaxa.rarefied.min <- list()

 for (ds in datasets.all){
   #convert counts to percentage 
   rarefied.min.prop[[ds]] <- transform_sample_counts(rarefied.min[[ds]], #input data (phyloseq object)
                                                function(x) x / sum(x)*100) #function to be used
   total[[ds]] <- sample_sums(rarefied.min.prop[[ds]]) #all should be approx 100 (since percentage)
   
   #get the name of remaining taxas
   subtaxa.rarefied.min[[ds]] <- taxa_names(rarefied.min.prop[[ds]])
   
   #using the result on the original table (keep counts not percentages)
   rarefied.min.int[[ds]] <- prune_taxa(subtaxa.rarefied.min[[ds]], rarefied.min[[ds]])
 }

saveRDS(rarefied.min.prop, "R_Objects/rarefied.min.prop.rds") 
saveRDS(rarefied.min.int, "R_Objects/rarefied.min.int.rds")

##Sanity Check
for (ds in datasets.all){
  print(setequal(
    taxa_names(rarefied.min.prop[[ds]]),
    taxa_names(rarefied.min.int[[ds]]))) #should be all TRUE
}



#########################################
##### 13.4. Human Contamination  ########
#########################################
contaminants <- c("Bacteroides","Bifidobacterium","Corynebacterium",
                  "Cutibacterium","Escherichia","Faecalibacterium", 
                  "Haemophilus", "Klebsiella", "Lactobacillus", "Listeria", 
                  "Moraxella", "Neisseria", "Porphyromonas", "Prevotella", 
                  "Propionibacterium", "Salmonella", "Shigella", 
                  "Staphylococcus", "Streptococcus", "Veillonella")

rare.genus <- list()
rare.wide.melt <- list()
contaminants.rarefied.glom <- list()

rare.ps.clean.wide.melt <- list()
contaminants.rarefied.raw <- list()

rare.min.no.contam <- list()

TAX <- list()
metadata <- list()
OTU <- list()

for (ds in datasets.all){
  
  # Aggregate at genus level since contaminants are genus-level
  rare.genus[[ds]] <- tax_glom(
    rarefied.min[[ds]],
    taxrank = "Genus"
  )
  
  # Make into long dataframe
  rare.wide.melt[[ds]] <- psmelt(
    rare.genus[[ds]]
  )
  
  # Percentage contaminants (genus-aggregated)
  contaminants.rarefied.glom[[ds]] <-
    sum(
      rare.wide.melt[[ds]]$Abundance[
        rare.wide.melt[[ds]]$Genus %in% contaminants
      ]
    ) /
    sum(rare.wide.melt[[ds]]$Abundance) * 100
  
  # Sanity check without tax_glom
  rare.ps.clean.wide.melt[[ds]] <- psmelt(
    rarefied.min[[ds]]
  )
  
  contaminants.rarefied.raw[[ds]] <-
    sum(
      rare.ps.clean.wide.melt[[ds]]$Abundance[
        rare.ps.clean.wide.melt[[ds]]$Genus %in% contaminants
      ]
    ) /
    sum(rare.ps.clean.wide.melt[[ds]]$Abundance) * 100
  
  # Remove contaminants
  rare.min.no.contam[[ds]] <-
    subset_taxa(
      rarefied.min[[ds]],
      !(Genus %in% contaminants)
    )
  
  # Final working data
  TAX[[ds]] <- as.matrix(
    as.data.frame(rare.min.no.contam[[ds]]@tax_table)
  )
  
  metadata[[ds]] <- as.matrix(
    as.data.frame(rare.min.no.contam[[ds]]@sam_data)
  )
  
  OTU[[ds]] <- as.matrix(
    as.data.frame(rare.min.no.contam[[ds]]@otu_table)
  )
}

saveRDS(TAX, "R_Objects/TAX.rds")
saveRDS(metadata, "R_Objects/metadata.rds")
saveRDS(OTU, "R_Objects/OTU.rds")

#Sanity Check
for (ds in datasets.all){
  
  cat("\nDataset:", ds, "\n")
  
  cat("Taxa before:",
      ntaxa(rarefied.min[[ds]]), "\n")
  
  cat("Taxa after:",
      ntaxa(rare.min.no.contam[[ds]]), "\n")
  
  cat("Contaminants remaining:",
      sum(
        tax_table(rare.min.no.contam[[ds]])[, "Genus"] %in%
          contaminants,
        na.rm = TRUE
      ),"\n")
}

for (ds in datasets.all){
  
  before <- sum(otu_table(rarefied.min[[ds]]))
  after  <- sum(otu_table(rare.min.no.contam[[ds]]))
  
  cat(
    "\n", ds, "\n",
    "Reads retained:",
    round(after/before*100, 2),
    "%\n"
  )
}

##### Here onwards: Used for testing of different models and plots #####
# For actual models and plot used - refer to atacama_project_figures #

###############################################
##### 14. Alpha Diversity Calculation  ########
###############################################

for (ds in datasets.all){
  
  alpha <- estimate_richness(
    rare.min.no.contam[[ds]],
    measures = c("Shannon")
  )
  
  sample_data(rare.min.no.contam[[ds]])$Shannon.Diversity <- alpha$Shannon
}

saveRDS(rare.min.no.contam, "R_Objects/rare.min.no.contam.rds")

##################################################################
##### 15. Plot Shannon Diversity vs. Moisture (Not Used)  ########
##################################################################
#15.1: Extract all Metadata from Datasets
all.alpha <- data.frame()

for (ds in datasets.all){
  temp <- data.frame(sample_data(rare.min.no.contam[[ds]]))
  temp$Dataset <- ds
  all.alpha <- rbind(all.alpha, temp)
}

saveRDS(all.alpha, "R_Objects/all.alpha.rds")

#15.2: Try Different Models to Plot 

#(1) Simply Plot all the Datasets Together 
library(ggplot2)
library(ggeffects)

pdf("Results/together_linear_plot.pdf", width = 15, height =8.5)

ggplot(all.alpha,
       aes(x = Moisture,
           y = Shannon.Diversity,
           colour = Dataset)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(x = "Moisture (%)",
       y = "Shannon Diversity",
       title = "Alpha Diversity vs Moisture")

dev.off()

m0 <- lm(Shannon.Diversity ~ Moisture, 
         data = all.alpha)

summary(m0)
anova(m0)
#not very useful model 
#Adjusted R Square = 0.008
#p-value for moisture coefficient = 0.247


#(2) Dataset-Adjusted Model (Linear Model)
#Using Dataset as a variable
m1 <- lm(Shannon.Diversity ~ Moisture + Dataset,
  data = all.alpha)

#Plot
pdf("Results/m1_plot.pdf", width = 15, height =8.5)

pred <- ggpredict(
  m1,
  terms = "Moisture [all]"
)

ggplot() +
  geom_point(data = all.alpha,
    aes(x = Moisture,
        y = Shannon.Diversity,
        colour = Dataset),
    size = 2) +
  geom_ribbon(
    data = pred,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2) +
  geom_line(data = pred,
    aes(x = x,
        y = predicted),
    linewidth = 1.2,
    colour = "black") +
  labs(
    x = "Moisture",
    y = "Shannon Diversity",
    colour = "Dataset"
  ) +
  theme_bw()
dev.off()

summary(m1)
anova(m1)
#good model 
#Adjusted R Square = 0.792
#Moisture p-value = 7.88e-07

#Trying without Dataset W
m1.noW <- lm(
  Shannon.Diversity ~ Moisture + Dataset,
  data = subset(all.alpha, Dataset != "W")
)

pdf("Results/m1.noW_plot.pdf", width = 15, height =8.5)

pred <- ggpredict(
  m1.noW,
  terms = "Moisture [all]"
)

ggplot() +
  geom_point(
    data = subset(all.alpha, Dataset != "W"),
    aes(
      x = Moisture,
      y = Shannon.Diversity,
      colour = Dataset
    ),
    size = 3,
    alpha = 0.8
  ) +
  geom_ribbon(
    data = pred,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    fill = "grey70",
    alpha = 0.3
  ) +
  geom_line(
    data = pred,
    aes(
      x = x,
      y = predicted
    ),
    colour = "black",
    linewidth = 1.2
  ) +
  labs(
    x = "Soil moisture",
    y = "Shannon diversity",
    colour = "Dataset"
  ) +
  theme_bw(base_size = 14)
dev.off()

summary(m1.noW)
#Adjusted R square = drops slightly 
#Moisture p-value = 4.54e-08, more significant 
#Proves that moisture remains significant even without the problematic dataset (conclusions not driven by problematic dataset)


#(3) Quadratic Data-Ajusted Model 
m2 <- lm(Shannon.Diversity ~ Moisture + I(Moisture^2) + Dataset,
         data = all.alpha)

pdf("Results/m2_plot.pdf", width = 15, height =8.5)
pred <- ggpredict(
  m2,
  terms = "Moisture [all]"
)

ggplot() +
  geom_point(
    data = all.alpha,
    aes(
      x = Moisture,
      y = Shannon.Diversity,
      colour = Dataset
    ),
    size = 2
  ) +
  geom_ribbon(
    data = pred,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.2
  ) +
  geom_line(
    data = pred,
    aes(
      x = x,
      y = predicted
    ),
    linewidth = 1.2,
    colour = "black"
  ) +
  labs(
    x = "Moisture",
    y = "Shannon Diversity",
    colour = "Dataset"
  ) +
  theme_bw()
dev.off()

summary(m2)
anova(m2)
#Even better model 
#Adjusted R Square = 0.844
#Moisture^2 p-value = 5.885e-09
#Moisture p-value = 0.005093

#Without W
m2.noW <- lm(
  Shannon.Diversity ~ Moisture + I(Moisture^2) + Dataset,
  data = subset(all.alpha, Dataset != "W")
)

pdf("Results/m2.noW_plot.pdf", width = 15, height =8.5)
pred <- ggpredict(
  m2.noW,
  terms = "Moisture [all]"
)

ggplot() +
  geom_point(
    data = subset(all.alpha, Dataset != "W"),
    aes(
      x = Moisture,
      y = Shannon.Diversity,
      colour = Dataset
    ),
    size = 2
  ) +
  geom_ribbon(
    data = pred,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.2
  ) +
  geom_line(
    data = pred,
    aes(
      x = x,
      y = predicted
    ),
    linewidth = 1.2,
    colour = "black"
  ) +
  labs(
    x = "Moisture",
    y = "Shannon Diversity",
    colour = "Dataset"
  ) +
  theme_bw()
dev.off()

summary(m2.noW)


### New with Number of Reads
m.new <- lm(Shannon.Diversity ~ Moisture + Reads,
         data = all.alpha)

pdf("Results/m.new_plot.pdf", width = 15, height =8.5)

pred <- ggpredict(
  m1,
  terms = "Moisture [all]"
)

ggplot() +
  geom_point(data = all.alpha,
             aes(x = Moisture,
                 y = Shannon.Diversity,
                 colour = Dataset),
             size = 2) +
  geom_ribbon(
    data = pred,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2) +
  geom_line(data = pred,
            aes(x = x,
                y = predicted),
            linewidth = 1.2,
            colour = "black") +
  labs(
    x = "Moisture",
    y = "Shannon Diversity",
    colour = "Dataset"
  ) +
  theme_bw()
dev.off()
#bubble plot by phylum and species (top 10) -> categorise into the low, medium and high 


#############################################################
##### 16. Plot Simpson Index Diversity vs. Moisture  ########
#############################################################
# Create Table to Build Model 
for (ds in datasets.all){
  
  alpha <- estimate_richness(rare.min.no.contam[[ds]],
    measures = c("Observed", "Simpson"))
  
  sample_data(rare.min.no.contam[[ds]])$Observed.Richness <- alpha$Observed
  sample_data(rare.min.no.contam[[ds]])$Simpson.Diversity <- alpha$Simpson
}


all.alpha <- data.frame()

for (ds in datasets.all){
  
  all.alpha <- rbind(
    all.alpha,
    data.frame(sample_data(rare.min.no.contam[[ds]]))
  )
}

all.alpha.noW <- subset(all.alpha, Dataset != "W")
saveRDS(all.alpha.noW, "R_Objects/all.alpha.noW.rds")

all.alpha.noW$logMoisture <- log(all.alpha.noW$Moisture)
saveRDS(all.alpha.noW, "R_Objects/all.alpha.noW.log.rds")

all.alpha.noW$logMoisture.centered <- all.alpha.noW$logMoisture - mean(all.alpha.noW$logMoisture)
saveRDS(all.alpha.noW, "R_Objects/all.alpha.noW.log.centre.rds")

# Model the Data: 
#(1) Simpson Diversity 
m.simpson <- lm(Simpson.Diversity ~ logMoisture.centered + Dataset,
  data = all.alpha.noW.log.centre)

summary(m.simpson)

#(2) Observed Richness
m.richness <- lm(Observed.Richness ~ logMoisture.centered + Dataset,
                data = all.alpha.noW.log.centre)

summary(m.richness)

# Plot the Models 
library(ggplot2)
library(ggeffects)

#(1) Simpson Diversity 
pdf("Results/m.simpson_plot.pdf", width = 8.5, height =15)
pred <- ggpredict(
  m.simpson,
  terms = "logMoisture.centered [all]"
)

ggplot() +
  geom_point(
    data = all.alpha.noW.log.centre,
    aes(
      x = logMoisture.centered,
      y = Simpson.Diversity,
      colour = Dataset
    ),
    size = 2
  ) +
  geom_ribbon(
    data = pred,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.2
  ) +
  geom_line(
    data = pred,
    aes(
      x = x,
      y = predicted
    ),
    colour = "black",
    linewidth = 1.2
  ) +
  labs(
    title = "Relationship Between Simpson Diversity and Soil Moisture",
    x = "Log-Transformed Soil Moisture (Centered)",
    y = "Simpson Diversity",
    colour = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(
    hjust = 0.5,
    face = "bold"))
dev.off()


#(2) Observed Richness 
pdf("Results/m.richness_plot.pdf", width = 8.5, height =15)
pred <- ggpredict(
  m.richness,
  terms = "logMoisture.centered [all]"
)

ggplot() +
  geom_point(
    data = all.alpha.noW,
    aes(
      x = logMoisture.centered,
      y = Observed.Richness,
      colour = Dataset
    ),
    size = 2
  ) +
  geom_ribbon(
    data = pred,
    aes(
      x = x,
      ymin = conf.low,
      ymax = conf.high
    ),
    alpha = 0.2
  ) +
  geom_line(
    data = pred,
    aes(
      x = x,
      y = predicted
    ),
    colour = "black",
    linewidth = 1.2
  ) +
  labs(
    title = "Relationshp Between Species Richness and Soil Moisture",
    x = "Log-Transformed Soil Moisture (centered)",
    y = "Species Richness",
    colour = "Dataset"
  ) +
  theme_bw() +
  theme(plot.title = element_text(
      hjust = 0.5,
      face = "bold"))
dev.off()




#######################################################
##### 17. Bubble/Bar Plots for Top 10 Taxonomy ########
#######################################################

#######################################################
##### 17.1 Top 10 Phyla Across Moisture Groups ########
#######################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(ggsci)

##### Prepare Data #####
rare.min.no.contam <- readRDS(
  "R_Objects/rare.min.no.contam.rds"
)

# Remove W
rare.min.no.contam.noW <- rare.min.no.contam[
  names(rare.min.no.contam) != "W"
]

# Merge phyloseq objects
combined.phy <- do.call(
  merge_phyloseq,
  rare.min.no.contam.noW
)

# Moisture groups
sample_data(combined.phy)$MoistureGroup <- cut(
  sample_data(combined.phy)$Moisture,
  breaks = c(-Inf, 2.5, 7.5, Inf),
  labels = c("Low", "Medium", "High"),
  right = FALSE
)

# Keep only Bacteria and Archaea
combined.phy <- subset_taxa(
  combined.phy,
  Kingdom %in% c("Bacteria", "Archaea")
)

##### Relative Abundance ######
phy.rel <- transform_sample_counts(
  combined.phy,
  function(x) x / sum(x)
)

##### Aggregate at Phylum Level ######
phy.phylum <- tax_glom(
  phy.rel,
  taxrank = "Phylum"
)

phy.phylum.df <- psmelt(
  phy.phylum
)

##### Handle Unassigned ######
phy.phylum.df$Phylum <- ifelse(
  is.na(phy.phylum.df$Phylum) |
    grepl("Unassigned", phy.phylum.df$Phylum),
  "Unassigned",
  as.character(phy.phylum.df$Phylum)
)

##### Top 10 Phyla Overall ######
top10.phyla <- phy.phylum.df %>%
  group_by(Phylum) %>%
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

##### Create Other Category ######
phy.phylum.df$Phylum <- ifelse(
  phy.phylum.df$Phylum %in%
    c(top10.phyla, "Unassigned"),
  phy.phylum.df$Phylum,
  "Other"
)

##### Mean Relative Abundance ######
plot.df <- phy.phylum.df %>%
  group_by(
    MoistureGroup,
    Phylum
  ) %>%
  summarise(
    Abundance = mean(Abundance),
    .groups = "drop"
  )

##### Force Bars to Sum to 100% ######
plot.df <- plot.df %>%
  group_by(MoistureGroup) %>%
  mutate(
    Abundance = Abundance / sum(Abundance)
  ) %>%
  ungroup()

##### Put Other and Unassigned at Bottom #####
phylum.order <- c(
  "Other",
  "Unassigned",
  setdiff(
    unique(plot.df$Phylum),
    c("Other", "Unassigned")
  )
)

plot.df$Phylum <- factor(
  plot.df$Phylum,
  levels = phylum.order
)

##### Plot Stacked Bar Plot ######
p.phylum <- ggplot(
  plot.df,
  aes(
    x = MoistureGroup,
    y = Abundance * 100,
    fill = Phylum
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  scale_fill_d3(palette = "category20") +
  theme_bw() +
  labs(
    title = "Top 10 Phyla Across Moisture Groups",
    x = "Moisture Group",
    y = "Mean Relative Abundance (%)"
  )

pdf(
  "Results/Top10Phylum_MoistureGroups.pdf",
  width = 8,
  height = 6
)

p.phylum

dev.off()

##### Create Bubble Plot Data Frame ######
bubble.df <- phy.phylum.df %>%
  filter(Phylum %in% top10.phyla) %>%
  group_by(
    MoistureGroup,
    Phylum
  ) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  )

##### Order Phyla by Overall Abundance ######
phylum.order <- bubble.df %>%
  group_by(Phylum) %>%
  summarise(
    OverallAbundance = sum(MeanAbundance)
  ) %>%
  arrange(desc(OverallAbundance)) %>%
  pull(Phylum)

bubble.df$Phylum <- factor(
  bubble.df$Phylum,
  levels = rev(phylum.order)
)

##### Plot Bubble Plot ######
pdf(
  "Results/Top10Phylum_bubble_plot.pdf",
  width = 8.5,
  height = 8
)

ggplot(
  bubble.df,
  aes(
    x = MoistureGroup,
    y = Phylum,
    size = MeanAbundance * 100
  )
) +
  geom_point(
    colour = "#1F4E79",
    alpha = 0.8) +
  scale_size_continuous(
    name = "Mean abundance (%)",
    range = c(1, 25)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5
    )
  )+
  labs(
    title = "Top 10 Phyla Within Each Moisture Group",
    x = "Moisture Group",
    y = "Phylum"
  )

dev.off()

#######################################################
##### 17.2 Top 10 Genus Across Moisture Groups ########
#######################################################
##### Aggregate at Genus Level #####
phy.genus <- tax_glom(
  phy.rel,
  taxrank = "Genus"
)

phy.genus.df <- psmelt(
  phy.genus
)

##### Handle Unassigned #####
phy.genus.df$Genus <- ifelse(
  is.na(phy.genus.df$Genus) |
    grepl("Unassigned", phy.genus.df$Genus),
  "Unassigned",
  as.character(phy.genus.df$Genus)
)

phy.genus.df$Phylum <- ifelse(
  is.na(phy.genus.df$Phylum) |
    grepl("Unassigned", phy.genus.df$Phylum),
  "Unassigned",
  as.character(phy.genus.df$Phylum)
)

#Keep the phylum information
genus.phylum <- phy.genus.df %>%
  select(
    Genus,
    Phylum
  ) %>%
  distinct()

##### Get Top 10 Actual Genera #####
top10.genus <- phy.genus.df %>%
  filter(
    Genus != "Unassigned"
  ) %>%
  group_by(Genus) %>%
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)


##### Create Other Category #####
phy.genus.df$Genus <- ifelse(
  phy.genus.df$Genus %in%
    c(top10.genus, "Unassigned"),
  phy.genus.df$Genus,
  "Other"
)

##### Mean Relative Abundance #####
plot.df.genus <- phy.genus.df %>%
  group_by(
    MoistureGroup,
    Genus
  ) %>%
  summarise(
    Abundance = mean(Abundance),
    .groups = "drop"
  )

##### Force Bars to Sum to 100% #####
plot.df.genus <- plot.df.genus %>%
  group_by(MoistureGroup) %>%
  mutate(
    Abundance = Abundance / sum(Abundance)
  ) %>%
  ungroup()

##### Put Other and Unassigned at Bottom #####
genus.order <- c(
  "Other",
  "Unassigned",
  setdiff(
    unique(plot.df.genus$Genus),
    c("Other", "Unassigned")
  )
)

plot.df.genus$Genus <- factor(
  plot.df.genus$Genus,
  levels = genus.order
)

##### Check Totals #####
plot.df.genus %>%
  group_by(MoistureGroup) %>%
  summarise(
    Total = sum(Abundance)
  )

##### Plot Stacked Bar Plot #####
p.genus <- ggplot(
  plot.df.genus,
  aes(
    x = MoistureGroup,
    y = Abundance * 100,
    fill = Genus
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  scale_fill_d3(palette = "category20") +
  theme_bw() +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Mean Relative Abundance (%)"
  )

pdf(
  "Results/Top10Genus_MoistureGroups.pdf",
  width = 8,
  height = 6
)

p.genus

dev.off()

##### Plot Faceted Bar Plot #####
plot.df.genus <- phy.genus.df %>%
  filter(
    Genus %in% top10.genus
  ) %>%
  group_by(
    MoistureGroup,
    Phylum,
    Genus
  ) %>%
  summarise(
    Abundance = mean(Abundance),
    .groups = "drop"
  )

#Order Genera Within Phyla 
genus.order <- plot.df.genus %>%
  group_by(
    Phylum,
    Genus
  ) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  ) %>%
  arrange(
    Phylum,
    desc(MeanAbundance)
  ) %>%
  pull(Genus)

plot.df.genus$Genus <- factor(
  plot.df.genus$Genus,
  levels = genus.order
)

#Order Phyla 
phylum.order <- plot.df.genus %>%
  group_by(Phylum) %>%
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  pull(Phylum)

plot.df.genus$Phylum <- factor(
  plot.df.genus$Phylum,
  levels = phylum.order
)

#Plot 
p.genus.faceted <- ggplot(
  plot.df.genus,
  aes(
    x = MoistureGroup,
    y = Abundance * 100,
    fill = Genus
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  facet_wrap(
    ~ Phylum,
    scales = "fixed",
    nrow = 1
  ) +
  scale_fill_d3(
    palette = "category20"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.title = element_text(
      face = "bold"
    )
  ) +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Mean Relative Abundance (%)",
    fill = "Genus"
  )

#Save Plot 
pdf(
  "Results/Top10Genus_FacetedByPhylum.pdf",
  width = 8,
  height = 6
)

p.genus.faceted

dev.off()


##### Plot Heatmap #####
p.heatmap <- ggplot(
  plot.df.genus,
  aes(
    x = MoistureGroup,
    y = Genus,
    fill = Abundance * 100
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  facet_grid(
    Phylum ~ .,
    scales = "free_y",
    space = "free_y",
  ) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Mean \nAbundance (%)"
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside", 
    legend.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.3), 
    legend.title = element_text(
      size = 9,
      face = "bold", 
      hjust = 0.5),
    legend.text = element_text(
      size = 9),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.3, "cm"),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold", 
      margin = margin(b = 10)
    ),
    strip.text.y = element_text(
      face = "bold",
      angle = 0
    ),
    panel.grid = element_blank(), 
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 15)),
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 15)), 
    axis.text.y = element_text(size = 10, face = "italic", colour = "black"), 
    axis.text.x = element_text(colour = "black")
    ) +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Genus")

#Save pdf
pdf(
  "Results/Top10Genus_Heatmap.pdf",
  width = 6,
  height = 6
)

p.heatmap

dev.off()

##### Bubble Plot #####
pdf(
  "Results/Top10Genus_bubble_plot.pdf",
  width = 8,
  height = 8
)

ggplot(
  plot.df.genus,
  aes(
    x = MoistureGroup,
    y = Genus,
    size = Abundance * 100,
    fill = Phylum
  )
) +
  geom_point(
    shape = 21,
    colour = "black",
    alpha = 0.8
  ) +
  scale_size_continuous(
    name = "Mean Abundance (%)", 
    range = c(2,20)
  ) +
  scale_fill_d3(palette = "category20") +
  theme_bw() +
  labs(
    title = "Top 10 Genera Across Moisture Groups",
    x = "Moisture Group",
    y = "Genus",
    fill = "Phylum"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.background = element_rect(
      colour = "black",
      fill = "white"
    ), 
    legend.title = element_text(
      face = "bold"
    ), 
    axis.text.y = element_text(
      face = "italic", 
      colour = "black"
    ), 
    axis.text.x = element_text(
      colour = "black"
    ), 
    axis.title.x = element_text(
      face = "bold", 
      margin = margin(t = 15)
    ),
    axis.title.y = element_text(
      face = "bold", 
      margin = margin(r = 15)
    )
  )

dev.off()
