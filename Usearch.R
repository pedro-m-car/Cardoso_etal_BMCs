Raquel 16S samples-Coral Pedro  V3-4 regions
author: "Marwan"
date: "8/9/2021"
output: html_document

This script includes the use of the new database DTGB and BLCA method. Also includes chloroplast removal and deseq normalisation method.
---

#Marwan script
install.packages("seqinr")
library(seqinr)
library(dplyr)
system("ln -s 00_DataNeeded/usearch/usearch11.0.667_i86osx64 usearch11.0.667 ; chmod +x usearch11.0.667") # MacOS
system('./usearch11.0.667') # Test

# Quality filtering
dir.create("02_PreProcessedData")
dir.create("03_MergedData")
dir.create("04_FilteredData")
dir.create("05_TrimmedData")
dir.create("06_DereplicatedData")

files = dir(path = '01_RawData/', pattern = '_R1_')
files
for (file1 in files)
{
  file2 = gsub(x = file1, pattern = '_R1_', replacement = '_R2_')
  filename = gsub(x = file1, pattern = "_L001_R1_001.fastq.*", replacement = '', perl = T)
  
  #	Quality trimming
  command <- paste('java -jar 00_DataNeeded/Trimmomatic-0.38/trimmomatic-0.38.jar PE 01_RawData/', file1, " 01_RawData/", file2, " 02_PreProcessedData/", filename, "_pF.fastq 02_PreProcessedData/", filename, "_upF.fastq 02_PreProcessedData/",filename, "_pR.fastq 02_PreProcessedData/",filename, "_upR.fastq  SLIDINGWINDOW:4:15 MINLEN:100", sep = "")
  system(command)
  
  #	Merge paired reads (range needs to be adjusted)
  command = paste("./usearch11.0.667 -fastq_mergepairs 02_PreProcessedData/", filename, "_pF.fastq -reverse 02_PreProcessedData/", filename, "_pR.fastq -fastqout 03_MergedData/",filename, ".fastq -relabel @ -fastq_maxdiffs 5 -fastq_pctid 80 -fastq_minmergelen 250 -fastq_maxmergelen 550 -sample ", filename, sep = "")
  system(command)
  
  # #	Quality filtering
  command = paste("./usearch11.0.667 -fastq_filter 03_MergedData/", filename, ".fastq -fastaout 04_FilteredData/", filename, ".fasta -fastq_maxns 1 -fastq_maxee 1", sep = "")
  system(command)
  
  # #	Check and remove primers
  FQ = read.fasta(file = paste0("04_FilteredData/", filename, ".fasta"), as.string = T)
  for_primer_present = grep(pattern = "^[ATGC]{0,2}CCTACGGG[ATGC]GGC[AT]GCAG", x = FQ, ignore.case = T)
  rev_primer_present = grep(pattern = "GGATTAGATACCC[CGT][AGT]GTAGTC[ATGC]{0,2}$", x = FQ, ignore.case = T)
  both_present = intersect(for_primer_present, rev_primer_present)
  
  FQ1 = FQ[both_present]
  IDs = names(FQ1)
  SEQs = as.character(FQ1)
  SEQs = gsub(pattern = "^[ATGC]{0,2}CCTACGGG[ATGC]GGC[AT]GCAG", replacement = "", SEQs, ignore.case = T)
  SEQs = gsub(pattern = "GGATTAGATACCC[CGT][AGT]GTAGTC[ATGC]{0,2}$", replacement = "", SEQs, ignore.case = T)
  OUT = file(description = paste0("05_TrimmedData/", filename, ".fasta"), open = "w")
  for(i in 1:length(IDs)) write(x = paste0(">", IDs[i], "\n", SEQs[i]), file = OUT)
  close(OUT)
  #	Dereplication
  command = paste("./usearch11.0.667  -fastx_uniques 05_TrimmedData/", filename, ".fasta -fastaout 06_DereplicatedData/", filename, ".fasta -sizeout", sep = "")
  system(command)
}





#below is an example of what the data looks like when merging and trimming.



rm(list=ls())

#Check merging range - checked and changed to 250-550 (this was the most appropriate for the current study)

seq_lengths = NULL
for(f in dir(path = '05_TrimmedData/', full.names = T))
{
  print(f)
  fasta_file = read.fasta(f)
  seq_lengths = c(seq_lengths, as.numeric(lengths(fasta_file)))
}
hist(seq_lengths) # USe in line 27/28 to adjust range
quantile(seq_lengths)

#0%    25%   50%   75%   100% 
# 213  402  404  424  512

#	Join
dir.create("07_JoinedData")

# Please decomment
system("cat 06_DereplicatedData/*.fasta > 07_JoinedData/AllSamples.fasta") # Linux and MacOS
#shell("type 06_DereplicatedData\\*.fasta > 07_JoinedData\\AllSamples.fasta") # Windows

FNA = readLines("07_JoinedData/AllSamples.fasta")
FNA[grep(pattern = ">", x = FNA, invert = T)] = toupper(FNA[grep(pattern = ">", x = FNA, invert = T)])
write(x = FNA, file = "07_JoinedData/AllSamples2.fasta")

#	Dereplication


dir.create("08_UniqueSequences")
system("./usearch11.0.667 -fastx_uniques 07_JoinedData/AllSamples2.fasta -fastaout 08_UniqueSequences/AllSamples_uniques.fasta -sizein -sizeout -strand both")


#	Generating unique sequences using UNOISE
dir.create("09_DenoisedSequences")

system("./usearch11.0.667 -unoise3 08_UniqueSequences/AllSamples_uniques.fasta -zotus 09_DenoisedSequences/AllSamples_denoised.fasta")

#	Chimera Removal
dir.create("10_UchimeReference")

system("./usearch11.0.667 -uchime2_ref 09_DenoisedSequences/AllSamples_denoised.fasta -db 00_DataNeeded/db_GTDB202/ssu_all_r202.fna -strand plus -mode high_confidence -notmatched 10_UchimeReference/AllSamples_unoise_nc.fasta")

#	OTU table generation
# An OTU table is made by the otutab command
dir.create("11_OtuTable")

system("./usearch11.0.667 -otutab 07_JoinedData/AllSamples2.fasta -otus 10_UchimeReference/AllSamples_unoise_nc.fasta -id 0.97 -otutabout 11_OtuTable/AllSamples_unoise_otu_table1.txt")

#instructions on how to install Blast can be accessed here (https://urldefense.com/v3/__https://github.com/mhahsler/rBLAST__;!!Nmw4Hv0!1PDyPa3KQeFk21XjzqkeuBG4E7erIuZki0r-r9xb4_t0RB82zTp7RELp7dcSRbLFaZBUdMrHQm0tcukcmtxST-KKNFUpsQ$ )
#	Mapping of OTUs on Reference Database

library(tidyverse)
library(reshape2)

##############################################################################################################
#################################################### BLCA ####################################################
##############################################################################################################

# add clustalo and muscle to system path (required by BLCA)
system("cp 00_DataNeeded/clustalo/clustal-omega-1.2.3-macosx 00_DataNeeded/clustalo/clustalo ; chmod +x 00_DataNeeded/clustalo/clustalo") # MacOS
system("cp 00_DataNeeded/muscle/muscle3.8.31_i86darwin64 00_DataNeeded/muscle/muscle ; chmod +x 00_DataNeeded/muscle/muscle") # MacOS
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/clustalo", sep = ":"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/muscle", sep = ":"))

dir.create("12_TaxAssignmentGTDB_BLCA")

# run BLCA agsinst GTDB SSU database

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

# r89
# system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.taxonomy -q 00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.fasta -i 10_UchimeReferenceGTDB/zOTU_nc.fasta -o 12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

# r95
system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/db_GTDB202/ssu_all_r202_BLCAparsed.taxonomy -q 00_DataNeeded/db_GTDB202/ssu_all_r202_BLCAparsed.fasta -i 10_UchimeReference/AllSamples_unoise_nc.fasta -o 12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')



# remove taxonomic ranks from BLCA's classification and put number in parenthesis (make it easy to read)
m = 1
for (each_line in readLines(file('12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt',open="r")) ){
  each_line_split = strsplit(each_line, '\t')
  OTU_ID = each_line_split[[1]][1]
  taxonomy = each_line_split[[1]][2]
  taxonomy_split = strsplit(taxonomy, ';')
  taxonomy_no_rank = ''
  n = 1
  for (taxon in taxonomy_split[[1]]){
    if (n%%2 == 1){
      taxon_split = strsplit(taxon, ':')
      if (length(taxon_split[[1]]) ==2)
      {taxon_no_rank = taxon_split[[1]][2]} 
      else 
      {taxon_no_rank = taxon_split[[1]][1]}
      taxonomy_no_rank = paste(taxonomy_no_rank, taxon_no_rank, sep = ");")} 
    else 
    {taxonomy_no_rank = paste(taxonomy_no_rank, taxon, sep = "(")}
    n = n + 1
  }
  taxonomy_no_rank = paste(taxonomy_no_rank, ')', sep = "")
  taxonomy_no_rank = substr(taxonomy_no_rank, 3, nchar(taxonomy_no_rank))
  if (taxonomy_no_rank == "Unclassified)"){taxonomy_no_rank = "Unclassified"}
  
  taxonomy_no_rank_with_OTU = paste(OTU_ID, taxonomy_no_rank, sep = "\t")
  
  if (m == 1)
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=FALSE)} 
  else 
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=TRUE)}
  m = m +1
}


# # Merge Table and Taxonomy
dir.create("13_FinalOtuTableGTDB_BLCA")
OTU = read.delim("11_OtuTable/AllSamples_unoise_otu_table1.txt", header = T)
#I manually removed the numbers in parentesis after the taxonomy e.g. (100.0000). Then I saved the file as "AllSamples_unoise_BLCA_out.2filtered2_split.txt". here the taxonomy is split.
TAX = read.delim("12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt", header = F)
names(TAX) = c("X.OTU.ID", "Taxonomy")
# 
OTU_TAX = merge(OTU, TAX, by = "X.OTU.ID")
write.table(OTU_TAX, "13_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# 
rm(list=ls(all=TRUE))

dir.create("12_TaxAssignmentSILVA_BLCA")

# run BLCA agsinst GTDB SSU database

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

# r89
# system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.taxonomy -q 00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.fasta -i 10_UchimeReferenceGTDB/zOTU_nc.fasta -o 12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

# r95
system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/db_SILVA_138/SILVA_138_SSURef_NR99_tax_silva_BLCAparsed.taxonomy -q 00_DataNeeded/db_SILVA_138/SILVA_138_SSURef_NR99_tax_silva_BLCAparsed.fasta -i 10_UchimeReference/AllSamples_unoise_nc.fasta -o 12_TaxAssignmentSILVA_BLCA/AllSamples_unoise_BLCA_out.1.txt')



# remove taxonomic ranks from BLCA's classification and put number in parenthesis (make it easy to read)
m = 1
for (each_line in readLines(file('12_TaxAssignmentSILVA_BLCA/AllSamples_unoise_BLCA_out.1.txt',open="r")) ){
  each_line_split = strsplit(each_line, '\t')
  OTU_ID = each_line_split[[1]][1]
  taxonomy = each_line_split[[1]][2]
  taxonomy_split = strsplit(taxonomy, ';')
  taxonomy_no_rank = ''
  n = 1
  for (taxon in taxonomy_split[[1]]){
    if (n%%2 == 1){
      taxon_split = strsplit(taxon, ':')
      if (length(taxon_split[[1]]) ==2)
      {taxon_no_rank = taxon_split[[1]][2]} 
      else 
      {taxon_no_rank = taxon_split[[1]][1]}
      taxonomy_no_rank = paste(taxonomy_no_rank, taxon_no_rank, sep = ");")} 
    else 
    {taxonomy_no_rank = paste(taxonomy_no_rank, taxon, sep = "(")}
    n = n + 1
  }
  taxonomy_no_rank = paste(taxonomy_no_rank, ')', sep = "")
  taxonomy_no_rank = substr(taxonomy_no_rank, 3, nchar(taxonomy_no_rank))
  if (taxonomy_no_rank == "Unclassified)"){taxonomy_no_rank = "Unclassified"}
  
  taxonomy_no_rank_with_OTU = paste(OTU_ID, taxonomy_no_rank, sep = "\t")
  
  if (m == 1)
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentSILVA_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=FALSE)} 
  else 
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentSILVA_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=TRUE)}
  m = m +1
}


# # Merge Table and Taxonomy
dir.create("13_FinalOtuTableSILVA_BLCA")
OTU = read.delim("11_OtuTable/AllSamples_unoise_otu_table1.txt", header = T)
#I manually removed the numbers in parentesis after the taxonomy e.g. (100.0000). Then I saved the file as "AllSamples_unoise_BLCA_out.2filtered2_split.txt". here the taxonomy is split.
TAX = read.delim("12_TaxAssignmentSILVA_BLCA/AllSamples_unoise_BLCA_out.2.txt", header = F)
names(TAX) = c("X.OTU.ID", "Taxonomy")
# 
OTU_TAX = merge(OTU, TAX, by = "X.OTU.ID")
write.table(OTU_TAX, "13_FinalOtuTableSILVA_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# 

###### now we need to remnove the chloroplast first
dir.create("14_Deseqnormalisation")
## copy the final OTU table to the new folder and proceed with the script below:


##now deseq normalisation

#install.packages("BiocManager")
#BiocManager::install('phyloseq')

library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape")
library("gridExtra")
library("tidyverse")
library("dplyr")
library(RColorBrewer)
library("Trial/R/functions.R")

AllSamples <- read.delim("14_Deseqnormalisation/AllSamples_unoise_otu_table_BLCA_filtered2.txt")
head(AllSamples)

#separating taxonomy
taxon_table <- select(AllSamples, one_of(c("X.OTU.ID", "Taxonomy")))
taxon_table <- separate(taxon_table, Taxonomy, into = (c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")), sep = ";", extra = "merge", fill = "right")
taxon_table <- remove_rownames(taxon_table)
taxon_table <- column_to_rownames(taxon_table, var = "X.OTU.ID")
taxon_table <- as.matrix(taxon_table)
#export taxonomy table if needed for later
write.csv(taxon_table, "14_Deseqnormalisation/taxonomy_table.csv")


# Data normalisation using DESeq2
# Two methods for DESeq2 normalisation are described here - using phyloseq package and using just DESeq2 functions. Both methods should produced excatly the same outputs


## if you have the new R version you will need to install DESeq2 as instructions below:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


library(phyloseq)
library(DESeq2)



# OTU table includes raw counts, displaying OTUs as row names and samples in columns
otu_table <- read.csv("14_Deseqnormalisation/AllSamples_unoise_otu_table_BLCA_filtered2.csv", row.names = 1, stringsAsFactors = FALSE)
head(otu_table)
otu_table[1:5, 1:5]
dim(otu_table) #  4257 OTUs and 36 samples 

# metadata object
# metadata object includes all relevant information about the study, and also displays samples as row names
fact <- read.csv("14_Deseqnormalisation/metafactors.csv",)
head(fact)
str(fact)
# include samples as row names
rownames(fact) <- fact$Diet

# make sure samples are displayed in the same order in otu_table and fact objects
# rearrange objects if necessary
colnames(otu_table) == rownames(fact)

# taxa object
# taxa object includes OTU_ID as row names and taxon as columns
tax_table <- read.csv("14_Deseqnormalisation/taxonomy_table.csv", row.names = 1, stringsAsFactors = FALSE)
head(tax_table)

# make sure OTUs are displayed in the same order in otu_table and tax_table objects
# rearrange objects if necessary
rownames(otu_table) == rownames(tax_table)
sum(rownames(otu_table) == rownames(tax_table))
sum(rownames(otu_table) != rownames(tax_table)) # 0 diff

# check differences in library size
# total number of sequences per sample
colSums(otu_table)
summary(colSums(otu_table))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#15563   23370   35064   33728   40944   61023 

# max difference in library size
max(colSums(otu_table))/min(colSums(otu_table)) # 3.921031
# definitely needs normalisation

#####################################################
#### DESeq2 normalisation using phyloseq package ####
#####################################################

# more info here:
# vignette("phyloseq-mixture-models")

# make phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax_table)),
               sample_data(fact),
               otu_table(as.matrix(otu_table), taxa_are_rows = TRUE))

# check if physeq object looks fine
ps
head(sample_data(ps), 10)
head(otu_table(ps),10)
head(tax_table(ps),10)

# make DESeq2 object
# DESeq2 doesn't work very well when using multiple factor models. So, instead of specifying individual factors (e.g. Time + Fertiliser), create an identifier that combines information from factors to be tested (in this case, column "TF" combines time + fertiliser)
head(sample_data(ps)$Diet, 10)
# time should appear as a factor with n levels

dds <- phyloseq_to_deseq2(ps, ~ Diet)

# Make a copy of your phyloseq object, which you will then modify with VST values
# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType = "local")

# visualise size factor
sizeFactors(dds)
barplot(sizeFactors(dds), las = 2, ylab = "Size factor normalisation", col = fact$treatment)

# normalised counts
deseq_otu_table <- otu_table(counts(dds, normalize = TRUE), taxa_are_rows = TRUE)
head(deseq_otu_table)

# round normalised counts
deseq_otu_table_round <- round(deseq_otu_table)

# save DESeq2-normalised OTU table
write.csv(deseq_otu_table_round, "deseq_otu_table_phyloseq.csv")



########THIS IS ANOTHER METHOD TO RUN DESEQ, you only need to use one and i have been using the first one above#####################
##########################################################
#### DESeq2 normalisation using only DESeq2 functions ####
##########################################################

# create DESeq2 object
count_deseq <- DESeqDataSetFromMatrix(countData = otu_table,
                                      colData = fact,
                                      design = ~ treatment)

# inspect DESeq object
count_deseq
attributes(count_deseq)

# positive counts normalisaton accounts for zero-inflation
count_deseq <- estimateSizeFactors(count_deseq, type = "poscounts")
sizeFactors(count_deseq)

# plot size factors
barplot(sizeFactors(count_deseq), las = 2, ylab = "Size factor normalisation", col = fact$treatment)

# calculation of size-factor corrected data
# size-factor corrected data are calculated by dividing the raw counts by the sample size factor and adding 0.5 to correct for the zeros
# so, zeros are just replaced with 0.5
count_deseq_norm <- sapply(row.names(count_deseq), function(x){
  plotCounts(count_deseq, x, "treatment", returnData = TRUE, normalized = TRUE)$count
})
count_deseq_norm <- data.frame(t(count_deseq_norm))
colnames(count_deseq_norm) <- colnames(count_deseq)

# note the difference between the normalised and raw counts
count_deseq_norm[1:5, 1:5]
otu_table[1:5, 1:5]

# remove the addition of 0.5 from all entries
count_deseq_norm <- count_deseq_norm - 0.5
count_deseq_norm[1:5, 1:5]

# round normalised counts
count_deseq_norm_round <- round(count_deseq_norm)

# check for zeros across samples - produced by normalisation
dim(count_deseq_norm_round[rowSums(count_deseq_norm_round) == 0, ])

# remove OTUs that have no counts (sum = 0) across all samples
count_deseq_norm_round <- count_deseq_norm_round[rowSums(count_deseq_norm_round) > 0, ]
dim(count_deseq_norm_round)

# check singletons produced by normalisation - keep these sequences
dim(count_deseq_norm_round[rowSums(count_deseq_norm_round) == 1, ])

# save DESeq2-normalised OTU table
write.csv(count_deseq_norm_round, "deseq_otu_table.csv")
