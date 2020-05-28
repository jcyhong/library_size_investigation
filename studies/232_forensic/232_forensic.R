library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/232_forensic/46811_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
#otu = otu[, colSums(otu)>=1000]
dim(otu) # 2223 taxa (rows) and 115 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/232_forensic/232_20170409-171325.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)

############ GLOBAL VARIABLE ##############
STUDY_ID = 232
SUMMARY = "forensic"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "sample_type", "TSS", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "GMPR2", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "CSS", "", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", FALSE, write = WRITE)

