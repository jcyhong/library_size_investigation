library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

otu <- otu_table(qa10394)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 9719 taxa x 1499 samples
taxa_are_rows(otu)
otu = as.data.frame(otu)


sample_df <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/10394_fecal_storage/sample_data_10394.txt', header=TRUE, sep='\t')
sample_df$name <- paste0('10394.', sample_df$description)
# Some of the observations seem to be ill-defined?
well_defined <- !grepl('BLANK', sample_df$host_subject_id) & 
  !grepl('FTA', sample_df$host_subject_id) &
  !grepl('mistake', sample_df$host_subject_id)

df_sample <- sample_df %>% filter(well_defined)

############ GLOBAL VARIABLE ##############
STUDY_ID = 10394
SUMMARY = "fecal storage"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "name", "sample_storage_temp_treatment", "TSS", "temp", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "sample_storage_temp_treatment", "GMPR", "temp", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "sample_storage_temp_treatment", "CSS", "temp", TRUE, write = WRITE)

show_plot(otu, df_sample, "name", "sample_preservation_method", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "sample_preservation_method", "GMPR", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "sample_preservation_method", "CSS", "chemicals", TRUE, write = WRITE)

show_plot(otu, df_sample, "name", "duration_of_storage", "TSS", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "duration_of_storage", "GMPR", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "duration_of_storage", "CSS", "time", TRUE, write = WRITE)


