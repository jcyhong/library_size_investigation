library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")

source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/353_storage/47627_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
dim(otu) # 8156 taxa (rows) and 144 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/353_storage/353_20180101-114721.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)


############ GLOBAL VARIABLE ##############
STUDY_ID = 353
SUMMARY = "storage"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "body_habitat", "TSS", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "body_habitat", "GMPR2", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "body_habitat", "CSS", "body_site", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "env_feature", "TSS", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_feature", "GMPR2", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_feature", "CSS", "", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "storage_days", "TSS","time",TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "storage_days", "GMPR2", "time",TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "storage_days", "CSS", "time",TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "storage_temp", "TSS", "temp", TRUE,write = WRITE)
show_plot(otu, df_sample, "sample_name", "storage_temp", "GMPR2", "temp", TRUE,write = WRITE)
show_plot(otu, df_sample, "sample_name", "storage_temp", "CSS", "temp", TRUE,write = WRITE)

