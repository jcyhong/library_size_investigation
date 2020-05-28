library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11815_bat/63793_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 20893 taxa and 1210 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11815_bat/11815_20181108-170610.txt', header=TRUE, sep='\t')
df_sample$country = sapply(strsplit(as.character(df_sample$geo_loc_name), ":"), `[[`, 1)
############ GLOBAL VARIABLE ##############
STUDY_ID = 11815
SUMMARY = "bat"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "sample_name", "geo_loc_name", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "GMPR2", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "country", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "GMPR", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "GMPR2", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "host_sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_sex", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_sex", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sample_type", "TSS", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "GMPR", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "CSS", "body_site", FALSE, write = WRITE)
