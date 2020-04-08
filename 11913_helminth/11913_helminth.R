library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11757_regional_variation/56281_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 13735 rows and 7009 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11757_regional_variation/regional_variation.txt', header=TRUE, sep='\t')
df_sample$name <- paste0('11757.', df_sample$host_subject_id)

############ GLOBAL VARIABLE ##############
STUDY_ID = 11757
SUMMARY = "regional variation"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "name", "districts", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "districts", "GMPR", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "name", "districts", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "name", "occupation", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "occupation", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "occupation", "CSS", "demographics", FALSE, write = WRITE)

# show_plot(otu, df_sample, "name", "gender", "TSS", "demographics", FALSE, write = WRITE)
# show_plot(otu, df_sample, "name", "gender", "GMPR", "demographics", FALSE, write = WRITE)
# show_plot(otu, df_sample, "name", "gender", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample, "name", "education", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "education", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "education", "CSS", "demographics", FALSE, write = WRITE)


show_plot(otu, df_sample, "name", "antibiotics", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "antibiotics", "GMPR", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "antibiotics", "CSS", "disease", FALSE, write = WRITE)

show_plot(otu, df_sample, "name", "diarrhea", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "diarrhea", "GMPR", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "name", "diarrhea", "CSS", "disease", FALSE, write = WRITE)

