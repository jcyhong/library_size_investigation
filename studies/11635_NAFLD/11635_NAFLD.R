library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11635_NAFLD/53596_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
dim(otu) # 5353 taxa and 220 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11635_NAFLD/11635_20190501-133407.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11635
SUMMARY = "NAFLD"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "adv_fibrosis", "TSS", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "adv_fibrosis", "GMPR2", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "adv_fibrosis", "CSS", "disease", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "diabetes_t2", "TSS", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diabetes_t2", "GMPR2", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diabetes_t2", "CSS", "disease", TRUE, write = WRITE)


show_plot(otu, df_sample, "sample_name", "cir_proband", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "cir_proband", "GMPR2", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "cir_proband", "CSS", "chemicals", TRUE, write = WRITE)


show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", TRUE, write = WRITE)
