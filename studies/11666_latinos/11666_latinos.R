library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11666_latinos/49921_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
dim(otu) # 5583 rows and 559 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11666_latinos/11666_20191210-120214.txt', header=TRUE, sep='\t')

well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)


############ GLOBAL VARIABLE ##############
STUDY_ID = 11666
SUMMARY = "latinos"
WRITE = TRUE
##########################################
show_plot(otu, df_sample, "sample_name", "hispanic_origin", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "hispanic_origin", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "hispanic_origin", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "diabetes_self_v2", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diabetes_self_v2", "GMPR2", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diabetes_self_v2", "CSS", "disease", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "us_born_v2", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "us_born_v2", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "us_born_v2", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "education_c3_v2", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "education_c3_v2", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "education_c3_v2", "CSS", "demographics", FALSE, write = WRITE)

