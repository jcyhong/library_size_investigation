library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library("biomformat")
library(GMPR)
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/232_forensic/46811_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 2223 taxa (rows) and 115 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/232_forensic/232_20170409-171325.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)

# Exploration
lib_sizes <- Matrix::colSums(otu)

sample_df <- data.frame(sample_name=colnames(otu), lib_sizes=lib_sizes,
                        stringsAsFactors=F)
sample_df <- sample_df %>% inner_join(df_sample, by='sample_name')
sample_df %>% 
  group_by(sample_type) %>% 
  summarize(count=n(), mean=mean(lib_sizes), sd=sd(lib_sizes))




############ GLOBAL VARIABLE ##############
STUDY_ID = 232
SUMMARY = "forensic"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "sample_name", "sample_type", "TSS", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "GMPR", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "CSS", "", TRUE, write = WRITE)

show_plot(otu, df_sample[df_sample!="Missing: Not provided"], "sample_name", "sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample!="Missing: Not provided"], "sample_name", "sex", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample!="Missing: Not provided"], "sample_name", "sex", "CSS", "demographics", FALSE, write = WRITE)

