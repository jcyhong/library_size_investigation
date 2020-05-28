library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")

source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")


path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/213_soil/47636_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
dim(otu) # 6099 taxa (rows) and 48 samples (cols)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/213_soil/213_20180101-114722.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)


############ GLOBAL VARIABLE ##############
STUDY_ID = 213
SUMMARY = "soil"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "carbon_substrate", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "carbon_substrate", "GMPR2", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "carbon_substrate", "CSS", "chemicals", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "env_biome", "TSS", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "GMPR2", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "CSS", "habitat", TRUE, write = WRITE)

