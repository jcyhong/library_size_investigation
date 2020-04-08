library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library("biomformat")
library(GMPR)
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")


path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/213_soil/47636_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 6099 taxa (rows) and 48 samples (cols)
# Since taxa are rows (so observations are columns), 
# we sum over the columns to get library sizes. 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/213_soil/213_20180101-114722.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)

# Exploration
lib_sizes <- Matrix::colSums(otu)
summary(lib_sizes)
sample_df <- data.frame(sample_name=colnames(otu), lib_sizes=lib_sizes,
                        stringsAsFactors=F)
sample_df <- sample_df %>% inner_join(df_sample, by='sample_name')
sample_df <- sample_df[sample_df$lib_sizes < 10000,]
sample_df %>% 
  group_by(sample_type) %>% 
  summarize(count=n(), mean=mean(lib_sizes), sd=sd(lib_sizes))



############ GLOBAL VARIABLE ##############
STUDY_ID = 213
SUMMARY = "soil"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "sample_name", "carbon_substrate", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "carbon_substrate", "GMPR", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "carbon_substrate", "CSS", "chemicals", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "env_biome", "TSS", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "GMPR", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "CSS", "habitat", TRUE, write = WRITE)

