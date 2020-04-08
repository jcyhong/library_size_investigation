library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library("biomformat")
library(GMPR)
library(metagenomeSeq)
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")


path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/317_hand_surface_bacteria/45546_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 5621 taxa (rows) and 179 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/317_hand_surface_bacteria/317_20180101-113932.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)

# # Exploration with TSS
# lib_sizes <- Matrix::colSums(otu)
# summary(lib_sizes)
# sample_df <- data.frame(sample_name=as.factor(colnames(otu)), 
#                         lib_sizes=lib_sizes)
# sample_df <- sample_df %>% inner_join(df_sample, by='sample_name')
# sample_df %>% 
#   group_by(sample_type) %>% 
#   summarize(count=n(), mean=mean(lib_sizes), sd=sd(lib_sizes))


############ GLOBAL VARIABLE ##############
STUDY_ID = 317
SUMMARY = "hand bacteria"
WRITE = FALSE
##########################################

show_plot(otu, df_sample[df_sample$anatomical_body_site != "FMA:Palm",], "sample_name", "anatomical_body_site", "TSS", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$anatomical_body_site != "FMA:Palm",], "sample_name", "anatomical_body_site", "GMPR", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$anatomical_body_site != "FMA:Palm",], "sample_name", "anatomical_body_site", "CSS", "body_site", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "gender", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR", "gender", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "gender", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "dominant_hand", "TSS", "habit", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "dominant_hand", "GMPR", "habit", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "dominant_hand", "CSS", "habit", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "time_since_last_wash", "TSS", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "time_since_last_wash", "GMPR", "time", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "time_since_last_wash", "CSS", "time", TRUE, write = WRITE)

