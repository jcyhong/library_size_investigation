library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library("biomformat")
library(GMPR)
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/349_atherosclerosis/44665_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 4623 taxa (rows) and 73 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/349_atherosclerosis/349_20180101-113736.txt', header=TRUE, sep='\t')
well_defined <- !grepl('BLANK', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)

# Exploration with TSS
lib_sizes <- Matrix::colSums(otu)
summary(lib_sizes)
sample_df <- data.frame(sample_name=as.factor(colnames(otu)), 
                        lib_sizes=lib_sizes)
sample_df <- sample_df %>% inner_join(df_sample, by='sample_name')

sample_df %>% 
  group_by(sample_type) %>% 
  summarize(count=n(), mean=mean(lib_sizes), sd=sd(lib_sizes))


############ GLOBAL VARIABLE ##############
STUDY_ID = 349
SUMMARY = "atherosclerosis"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "sample_name", "atherosclerosis", "TSS", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "atherosclerosis", "GMPR", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "atherosclerosis", "CSS", "disease", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sample_type", "TSS", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "GMPR", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "CSS", "body_site", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "chronic_condition", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "chronic_condition", "GMPR", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "chronic_condition", "CSS", "disease", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "diabetes", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diabetes", "GMPR", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diabetes", "CSS", "disease", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "hypertension", "TSS", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "hypertension", "GMPR", "disease", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "hypertension", "CSS", "disease", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "gender", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR", "gender", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "gender", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "smoker", "TSS", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "smoker", "GMPR", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "smoker", "CSS", "habit", FALSE, write = WRITE)





