library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11566_ketogenic_diet/57715_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 2272 taxa and 24 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11566_ketogenic_diet/11566_20180927-085224.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11566
SUMMARY = "ketogenic_diet"
WRITE = TRUE
##########################################

show_plot(otu, df_sample[df_sample$kdday==14, ], "sample_name", "kdtreatment1", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdday==14, ], "sample_name", "kdtreatment1", "GMPR", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdday==14, ], "sample_name", "kdtreatment1", "GMPR2", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdday==14, ], "sample_name", "kdtreatment1", "CSS", "chemicals", TRUE, write = WRITE)

show_plot(otu, df_sample[df_sample$kdtreatment1=="CD",], "sample_name", "kdday", "TSS", "time", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdtreatment1=="CD",], "sample_name", "kdday", "GMPR", "time", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdtreatment1=="CD",], "sample_name", "kdday", "GMPR2", "time", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdtreatment1=="CD",], "sample_name", "kdday", "CSS", "time", TRUE, write = WRITE)

show_plot(otu, df_sample[df_sample$kdtreatment1=="KD",], "sample_name", "kdday", "TSS", "time", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdtreatment1=="KD",], "sample_name", "kdday", "GMPR", "time", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdtreatment1=="KD",], "sample_name", "kdday", "GMPR2", "time", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$kdtreatment1=="KD",], "sample_name", "kdday", "CSS", "time", TRUE, write = WRITE)


show_plot(otu, df_sample, "sample_name", "kdtreatment1", "TSS", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "kdtreatment1", "GMPR", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "kdtreatment1", "GMPR2", "chemicals", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "kdtreatment1", "CSS", "chemicals", TRUE, write = WRITE)
