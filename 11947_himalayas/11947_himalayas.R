library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11947_himalayas/65668_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 3663 rows and 78 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11947_himalayas/11947_20181001-151811.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 11947
SUMMARY = "himalayas"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "country", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "CSS", "demographics", FALSE, write = WRITE)


show_plot(otu, df_sample, "sample_name", "population", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "population", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "population", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample[df_sample$population != "Europeans",], "sample_name", "population", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$population != "Europeans",], "sample_name", "population", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$population != "Europeans",], "sample_name", "population", "CSS", "demographics", FALSE, write = WRITE)


show_plot(otu, df_sample, "sample_name", "smok", "TSS", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "smok", "GMPR", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "smok", "CSS", "habit", FALSE, write = WRITE)
