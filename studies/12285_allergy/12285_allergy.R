library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")


path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12285_allergy/71476_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 7375 rows and 65 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12285_allergy/12285_20190425-073701.txt', header=TRUE, sep='\t')
df_sample$sample_name <- as.character(df_sample$sample_name)

############ GLOBAL VARIABLE ##############
STUDY_ID = 12285
SUMMARY = "allergy"
WRITE = TRUE
##########################################
# show_plot(otu, df_sample, "sample_name", "airconditioning", "TSS", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "airconditioning", "GMPR", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "airconditioning", "CSS", FALSE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "carpets", "TSS", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "carpets", "GMPR", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "carpets", "CSS", FALSE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "dustmitecovers", "TSS", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "dustmitecovers", "GMPR", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "dustmitecovers", "CSS", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "fan", "TSS", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "fan", "GMPR2", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "fan", "CSS", "habit", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "neighborhood", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "neighborhood", "GMPR2", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "neighborhood", "CSS", "geo_loc", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "floweringplants", "TSS", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "floweringplants", "GMPR2", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "floweringplants", "CSS", "habit", FALSE, write = WRITE)

