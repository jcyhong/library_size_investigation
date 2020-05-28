library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")


path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12080_immigration/63780_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 6052 rows and 647 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12080_immigration/12080_20181107-185650.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 12080
SUMMARY = "immigration"
WRITE = TRUE
##########################################
show_plot(otu, df_sample, "sample_name", "type_location_before_us", "TSS", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "type_location_before_us", "GMPR2", "geo_loc", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "type_location_before_us", "CSS", "geo_loc", TRUE,write = WRITE)

show_plot(otu, df_sample, "sample_name", "sample_group", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_group", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_group", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "religion", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "religion", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "religion", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "highest_education", "TSS","demographics", TRUE,  write = WRITE)
show_plot(otu, df_sample, "sample_name", "highest_education", "GMPR2", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "highest_education", "CSS", "demographics", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "ethnicity", "TSS", "demographics", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "ethnicity", "GMPR2", "demographics",  TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "ethnicity", "CSS", "demographics", TRUE,  write = WRITE)

show_plot(otu, df_sample, "sample_name", "alcohol_use", "TSS", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "alcohol_use", "GMPR2", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "alcohol_use", "CSS", "habit", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "bmi_class", "TSS", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "bmi_class", "GMPR2", "", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "bmi_class", "CSS", "", TRUE, write = WRITE)

