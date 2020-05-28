library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")

source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

data("RISK_CCFA")

otu <- otu_table(RISK_CCFA)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 9511 taxa x 1359 samples
taxa_are_rows(otu)
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/1939_CD/CD_sample_info.txt', header=TRUE, sep='\t')

############ GLOBAL VARIABLE ##############
STUDY_ID = 1939
SUMMARY = "CD"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "type_sample", "TSS", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "type_sample", "GMPR", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "type_sample", "CSS", "body_site", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "body_site", "TSS", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "body_site", "GMPR", "body_site", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "body_site", "CSS", "body_site", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "antibiotics", "TSS", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "antibiotics", "GMPR", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "antibiotics", "CSS", "disease", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "diagnosis", "TSS", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diagnosis", "GMPR", "disease", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "diagnosis", "CSS", "disease", TRUE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "race", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "race", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "race", "CSS", "demographics", FALSE, write = WRITE)
