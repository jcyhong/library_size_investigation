library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")

source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")
load("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/soil/soil.RData")

otu <- otu_table(phyloObject)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
dim(otu) # 11352 taxa x 92 samples
taxa_are_rows(otu)
otu = as.data.frame(otu)

df_sample <- as.data.frame(as.matrix(sample_data(phyloObject)))


############ GLOBAL VARIABLE ##############
STUDY_ID = NA
SUMMARY = "soil from Ulas"
WRITE = F
##########################################

show_plot(otu, df_sample, "id", "aggsize", "TSS", "attribute", TRUE, write = WRITE)
show_plot(otu, df_sample, "id", "aggsize", "GMPR2", "attribute", TRUE, write = WRITE)
show_plot(otu, df_sample, "id", "aggsize", "CSS", "attribute", TRUE, write = WRITE)


