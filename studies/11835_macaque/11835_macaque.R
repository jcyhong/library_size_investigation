library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11835_macaque/59993_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 5182 taxa and 30 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11835_macaque/11835_20181024-101425.txt', header=TRUE, sep='\t')
df_sample$sample_name = as.factor(df_sample$sample_name)

############ GLOBAL VARIABLE ##############
STUDY_ID = 11835
SUMMARY = "macaque"
WRITE = TRUE
##########################################

show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "env_biome", "TSS", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "env_biome", "GMPR2", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "env_biome", "CSS", "habitat", TRUE, write = WRITE)

show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "life_stage", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "life_stage", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "life_stage", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "sample_source", "TSS", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "sample_source", "GMPR2", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "sample_source", "CSS", "geo_loc", FALSE, write = WRITE)

show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "sample_source", "TSS", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "sample_source", "GMPR2", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "sample_source", "CSS", "geo_loc", FALSE, write = WRITE)

show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "sex", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name != "human", ], "sample_name", "sex", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "sex", "GMPR2", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample[df_sample$host_common_name == "human", ], "sample_name", "sex", "CSS", "demographics", FALSE, write = WRITE)


