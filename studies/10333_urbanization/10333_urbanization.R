library(MicrobeDS)
library(dplyr)
library(phyloseq)
library(metagenomeSeq)
library("biomformat")
source("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/functions.R")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/10333_urbanization/44761_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]
otu = otu[, colSums(otu)>=1000]
dim(otu) # 28423 rows and 750 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/10333_urbanization/10333_20190808-130957.txt', header=TRUE, sep='\t')

well_defined <- !grepl('blank', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)

############ GLOBAL VARIABLE ##############
STUDY_ID = 10333
SUMMARY = "urbanization"
WRITE = F
##########################################

show_plot(otu, df_sample, "sample_name", "env_biome", "TSS", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "GMPR", "habitat", TRUE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "CSS", "habitat", TRUE, write = WRITE)

# 
# show_plot(otu, df_sample, "sample_name", "env_material", "TSS", "", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "env_material", "GMPR", "", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "env_material", "CSS", "", FALSE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "env_package", "TSS", "", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "env_package", "GMPR", "", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "env_package", "CSS", "", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "cleaning_frequency", "TSS", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "cleaning_frequency", "GMPR", "habit", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "cleaning_frequency", "CSS", "habit", FALSE, write = WRITE)


show_plot(otu, df_sample, "sample_name", "country", "TSS", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "GMPR", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "CSS", "geo_loc", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "ethnicity", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "ethnicity", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "ethnicity", "CSS", "demographics", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "geo_loc_name", "TSS", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "GMPR", "geo_loc", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "CSS", "geo_loc", FALSE, write = WRITE)

# show_plot(otu, df_sample, "sample_name", "host_body_habitat", "TSS", "body_site", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "host_body_habitat", "GMPR", "body_site", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "host_body_habitat", "CSS", "body_site", FALSE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "host_body_site", "TSS", "body_site", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "host_body_site", "GMPR", "body_site", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "host_body_site", "CSS", "body_site", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sample_type", "TSS", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "GMPR", "body_site", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "CSS", "body_site", FALSE, write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR", "demographics", FALSE, write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", "demographics", FALSE, write = WRITE)

# 
# show_plot(otu, df_sample, "sample_name", "fans", "TSS", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "fans", "GMPR", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "fans", "CSS", "habit", FALSE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "washer", "TSS", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "washer", "GMPR", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "washer", "CSS", "habit", FALSE,  write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "cats", "TSS", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "cats", "GMPR", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "cats", "CSS", "habit", FALSE, write = WRITE)
# 
# show_plot(otu, df_sample, "sample_name", "dogs", "TSS", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "dogs", "GMPR", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "dogs", "CSS", "habit", FALSE, write = WRITE)

# show_plot(otu, df_sample, "sample_name", "wireless", "TSS", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "wireless", "GMPR", "habit", FALSE, write = WRITE)
# show_plot(otu, df_sample, "sample_name", "wireless", "CSS", "habit", FALSE, write = WRITE)



