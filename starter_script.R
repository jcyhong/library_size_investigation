# Starter script for processed microbiome data
#
# Remark: This is merely a starter script. Feel free to 
# do as much exploration as you would like; wite your own
# code, etc.
#
# The repository https://github.com/twbattaglia/MicrobeDS
# contains several microbiome datasets in the form of
# a phyloseq object.
#
# To install the package, run
# devtools::install_github("twbattaglia/MicrobeDS")

library(MicrobeDS)
library(dplyr)
library(ggplot2)
data("qa10394")

otu <- otu_table(qa10394)
dim(otu)
taxa_are_rows(otu)

# Be careful! Sometimes the OTU table has taxa as rows;
# sometimes it has taxa as columns. 
# For this data set, we have taxa as rows

# To identify the grouping (or in general, the features) of the observations,
# we need the sample table.
# The sample table can be downloaded from Qiita:
# https://qiita.ucsd.edu/study/description/10394
# Click Sample Information -> Sample Information to download the sample data.
# For simplicity, I had uploaded the sample data to the repo.

# Be careful! Sometimes the order of the observations in sample data is 
# different from the order in the OTU tables.
# (i.e.: the first observation in the sample data might not correspond to 
# the first observation in OTU table)

sample_data_10394 <- read.table('sample_data_10394.txt', header=TRUE, sep='\t')
sample_data_10394$name <- paste0('10394.', sample_data_10394$description)
sample_df <- data.frame(name=colnames(otu_table(qa10394)), 
                        stringsAsFactors=F)
sample_df <- sample_df %>% left_join(sample_data_10394, by='name')
# Some of the observations seem to be ill-defined?
well_defined <- !grepl('BLANK', sample_df$host_subject_id) & 
  !grepl('FTA', sample_df$host_subject_id) &
  !grepl('mistake', sample_df$host_subject_id)

# Suppose we group the observations by subjects.
# There are 15 subjects: 10 humans and 5 dogs.
unique(sample_df[well_defined, ]$host_subject_id)

# Since taxa are rows (so observations are columns), 
# we sum over the columns to get library sizes. 
lib_sizes <- colSums(otu)
summary(lib_sizes)
# We can see that library sizes vary a lot.

sample_df_well_defined <- cbind(sample_df, lib_sizes) %>% filter(well_defined)
sample_df_well_defined %>%
  group_by(host_subject_id) %>%
  summarize(mean=mean(lib_sizes),
            sd=sd(lib_sizes))
lib_sizes_aov <- aov(lib_sizes ~ host_subject_id, 
                     data = sample_df_well_defined)
summary(lib_sizes_aov)

# We reject the null that all groups have the same mean library size.
# However, normality assumptions might not be reasonable.
# We proceed with a non-parametric test.

kruskal.test(lib_sizes ~ host_subject_id, 
             data = sample_df_well_defined)
# We reject the null that all groups have library sizes distribution.
# We might want to do some visualization instead, since our 
# project is exploratory.
# For example, we can plot the distribution of library sizes by groups.

# Original scale
ggplot(sample_df_well_defined, aes(x=lib_sizes, color=host_subject_id)) +
  geom_density() +
  scale_color_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(15))

# Log scale
ggplot(sample_df_well_defined, aes(x=log(lib_sizes), color=host_subject_id)) +
  geom_density() +
  scale_color_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(15))

# The library size distributions do not seem too different when grouped by subjects?
# What if we instead group by duration_of_storage AND/OR
# sample_storage_temp_treatment AND/OR 
# sample_preservation_method?
