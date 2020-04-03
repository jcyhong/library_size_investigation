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
library(phyloseq)
library(RColorBrewer)
data("qa10394")

otu <- otu_table(qa10394)
dim(otu) # 9719 taxa x 1522 samples
taxa_are_rows(otu)
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]


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

sample_df <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/10394_fecal_storage/sample_data_10394.txt', header=TRUE, sep='\t')
sample_df$name <- paste0('10394.', sample_df$description)
# Some of the observations seem to be ill-defined?
well_defined <- !grepl('BLANK', sample_df$host_subject_id) & 
  !grepl('FTA', sample_df$host_subject_id) &
  !grepl('mistake', sample_df$host_subject_id)

df_sample <- sample_df %>% filter(well_defined)


# Adapted from https://github.com/lichen-lab/GMPR/blob/master/R/GMPR.R
GMPR<-function(OTUmatrix, min_ct=2, intersect_no=4){
  #datatype check
  if(!(class(OTUmatrix) %in% c("data.frame","matrix")))
    stop("Unknown datatype of object \"OTUmatrix\".")
  rownames(OTUmatrix)->SampleName
  if(ncol(OTUmatrix)<nrow(OTUmatrix))
    warning("Sample size is larger than OTU number. Samples should be rows.")
  if(length(OTUmatrix[OTUmatrix<1 & OTUmatrix>0])>=0.1*ncol(OTUmatrix)*nrow(OTUmatrix))
    stop("More than 10% values are fractional, please check.")
  apply(OTUmatrix, MARGIN = 2, as.integer)->OTUmatrix
  
  gmpr(OTUmatrix,min_ct,intersect_no)->size.factor
  names(size.factor)<-SampleName
  size.factor[abs(size.factor-1)<1e-10]<-NA
  size.factor
}

# Adapted from http://www.metagenomics.wiki/tools/16s/norm/css
CSS <- function(OTUmatrix) {
  metaSeqObj = newMRexperiment(OTUmatrix)
  metaSeqObj_css = cumNorm(metaSeqObj, p = cumNormStatFast(metaSeqObj))
  df_otu_css = data.frame(MRcounts(metaSeqObj_css, norm=TRUE, log=TRUE))
  Matrix::colSums(df_otu_css)
}

get_lib_sizes <- function(otu, method = "TSS") {
  if (method == "TSS") {
    lib_sizes <- Matrix::colSums(otu)
  } else if (method == "GMPR") {
    lib_sizes <- GMPR(t(otu))
  } else if (method == "CSS") {
    lib_sizes <- CSS(otu)
  } 
  lib_sizes
}


# Create a function to produce kruskal-wallis test result and ggplot
show_plot <- function(otu, sample_table, sample_id_col, col, method = "TSS", log = FALSE, write = FALSE) {
  lib_sizes <- get_lib_sizes(otu, method)
  sample_df <- data.frame(lib_sizes=lib_sizes)
  sample_df[sample_id_col]=as.factor(colnames(otu))
  sample_df <- sample_df %>% inner_join(sample_table, by=sample_id_col)
  
  sample_df <- sample_df[sample_df[col]!="",]
  
  krus_pval = kruskal.test(sample_df$lib_sizes ~ sample_df[,col])$p.value
  if (log) {
    ggplot(sample_df, aes(x=log(lib_sizes), color=factor(sample_df[,col]))) +
      geom_density() +
      labs(title = sprintf("Log library size(%s) distribution by %s", method, col),
           caption = sprintf("p value is: %E", krus_pval),
           color = col)
  } else {
    ggplot(sample_df, aes(x=lib_sizes, color=factor(sample_df[,col]))) +
      geom_density() +
      labs(title = sprintf("library size(%s) distribution by %s", method, col),
           caption = sprintf("p value is: %E", krus_pval),
           color = col)
  }
  if (write) {
    result <- data.frame(study_id = STUDY_ID, summary=SUMMARY, method = method, group = col, p_val = krus_pval)
    write.table(result, file = '~/Documents/UC Berkeley/spring2020/Researches/librarysize/compiled.txt', row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

############ GLOBAL VARIABLE ##############
STUDY_ID = 10394
SUMMARY = "fecal storage"
WRITE = FALSE
##########################################

show_plot(otu, df_sample, "name", "host_subject_id", "TSS", write = WRITE)
show_plot(otu, df_sample, "name", "host_subject_id", "GMPR", write = WRITE)
show_plot(otu, df_sample, "name", "host_subject_id", "CSS", write = WRITE)


show_plot(otu, df_sample, "name", "sample_storage_temp_treatment", "TSS", write = WRITE)
show_plot(otu, df_sample, "name", "sample_storage_temp_treatment", "GMPR", write = WRITE)
show_plot(otu, df_sample, "name", "sample_storage_temp_treatment", "CSS", write = WRITE)

show_plot(otu, df_sample, "name", "sample_preservation_method", "TSS", write = WRITE)
show_plot(otu, df_sample, "name", "sample_preservation_method", "GMPR", write = WRITE)
show_plot(otu, df_sample, "name", "sample_preservation_method", "CSS", write = WRITE)

show_plot(otu, df_sample, "name", "duration_of_storage", "TSS", write = WRITE)
show_plot(otu, df_sample, "name", "duration_of_storage", "GMPR", write = WRITE)
show_plot(otu, df_sample, "name", "duration_of_storage", "CSS", write = WRITE)


