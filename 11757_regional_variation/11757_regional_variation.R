library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11757_regional_variation/56281_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 13735 rows and 7009 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]


df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/11757_regional_variation/regional_variation.txt', header=TRUE, sep='\t')
df_sample$name <- paste0('11757.', df_sample$host_subject_id)


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
STUDY_ID = 11757
SUMMARY = "regional variation"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "name", "districts", "TSS", write = WRITE)
show_plot(otu, df_sample, "name", "districts", "GMPR", write = WRITE)
show_plot(otu, df_sample, "name", "districts", "CSS", write = WRITE)

show_plot(otu, df_sample[df_sample$occupation!="NA",], "name", "occupation", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$occupation!="NA",], "name", "occupation", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$occupation!="NA",], "name", "occupation", "CSS", write = WRITE)

show_plot(otu, df_sample[df_sample$education!="NA",], "name", "education", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$education!="NA",], "name", "education", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$education!="NA",], "name", "education", "CSS", write = WRITE)


show_plot(otu, df_sample[df_sample$antibiotics!="NA",], "name", "antibiotics", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$antibiotics!="NA",], "name", "antibiotics", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$antibiotics!="NA",], "name", "antibiotics", "CSS", write = WRITE)

show_plot(otu, df_sample[df_sample$diarrhea!="NA",], "name", "diarrhea", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$diarrhea!="NA",], "name", "diarrhea", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$diarrhea!="NA",], "name", "diarrhea", "CSS", write = WRITE)

