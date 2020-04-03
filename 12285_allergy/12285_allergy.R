library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library(GMPR)
library(metagenomeSeq)
library("biomformat")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12285_allergy/71476_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 7375 rows and 65 samples (cols)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu[, which(lapply(1:ncol(otu), function(i) {length(which(otu[, i]>0))})>1)]

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/12285_allergy/12285_20190425-073701.txt', header=TRUE, sep='\t')
df_sample$sample_name <- as.character(df_sample$sample_name)

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
  sample_df <- sample_df[sample_df[col]!="not provided",]
  
  krus_pval = kruskal.test(sample_df$lib_sizes ~ sample_df[,col])$p.value
  if (write) {
    result <- data.frame(study_id = STUDY_ID, summary=SUMMARY, method = method, group = col, p_val = krus_pval)
    write.table(result, file = '~/Documents/UC Berkeley/spring2020/Researches/librarysize/compiled.txt', row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  if (log) {
    plot = ggplot(sample_df, aes(x=log(lib_sizes), color=factor(sample_df[,col]))) +
      geom_density() +
      labs(title = sprintf("Log library size(%s) distribution by %s", method, col),
           caption = sprintf("p value is: %E", krus_pval),
           color = col)
  } else {
    plot = ggplot(sample_df, aes(x=lib_sizes, color=factor(sample_df[,col]))) +
      geom_density() +
      labs(title = sprintf("library size(%s) distribution by %s", method, col),
           caption = sprintf("p value is: %E", krus_pval),
           color = col)
  }
  plot
}

############ GLOBAL VARIABLE ##############
STUDY_ID = 12285
SUMMARY = "allergy"
WRITE = FALSE
##########################################
show_plot(otu, df_sample, "sample_name", "airconditioning", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "airconditioning", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "airconditioning", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "carpets", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "carpets", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "carpets", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "dustmitecovers", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "dustmitecovers", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "dustmitecovers", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "fan", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "fan", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "fan", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "neighborhood", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "neighborhood", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "neighborhood", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "floweringplants", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "floweringplants", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "floweringplants", "CSS", write = WRITE)

