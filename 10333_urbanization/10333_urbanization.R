library(MicrobeDS)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)
library("biomformat")

path_to_data <- "~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/10333_urbanization/44761_otu_table.biom"
otu <- otu_table(import_biom(path_to_data))
dim(otu) # 28423 rows and 754 samples (rows)
taxa_are_rows(otu) # taxa are rows (so observations are columns), 
otu = as.data.frame(otu)
otu = otu %>% select(which(colSums(otu)>2)) # https://rdrr.io/bioc/metagenomeSeq/src/R/cumNormStatFast.R

lib_sizes <- Matrix::colSums(otu)
summary(lib_sizes)

df_sample <- read.table('~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/10333_urbanization/10333_20190808-130957.txt', header=TRUE, sep='\t')

well_defined <- !grepl('blank', df_sample$sample_name)
df_sample <- df_sample %>% filter(well_defined)


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
STUDY_ID = 10333
SUMMARY = "urbanization"
WRITE = TRUE
##########################################

show_plot(otu, df_sample, "sample_name", "env_biome", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_biome", "CSS", write = WRITE)


show_plot(otu, df_sample, "sample_name", "env_material", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_material", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_material", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "env_package", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_package", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "env_package", "CSS", write = WRITE)

show_plot(otu, df_sample[df_sample$cleaning_frequency!="not provided",], "sample_name", "cleaning_frequency", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$cleaning_frequency!="not provided",], "sample_name", "cleaning_frequency", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$cleaning_frequency!="not provided",], "sample_name", "cleaning_frequency", "CSS", write = WRITE)


show_plot(otu, df_sample, "sample_name", "country", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "country", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "fans", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "fans", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "fans", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "washer", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "washer", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "washer", "CSS", write = WRITE)

show_plot(otu, df_sample[df_sample$cats!="not provided",], "sample_name", "cats", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$cats!="not provided",], "sample_name", "cats", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$cats!="not provided",], "sample_name", "cats", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "dogs", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "dogs", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "dogs", "CSS", write = WRITE)

show_plot(otu, df_sample[df_sample$ethnicity!="not applicable",], "sample_name", "ethnicity", "TSS", write = WRITE)
show_plot(otu, df_sample[df_sample$ethnicity!="not applicable",], "sample_name", "ethnicity", "GMPR", write = WRITE)
show_plot(otu, df_sample[df_sample$ethnicity!="not applicable",], "sample_name", "ethnicity", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "geo_loc_name", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "geo_loc_name", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "host_body_habitat", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_body_habitat", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_body_habitat", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "host_body_site", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_body_site", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "host_body_site", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "sample_type", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "sample_type", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "sex", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "sex", "CSS", write = WRITE)

show_plot(otu, df_sample, "sample_name", "wireless", "TSS", write = WRITE)
show_plot(otu, df_sample, "sample_name", "wireless", "GMPR", write = WRITE)
show_plot(otu, df_sample, "sample_name", "wireless", "CSS", write = WRITE)



