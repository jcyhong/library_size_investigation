# Library size distributions for grouped microbiome data

## Abstract
In microbiome data analysis, uneven library sizes among observations present great challenges for proper statistical inference. The library size of an observation is the total number of individuals present in the observation. In practice, it is not uncommon for library sizes to vary across several orders of magnitude. On the other hand, observations in microbiome data are often grouped (for example, fecal samples from patients in different treatment groups). The goal of the project is to investigate whether library sizes are typically dependent on the group membership in grouped microbiome data. 

## Descriptions of folders and files:
- **compiled.txt**:
	- Meta-data frame of results with header: study_id, summary, method, group, p_val, effsize, category, interest
- **analysis.R**:
	- contains analysis of the metadata.
- **studies**:
	- Contains 28 folders of studies analyzed in this project. Each folder is named in the format "QIITAid_DescriptionOfTheStudy"
	- Each folder contains 3 files: a txt file for sample information table, an R file with scripts to calculate results and generate plots, and a .biom file for accessing the OTU table.
- **functions.R**:
	- Contains all the important functions used in the analysis of individual studies.
- **presentation.pdf**:
	- Slides for the project presentation.