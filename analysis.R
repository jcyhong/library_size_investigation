library(dplyr)
library(knitr)
library(kableExtra)

data = read.csv("~/Documents/UC Berkeley/spring2020/Researches/librarysize/DATA/compiled.txt", sep = "")

head(data) %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),full_width = F,position = "center")

# Combine GMPR and GMPR2
data$method2 = replace(data$method, data$method=="GMPR", "GMPR2")

data_int = data[data$interest==1,]

# Number of unique studies: 28
length(unique(data$study_id))

# Check all indices
for (i in unique(data$study_id)) {
  cat(i,"\n")
  print(which(data$study_id == i))
}

# How many distinct studies have GMPR2 only: 23
length(unique(data[data$method =="GMPR2",]$study_id))

# Percentage of significant p_val for methods, inyterest, category
data_int %>%
  group_by(method2) %>%
  summarise(num_group = n(), 
            pct_signif_p=paste0(round(100*sum(p_val<0.05)/n(),2),'%'), 
            pct_large_effsize=paste0(round(100*sum(effsize>0.14)/n(),2),'%')) %>% kable() %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),full_width = F) %>% row_spec(2,bold=T, color="red")

data %>%
  group_by(interest) %>%
  summarise(num_group = n()/3, per=paste0(round(100*sum(p_val<0.05)/n(),2),'%'))

# popular category: demographics, disease, geo_loc
data_int[!data_int$category %in% c("", "attribute"),] %>% filter(method2=="GMPR2") %>% 
  group_by(category) %>%
  summarise(num_group = n(), pct_signif_p=paste0(round(100*sum(p_val<0.05)/n(),2),'%'), pct_large_effsize=paste0(round(100*sum(effsize>=0.14)/n(),2),'%')) %>% kable() %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),full_width = F) %>% row_spec(c(1, 9),bold=T, color = "red")

data_int %>%
  group_by(category) %>%
  summarise(num_group = n()/3, per=paste0(round(100*sum(p_val<0.05)/n(),2),'%'))

data_int[data_int$method=="TSS",] %>%
  group_by(category) %>%
  summarise(num_group = n(), per=paste0(round(100*sum(p_val<0.05)/n(),2),'%'))

data_int[data_int$method2=="GMPR2",] %>%
  group_by(category) %>%
  summarise(num_group = n(), per=paste0(round(100*sum(p_val<0.05)/n(),2),'%'))


data_int[data_int$method=="CSS",] %>%
  group_by(category) %>%
  summarise(num_group = n(), per=paste0(round(100*sum(p_val<0.05)/n(),2),'%'))


# Percentage of large effect size for methods, inyterest, category
data %>%
  group_by(method2) %>%
  summarise(per=paste0(round(100*sum(effsize>0.5)/n(),2),'%'))





