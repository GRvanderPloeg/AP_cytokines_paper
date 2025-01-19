library(tidyverse)

wd <- "."
setwd(wd)

path = "./Tax4Fun2_output/Athina_Tax4Fun2_rarefied_6k/functional_prediction.txt"
input_df = readRDS("Athina_Tax4Fun2_rarefied_annotated_6k.RDS") %>% as_tibble()

input_df_indexed = input_df %>% mutate(consistent_index=1:nrow(.))
output_df = read.csv(path, sep="\t", header=T) %>% as_tibble()

#microbiome.raw = read.csv("../0. Raw data input/20221005_wp2/count-table.tsv", sep="\t")
#taxa = read.csv("../0. Raw data input/20221005_wp2/taxonomic-classification.tsv", sep="\t")
#header = colnames(df)
#fixedHeader = header %>% str_split_fixed("X", 2)
#newHeader = c(header[1], fixedHeader[2:nrow(fixedHeader),2])
#colnames(df) = newHeader

ko2name = output_df[,c(1,ncol(output_df))]
colnames(ko2name) = c("KO", "description")
ko2name = ko2name %>% as_tibble()

output_df = output_df %>% select(-KO, -description)

dataFunc = output_df %>% as_tibble()

wrangle2 = function(df, path){
  dataOutput = t(df) %>% as_tibble()
  subjectOutput = input_df
  featureOutput = ko2name
  
  write.table(dataOutput, paste0(path, "_functional_analysis_raw.csv"), row.names=FALSE, col.names=FALSE)
  write.table(featureOutput, paste0(path, "_functional_analysis_raw_features.csv"), row.names=FALSE, col.names=FALSE)
  write.table(subjectOutput, paste0(path, "_functional_analysis_raw_subjects.csv"), row.names=FALSE, col.names=FALSE)
}

wrangle2(dataFunc, "./Tax4fun2_output/Athina")
# 
# ```{r wrangle into correct format}
# wrangle = function(df, samplingSite, path){
#   df2 = t(df) %>% as_tibble()
#   colnames(df2) = ko2name$KO
#   metadata = input_df %>% filter(niche==samplingSite)
#   df2 = df2 %>% mutate(subject=metadata$subject, visit=metadata$visit)
#   
#   result = df2 %>% pivot_longer(-c(subject, visit)) %>% mutate(newname=paste0(name,"_t",visit)) %>% select(-visit,-name) %>% pivot_wider(names_from=newname) %>% arrange(subject)
#   
#   t1 = paste0(ko2name$KO, "_t1")
#   t2 = paste0(ko2name$KO, "_t2")
#   t3 = paste0(ko2name$KO, "_t3")
#   t4 = paste0(ko2name$KO, "_t4")
#   t5 = paste0(ko2name$KO, "_t5")
#   t6 = paste0(ko2name$KO, "_t6")
#   t7 = paste0(ko2name$KO, "_t7")
#   sortedColumns = c(t1,t2,t3,t4,t5,t6,t7)
#   
#   dataOutput = result %>% select(sortedColumns)
#   featureOutput = ko2name
#   subjectOutput = result %>% select(subject) %>% left_join(rf)
#   
#   write.table(dataOutput, paste0(path, "_functional_analysis_raw.csv"), row.names=FALSE, col.names=FALSE)
#   write.table(featureOutput, paste0(path, "_functional_analysis_raw_features.csv"), row.names=FALSE, col.names=FALSE)
#   write.table(subjectOutput, paste0(path, "_functional_analysis_raw_subjects.csv"), row.names=FALSE, col.names=FALSE)
# }
# 
# wrangle(tongueFunc, "tongue", "./Tongue")
# wrangle(lowlingFunc, "lower jaw, lingual", "./Lowling")
# wrangle(lowinterFunc, "lower jaw, interproximal", "./Lowinter")
# wrangle(uplingFunc, "upper jaw, lingual", "./Upling")
# wrangle(upinterFunc, "upper jaw, interproximal", "./Upinter")
# wrangle(salivaFunc, "saliva", "./Saliva")
# ```