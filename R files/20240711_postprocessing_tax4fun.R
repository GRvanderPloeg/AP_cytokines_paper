# postprocessing of Tax4Fun results

# Libraries
library(tidyverse)

# Load input microbiome count data
microbiome.raw = read.delim("./Data/20240429_microbiome_counts.csv", sep=" ", header=FALSE)
taxa = read.delim("./Data/20240429_taxonomy.csv", sep=" ", header=FALSE)

# Load Tax4Fun2 input dataset
input_df = readRDS("./Athina_Tax4Fun2_rarefied_annotated_6k.RDS")
input_df_numeric = readRDS("./Athina_Tax4Fun2_rarefied_6k.RDS")

# Load log file with FTUs and FSUs
offset = 7
raw_data = read.csv("./20240724_Tax4Fun2_output/Athina_Tax4Fun2_rarefied_6k/logfile2.txt", sep=" ", header=FALSE, skip=3)
FTUs = raw_data[1:nrow(input_df),] %>% as_tibble()
FSUs = raw_data[(nrow(input_df)+offset):nrow(raw_data),] %>% as_tibble()

# Make plots
FTUs %>% mutate(FTU=as.numeric(V2)) %>% ggplot(aes(x=FTU)) + geom_histogram(bins=20, col="black") + xlab("Fraction of unused taxonomic units")

FSUs %>%
  mutate(FSU=as.numeric(V2)) %>%
  ggplot(aes(x=1,y=FSU)) + 
  geom_violin() +
  geom_boxplot(width=0.15) +
  xlab("Sampling location") +
  ylab("Fraction of sequences unused (FSU)") +
  ggtitle("FSUs - rarefied to 6k")