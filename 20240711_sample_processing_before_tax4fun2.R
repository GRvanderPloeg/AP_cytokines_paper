wd <- "."
setwd(wd)

library(vegan)
library(Biostrings)
library(tidyverse)

counts <- read.delim("./Data/20240429_microbiome_counts.csv", sep=" ", header=FALSE)
taxonomy <- read.delim("./Data/20240429_taxonomy.csv", sep=" ", header=FALSE)
sampleMeta <- read.delim("./Data/20240429_microbiome_sampleMeta.csv", sep=" ", header=FALSE)
rare_threshold = 6500 # also consider 15k,  25k, 40k

seqs = readDNAStringSet("./Data/ACTA094_Root/ACTA094_Root/Root_zotus.fasta")

# 20240724 REORDER Taxonomy by ZOTU to make it the same as seqs (ordered by zOTU number)
taxonomy = taxonomy %>% as_tibble() %>% mutate(Zotu_number = as.numeric(str_split_fixed(taxonomy$V8, "Zotu", 2)[,2])) %>% arrange(Zotu_number)
taxonomy$representative_sequence = paste(seqs)
otus <- DNAStringSet(taxonomy$representative_sequence)
names(otus) <- taxonomy$V8
colnames(counts) = taxonomy$V8

# Diagnostics
gm <- apply(counts,1,function(x) exp(mean(log(x[x>0]))))
gm2 <- apply(counts,1,function(x) exp(mean(log(c(x[x>0],rep(0.5,length(which(x==0))))))))

hist(rowSums(counts))
summary(rowSums(counts))

# Boxplot of all sequencing depths
rowSums(counts) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_boxplot() + geom_vline(xintercept=rare_threshold,col="red")

# Rarification curves
rarecurve(counts, step=100, label=FALSE)

# IF we rarefy to n reads, what do we remove?
mask = which(rowSums(counts)>rare_threshold)
counts_filtered <- counts[mask,]
sampleMeta_filtered = sampleMeta[mask,]
#tab10 <- counts_control[which(rowSums(counts_control[,6:ncol(counts_control)])>rare_threshold),]
#sn10 <- table(tab10$subject,tab10$niche)
#rm10 <- cbind(rep(levels(as.factor(tab10$subject)),times=length(levels(as.factor(tab10$niche)))),
#              rep(levels(as.factor(tab10$niche)),each=length(levels(as.factor(tab10$subject)))))[which(!sn10 %in% c(0,7)),]

#counts_filtered <- tab10[apply(tab10[,c("subject","niche")],1,function(x) length(which(rm10[,1]==x[1]&rm10[,2]==x[2]))==0),]
counts_filtered <- counts_filtered[,which(colSums(counts_filtered)>0)]
taxonomy_filtered <- taxonomy[which(taxonomy$V8 %in% colnames(counts_filtered)),]

# Rarefication
set.seed(607)
counts_rarefied <- data.frame(rrarefy(counts_filtered,rare_threshold))
counts_rarefied <- counts_rarefied[,which(colSums(counts_rarefied)>0)]
taxonomy_rarefied <- taxonomy_filtered[which(taxonomy_filtered$V8 %in% colnames(counts_rarefied)),]

otu_rarefied <- DNAStringSet(taxonomy_rarefied$representative_sequence)
names(otu_rarefied) <- taxonomy_rarefied$V8
counts_rarefied_numeric <- data.frame(t(counts_rarefied[,5:ncol(counts_rarefied)]))
colnames(counts_rarefied_numeric) <- sampleMeta_filtered$V3
counts_rarefied_numeric$OTU <- rownames(counts_rarefied_numeric)

# Diagnostics
#gm.10 <- apply(tabr10[,6:ncol(tabr10)],1,function(x) exp(mean(log(c(x[x>0])))))
#gm2.10 <- apply(tabr10[,6:ncol(tabr10)],1,function(x) exp(mean(log(c(x[x>0],rep(0.5,length(which(x==0))))))))

counts_rarefied = cbind(sampleMeta_filtered, counts_rarefied)

# Save prepared input data for Tax4fun2
saveRDS(counts_rarefied,"Athina_Tax4Fun2_rarefied_annotated_6k.RDS")
saveRDS(counts_rarefied_numeric,"Athina_Tax4Fun2_rarefied_6k.RDS")
writeXStringSet(otu_rarefied,"Athina_Tax4Fun2_rarefied_6k_OTUs.fa")

