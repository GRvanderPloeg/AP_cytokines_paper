# exploration
library(tidyverse)
library(ggpubr)

df = read.csv("./Data/export_txt_6_10_20_results_combi_merged.txt", sep="\t")
raw_df = read.csv("./Data/Working data_280620.csv")
colnames(raw_df) = c("plate", "sample", "sample_id", "subject_id", "gender", "age", "DMFT", "duplicate", "visit", "case_control",
                     "VEGF", "VEGF_mean", "VEGF_cv", "CRP", "CRP_mean", "CRP_cv", "GM_CSF", "GM_CSF_mean", "GM_CSF_cv", "IL1alpha",
                     "IL1alpha_mean", "IL1alpha_cv", "IL1beta", "IL1beta_mean", "IL1beta_cv", "IL4", "IL4_mean", "IL4_cv", "IL6",
                     "IL6_mean", "IL6_cv", "IL8", "IL8_mean", "IL8_cv", "IL10", "IL10_mean", "IL10_cv", "IL12p70", "IL12p70_mean",
                     "IL12p70_cv", "IL17A", "IL17A_mean", "IL17A_cv", "IFN_gamma", "IFN_gamma_mean", "IFN_gamma_cv", "MIP1alpha",
                     "MIP1alpha_mean", "MIP1alpha_cv", "OPG", "OPG_mean", "OPG_cv", "TNFalpha", "TNFalpha_mean", "TNFalpha_cv")

metadata = raw_df %>% select(subject_id, gender, age, DMFT, visit, case_control) %>% unique() %>% as_tibble()
#raw_df %>% as_tibble() %>% mutate(id=paste0(subject_id,visit)) %>% ggplot(aes(x=as.factor(duplicate),y=OPG) ) + geom_boxplot() + geom_point() + geom_line(aes(group=id))

# Mean is okay for: VEGF, CRP, IL4, IL6, IL8, IL10, IL12p70, IFN_gamma, MIP1alpha, TNFalpha, OPG (only slightly problematic)
# Only replicate 1 is okay for: IL1beta
# Only replicate 2 is okay for: GM_CSF, IL17A
# Problematic: IL1alpha

# Magic line that fixes everything
new_df = raw_df %>% select(subject_id, visit, duplicate, VEGF, CRP, GM_CSF, IL1alpha, IL1beta, IL4, IL6, IL8, IL10, IL12p70, IL17A, IFN_gamma, MIP1alpha, OPG, TNFalpha) %>% pivot_longer(-c(subject_id,visit,duplicate)) %>% mutate(newname=paste0(name,"_",duplicate)) %>% select(-duplicate,-name) %>% pivot_wider(names_from=newname) %>% mutate(VEGF_final = rowMeans(select(., starts_with("VEGF_")), na.rm = TRUE), CRP_final = rowMeans(select(., starts_with("CRP_")), na.rm = TRUE), GM_CSF_final = GM_CSF_2, IL1alpha_final = rowMeans(select(., starts_with("IL1alpha_")), na.rm = TRUE), IL1beta_final = IL1beta_1, IL4_final = rowMeans(select(., starts_with("IL4_")), na.rm = TRUE), IL6_final = rowMeans(select(., starts_with("IL6_")), na.rm = TRUE), IL8_final = rowMeans(select(., starts_with("IL8_")), na.rm = TRUE), IL10_final = rowMeans(select(., starts_with("IL10_")), na.rm = TRUE), IL12p70_final = rowMeans(select(., starts_with("IL12p70_")), na.rm = TRUE), IL17A_final = IL17A_2, IFN_gamma_final = rowMeans(select(., starts_with("IFN_gamma_")), na.rm = TRUE), MIP1alpha_final = rowMeans(select(., starts_with("MIP1alpha_")), na.rm = TRUE), OPG_final = rowMeans(select(., starts_with("OPG_")), na.rm = TRUE), TNFalpha_final = rowMeans(select(., starts_with("TNFalpha_")), na.rm = TRUE))
new_df = new_df %>% select(subject_id, visit, VEGF_final, CRP_final, GM_CSF_final, IL1alpha_final, IL1beta_final, IL4_final, IL6_final, IL8_final, IL10_final, IL12p70_final, IL17A_final, IFN_gamma_final, MIP1alpha_final, OPG_final, TNFalpha_final)
new_df = new_df %>% left_join(metadata)
new_df$visit = as.numeric(as.factor(new_df$visit))
new_df = new_df %>% mutate(SubjectID = subject_id, Visit = visit) %>% left_join(df %>% select(SubjectID, Visit, RANKL_Conc_Mean)) %>% select(-SubjectID,-Visit)

new_df_numeric = new_df %>% select(-subject_id, -visit, -gender, -age, -DMFT, -case_control)
new_df_metadata = new_df %>% select(subject_id, visit, gender, age, DMFT, case_control)

a = df %>% ggplot(aes(x=VGEF_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
b=df %>% ggplot(aes(x=CRP_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
c=df %>% ggplot(aes(x=GM_CSF_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
d=df %>% ggplot(aes(x=IL_1alpha_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
e=df %>% ggplot(aes(x=IL1beta_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
f=df %>% ggplot(aes(x=IL4_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
g=df %>% ggplot(aes(x=IL6_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
h=df %>% ggplot(aes(x=IL8_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
i=df %>% ggplot(aes(x=IL10_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
j=df %>% ggplot(aes(x=IL_12p70_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
k=df %>% ggplot(aes(x=IL_17A_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
l=df %>% ggplot(aes(x=IFN_gama_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
m=df %>% ggplot(aes(x=MIP1_alpha_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
n=df %>% ggplot(aes(x=OPG_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
o=df %>% ggplot(aes(x=TNF_alpha_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
p=df %>% ggplot(aes(x=RANKL_Conc_Mean, fill=as.factor(case_control))) + geom_histogram()
ggarrange(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p, common.legend=TRUE)

df.meta = df[,1:13] %>% as_tibble()
df.numeric = df[,14:ncol(df)] %>% as_tibble()

# Autoscale for PCA
df.numeric.norm = apply(df.numeric, 2, function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)})

# Remove missing sample
df.numeric.norm = df.numeric.norm[-147,]

# PCA
model = prcomp(df.numeric.norm)
model$rotation %>% as_tibble() %>% ggplot(aes(x=PC1,y=PC2,label=colnames(df.numeric.norm))) + geom_point() + geom_text()
biplot(model, choices=1:2)