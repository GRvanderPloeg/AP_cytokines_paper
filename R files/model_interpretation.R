# parafac work
library(tidyverse)
library(ggpubr)
library(paramGUI)
library(pracma)
library(ggrepel)
library(parafac4microbiome)
# 
# importPARAFAC = function(path, featureNames, dayVector, featureColumnOfInterest){
#   #  Output:
#   #    [[1]]: Loadings in ID mode
#   #    [[2]]: Loadings in feature mode
#   #    [[3]]: Loadings in time mode
#   #    [[4]]: Data as modeled, wide format
#   #    [[5]]: Input data, wide format
#   #    [[6]]: Data as modeled, long format
#   #    [[7]]: Input data, long format
#   #    [[8]]: Data as modeled, component 1, wide format
#   #    [[9]]: Data as modeled, component 1, long format
#   #   [[**]]: and so on until all components have a list item.
#   
#   id_mode = read.csv(paste0(path,"_individual_mode.csv"), header=FALSE)
#   feature_mode = read.csv(paste0(path,"_feature_mode.csv"), header=FALSE)
#   time_mode = read.csv(paste0(path,"_time_mode.csv"), header=FALSE)
#   numComponents = ncol(time_mode)
#   
#   if(dim(id_mode)[2] == (numComponents+1)){
#     colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject")
#     id_mode = id_mode %>% left_join(rf)
#   }
#   colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject", "RFgroup")
#   id_mode = as_tibble(id_mode)
#   
#   colnames(feature_mode) = c(paste0("Component_", 1:numComponents), featureNames)
#   feature_mode = as_tibble(feature_mode)
#   
#   colnames(time_mode) = c(paste0("Component_", 1:numComponents))
#   time_mode = as_tibble(time_mode) %>% mutate(days=dayVector)
#   
#   rawmodel = scan(paste0(path,"_model.csv"), sep=",")
#   model_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
#   
#   for(i in 1:nrow(id_mode)){
#     model_wide[i,] = rawmodel[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
#   }
#   column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
#   
#   for(i in 2:length(dayVector)){
#     column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
#   }
#   colnames(model_wide) = column_names
#   
#   model_wide = as_tibble(model_wide)
#   
#   rawinput = scan(paste0(path,"_input.csv"), sep=",")
#   input_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
#   
#   for(i in 1:nrow(id_mode)){
#     input_wide[i,] = rawinput[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
#   }
#   column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
#   
#   for(i in 2:length(dayVector)){
#     column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
#   }
#   colnames(input_wide) = column_names
#   
#   input_wide = as_tibble(input_wide)
#   model_long = make_longer(model_wide, id_mode$subject)
#   input_long = make_longer(input_wide, id_mode$subject)
#   
#   result = list(id_mode, feature_mode, time_mode, model_wide, input_wide, model_long, input_long)
#   listIterator = 8
#   
#   for(i in 1:numComponents){
#     modeled_component_raw = scan(paste0(path,"_component_", i, ".csv"), sep=",")
#     modeled_component_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
#     
#     for(i in 1:nrow(id_mode)){
#       modeled_component_wide[i,] = modeled_component_raw[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
#     }
#     column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
#     
#     for(i in 2:length(dayVector)){
#       column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
#     }
#     colnames(modeled_component_wide) = column_names
#     
#     modeled_component_wide = as_tibble(modeled_component_wide)
#     modeled_component_long = make_longer(modeled_component_wide, id_mode$subject)
#     
#     result[[listIterator]] = modeled_component_wide
#     listIterator = listIterator + 1
#     result[[listIterator]] = modeled_component_long
#     listIterator = listIterator + 1
#   }
#   
#   return(result)
# }
# 
# make_longer = function(data, subjects){
#   df = data %>% mutate(subject=subjects) %>% pivot_longer(-subject, names_to=c("asv", "visit"), names_sep="_t") %>% as.data.frame()
#   df$visit = as.integer(df$visit)
#   df$value = as.double(df$value)
#   df = df %>% as_tibble() %>% pivot_wider(id_cols=c("subject","visit"), names_from=asv, values_from=value) %>% select(-subject,-visit)
#   return(df)
# }
# 
# correctPARAFACloadings = function(model, modeToCorrect){
#   A = model[[1]] %>% select(starts_with("Component")) %>% as.matrix()
#   B = model[[2]] %>% select(starts_with("Component")) %>% as.matrix()
#   C = model[[3]] %>% select(starts_with("Component")) %>% as.matrix()
#   
#   if(modeToCorrect == 1){
#     F = kroneckercol(C, B) %>% as.matrix()
#     Ftilde = gramSchmidt(F)$Q
#     T = pinv(F) %*% Ftilde
#     Atilde = A %*% T
#     result = Atilde
#   }
#   else if(modeToCorrect == 2){
#     F = kroneckercol(A, C) %>% as.matrix()
#     Ftilde = gramSchmidt(F)$Q
#     T = pinv(F) %*% Ftilde
#     Btilde = B %*% T
#     result = Btilde
#   }
#   else if(modeToCorrect == 3){
#     F = kroneckercol(B, A) %>% as.matrix()
#     Ftilde = gramSchmidt(F)$Q
#     T = pinv(F) %*% Ftilde
#     Ctilde = C %*% T
#     result = Ctilde
#   }
#   
#   return(result)
# }
# 
# biplotAcrossTime = function(model, filteredFeatures, titleString){
#   A = model[[1]] %>% select(starts_with("Component")) %>% as.matrix()
#   B = model[[2]] %>% select(starts_with("Component")) %>% as.matrix()
#   C = model[[3]] %>% select(starts_with("Component")) %>% as.matrix()
#   
#   F = kroneckercol(A,C)
#   Btilde = gramSchmidt(B)$Q
#   T = pinv(B) %*% Btilde
#   Ftilde = F %*% T
#   
#   numComponents = ncol(model[[3]]) - 1
#   
#   if(numComponents == 2){
#     enlargingFactor = mean(apply(Ftilde, 2, max) / apply(Btilde, 2, max))
#     Btilde = Btilde * enlargingFactor
#     
#     Ftilde = Ftilde %>% as_tibble() %>% mutate(subject = rep(model[[1]]$subject,each=nrow(model[[3]])), timepoint = rep(model[[3]]$days, nrow(model[[1]])))
#     colnames(Ftilde) = c("Component_1", "Component_2", "subject", "timepoint")
#     featureLabels = model[[2]] %>% left_join(filteredFeatures) %>% filter(keep==TRUE)
#     featureLabels = featureLabels %>% select(featureNames) %>% pull()
#     Btilde = Btilde %>% as_tibble() %>% mutate(featureNames = model[[2]]$featureNames) %>% left_join(filteredFeatures) %>% filter(keep==TRUE) %>% select(V1,V2)
#     colnames(Btilde) = c("Component_1", "Component_2")
#     Ftilde = Ftilde %>% left_join(model[[1]] %>% select(subject,RFgroup))
#     
#     ggplot() + 
#       geom_path(data = Ftilde, aes(x=Component_1, y=Component_2, group=as.factor(subject), col=as.factor(RFgroup))) +
#       geom_segment(data = Btilde, aes(x=0,y=0,xend=Component_1,yend=Component_2), arrow=arrow(), color="red") +
#       geom_label_repel(data = Btilde, aes(x=Component_1, y=Component_2, label=featureLabels), max.overlaps=5, segment.color=NA) +
#       scale_x_continuous(name="Component 1",sec.axis=sec_axis(name="Component 1",~./enlargingFactor)) +
#       scale_y_continuous(name="Component 2",sec.axis=sec_axis(name="Component 2",~./enlargingFactor)) +
#       theme_classic() +
#       theme(axis.text.y.right=element_text(colour="red"), axis.title.y.right=element_text(colour="red"), axis.title.x.top=element_text(colour="red"), axis.text.x.top=element_text(colour="red")) +
#       ggtitle(titleString)
#   }
#   else{
#     enlargingFactor = apply(Ftilde, 2, max) / apply(Btilde, 2, max)
#     Btilde = Btilde * enlargingFactor
#     
#     Ftilde = Ftilde %>% as_tibble() %>% mutate(subject = rep(model[[1]]$subject,each=nrow(model[[3]])), timepoint = rep(model[[3]]$days, nrow(model[[1]])))
#     colnames(Ftilde) = c("Component_1", "subject", "timepoint")
#     featureLabels = model[[2]] %>% left_join(filteredFeatures) %>% filter(keep==TRUE)
#     featureLabels = featureLabels %>% select(featureNames) %>% pull()
#     Btilde = Btilde %>% as_tibble() %>% mutate(featureNames = model[[2]]$featureNames) %>% left_join(filteredFeatures) %>% filter(keep==TRUE) %>% select(V1)
#     colnames(Btilde) = c("Component_1")
#     Ftilde = Ftilde %>% left_join(model[[1]], by="subject")
#     
#     set.seed(1)
#     jitterpoints = rep(jitter(rep(0, nrow(model[[1]]))), each=nrow(model[[3]]))
#     Ftilde = Ftilde %>% mutate(Component_2=jitterpoints, Component_1=Component_1.x) %>% select(-Component_1.x, -Component_1.y)
#     
#     ggplot() + 
#       geom_path(data = Ftilde, aes(x=Component_1, y=Component_2, group=as.factor(subject), col=as.factor(RFgroup))) +
#       geom_segment(data = Btilde, aes(x=0,y=0,xend=Component_1,yend=0), arrow=arrow(), color="red") +
#       geom_label_repel(data = Btilde, aes(x=Component_1, y=0, label=featureLabels), max.overlaps=5, segment.color=NA) +
#       scale_x_continuous(name="Component 1",sec.axis=sec_axis(name="Component 1",~./enlargingFactor)) +
#       scale_y_continuous(name="Component 2",sec.axis=sec_axis(name="Component 2",~./enlargingFactor)) +
#       theme_classic() +
#       theme(axis.text.y.right=element_text(colour="red"), axis.title.y.right=element_text(colour="red"), axis.title.x.top=element_text(colour="red"), axis.text.x.top=element_text(colour="red")) +
#       ggtitle(titleString)
#   }
#   
# }

#timepoints = c(7, 14, 21, 28, 58)
timepoints = c(-14, -7, 0, 7, 31)

# Import the models
onlyCase = importPARAFAC("./20230913_run_onlyCase/PARAFAC models/", "Athina", numComponents=2, loadRDS=FALSE, header=FALSE)
onlyControl = importPARAFAC("./20230913_run_onlyControl/PARAFAC models/", "Athina", numComponents=3, loadRDS=FALSE, header=FALSE)

colnames(onlyCase$dataset$mode1) = c("subject", "Symptoms")
colnames(onlyCase$dataset$mode2) = c("cytokine")
colnames(onlyCase$dataset$mode3) = c("extraction")
onlyCase$dataset$mode3 = onlyCase$dataset$mode3 %>% mutate(timepoint = timepoints)

colnames(onlyControl$dataset$mode1) = c("subject", "control")
colnames(onlyControl$dataset$mode2) = c("cytokine")
colnames(onlyControl$dataset$mode3) = c("extraction")
onlyControl$dataset$mode3 = onlyControl$dataset$mode3 %>% mutate(timepoint = timepoints)

# Plot the models
colourCols = c("GenderID", "Phylum", "")
legendTitles = c("Gender identity", "Phylum", "")
xLabels = c("Subject index", "Feature index", "Time index")
legendColNums = c(2,4,0)
arrangeModes = c(TRUE, TRUE, FALSE)
continuousModes = c(FALSE,FALSE,TRUE)

varExp_t3 = round(multiway::sumsq(tongue_t3$reconstructedData) / multiway::sumsq(tongue_t3$dataset$data, na.rm=TRUE) * 100,2)
varExp_t12 = round(multiway::sumsq(tongue_t12$reconstructedData) / multiway::sumsq(tongue_t12$dataset$data, na.rm=TRUE) * 100,2)
varExp_t3_t12 = round(multiway::sumsq(tongue_t12$reconstructedData) / multiway::sumsq(tongue_t12$dataset$data, na.rm=TRUE) * 100,2)
varExp_t = round(multiway::sumsq(tongue_t$reconstructedData) / multiway::sumsq(tongue_t$dataset$data, na.rm=TRUE) * 100,2)

plotPARAFACmodel(list("A"=tongue_t3$model[[1]], "B"=tongue_t3$model[[2]], "C"=tongue_t3$model[[3]]), tongue_t3$dataset, 1, colourCols, legendTitles, xLabels, legendColNums, arrangeModes, continuousModes, overallTitle=paste0("NPLS Tongue T3 (", varExp_t3, " %)"))
plotPARAFACmodel(list("A"=tongue_t12$model[[1]], "B"=tongue_t12$model[[2]], "C"=tongue_t12$model[[3]]), tongue_t12$dataset, 1, colourCols, legendTitles, xLabels, legendColNums, arrangeModes, continuousModes, overallTitle=paste0("NPLS Tongue T12 (", varExp_t12, " %)"))


# cytokinesModel = importPARAFAC("./20230903_run_onlyCase/PARAFAC models/Cytokines", "featureNames", timepoints, "featureNames")
# 
# metadata = read.csv("./input_deduplicated_metadata_RvdP.csv", header=FALSE, sep=" ") %>% as_tibble()
# colnames(metadata) = c("subject", "visit", "sex", "age", "DMFT", "case_control")
# 
# raw_data = read.csv("./input_deduplicated_RvdP.csv", header=FALSE, sep=" ") %>% as_tibble()
# colnames(raw_data) = cytokinesModel[[2]]$featureNames
# raw_data = cbind(metadata, raw_data) %>% as_tibble()
# 
# subject_congruences = read.csv("./20230815_run_onlyCase_t6_removed/Cytokines_individual_congruence_loadings.csv", header=FALSE) %>% as_tibble()
# feature_congruences = read.csv("./20230815_run_onlyCase_t6_removed/Cytokines_feature_congruence_loadings.csv", header=FALSE) %>% as_tibble()
# 
# feature_congruences %>% ggplot(aes(x=V1,y=0)) + geom_point() + geom_vline(xintercept=-0.5) + geom_vline(xintercept=0.5)
# (colSums(cytokinesModel[[6]]^2) / colSums(cytokinesModel[[7]]^2, na.rm=TRUE)) %>% as_tibble() %>% mutate(name=cytokinesModel[[2]]$featureNames) %>% ggplot(aes(x=name,y=value)) + geom_bar(stat="identity")
# 
# temp = kmeans(cytokinesModel[[1]]$Component_1, 2)$cluster
# subjectClustering = cbind(cytokinesModel[[1]]$subject, temp) %>% as_tibble()
# colnames(subjectClustering) = c("subject", "cluster")
# 
# filteredFeatures = cytokinesModel[[2]] %>% mutate(keep=TRUE)
# biplotAcrossTime(cytokinesModel, filteredFeatures, "temp")
# 
# a = model[[1]] %>% ggplot(aes(x=as.factor(subject),y=Component_1,fill=as.factor(RFgroup))) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90)) + ggtitle("Subject scores")
# b = model[[2]] %>% ggplot(aes(x=as.factor(featureNames),y=Component_1)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90)) + ggtitle("Feature loadings")
# c = model[[3]] %>% ggplot(aes(x=timepoints,y=Component_1,label=timepoints)) + geom_line() + geom_point() + ggtitle("Time loadings")
# 
# d = model[[1]] %>% ggplot(aes(x=as.factor(subject),y=Component_2,fill=as.factor(RFgroup))) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90)) + ggtitle("Subject scores")
# e = model[[2]] %>% ggplot(aes(x=as.factor(featureNames),y=Component_2)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90)) + ggtitle("Feature loadings")
# f = model[[3]] %>% ggplot(aes(x=timepoints,y=Component_2,label=timepoints)) + geom_line() + geom_point() + ggtitle("Time loadings")