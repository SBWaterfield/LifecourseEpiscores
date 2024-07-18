library(tidyverse)
library(ggplot2)
library(meffil)
library(reshape)

#==============================================================================#
#------------------------------------------------------------------------------#
#---------------------Comparing Age 24 Model types-----------------------------#
#------------------------------------------------------------------------------#
#==============================================================================#
# Penlasied regression results
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PenRegAge24ResidProtNoSexChr_ForcedModels.Rdata")
Age24PR <- Age24PenReg[[2]]
colnames(Age24PR) <- c('Cor_PR', 'P_PR', 'Episcore')
rm(Age24PenReg)

Age24Models <- Age24PR

# Neural network results
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/NeuralNetAge24500CGResidProt.Rdata")

Age24NN <- Age24NN[[2]]
colnames(Age24NN) <- c('Cor_NN', 'P_NN', 'Episcore')

Age24Models <- merge(Age24Models, Age24NN, by='Episcore')


# Random Forest models 
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/RanForAge24500CGResidProt.Rdata")

colnames(Age24RF) <- c('Cor_RF', 'P_RF', 'Episcore')

Age24Models <- merge(Age24Models, Age24RF, by='Episcore')


#XG boost models
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/XGBoostAge24500CG500CGResidProt.Rdata")


colnames(Age24XGBoost) <- c('Cor_XG', 'P_XG', 'Episcore')

Age24Models <- merge(Age24Models, Age24XGBoost, by='Episcore')


# PCA model # 500 CpG
# load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PCAReg500_Age24_ResidProt.Rdata")
# Age24PCA <- Age24PenRegDNAmReduced[[2]]
# colnames(Age24PCA) <- c('Cor_PCA', 'P_PCA', 'Episcore')
# 
# Age24Models <- merge(Age24Models, Age24PCA, by='Episcore')

# PCA Model Epigenome wide
PCAmod1 <- read.csv(file='Age24PCA_Cor_Results1.csv')
PCAmod2 <- read.csv(file='Age24PCA_Cor_Results2.csv')
PCAmod <- rbind(PCAmod1, PCAmod2)
PCAmod <- PCAmod[2:4]
colnames(PCAmod) <- c('Cor_PCA', 'P_PCA', 'Episcore')
Age24Models <- merge(Age24Models, PCAmod, by='Episcore')
rm(PCAmod1, PCAmod2)


# SVM models
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/SVM500_Age24_ResidProt.Rdata")
Age24SVM <- Age24PenRegDNAmReduced[[2]]
colnames(Age24SVM) <- c('Cor_SVM', 'P_SVM', 'Episcore')


Age24Models <- merge(Age24Models, Age24SVM, by='Episcore')

#  KNN models
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/KNN500_Age24_ResidProt.Rdata")
Age24KNN <- Age24PenRegDNAmReduced[[2]]
colnames(Age24KNN) <- c('Cor_KNN', 'P_KNN', 'Episcore')

Age24Models <- merge(Age24Models, Age24KNN, by='Episcore')


# Bayes net anaylsis
library(readr)
Age24BN <- read_csv("Age24ResidProtCor_BN500.csv")
Age24BN <- Age24BN[,2:4]

colnames(Age24BN) <- c('Cor_BN', 'P_BN', 'Episcore')

Age24Models <- merge(Age24Models, Age24BN)


# Reuced DNAm set penalised regression
Age24PR500 <- read_csv("Age24ResidProtCor_DNAmReduced_tabletest.csv")
Age24PR500 <- Age24PR500[,2:4]

colnames(Age24PR500) <- c('Cor_PR500', 'P_PR500', 'Episcore')

Age24Models <- merge(Age24Models, Age24PR500)

#------------------------------------------------------------------------------#
# Load in EWAS-PCA Models 

#ELNET
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCAELNETAge24_EPiscores.Rdata")
EWASElnet <- Age24EWASPCAPenReg[[2]]
colnames(EWASElnet) <- c('Cor_EWAS_PCA_Elnet', 'P_EWAS_PCA_Elnet', 'Episcore')
EWASElnet$Episcore <- paste0(EWASElnet$Episcore, '_F24', sep="")


#PCA regression
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCAReg24_PenReg_Episcores.Rdata")
Age24EWASPCReg <- Age24EWASPCAPenReg[[2]]
colnames(Age24EWASPCReg) <- c('Cor_EWAS1K_PCReg', 'P_EWAS1K_PCReg', 'Episcore')
Age24EWASPCReg$Episcore <- paste0(Age24EWASPCReg$Episcore, '_F24', sep="")
Age24EWASModels <- merge(EWASElnet, Age24EWASPCReg)



# XGBoost
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCAXGBoostAge24_EPiscores.Rdata")
EWASXGBoost <- Age24EWASPCAPenReg[[2]]
colnames(EWASXGBoost) <- c('Cor_EWAS_PCA_XGBoost', 'P_EWAS_PCA_XGBoost', 'Episcore')
EWASXGBoost$Episcore <- paste0(EWASXGBoost$Episcore, '_F24', sep="")
Age24EWASModels <- merge(Age24EWASModels, EWASXGBoost)


  
#SVM
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCASVMAge24_EPiscores.Rdata")
EWASSVM <- Age24EWASPCAPenReg[[2]]
colnames(EWASSVM) <- c('Cor_EWAS_PCA_SVM', 'P_EWAS_PCA_SVM', 'Episcore')
EWASSVM$Episcore <- paste0(EWASSVM$Episcore, '_F24', sep="")
Age24EWASModels <- merge(Age24EWASModels, EWASSVM)

# Random Forest
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCARFAge24_EPiscores.Rdata")
EWASRF <- Age24EWASPCAPenReg[[2]]
colnames(EWASRF) <- c('Cor_EWAS_PCA_RF', 'P_EWAS_PCA_RF', 'Episcore')
EWASRF$Episcore <- paste0(EWASRF$Episcore, '_F24', sep="")
Age24EWASModels <- merge(Age24EWASModels, EWASRF)


#Neural Net
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCANeuralNetAge24_EPiscores.Rdata")
EWASNeuralNet <- Age24EWASPCAPenReg[[2]]
colnames(EWASNeuralNet) <- c('Cor_EWAS_PCA_NeuralNet', 'P_EWAS_PCA_NeuralNet', 'Episcore')
EWASNeuralNet$Episcore <- paste0(EWASNeuralNet$Episcore, '_F24', sep="")
Age24EWASModels <- merge(Age24EWASModels, EWASNeuralNet)


#KNN
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCAKNNAge24_EPiscores.Rdata")
EWASKNN <- Age24EWASPCAPenReg[[2]]
colnames(EWASKNN) <- c('Cor_EWAS_PCA_KNN', 'P_EWAS_PCA_KNN', 'Episcore')
EWASKNN$Episcore <- paste0(EWASKNN$Episcore, '_F24', sep="")
Age24EWASModels <- merge(Age24EWASModels, EWASKNN)


#Bayes Net
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EWASMothers_PCASBNAge24_EPiscores.Rdata")
EWASBN <- Age24EWASPCAPenReg[[2]]
colnames(EWASBN) <- c('Cor_EWAS_PCA_BN', 'P_EWAS_PCA_BN', 'Episcore')
EWASBN$Episcore <- paste0(EWASBN$Episcore, '_F24', sep="")
Age24EWASModels <- merge(Age24EWASModels, EWASBN)


# Combine All models togeher 
Age24Models <- merge(Age24Models, Age24EWASModels)

# Quick descriptive stats of epigenome wide vs 5000cpg vs EWASPCA models
Epigenomewide <- Age24Models
Epigenomewide$SigRes <- ifelse(
  (Epigenomewide$P_PR <0.05 | Epigenomewide$P_PCA <0.05),
  T, F)
Epigenomewide <- Epigenomewide[Epigenomewide$SigRes==T,]
Epigenomewide$SigRes <- ifelse(
  (Epigenomewide$Cor_PR >0.1 | Epigenomewide$Cor_PCA >0.1),
  T, F)
Epigenomewide <- Epigenomewide[Epigenomewide$SigRes==T,]
Epigenomewide <- Epigenomewide[rowSums(is.na(Epigenomewide)) != ncol(Epigenomewide), ]

Reduced500 <- Age24Models
Reduced500$SigRes <- ifelse(
  (Reduced500$P_NN <0.05 | Reduced500$P_RF <0.05 |
     Reduced500$P_XG <0.05 | Reduced500$P_SVM <0.05 |
     Reduced500$P_KNN <0.05 |Reduced500$P_BN<0.05 | Reduced500$P_PR500<0.05),
  T, F)
Reduced500 <- Reduced500[Reduced500$SigRes==T,]
Reduced500$SigRes <- ifelse(
  (Reduced500$Cor_NN >0.1 | Reduced500$Cor_RF >0.1 |
     Reduced500$Cor_XG >0.1 | Reduced500$Cor_SVM >0.1 |
     Reduced500$Cor_KNN >0.1 |Reduced500$Cor_BN>0.1 | Reduced500$Cor_PR500>0.1),
  T, F)
Reduced500 <- Reduced500[Reduced500$SigRes==T,]
Reduced500 <- Reduced500[rowSums(is.na(Reduced500)) != ncol(Reduced500), ]


EWASPCA <- Age24Models
EWASPCA$SigRes <- ifelse(
  (EWASPCA$P_EWAS1K_PCReg <0.05 |
     EWASPCA$P_EWAS_PCA_Elnet < 0.05 | EWASPCA$P_EWAS_PCA_NeuralNet < 0.05 | 
     EWASPCA$P_EWAS_PCA_XGBoost < 0.05 | EWASPCA$P_EWAS_PCA_RF < 0.05 | 
     EWASPCA$P_EWAS_PCA_SVM < 0.05 | EWASPCA$P_EWAS_PCA_KNN < 0.05 | 
     EWASPCA$P_EWAS_PCA_BN < 0.05),
  T, F)
EWASPCA <- EWASPCA[EWASPCA$SigRes==T,]
EWASPCA$SigRes <- ifelse(
  (EWASPCA$Cor_EWAS1K_PCReg >0.1 |
     EWASPCA$Cor_EWAS_PCA_Elnet >0.1 | EWASPCA$Cor_EWAS_PCA_NeuralNet >0.1 | 
     EWASPCA$Cor_EWAS_PCA_XGBoost >0.1 | EWASPCA$Cor_EWAS_PCA_RF >0.1 | 
     EWASPCA$Cor_EWAS_PCA_SVM >0.1 | EWASPCA$Cor_EWAS_PCA_KNN >0.1 | 
     EWASPCA$Cor_EWAS_PCA_BN >0.1),
  T, F)
EWASPCA <- EWASPCA[EWASPCA$SigRes==T,]
EWASPCA <- EWASPCA[rowSums(is.na(EWASPCA)) != ncol(EWASPCA), ]


colnames(Age24PR) <- c('Cor', 'P', 'Episcore')
Age24PR$Model <- 'Elnet'
Age24PR$FDR <- p.adjust(Age24PR$P, method = 'fdr')
FDRModels <- Age24PR


colnames(Age24NN) <- c('Cor', 'P', 'Episcore')
Age24NN$Model <- 'NN'
Age24NN$FDR <- p.adjust(Age24NN$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24NN)

colnames(Age24RF) <- c('Cor', 'P', 'Episcore')
Age24RF$Model <- 'RF'
Age24RF$FDR <- p.adjust(Age24RF$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24RF)


colnames(Age24XGBoost) <- c('Cor', 'P', 'Episcore')
Age24XGBoost$Model <- 'XGBoost'
Age24XGBoost$FDR <- p.adjust(Age24XGBoost$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24XGBoost)


colnames(PCAmod) <- c('Cor', 'P', 'Episcore')
PCAmod$Model <- 'PCReg'
PCAmod$FDR <- p.adjust(PCAmod$P, method = 'fdr')
FDRModels <- rbind(FDRModels, PCAmod)


colnames(Age24SVM) <- c('Cor', 'P', 'Episcore')
Age24SVM$Model <- 'SVM'
Age24SVM$FDR <- p.adjust(Age24SVM$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24SVM)


colnames(Age24KNN) <- c('Cor', 'P', 'Episcore')
Age24KNN$Model <- 'KNN'
Age24KNN$FDR <- p.adjust(Age24KNN$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24KNN)


colnames(Age24BN) <- c('Cor', 'P', 'Episcore')
Age24BN$Model <- 'BN'
Age24BN$FDR <- p.adjust(Age24BN$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24BN)


colnames(Age24PR500) <- c('Cor', 'P', 'Episcore')
Age24PR500$Model <- 'Elnet500'
Age24PR500$FDR <- p.adjust(Age24PR500$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24PR500)

colnames(EWASElnet) <- c('Cor', 'P', 'Episcore')
EWASElnet$Model <- 'EWAS_Elnet'
EWASElnet$FDR <- p.adjust(EWASElnet$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASElnet)


colnames(Age24EWASPCReg) <- c('Cor', 'P', 'Episcore')
Age24EWASPCReg$Model <- 'EWAS_PCReg'
Age24EWASPCReg$FDR <- p.adjust(Age24EWASPCReg$P, method = 'fdr')
FDRModels <- rbind(FDRModels, Age24EWASPCReg)


colnames(EWASXGBoost) <- c('Cor', 'P', 'Episcore')
EWASXGBoost$Model <- 'EWAS_XGBoost'
EWASXGBoost$FDR <- p.adjust(EWASXGBoost$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASXGBoost)

colnames(EWASSVM) <- c('Cor', 'P', 'Episcore')
EWASSVM$Model <- 'EWAS_SVM'
EWASSVM$FDR <- p.adjust(EWASSVM$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASSVM)


colnames(EWASRF) <- c('Cor', 'P', 'Episcore')
EWASRF$Model <- 'EWAS_RF'
EWASRF$FDR <- p.adjust(EWASRF$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASRF)


colnames(EWASNeuralNet) <- c('Cor', 'P', 'Episcore')
EWASNeuralNet$Model <- 'EWAS_NN'
EWASNeuralNet$FDR <- p.adjust(EWASNeuralNet$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASNeuralNet)


colnames(EWASKNN) <- c('Cor', 'P', 'Episcore')
EWASKNN$Model <- 'EWAS_KNN'
EWASKNN$FDR <- p.adjust(EWASKNN$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASKNN)


colnames(EWASBN) <- c('Cor', 'P', 'Episcore')
EWASBN$Model <- 'EWAS_BN'
EWASBN$FDR <- p.adjust(EWASBN$P, method = 'fdr')
FDRModels <- rbind(FDRModels, EWASBN)

# Super stringent P values
FDRModels$StringentFDR <- p.adjust(FDRModels$P, method = 'fdr')
FDRModelsStringent <- FDRModels[FDRModels$StringentFDR<=0.05,]
table(FDRModelsStringent$Model)

# Just Raw P values
FDRModels <- FDRModels[FDRModels$P<=0.05,]
FDRModels <- FDRModels[FDRModels$Cor>=0.1,]
table(FDRModels$Model)

# Method based FDR
FDRModels <- FDRModels[FDRModels$FDR<=0.05,]
table(FDRModels$Model)


# Comparison with Gadd results
library(readr)
Gadd_SomaResults <- read_csv("Gadd_SomaResults.csv")
Gadd_OlinkResults <- read_csv("Gadd_OlinkResults.csv")



#------------------------------------------------------------------------------#
#Plot results for 500 Var models
# Singificant figures for data
Age24Models[,c(2:ncol(Age24Models))] <- signif(Age24Models[,c(2:ncol(Age24Models))], 2)

# Add in a column saying if there is at least one significant result
Age24Models$SigRes <- ifelse(
                      (Age24Models$P_PR <0.05 | Age24Models$P_NN <0.05 | Age24Models$P_RF <0.05 |
                         Age24Models$P_XG <0.05 | Age24Models$P_PCA <0.05 | Age24Models$P_SVM <0.05 |
                         Age24Models$P_KNN <0.05 |Age24Models$P_BN<0.05 | Age24Models$P_PR500<0.05 |
                         Age24Models$P_EWAS1K_PCReg <0.05 |
                         Age24Models$P_EWAS_PCA_Elnet < 0.05 | Age24Models$P_EWAS_PCA_NeuralNet < 0.05 | 
                         Age24Models$P_EWAS_PCA_XGBoost < 0.05 | Age24Models$P_EWAS_PCA_RF < 0.05 | 
                         Age24Models$P_EWAS_PCA_SVM < 0.05 | Age24Models$P_EWAS_PCA_KNN < 0.05 | 
                         Age24Models$P_EWAS_PCA_BN < 0.05),
                      T, F)


Age24ModelsFull <- Age24Models
Age24Models <- Age24Models[Age24Models$SigRes == T,]
sum(Age24Models$SigRes==T, na.rm = T)

#Sig result AND positive correlation 
Age24Models$SigRes <- ifelse(
  (Age24Models$Cor_PR >=0.1 | Age24Models$Cor_NN >=0.1 | Age24Models$Cor_RF >=0.1 |
     Age24Models$Cor_XG >=0.1 | Age24Models$Cor_PCA >=0.1 | Age24Models$Cor_SVM >=0.1 |
     Age24Models$Cor_KNN >=0.1 |Age24Models$Cor_BN>=0.1 | Age24Models$Cor_PR500>=0.1 |
     Age24Models$Cor_EWAS1K_PCReg>=0.1 |
     Age24Models$Cor_EWAS_PCA_Elnet >=0.1 | Age24Models$Cor_EWAS_PCA_NeuralNet >=0.1 | 
     Age24Models$Cor_EWAS_PCA_XGBoost >=0.1 | Age24Models$Cor_EWAS_PCA_RF >=0.1 | 
     Age24Models$Cor_EWAS_PCA_SVM >=0.1 | Age24Models$Cor_EWAS_PCA_KNN >=0.1 | 
     Age24Models$Cor_EWAS_PCA_BN >=0.1),
  T, F)
Age24Models <- Age24Models[Age24Models$SigRes == T,]



Age24Models <- Age24Models[!is.na(Age24Models$Episcore), ] 

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

Age24ModelsHeatmap <- data.frame(t(Age24Models))
Age24ModelsHeatmap <- header.true(Age24ModelsHeatmap)
Heatmapnames <- c('ElasticNet', 'NeuralNet', 'RandomForest', 'XGBoost', 'PCARegression(10PC)', 'SVM', 'KNN', 'BayesNet', 'ElNet500', 'EWAS_PCARegression(10PC)',
                  'EWAS_PCA_ElasticNet', 'EWAS_PCA_XGBoost','EWAS_PCA_SVM','EWAS_PCA_RF','EWAS_PCA_NeuralNet',
                  'EWAS_PCA_KNN','EWAS_PCA_BayesNet')
Age24ModelsHeatmap <- Age24ModelsHeatmap[c(1,3,5,7,9,11,13,15,17,19,
                                           21,23,25,27,29,31,33),]

colnames(Age24ModelsHeatmap) <- gsub('_F24', '', colnames(Age24ModelsHeatmap))


# Make NAs zero
Age24ModelsHeatmap[] <- lapply(Age24ModelsHeatmap, function(x) as.numeric(replace(x, is.na(x), 0)))
rownames(Age24ModelsHeatmap) <- Heatmapnames
Age24ModelsHeatmap <- as.matrix(Age24ModelsHeatmap)
rowdim <- row.names(Age24ModelsHeatmap)
coldim <- colnames(Age24ModelsHeatmap)
Age24ModelsHeatmap <- matrix(as.numeric(Age24ModelsHeatmap), ncol = ncol(Age24ModelsHeatmap))
rownames(Age24ModelsHeatmap) <- rowdim
colnames(Age24ModelsHeatmap) <- coldim



library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-max(Age24ModelsHeatmap), 0, max(Age24ModelsHeatmap)), c("red", "white", "blue"))
Heatmap(Age24ModelsHeatmap, col = col_fun, name = 'Correlation (r)',
        layer_fun = function(j, i, x, y, width, height, fill) {
          v = pindex(Age24ModelsHeatmap, i, j)
          l = v == 0
          s = v >=0.1
          grid.text(sprintf("%s", '-'), x[l], y[l], gp = gpar(fontsize = 10))
          grid.text(sprintf("%s", '*'), x[s], y[s], gp = gpar(fontsize = 15))
        })

# Results per Method
FullResults <-t(Age24ModelsHeatmap)

FullModelResults <- data.frame()
for (col in colnames(FullResults)){
  
  df <- as.data.frame(FullResults[,col])
  df$Protein <- rownames(df)
  df <- df[df$`FullResults[, col]` >=0.1,]
  if (nrow(df >=1)){
  df$Method <- col
  FullModelResults <- rbind(FullModelResults, df)
  }
}

table(FullModelResults$Method)
colnames(FullModelResults) <- c('Cor', 'Protein', 'Model')

ggplot(FullModelResults, aes(Model, Cor)) +   
  geom_bar(aes(fill = Protein), position = "dodge", stat="identity")

# Plot these results 
Episcores <- Age24Models$Episcore
minx <- min(c(Age24Models$P_PR, Age24Models$P_NN, Age24Models$P_RF, Age24Models$P_XG, 
              Age24Models$P_PCA, Age24Models$P_SVM, Age24Models$P_KNN,Age24Models$P_BN, Age24Models$P_PR500,
              Age24Models$P_EWAS1K_PCReg), na.rm = T)
maxx <- max(c(Age24Models$P_PR, Age24Models$P_NN, Age24Models$P_RF, Age24Models$P_XG, 
              Age24Models$P_PCA, Age24Models$P_SVM, Age24Models$P_KNN,  Age24Models$P_PR500,
              Age24Models$P_EWAS1K_PCReg), na.rm = T)
miny <- min(c(Age24Models$Cor_PR, Age24Models$Cor_NN, Age24Models$Cor_RF, Age24Models$Cor_XG, 
              Age24Models$Cor_PCA, Age24Models$Cor_SVM, Age24Models$Cor_KNN, Age24Models$Cor_BN, Age24Models$Cor_PR500,
              Age24Models$Cor_EWAS1K_PCReg), na.rm = T)
maxy <- max(c(Age24Models$Cor_PR, Age24Models$Cor_NN, Age24Models$Cor_RF, Age24Models$Cor_XG, 
              Age24Models$Cor_PCA, Age24Models$Cor_SVM, Age24Models$Cor_KNN, Age24Models$Cor_BN, Age24Models$Cor_PR500,
              Age24Models$Cor_EWAS1K_PCReg), na.rm = T)

library(ggplot2)
prot <- 'IL6_F24'
for (prot in Episcores){

  df <- Age24Models[Age24Models$Episcore == prot,]
  dfcor <- df[, c(2,4,6,8,10,12,14,16,18,20)]
  dfpval <- df[,c(3,5,7,9,11,13,15,17,19,21)]
  df <- data.frame(cbind(t(dfcor), t(dfpval), c('PenReg', 'NeuralNet', 'RandomForest', 'XGBoost', 'PCA', 'SVM', 'KNN', 'BayesNet', 'PenReg500', 'PCAMother')))
  colnames(df) <- c('Correlation', 'Pval', 'Model')
  df$Model <- factor(df$Model, levels = c('PenReg', 'PenReg500', 'PCA', 'RandomForest', 'XGBoost', 'SVM', 'KNN', 'BayesNet', 'NeuralNet', 'PCAMother'))
  
  
  df$Correlation <-as.numeric(df$Correlation)
  df$Pval <-as.numeric(df$Pval)
  df$Pval <- log(df$Pval)
  colnames(df) <- c('Correlation', 'Pval', 'Model')
  
  # # Scatter plot by group
  # group.colors <- c(PenReg = "Red", PenReg500 = "Red", RandomForest ="orange", XGBoost = "Purple",
  #                   PCA = 'black', SVM = 'Green', KNN = 'Cyan',  NeuralNet = "Blue")
  # group.shapes <- c(PenReg = "circle", PenReg500 = "circle cross", RandomForest ="circle", 
  #                   XGBoost = "circle", PCA = 'circle', SVM = 'circle', KNN = 'circle', NeuralNet = "circle")
  # 
  # 
  # # tidied up plot
  # print(
  #   ggplot(data=subset(df, !is.na(Correlation)), aes(x = Pval, y = Correlation, color = Model, shape = Model)) +
  #     geom_point(size=3) + xlab('Log P Value') + ylab('Correlation')  +
  #     xlim(log(minx), log(maxx)) +ylim(miny-0.05, maxy+0.05) +
  #     geom_vline(xintercept = log(0.05),linetype = 2) + geom_hline(yintercept = 0.1, linetype = 2) +
  #     geom_hline(yintercept = 0.0) +
  #     ggtitle(paste('Correlation between Episcore and Protein expression For: ',prot )) +
  #     scale_colour_manual(values=group.colors) + scale_shape_manual(values=group.shapes) +
  #     theme_classic() +  
  #     geom_rect(aes(xmin = -Inf, xmax = log(0.05), ymin = 0.1, ymax = Inf), alpha = 0.015,fill = "green", inherit.aes = F)
  # )
  
  sigtest <- ifelse(dfpval <0.05, '*', '')
  sigtest <- sigtest[!is.na(sigtest)]
  prot <- gsub('_F24', '', prot)
  
  print(
    ggplot(data=subset(df, !is.na(Correlation)), aes(x=Model, y=Correlation)) +
      geom_bar(stat='identity', fill='steelblue') + ylim(miny-0.05, maxy+0.05) +
      geom_hline(yintercept = 0.1,linetype = 2) + geom_hline(yintercept = -0.1,linetype = 2) +
      geom_hline(yintercept = 0,linetype = 1) +
      ggtitle(paste('Correlation between Episcore and Protein expression For: ',prot )) +
      ylab('Correlation (r)') + xlab('Model') +
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
      geom_text(aes(label = sigtest),
                hjust = 0.5,
                color = "black",
                size = 10)
  )
 Sys.sleep(5)
}


#==============================================================================#
#------------------------------------------------------------------------------#
#-----------------------Comparing FOM1 Model types-----------------------------#
#------------------------------------------------------------------------------#
#==============================================================================#
# # Penlasied regression results
# load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PenRegFOM1ResidProt.Rdata")
# FOM1PR <- FOM1PenReg[[2]]
# colnames(FOM1PR) <- c('Cor_PR', 'P_PR', 'Episcore')
# rm(FOM1PenReg)
# 
# FOM1Models <- FOM1PR
# 
# # Neural network results
# load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/NeuralNetFOM1_500CG.Rdata")
# 
# FOM1NN <- FOM1NN[[2]]
# colnames(FOM1NN) <- c('Cor_NN', 'P_NN', 'Episcore')
# 
# FOM1Models <- merge(FOM1Models, FOM1NN, by='Episcore')
# 
# # Random Forest models 
# load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/RanForFOM1_500CGResidProt.Rdata")
# 
# colnames(FOM1RF) <- c('Cor_RF', 'P_RF', 'Episcore')
# 
# FOM1Models <- merge(FOM1Models, FOM1RF, by='Episcore')
# 
# 
# #XG boost models
# load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/XGBoostAgeFOM1_500CGResidProt.Rdata")
# 
# 
# colnames(FOM1XGBoost) <- c('Cor_XG', 'P_XG', 'Episcore')
# 
# FOM1Models <- merge(FOM1Models, FOM1XGBoost, by='Episcore')
# 
# # Penalised regression with reduced DNAm 
# load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PenRegFOM1_500CG.Rdata")
# 
# FOM1PR500 <- FOM1PenReg[[2]]
# colnames(FOM1PR500) <- c('Cor_PR500', 'P_PR500', 'Episcore')
# rm(FOM1PenReg)
# 
# 
# FOM1Models <- merge(FOM1Models, FOM1PR500)
# #------------------------------------------------------------------------------#
# #Plot results 
# 
# # Singificant figures for data
# FOM1Models[,c(2:ncol(FOM1Models))] <- signif(FOM1Models[,c(2:ncol(FOM1Models))], 2)
# 
# # Add in a column saying if there is at least one significant result
# FOM1Models$SigRes <- ifelse(
#   (FOM1Models$P_PR <0.05 | FOM1Models$P_NN <0.05 | FOM1Models$P_RF <0.05 |
#      FOM1Models$P_XG <0.05 | FOM1Models$P_PR500<0.05),
#   T, F)
# 
# 
# FOM1Models <- FOM1Models[FOM1Models$SigRes == T,]
# FOM1Models <- FOM1Models[!is.na(FOM1Models$Episcore), ] 
# 
# # Plot these results 
# Episcores <- FOM1Models$Episcore
# minx <- min(c(FOM1Models$P_PR, FOM1Models$P_NN, FOM1Models$P_RF, FOM1Models$P_XG, FOM1Models$P_PR500), na.rm = T)
# maxx <- max(c(FOM1Models$P_PR, FOM1Models$P_NN, FOM1Models$P_RF, FOM1Models$P_XG, FOM1Models$P_PR500), na.rm = T)
# miny <- min(c(FOM1Models$Cor_R, FOM1Models$Cor_NN, FOM1Models$Cor_RF, FOM1Models$Cor_XG, FOM1Models$Cor_PR500), na.rm = T)
# maxy <- max(c(FOM1Models$Cor_R, FOM1Models$Cor_NN, FOM1Models$Cor_RF, FOM1Models$Cor_XG, FOM1Models$Cor_PR500), na.rm = T)
# 
# for (prot in Episcores){
#   
#   df <- FOM1Models[FOM1Models$Episcore == prot,]
#   dfcor <- df[, c(2,4,6,8,10)]
#   dfpval <- df[,c(3,5,7,9,11)]
#   df <- data.frame(cbind(t(dfcor), t(dfpval), c('PenReg', 'NeuralNet', 'RandomForest', 'XGBoost', 'PenReg500')))
#   colnames(df) <- c('Correlation', 'Pval', 'Model')
#   
#   df$Correlation <-as.numeric(df$Correlation)
#   df$Pval <-as.numeric(df$Pval)
#   df$Pval <- log(df$Pval)
#   colnames(df) <- c('Correlation', 'Pval', 'Model')
#   
#   # Scatter plot by group
#   group.colors <- c(PenReg = "Red", PenReg500 = "Red", RandomForest ="orange", XGBoost = "Purple", NeuralNet = "Blue")
#   group.shapes <- c(PenReg = "circle", PenReg500 = "circle cross", RandomForest ="circle", XGBoost = "circle", NeuralNet = "circle")
#   print(
#     ggplot(data=subset(df, !is.na(Correlation)), aes(x = Pval, y = Correlation, color = Model, shape = Model)) +
#       geom_point(size=4) + xlab('Log P Value') + ylab('Correlation')  +
#       xlim(log(minx), log(maxx)) +ylim(miny-0.05, maxy+0.05) +
#       geom_vline(xintercept = log(0.05)) + geom_hline(yintercept = 0.1) +
#       geom_hline(yintercept = 0.0, linetype = 2 ) +
#       ggtitle(paste('Correlation between Episcore and Protein expression For: ',prot )) +
#       scale_colour_manual(values=group.colors) + scale_shape_manual(values=group.shapes)
#   )
#   
# }
# 
# 
# 
# 









