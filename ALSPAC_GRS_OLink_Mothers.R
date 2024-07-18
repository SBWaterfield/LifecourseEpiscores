

library(dplyr)

# Load in mothers data
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EpiscoresFOM1.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/MothersVars.RData")

library(labelled)
MothersVars <- remove_labels(MothersVars)
MothersVars[MothersVars<0] <- NA
 
# Function to deal with extracted genotype data
meffil.extract.genotypes <- function(filenames, verbose=F) {
  stopifnot(all(sapply(filenames, file.exists)))
  
  table.list <- lapply(filenames, function(filename) {
    read.table(filename, header=T)
  })
  
  sample.names <- table.list[[1]][,1] ## family id
  genotypes <- lapply(table.list, function(genotype.table) {
    genotype.table[match(sample.names, genotype.table[,1]), ## match family ids
                   -(1:6),
                   drop=F]
  })
  genotypes <- do.call(cbind, genotypes)
  colnames(genotypes) <- sub("_.*", "", colnames(genotypes)) 
  rownames(genotypes) <- sample.names
  
  t(as.matrix(genotypes))
}

#------------------------------------------------------------------------------#
# Load in data 
GRSbuilder <- function(files, modname){
  
  
  table.list <- lapply(files, function(filename) {
    read.table(filename, header=T)
  })
  
  sample.names <- table.list[[1]][,1] ## family id
  genotypes <- lapply(table.list, function(genotype.table) {
    genotype.table[match(sample.names, genotype.table[,1]), ## match family ids
                   -(1:6),
                   drop=F]
  })
  genotypes <- do.call(cbind, genotypes)
  rownames(genotypes) <- sample.names
  
  
  filename <- paste("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/Xu_GRS_Mods/", modname, '.txt', sep = '')
  Mod <- read.csv(filename, sep = '\t')
  Mod$RSID <- paste0(Mod$rsid, '_', Mod$effect_allele)
  Mod$FlipRSID <- paste0(Mod$rsid, '_', Mod$other_allele)
  
  # If ALSPAC reference SNP == model effect SNP, just make a new SNP column
  for (rsid in Mod$RSID){
    if (rsid %in% names(genotypes)){
      newname <- sub("_.*", "", rsid)
      genotypes[[newname]] <- genotypes[[rsid]] # Make SNP column with no _A/C/T/G 
    }
    
  }
  
  # if SNPs are flipped, flip
  for (rsid in Mod$FlipRSID){
    if (rsid %in% names(genotypes)){
      newname <- sub("_.*", "", rsid)
      
      genotypes[[rsid]] <- case_match(genotypes[[rsid]], 2~0, 1~1, 0~2)
      
      genotypes[[newname]] <- genotypes[[rsid]] 
      
    }
  }
  
  
  
  
  data <- t(as.matrix(genotypes))
  Mod <- Mod[Mod$rsid %in% rownames(data),]
  data <- data[row.names(data) %in% Mod$rsid,]
  data <- as.data.frame(data)
  if (ncol(data) > 1){
    Mod <- Mod[match(rownames(data), Mod$rsid),]
  }
  
  GRS <- data * Mod$effect
  if (ncol(GRS) >1){
    GRS <- colSums(GRS, na.rm = T)
  }
  GRS <- as.data.frame(GRS)
  GRS$ID <- rownames(GRS)
  GRS$ID <- gsub('X', '', GRS$ID)
  # Get IDs that end in a A or B (kids only) 
  GRS <- subset(GRS,grepl("^.+(M)$",ID))
  
  return(GRS)
}


# All SNPs of interest for all GRS
SNPs <- c('data_chr01.raw', 'data_chr02.raw', 'data_chr03.raw', 'data_chr04.raw', 
          'data_chr05.raw', 'data_chr06.raw', 'data_chr07.raw', 'data_chr08.raw', 
          'data_chr09.raw', 'data_chr10.raw', 'data_chr11.raw', 'data_chr12.raw',
          'data_chr14.raw', 'data_chr17.raw', 'data_chr19.raw', 'data_chr20.raw',
          'data_chr21.raw', 'data_chr22.raw')

# Build GRS
ADAOlink <- GRSbuilder(SNPs, 'ADA_OlinkMod')
AXIN1Olink <- GRSbuilder(SNPs, 'AXIN1_OlinkMod')
CASP8Olink <- GRSbuilder(SNPs,'CASP8_OlinkMod')
CCL11Soma <- GRSbuilder(SNPs,'CCL11_SomaMod')
CCL19Olink <- GRSbuilder(SNPs,'CCL19_OlinkMod')
CCL19Soma <- GRSbuilder(SNPs,'CCL19_SomaMod')
CCL20Olink <- GRSbuilder(SNPs, 'CCL20_OlinkMod')
CCL23Olink <- GRSbuilder(SNPs, 'CCL23_OlinkMod')
CCL23Soma <- GRSbuilder(SNPs, 'CCL23_SomaMod')
CCL25Olink <- GRSbuilder(SNPs, 'CCL25_OlinkMod')
CCL25Soma1 <- GRSbuilder(SNPs, 'CCL25_SomaMod1')
CCL25Soma2 <- GRSbuilder(SNPs, 'CCL25_SomaMod2')
CD6Olink <- GRSbuilder(SNPs, 'CD6_OlinkMod')
CD244Olink <- GRSbuilder(SNPs, 'CD244_OlinkMod')
CDCP1Olink <- GRSbuilder(SNPs, 'CDCP1_OlinkMod')
CSF1Olink <- GRSbuilder(SNPs, 'CSF1_OlinkMod')
CXCL1Soma <- GRSbuilder(SNPs, 'CXCL1_SomaMod')
CXCL5Olink <- GRSbuilder(SNPs, 'CXCL5_OlinkMod')
CXCL5Soma <- GRSbuilder(SNPs, 'CXCL5_SomaMod')
CXCL1Olink <- GRSbuilder(SNPs, 'CXCL1_OlinkMod')
CXCL9Olink <- GRSbuilder(SNPs, 'CXCL9_OlinkMod')
CXCL10Olink <- GRSbuilder(SNPs, 'CXCL10_OlinkMod')
CXCL10Soma <- GRSbuilder(SNPs, 'CXCL10_SomaMod')
CXCL11Soma <- GRSbuilder(SNPs, 'CXCL11_SomaMod')
FGF5Olink <- GRSbuilder(SNPs, 'FGF5_OlinkMod')
FGF19Olink <- GRSbuilder(SNPs, 'FGF19_OlinkMod')
FGF19Soma <- GRSbuilder(SNPs, 'FGF19_SomaMod')
FGF21Olink <- GRSbuilder(SNPs, 'FGF21_OlinkMod')
Flt3Soma <- GRSbuilder(SNPs, 'Flt3_SomaMod')
GDNFOlink <- GRSbuilder(SNPs, 'GDNF_OlinkMod')
HGFSoma <- GRSbuilder(SNPs, 'HGF_SomaMod')
IFNGammaSoma <- GRSbuilder(SNPs, 'IFNgamma_SomaMod')
IL5Soma <- GRSbuilder(SNPs, 'IL5_SomaMod')
IL6Olink <- GRSbuilder(SNPs, 'IL6_OlinkMod')
IL6Soma <- GRSbuilder(SNPs, 'IL6_SomaMod')
IL10Olink <- GRSbuilder(SNPs, 'IL10_OlinkMod')
IL10Soma <- GRSbuilder(SNPs, 'IL10_SomaMod')
IL10RBOlink <- GRSbuilder(SNPs, 'IL10RB_OlinkMod')
IL10RBSoma <- GRSbuilder(SNPs, 'IL10RB_SomaMod')
IL12BOlink <- GRSbuilder(SNPs, 'IL12B_OlinkMod')
IL17ASoma <- GRSbuilder(SNPs, 'IL17A_SomaMod')
IL18R1Olink <- GRSbuilder(SNPs, 'IL18R1_OlinkMod')
IL18R1Soma1 <- GRSbuilder(SNPs, 'IL18R1_SomaMod1')
IL18R1Soma2 <- GRSbuilder(SNPs, 'IL18R1_SomaMod2')
IL22RA1Soma <- GRSbuilder(SNPs, 'IL22RA1_SomaMod')
LIFSoma <- GRSbuilder(SNPs, 'LIF_SomaMod')
LIFROlink <- GRSbuilder(SNPs, 'LIFR_OlinkMod')
LIFRSoma <- GRSbuilder(SNPs, 'LIFR_SomaMod')
MMP1Soma <- GRSbuilder(SNPs, 'MMP1_SomaMod')
MMP10Olink <- GRSbuilder(SNPs, 'MMP10_OlinkMod')
MMP10Soma <- GRSbuilder(SNPs, 'MMP10_SomaMod')
OSMSoma <- GRSbuilder(SNPs, 'OSM_SomaMod')
SIRT2Olink <- GRSbuilder(SNPs, 'SIRT2_OlinkMod')
SIRT2Soma <- GRSbuilder(SNPs, 'SIRT2_SomaMod')
SLAMF1Olink <- GRSbuilder(SNPs, 'SLAMF1_OlinkMod')
TGFaOlink <- GRSbuilder(SNPs, 'TGFa_OlinkMod')
TNFSF14Olink <- GRSbuilder(SNPs, 'TNFSF14_OlinkMod')
TNFSF14Soma <- GRSbuilder(SNPs, 'TNFSF14_SomaMod')
VEGFASoma1 <- GRSbuilder(SNPs, 'VEGFA_SomaMod1')
VEGFASoma2 <- GRSbuilder(SNPs, 'VEGFA_SomaMod2')
VEGFASoma3 <- GRSbuilder(SNPs, 'VEGFA_SomaMod3')

GRSdf <- cbind(ADAOlink, AXIN1Olink, CASP8Olink, CCL11Soma, CCL19Olink, CCL19Soma, 
               CCL20Olink, CCL23Olink, CCL23Soma, CCL25Olink, CCL25Soma1, CCL25Soma2, 
               CD6Olink, CD244Olink, CDCP1Olink, CSF1Olink, CXCL1Soma, CXCL5Olink, 
               CXCL5Soma, CXCL1Olink, CXCL9Olink, CXCL10Olink, CXCL10Soma,CXCL11Soma,
               FGF5Olink, FGF19Olink, FGF19Soma, FGF21Olink, Flt3Soma, GDNFOlink,  
               HGFSoma, IFNGammaSoma, IL5Soma, IL6Olink, IL6Soma, IL10Olink, IL10Soma, 
               IL10RBOlink, IL10RBSoma, IL12BOlink, IL17ASoma, IL18R1Olink, IL18R1Soma1, 
               IL18R1Soma2, IL22RA1Soma, LIFSoma, LIFROlink, LIFRSoma, MMP1Soma, 
               MMP10Olink, MMP10Soma, OSMSoma, SIRT2Olink, SIRT2Soma, SLAMF1Olink, 
               TGFaOlink, TNFSF14Olink, TNFSF14Soma, VEGFASoma1, VEGFASoma2, VEGFASoma3) 
GRSdf <- GRSdf[,c(2, seq(1, ncol(GRSdf), 2))]
colnames(GRSdf) <- c('ID','ADAOlink', 'AXIN1Olink', 'CASP8Olink', 'CCL11Soma', 'CCL19Olink', 'CCL19Soma', 
                     'CCL20Olink', 'CCL23Olink', 'CCL23Soma', 'CCL25Olink', 'CCL25Soma1', 'CCL25Soma2', 
                     'CD6Olink', 'CD244Olink', 'CDCP1Olink', 'CSF1Olink', 'CXCL1Soma', 'CXCL5Olink', 
                     'CXCL5Soma', 'CXCL1Olink', 'CXCL9Olink', 'CXCL10Olink', 'CXCL10Soma','CXCL11Soma',
                     'FGF5Olink', 'FGF19Olink', 'FGF19Soma', 'FGF21Olink', 'Flt3Soma', 'GDNFOlink',  
                     'HGFSoma', 'IFNGammaSoma', 'IL5Soma', 'IL6Olink', 'IL6Soma', 'IL10Olink', 'IL10Soma', 
                     'IL10RBOlink', 'IL10RBSoma', 'IL12BOlink', 'IL17ASoma', 'IL18R1Olink', 'IL18R1Soma1', 
                     'IL18R1Soma2', 'IL22RA1Soma', 'LIFSoma', 'LIFROlink', 'LIFRSoma', 'MMP1Soma', 
                     'MMP10Olink', 'MMP10Soma', 'OSMSoma', 'SIRT2Olink', 'SIRT2Soma', 'SLAMF1Olink', 
                     'TGFaOlink', 'TNFSF14Olink', 'TNFSF14Soma', 'VEGFASoma1', 'VEGFASoma2', 'VEGFASoma3')
rm(ADAOlink, AXIN1Olink, CASP8Olink, CCL11Soma, CCL19Olink, CCL19Soma, 
   CCL20Olink, CCL23Olink, CCL23Soma, CCL25Olink, CCL25Soma1, CCL25Soma2, 
   CD6Olink, CD244Olink, CDCP1Olink, CSF1Olink, CXCL1Soma, CXCL5Olink, 
   CXCL5Soma, CXCL1Olink, CXCL9Olink, CXCL10Olink, CXCL10Soma,CXCL11Soma,
   FGF5Olink, FGF19Olink, FGF19Soma, FGF21Olink, Flt3Soma, GDNFOlink,  
   HGFSoma, IFNGammaSoma, IL5Soma, IL6Olink, IL6Soma, IL10Olink, IL10Soma, 
   IL10RBOlink, IL10RBSoma, IL12BOlink, IL17ASoma, IL18R1Olink, IL18R1Soma1, 
   IL18R1Soma2, IL22RA1Soma, LIFSoma, LIFROlink, LIFRSoma, MMP1Soma, 
   MMP10Olink, MMP10Soma, OSMSoma, SIRT2Olink, SIRT2Soma, SLAMF1Olink, 
   TGFaOlink, TNFSF14Olink, TNFSF14Soma, VEGFASoma1, VEGFASoma2, VEGFASoma3)
GRSdf$ID <- gsub('M', '', GRSdf$ID)


# Merge data together
Vars_GRS <- merge(GRSdf, MothersVars, by.x = 'ID', by.y = 'ALN')


# Load in ARIES sample sheets 
#load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Samples.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/FOM1DNAm.RData")


# Match samples to episcores
VarFOM1_FOM1DNAm <- Vars_GRS[match(FOM1DNAm$ALN, Vars_GRS$ID),]
Vars_All <- cbind(VarFOM1_FOM1DNAm, t(EpiscoresFOM1))


NonGRSEffects <- function(Protein, MethScore, GRS, VarExp = F){
  
  if (VarExp == F){
  fit <- lm(Vars_All[[Protein]] ~ Vars_All[[MethScore]]) #CV training
  MethCor <- sqrt(summary(fit)$adj.r.squared)
  if (is.na(MethCor)){
    MethCor <- sqrt(abs(summary(fit)$adj.r.squared))
    MethCor <- MethCor * -1
  }
  sumfit <- summary(fit)
  MethVals <- c(MethCor, pf(sumfit$fstatistic[1], sumfit$fstatistic[2], sumfit$fstatistic[3], lower.tail = F))
  
  fit <- lm(Vars_All[[Protein]] ~ Vars_All[[GRS]]) #GRS
  GRSCor <- sqrt(summary(fit)$adj.r.squared)
  if (is.na(GRSCor)){
    GRSCor <- sqrt(abs(summary(fit)$adj.r.squared))
    GRSCor <- GRSCor * -1
  }
  sumfit <- summary(fit)
  GRSVals <- c(GRSCor, pf(sumfit$fstatistic[1], sumfit$fstatistic[2], sumfit$fstatistic[3], lower.tail = F))
  
  fit <- lm(Vars_All[[Protein]] ~ Vars_All[[MethScore]] + Vars_All[[GRS]]) #both
  ComCor <- sqrt(summary(fit)$adj.r.squared)
  if (is.na(ComCor)){
    ComCor <- sqrt(abs(summary(fit)$adj.r.squared))
    ComCor <- GRSCor * -1
  }  
  sumfit <- summary(fit)
  ComboVals <- c(ComCor, pf(sumfit$fstatistic[1], sumfit$fstatistic[2], sumfit$fstatistic[3], lower.tail = F))
  
  AOV <- anova(lm(Vars_All[[Protein]] ~ Vars_All[[GRS]]),
               lm(Vars_All[[Protein]] ~ Vars_All[[MethScore]] + Vars_All[[GRS]]))
  
  AOVVals <- AOV[2,5:6]
  
  N <- fit[["df.residual"]]+3
  
  return(as.data.frame(c(Protein, MethVals, GRSVals, ComboVals, AOVVals,N), col.names = c('Protein',
                                                                                        'MethCor', 'MethP',
                                                                                        'GRSCor', 'GRSP',
                                                                                        'ComCor', 'ComP',
                                                                                        'AOVF', 'AOVP',
                                                                                        'N')))
  }else{
    
    
    
    vp <- var(Vars_All[[Protein]], na.rm = T)
    vp_pgs <- vp - var(residuals(lm(Vars_All[[Protein]] ~ Vars_All[[GRS]])))
    vp_episcore <- vp - var(residuals(lm(Vars_All[[Protein]] ~ Vars_All[[MethScore]])))
    vp_both <- vp - var(residuals(lm(Vars_All[[Protein]] ~ Vars_All[[GRS]] + Vars_All[[MethScore]])))
    
    ve_pgs <- (vp_both - vp_episcore) / vp
    ve_episcore <- (vp_both - vp_pgs) / vp
    ve_both <- vp_both / vp
    
    
    
    
    return(data.frame('Protein' = Protein, 'VE_Episcore' = ve_episcore, 'VE_PGS' = ve_pgs, 'VE_Combined' = ve_both))
    
  }
}

# Gadd Episcores and Xu GRS effects in ALSPAC
FOM1_Gadd_GRS <- NonGRSEffects('CCL11_FOM1', 'CCL11_Olink', 'CCL11Soma')
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CCL25_FOM1', 'CCL25', 'CCL25Olink'))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CD6_FOM1', 'CD6_Olink', 'CD6Olink'))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CXCL10_FOM1', 'CXCL10', 'CXCL10Soma'))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CXCL11_FOM1', 'CXCL11', 'CXCL11Soma'))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('MMP1_FOM1', 'MMP1', 'MMP1Soma'))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('OSM_FOM1', 'OSM_Olink', 'OSMSoma'))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('VEGFA_FOM1', 'VEGFA_Olink', 'VEGFASoma3'))

write.csv(FOM1_Gadd_GRS, "Gadd_Xu_Correlations_FOM1.csv")

# Variance explained Gadd Episcores and Xu GRS effects in ALSPAC
FOM1_Gadd_GRS <- NonGRSEffects('CCL11_FOM1', 'CCL11_Olink', 'CCL11Soma',T)
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CCL25_FOM1', 'CCL25', 'CCL25Olink',T))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CD6_FOM1', 'CD6_Olink', 'CD6Olink',T))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CXCL10_FOM1', 'CXCL10', 'CXCL10Soma',T))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('CXCL11_FOM1', 'CXCL11', 'CXCL11Soma',T))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('MMP1_FOM1', 'MMP1', 'MMP1Soma',T))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('OSM_FOM1', 'OSM_Olink', 'OSMSoma',T))
FOM1_Gadd_GRS <- rbind(FOM1_Gadd_GRS, NonGRSEffects('VEGFA_FOM1', 'VEGFA_Olink', 'VEGFASoma3',T))
FOM1_Gadd_GRS[, 2:4] <- FOM1_Gadd_GRS[, 2:4] * 100
write.csv(FOM1_Gadd_GRS, 'Gadd_XU_PercentVarianceFOM1.csv')

#==============================================================================#
# Plot results
library(readr)
library(tidyr)

Gadd_XuComparisons <- read_csv("Gadd_Xu_Correlations_FOM1.csv")
GaddXU_VarExp <- read.csv('Gadd_XU_PercentVarianceFOM1.csv')
Gadd_XuComparisons$Protein <- as.factor(Gadd_XuComparisons$Protein)
test1 <- gather(Gadd_XuComparisons,  Model, Cor, MethCor, GRSCor, ComCor)
test2 <- gather(Gadd_XuComparisons,  Model, P, MethP, GRSP, ComP)
test3 <- gather(GaddXU_VarExp, Model, Variance, VE_Episcore, VE_PGS, VE_Combined)
resdiff <- test2
resdiff$Diff <- resdiff$ComCor - resdiff$GRSCor

Moddetails <- cbind(test1, test2, test3)
Moddetails$ModelName <- gsub('Cor', '', Moddetails$Model) 
#Moddetails <- Moddetails[c('Protein', 'ModelName', 'Cor', 'P', 'AOVP')]
Moddetails$ModelName <- replace(Moddetails$ModelName, Moddetails$ModelName=='Com', 'Combined')
Moddetails$ModelName <- as.factor(Moddetails$ModelName)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels=c('Meth', 'GRS', 'Combined'))

signames <- as.character(unique(Moddetails[Moddetails$AOVP<0.05,]$Protein))
sigtest <- ifelse( (Moddetails$ModelName=='Combined' & Moddetails$AOVP<0.05), '*', '')
Moddetails$Protein <- gsub('_FOM1', '', Moddetails$Protein)

Moddetails<- subset(Moddetails, select = !duplicated(names(Moddetails)))
library(ggplot2)
ggplot(Moddetails, aes(x=Protein, y=Variance, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Variance Explained (%)') + xlab('Protein') + labs(fill='Model') + 
  ggtitle('Middle Age') + 
  scale_fill_discrete(breaks = c('Meth', 'GRS', 'Combined'), labels = c('Gadd Episcore', 'Xu PGS', 'Combined')) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +
  geom_text(aes(label = sigtest),
            hjust = -2,
            color = "black",
            size = 15)

#------------------------------------------------------------------------------#
#Within CV training results

NonGRSEffects <- function(Protein, MethScore, GRS, VarExp=F){
  
  HoldDF <- Vars_All
  HoldDF <- HoldDF[c(Protein, MethScore, GRS)]
  HoldDF <- HoldDF[complete.cases(HoldDF),]
  
  if(!VarExp){
    fit <- lm(HoldDF[[Protein]] ~ HoldDF[[MethScore]]) #CV training
    MethCor <- sqrt(summary(fit)$adj.r.squared)
    if (is.na(MethCor)){
      MethCor <- sqrt(abs(summary(fit)$adj.r.squared))
      MethCor <- MethCor * -1
    }
    sumfit <- summary(fit)
    MethVals <- c(MethCor, pf(sumfit$fstatistic[1], sumfit$fstatistic[2], sumfit$fstatistic[3], lower.tail = F))
    
    fit <- lm(HoldDF[[Protein]] ~ HoldDF[[GRS]]) #GRS
    GRSCor <- sqrt(summary(fit)$adj.r.squared)
    if (is.na(GRSCor)){
      GRSCor <- sqrt(abs(summary(fit)$adj.r.squared))
      GRSCor <- GRSCor * -1
    }
    sumfit <- summary(fit)
    GRSVals <- c(GRSCor, pf(sumfit$fstatistic[1], sumfit$fstatistic[2], sumfit$fstatistic[3], lower.tail = F))
    
    fit <- lm(HoldDF[[Protein]] ~ HoldDF[[MethScore]] + HoldDF[[GRS]]) #both
    ComCor <- sqrt(summary(fit)$adj.r.squared)
    if (is.na(ComCor)){
      ComCor <- sqrt(abs(summary(fit)$adj.r.squared))
      ComCor <- GRSCor * -1
    }  
    sumfit <- summary(fit)
    ComboVals <- c(ComCor, pf(sumfit$fstatistic[1], sumfit$fstatistic[2], sumfit$fstatistic[3], lower.tail = F))
    
    AOV <- anova(lm(HoldDF[[Protein]] ~ HoldDF[[GRS]]),
                 lm(HoldDF[[Protein]] ~ HoldDF[[MethScore]] + HoldDF[[GRS]]))
    
    AOVVals <- AOV[2,5:6]
    
    N <- fit[["df.residual"]]+3
    
    return(as.data.frame(c(Protein, MethVals, GRSVals, ComboVals, AOVVals,N), col.names = c('Protein',
                                                                                            'MethCor', 'MethP',
                                                                                            'GRSCor', 'GRSP',
                                                                                            'ComCor', 'ComP',
                                                                                            'AOVF', 'AOVP',
                                                                                            'N')))
    
  }else{
    
    
    vp <- var(HoldDF[[Protein]], na.rm = T)
    vp_pgs <- vp - var(residuals(lm(HoldDF[[Protein]] ~ HoldDF[[GRS]])))
    vp_episcore <- vp - var(residuals(lm(HoldDF[[Protein]] ~ HoldDF[[MethScore]])))
    vp_both <- vp - var(residuals(lm(HoldDF[[Protein]] ~ HoldDF[[GRS]] + HoldDF[[MethScore]])))
    
    ve_pgs <- (vp_both - vp_episcore) / vp
    ve_episcore <- (vp_both - vp_pgs) / vp
    ve_both <- vp_both / vp
    
    
    
    
    return(data.frame('Protein' = Protein, 'VE_Episcore' = ve_episcore, 'VE_PGS' = ve_pgs, 'VE_Combined' = ve_both))
    
  }
  
  
}


library(dplyr)

load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PenRegFOM1ResidProt.Rdata")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/MothersARIES.RData")
load("MothersVars.RData")

FOM1DNAm <- FOM1DNAm$samples
MothersVars[MothersVars<0] <- NA

# Grab CV training predictions
FOM1CV <- as.data.frame(do.call(bind_rows, FOM1PenReg[[2]]))
rownames(FOM1CV) <- names(FOM1PenReg[[2]])

# CV training results
FOM1CVResults <- read_csv("FOM1TrainResults.csv")
colnames(FOM1CVResults) <- c('index', 'Cor', 'P', 'Protein')
FOM1CVResults <- FOM1CVResults[FOM1CVResults$Cor>=0.1 & FOM1CVResults$P<=0.05,]

# Make ARIES file same order as CV training
FOM1DNAm <- FOM1DNAm[FOM1DNAm$Sample_Name %in% colnames(FOM1CV),]
FOM1CV <- FOM1CV[, FOM1DNAm$Sample_Name]

# Sub sample names for alnqlet
colnames(FOM1CV) <- FOM1DNAm$ALN
rownames(FOM1CV) <- paste0(rownames(FOM1CV), '_CV')

# Match samples to episcores
VarFOM1_FOM1DNAm <- MothersVars[match(FOM1DNAm$ALN, MothersVars$aln),]
Vars_Gadd <- cbind(VarFOM1_FOM1DNAm, t(FOM1CV))

# Merge GRS data with Gadd data
Vars_All <- merge(Vars_Gadd, GRSdf, by.x = 'ALN', by.y = 'ID')

library(labelled)
Vars_All <- remove_labels(Vars_All)


FOM1_CV_GRS <- NonGRSEffects('CCL11_FOM1', 'CCL11_FOM1_CV', 'CCL11Soma')
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CCL19_FOM1', 'CCL19_FOM1_CV', 'CCL19Olink'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CCL20_FOM1', 'CCL20_FOM1_CV', 'CCL20Olink'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CD6_FOM1', 'CD6_FOM1_CV', 'CD6Olink'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CXCL10_FOM1', 'CXCL10_FOM1_CV', 'CXCL10Soma'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CXCL11_FOM1', 'CXCL11_FOM1_CV', 'CXCL11Soma'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CXCL9_FOM1', 'CXCL9_FOM1_CV', 'CXCL9Olink'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('Flt3L_FOM1', 'Flt3L_FOM1_CV', 'Flt3Soma'))
# FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('IFNgamma_FOM1', 'IFNgamma_FOM1_CV', 'IFNGammaSoma')) <- SNP missing
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('IL18R1_FOM1', 'IL18R1_FOM1_CV', 'IL18R1Soma1'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('IL6_FOM1', 'IL6_FOM1_CV', 'IL6Olink'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('MMP1_FOM1', 'MMP1_FOM1_CV', 'MMP1Soma'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('MMP10_FOM1', 'MMP10_FOM1_CV', 'MMP10Olink'))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('OSM_FOM1', 'OSM_FOM1_CV', 'OSMSoma'))

FOM1_CV_GRS$FDR <- p.adjust(FOM1_CV_GRS$AOVP, method = 'fdr')
FOM1_CV_GRS$Bonferroni <- p.adjust(FOM1_CV_GRS$AOVP, method = 'bonferroni')

FOM1_CV_GRS$Protein <- as.factor(FOM1_CV_GRS$Protein)
test1 <- gather(FOM1_CV_GRS,  Model, Cor, MethCor, GRSCor, ComCor)
test2 <- gather(FOM1_CV_GRS,  Model, P, MethP, GRSP, ComP)
resdiff <- test2
resdiff$Diff <- resdiff$ComCor - resdiff$GRSCor
resdiff <- resdiff[resdiff$Model == 'MethP',]

FOM1_CV_GRS <- NonGRSEffects('CCL11_FOM1', 'CCL11_FOM1_CV', 'CCL11Soma',T)
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CCL19_FOM1', 'CCL19_FOM1_CV', 'CCL19Olink',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CCL20_FOM1', 'CCL20_FOM1_CV', 'CCL20Olink',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CD6_FOM1', 'CD6_FOM1_CV', 'CD6Olink',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CXCL10_FOM1', 'CXCL10_FOM1_CV', 'CXCL10Soma',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CXCL11_FOM1', 'CXCL11_FOM1_CV', 'CXCL11Soma',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('CXCL9_FOM1', 'CXCL9_FOM1_CV', 'CXCL9Olink',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('Flt3L_FOM1', 'Flt3L_FOM1_CV', 'Flt3Soma',T))
# FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('IFNgamma_FOM1', 'IFNgamma_FOM1_CV', 'IFNGammaSoma')) <- SNP missing
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('IL18R1_FOM1', 'IL18R1_FOM1_CV', 'IL18R1Soma1',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('IL6_FOM1', 'IL6_FOM1_CV', 'IL6Olink',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('MMP1_FOM1', 'MMP1_FOM1_CV', 'MMP1Soma',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('MMP10_FOM1', 'MMP10_FOM1_CV', 'MMP10Olink',T))
FOM1_CV_GRS <- rbind(FOM1_CV_GRS, NonGRSEffects('OSM_FOM1', 'OSM_FOM1_CV', 'OSMSoma',T))

test3 <- gather(FOM1_CV_GRS, Model, Variance, VE_Episcore, VE_PGS, VE_Combined)

Moddetails <- cbind(test1, test2, test3)
Moddetails$Variance <- Moddetails$Variance *100
Moddetails$ModelName <- gsub('Cor', '', Moddetails$Model) 
#Moddetails <- Moddetails[c('Protein', 'ModelName', 'Cor', 'P', 'AOVP')]
Moddetails$ModelName <- replace(Moddetails$ModelName, Moddetails$ModelName=='Com', 'Combined')
Moddetails$ModelName <- as.factor(Moddetails$ModelName)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels=c('Meth', 'GRS', 'Combined'))

signames <- as.character(unique(Moddetails[Moddetails$AOVP<0.05,]$Protein))
sigtest <- ifelse( (Moddetails$ModelName=='Combined' & Moddetails$AOVP<0.05), '*', '')

Moddetails<- subset(Moddetails, select = !duplicated(names(Moddetails)))
library(ggplot2)
ggplot(Moddetails, aes(x=gsub('_FOM1', '', Protein), y=Variance, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Variance Explained (%)') + xlab('Protein') + labs(fill='Model') + 
  ggtitle('Middle Age') + 
  scale_fill_discrete(breaks = c('Meth', 'GRS', 'Combined'), labels = c('ALSPAC Episcore', 'Xu PGS', 'Combined')) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +
  geom_text(aes(label = sigtest),
            hjust = -1,
            color = "black",
            size = 15)


FivePercenters <- Moddetails[Moddetails$ModelName == 'Meth',]
FivePercenters <- FivePercenters[FivePercenters$AOVP<0.05 & FivePercenters$Variance>=5,]
as.character(FivePercenters$Protein)
