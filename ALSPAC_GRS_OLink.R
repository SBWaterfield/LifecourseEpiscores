

library(dplyr)
getwd()
setwd("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS")


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
  GRS <- subset(GRS,grepl("^.+(A|B)$",ID))
  
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
IL2RBSoma <- GRSbuilder(SNPs, 'IL2RB_SomaMod')
IL2Soma <- GRSbuilder(SNPs, 'IL2_SomaMod')
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
               HGFSoma, IFNGammaSoma, IL2Soma, IL2RBSoma, IL5Soma, IL6Olink, IL6Soma, IL10Olink, IL10Soma, 
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
                     'HGFSoma', 'IFNGammaSoma', 'IL2Soma', 'IL2RBSoma','IL5Soma', 'IL6Olink', 'IL6Soma', 'IL10Olink', 'IL10Soma', 
                     'IL10RBOlink', 'IL10RBSoma', 'IL12BOlink', 'IL17ASoma', 'IL18R1Olink', 'IL18R1Soma1', 
                     'IL18R1Soma2', 'IL22RA1Soma', 'LIFSoma', 'LIFROlink', 'LIFRSoma', 'MMP1Soma', 
                     'MMP10Olink', 'MMP10Soma', 'OSMSoma', 'SIRT2Olink', 'SIRT2Soma', 'SLAMF1Olink', 
                     'TGFaOlink', 'TNFSF14Olink', 'TNFSF14Soma', 'VEGFASoma1', 'VEGFASoma2', 'VEGFASoma3')
rm(ADAOlink, AXIN1Olink, CASP8Olink, CCL11Soma, CCL19Olink, CCL19Soma, 
   CCL20Olink, CCL23Olink, CCL23Soma, CCL25Olink, CCL25Soma1, CCL25Soma2, 
   CD6Olink, CD244Olink, CDCP1Olink, CSF1Olink, CXCL1Soma, CXCL5Olink, 
   CXCL5Soma, CXCL1Olink, CXCL9Olink, CXCL10Olink, CXCL10Soma,CXCL11Soma,
   FGF5Olink, FGF19Olink, FGF19Soma, FGF21Olink, Flt3Soma, GDNFOlink,  
   HGFSoma, IFNGammaSoma, IL2Soma, IL2RBSoma, IL5Soma, IL6Olink, IL6Soma, IL10Olink, IL10Soma, 
   IL10RBOlink, IL10RBSoma, IL12BOlink, IL17ASoma, IL18R1Olink, IL18R1Soma1, 
   IL18R1Soma2, IL22RA1Soma, LIFSoma, LIFROlink, LIFRSoma, MMP1Soma, 
   MMP10Olink, MMP10Soma, OSMSoma, SIRT2Olink, SIRT2Soma, SLAMF1Olink, 
   TGFaOlink, TNFSF14Olink, TNFSF14Soma, VEGFASoma1, VEGFASoma2, VEGFASoma3)
# Load in Olink data
load("Data/OlinkVars.RData")

# Merge data together
Vars_GRS <- merge(GRSdf, Vars, by.x = 'ID', by.y = 'aln_qlet')

# Check best performer for proteins with mutiple scores
cor.test(Vars_GRS$CCL19Olink, Vars_GRS$CCL19_F24) #0.13
cor.test(Vars_GRS$CCL19Soma, Vars_GRS$CCL19_F24) #0.06

cor.test(Vars_GRS$CCL23Olink, Vars_GRS$CCL23_F24) #0.36
cor.test(Vars_GRS$CCL23Soma, Vars_GRS$CCL23_F24) #0.33

cor.test(Vars_GRS$CCL25Olink, Vars_GRS$CCL25_F24) #0.59
cor.test(Vars_GRS$CCL25Soma1, Vars_GRS$CCL25_F24) #0.21
cor.test(Vars_GRS$CCL25Soma2, Vars_GRS$CCL25_F24) #0.49

cor.test(Vars_GRS$CXCL5Olink, Vars_GRS$CXCL5_F24) #0.37
cor.test(Vars_GRS$CXCL5Soma, Vars_GRS$CXCL5_F24) #0.21

cor.test(Vars_GRS$CXCL10Olink, Vars_GRS$CXCL10_F24) #0.32
cor.test(Vars_GRS$CXCL10Soma, Vars_GRS$CXCL10_F24) #0.33

cor.test(Vars_GRS$FGF19Olink, Vars_GRS$FGF19_F24) #0.11
cor.test(Vars_GRS$FGF19Soma, Vars_GRS$FGF19_F24) #0.10

cor.test(Vars_GRS$IL6Olink, Vars_GRS$IL6_F24) #0.09
cor.test(Vars_GRS$IL6Soma, Vars_GRS$IL6_F24) #0.002

cor.test(Vars_GRS$IL10Olink, Vars_GRS$IL10_F24) #0.08
cor.test(Vars_GRS$IL10Soma, Vars_GRS$IL10_F24) #0.02

cor.test(Vars_GRS$IL10RBOlink, Vars_GRS$IL10RB_F24) #0.32
cor.test(Vars_GRS$IL10RBSoma, Vars_GRS$IL10RB_F24) #0.02

cor.test(Vars_GRS$IL18R1Olink, Vars_GRS$IL18R1_F24) #0.53
cor.test(Vars_GRS$IL18R1Soma1, Vars_GRS$IL18R1_F24) #0.54
cor.test(Vars_GRS$IL18R1Soma2, Vars_GRS$IL18R1_F24) #0.53

cor.test(Vars_GRS$LIFROlink, Vars_GRS$LIFR_F24) #0.16
cor.test(Vars_GRS$LIFRSoma, Vars_GRS$LIFR_F24) #0.14

cor.test(Vars_GRS$MMP10Olink, Vars_GRS$MMP10_F24) #0.32
cor.test(Vars_GRS$MMP10Soma, Vars_GRS$MMP10_F24) #0.27

cor.test(Vars_GRS$SIRT2Olink, Vars_GRS$SIRT2_F24) #0.05
cor.test(Vars_GRS$SIRT2Soma, Vars_GRS$SIRT2_F24) #0.05

cor.test(Vars_GRS$TNFSF14Olink, Vars_GRS$TNFSF14_F24) #0.11
cor.test(Vars_GRS$TNFSF14Soma, Vars_GRS$TNFSF14_F24) #-0.02

cor.test(Vars_GRS$VEGFASoma1, Vars_GRS$VEGFA_F24) #0.05
cor.test(Vars_GRS$VEGFASoma2, Vars_GRS$VEGFA_F24) #0.05
cor.test(Vars_GRS$VEGFASoma3, Vars_GRS$VEGFA_F24) #0.06

# Load in ALSPAC Olink data
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/OlinkVars.RData")

# Load in ARIES Gadd episcores
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Episcores.Rdata")

# Load in ARIES sample sheets 
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Samples.RData")


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

# Match samples to episcores
Var24_TwentyFourDNAm <- Vars[match(TwentyFourDNAm$alnqlet, Vars$aln_qlet),]
EpiscoreTwentyFour <- data.frame(t(EpiscoreTwentyFour))
Vars_Gadd <- cbind(Var24_TwentyFourDNAm, EpiscoreTwentyFour)
# Merge GRS data with Gadd data
Vars_All <- merge(Vars_Gadd, GRSdf, by.x = 'aln_qlet', by.y = 'ID')
library(labelled)
Vars_All <- remove_labels(Vars_All)


# Gadd Episcores and Xu GRS effects in ALSPAC
F24_Gadd_GRS <- NonGRSEffects('CCL11_F24', 'CCL11_Olink', 'CCL11Soma')
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CCL25_F24', 'CCL25', 'CCL25Olink'))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CD6_F24', 'CD6_Olink', 'CD6Olink'))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CXCL10_F24', 'CXCL10', 'CXCL10Soma'))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CXCL11_F24', 'CXCL11', 'CXCL11Soma'))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('MMP1_F24', 'MMP1', 'MMP1Soma'))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('OSM_F24', 'OSM_Olink', 'OSMSoma'))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('VEGFA_F24', 'VEGFA_Olink', 'VEGFASoma3'))

write.csv(F24_Gadd_GRS, "Gadd_Xu_Correlations24.csv")

# # Variance explained Gadd Episcores and Xu GRS effects in ALSPAC
F24_Gadd_GRS <- NonGRSEffects('CCL11_F24', 'CCL11_Olink', 'CCL11Soma', T)
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CCL25_F24', 'CCL25', 'CCL25Olink',T))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CD6_F24', 'CD6_Olink', 'CD6Olink',T))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CXCL10_F24', 'CXCL10', 'CXCL10Soma',T))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('CXCL11_F24', 'CXCL11', 'CXCL11Soma',T))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('MMP1_F24', 'MMP1', 'MMP1Soma',T))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('OSM_F24', 'OSM_Olink', 'OSMSoma',T))
F24_Gadd_GRS <- rbind(F24_Gadd_GRS, NonGRSEffects('VEGFA_F24', 'VEGFA_Olink', 'VEGFASoma3',T))
F24_Gadd_GRS[, 2:4] <- F24_Gadd_GRS[, 2:4] * 100
write.csv(F24_Gadd_GRS, 'Gadd_XU_PercentVariance24.csv')
#==============================================================================#
#Age 9 results

# Match samples to episcores
Var9_NineDNAm <- Vars[match(NineDNAm$alnqlet, Vars$aln_qlet),]
EpiscoreNine <- data.frame(t(EpiscoreNine))
Vars_Gadd <- cbind(Var9_NineDNAm, EpiscoreNine)

# Merge GRS data with Gadd data
Vars_All <- merge(Vars_Gadd, GRSdf, by.x = 'aln_qlet', by.y = 'ID')

# Gadd Episcores and Xu GRS effects in ALSPAC
F9_Gadd_GRS <- NonGRSEffects('CCL11_F9', 'CCL11_Olink', 'CCL11Soma')
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('CCL25_F9', 'CCL25', 'CCL25Olink'))
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('CD6_F9', 'CD6_Olink', 'CD6Olink'))
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('CXCL10_F9', 'CXCL10', 'CXCL10Soma'))
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('CXCL11_F9', 'CXCL11', 'CXCL11Soma'))
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('MMP1_F9', 'MMP1', 'MMP1Soma'))
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('OSM_F9', 'OSM_Olink', 'OSMSoma'))
F9_Gadd_GRS <- rbind(F9_Gadd_GRS, NonGRSEffects('VEGFA_F9', 'VEGFA_Olink', 'VEGFASoma3'))

write.csv(F9_Gadd_GRS, "Gadd_Xu_CorrelationsF9.csv")
#------------------------------------------------------------------------------#
# plot results
library(tidyr)
library(ggplot2)
library(readr)
library(tidyverse)

Gadd_XuComparisons <- read_csv("Gadd_Xu_Correlations24.csv")
GaddXU_VarExp <- read.csv('Gadd_XU_PercentVariance24.csv')
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
Moddetails$Protein <- gsub('_F24', '', Moddetails$Protein)
signames <- as.character(unique(Moddetails[Moddetails$AOVP<0.05,]$Protein))
sigtest <- ifelse( (Moddetails$ModelName=='Combined' & Moddetails$AOVP<0.05), '*', '')

Moddetails<- subset(Moddetails, select = !duplicated(names(Moddetails)))
library(ggplot2)
ggplot(Moddetails, aes(x=Protein, y=Variance, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Variance Explained (%)') + xlab('Protein') + labs(fill='Model') + 
  ggtitle('Age 24') + 
  scale_fill_discrete(breaks = c('Meth', 'GRS', 'Combined'), labels = c('Gadd Episcore', 'Xu PGS', 'Combined')) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +
  geom_text(aes(label = sigtest),
            hjust = -2,
            color = "black",
            size = 15)

#==============================================================================#
Gadd_XuComparisons <- read_csv("Gadd_Xu_CorrelationsF9.csv")

Gadd_XuComparisons$Protein <- as.factor(Gadd_XuComparisons$Protein)
test1 <- gather(Gadd_XuComparisons,  Model, Cor, MethCor, GRSCor, ComCor)
test2 <- gather(Gadd_XuComparisons,  Model, P, MethP, GRSP, ComP)
Moddetails <- cbind(test1, test2)
Moddetails$ModelName <- gsub('Cor', '', Moddetails$Model) 
Moddetails <- Moddetails[c('Protein', 'ModelName', 'Cor', 'P', 'AOVP')]
Moddetails$ModelName <- replace(Moddetails$ModelName, Moddetails$ModelName=='Com', 'Combined')
Moddetails$ModelName <- as.factor(Moddetails$ModelName)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels=c('Meth', 'GRS', 'Combined'))

signames <- as.character(unique(Moddetails[Moddetails$AOVP<0.05,]$Protein))
sigtest <- ifelse( (Moddetails$ModelName=='Combined' & Moddetails$AOVP<0.05), '*', '')

library(ggplot2)
ggplot(Moddetails, aes(x=Protein, y=Cor, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') + labs(fill='Model') + 
  ylim(c(-0.1, 0.6)) + 
  ggtitle('Age 9') + 
  scale_fill_discrete(breaks = c('Meth', 'GRS', 'Combined'), labels = c('Gadd Episcore', 'Xu PGS', 'Combined')) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +
  geom_text(aes(label = sigtest),
            hjust = -0.5,
            color = "black",
            size = 15)

#==============================================================================#
#Within CV training models 
library(dplyr)

load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PenRegAge9ResidProtNoSexChr_FullSamples.Rdata")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Samples.RData")
load("Data/OlinkVars.RData")


# Grab CV training predictions
Age9CV <- as.data.frame(do.call(bind_rows, Age9PenReg[[3]]))
rownames(Age9CV) <- names(Age9PenReg[[3]])

# CV training results
library(readr)
Age9TrainResults <- read_csv("Age9CVTrainResults.csv")
Age9TrainResults <- Age9TrainResults[,-1]
colnames(Age9TrainResults) <- c('Cor', 'SD', 'Protein')
Age9TrainResults$Protein <- gsub('_F9', '', as.character(Age9TrainResults$Protein))
Age9TrainResults$Transfer <- 'CVTraining'
Age9TrainResults$P <- NA
Age9CVResults <- Age9TrainResults[Age9TrainResults$Cor>=0.1,]

# Make ARIES file same order as CV training
NineDNAm <- NineDNAm[NineDNAm$Sample_Name %in% colnames(Age9CV),]
Age9CV <- Age9CV[, NineDNAm$Sample_Name]

# Sub sample names for alnqlet
colnames(Age9CV) <- NineDNAm$alnqlet
rownames(Age9CV) <- paste0(rownames(Age9CV), '_CV')

# Match samples to episcores
Var9_NineDNAm <- Vars[match(NineDNAm$alnqlet, Vars$aln_qlet),]
Vars_Gadd <- cbind(Var9_NineDNAm, t(Age9CV))

# Merge GRS data with Gadd data
Vars_All <- merge(Vars_Gadd, GRSdf, by.x = 'aln_qlet', by.y = 'ID')

library(labelled)
Vars_All <- remove_labels(Vars_All)

# Check how genetic addition effects CV training score
Protein <- 'IL6_F9'
MethScore<- 'IL6_F9_CV' 
GRS <- 'IL6Olink'
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

Age9CV_GRS <- NonGRSEffects('IL6_F9', 'IL6_F9_CV', 'IL6Olink')
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL2_F9', 'IL2_F9_CV', 'IL2Soma'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL2RB_F9', 'IL2RB_F9_CV', 'IL2RBSoma'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL5_F9', 'IL5_F9_CV', 'IL5Soma'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('CXCL9_F9', 'CXCL9_F9_CV', 'CXCL9Olink'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('LIFR_F9', 'LIFR_F9_CV', 'LIFROlink'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL18R1_F9', 'IL18R1_F9_CV', 'IL18R1Soma1'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('MMP10_F9', 'MMP10_F9_CV', 'MMP10Olink'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('CCL23_F9', 'CCL23_F9_CV', 'CCL23Olink'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('LIF_F9', 'LIF_F9_CV', 'LIFSoma'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('TGFalpha_F9', 'TGFalpha_F9_CV', 'TGFaOlink'))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('ADA_F9', 'ADA_F9_CV', 'ADAOlink'))


Age9CV_GRS$Protein <- as.factor(Age9CV_GRS$Protein)
test1 <- gather(Age9CV_GRS,  Model, Cor, MethCor, GRSCor, ComCor)
test2 <- gather(Age9CV_GRS,  Model, P, MethP, GRSP, ComP)

Age9CV_GRS <- NonGRSEffects('IL6_F9', 'IL6_F9_CV', 'IL6Olink',T)
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL2_F9', 'IL2_F9_CV', 'IL2Soma',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL2RB_F9', 'IL2RB_F9_CV', 'IL2RBSoma',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL5_F9', 'IL5_F9_CV', 'IL5Soma',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('CXCL9_F9', 'CXCL9_F9_CV', 'CXCL9Olink',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('LIFR_F9', 'LIFR_F9_CV', 'LIFROlink',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('IL18R1_F9', 'IL18R1_F9_CV', 'IL18R1Soma1',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('MMP10_F9', 'MMP10_F9_CV', 'MMP10Olink',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('CCL23_F9', 'CCL23_F9_CV', 'CCL23Olink',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('LIF_F9', 'LIF_F9_CV', 'LIFSoma',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('TGFalpha_F9', 'TGFalpha_F9_CV', 'TGFaOlink',T))
Age9CV_GRS <- rbind(Age9CV_GRS, NonGRSEffects('ADA_F9', 'ADA_F9_CV', 'ADAOlink',T))
test3 <- gather(Age9CV_GRS, Model, Variance, VE_Episcore, VE_PGS, VE_Combined)

resdiff <- test2
resdiff$Diff <- resdiff$ComCor - resdiff$GRSCor

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
ggplot(Moddetails, aes(x=gsub('_F9', '', Protein), y=Variance, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Variance Explained (%)') + xlab('Protein') + labs(fill='Model') + 
  ggtitle('Age 9') + 
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


# Age 9 CV train explaining  Age 24 Proteins
Age9CVTranssfer_GRS <- NonGRSEffects('IL6_F24', 'IL6_F9_CV', 'IL6Olink')
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('IL2_F24', 'IL2_F9_CV', 'IL2Soma'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('IL2RB_F24', 'IL2RB_F9_CV', 'IL2RBSoma'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('IL5_F24', 'IL5_F9_CV', 'IL5Soma'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('CXCL9_F24', 'CXCL9_F9_CV', 'CXCL9Olink'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('LIFR_F24', 'LIFR_F9_CV', 'LIFROlink'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('IL18R1_F24', 'IL18R1_F9_CV', 'IL18R1Soma1'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('MMP10_F24', 'MMP10_F9_CV', 'MMP10Olink'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('CCL23_F24', 'CCL23_F9_CV', 'CCL23Olink'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('LIF_F24', 'LIF_F9_CV', 'LIFSoma'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('TGFalpha_F24', 'TGFalpha_F9_CV', 'TGFaOlink'))
Age9CVTranssfer_GRS <- rbind(Age9CVTranssfer_GRS, NonGRSEffects('ADA_F24', 'ADA_F9_CV', 'ADAOlink'))
#------------------------------------------------------------------------------#
# Grab CV training predictions
library(readr)
library(tidyverse)
Age24TrainResults <- read_csv("Age24TrainResults.csv")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/PenRegAge24ResidProtNoSexChr_FullSamples.Rdata")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Samples.RData")
load("Data/OlinkVars.RData")

Age24CV <- as.data.frame(do.call(bind_rows, Age24PenReg[[3]]))
rownames(Age24CV) <- names(Age24PenReg[[3]])

# CV training results
Age24CVResults <- Age24PenReg[[2]]
colnames(Age24CVResults) <- c('Cor', 'P', 'Protein')
Age24CVResults$P <- as.numeric(Age24CVResults$P)
Age24CVResults$Cor <- as.numeric(Age24CVResults$Cor)
Age24CVResults <- Age24CVResults[Age24CVResults$Cor>=0.1 & Age24CVResults$P<=0.05,]

# Make ARIES file same order as CV training
TwentyFourDNAm <- TwentyFourDNAm[TwentyFourDNAm$Sample_Name %in% colnames(Age24CV),]
Age24CV <- Age24CV[, TwentyFourDNAm$Sample_Name]

# Sub sample names for alnqlet
colnames(Age24CV) <- TwentyFourDNAm$alnqlet
rownames(Age24CV) <- paste0(rownames(Age24CV), '_CV')

# Match samples to episcores
Var24_TwentyFourDNAm <- Vars[match(TwentyFourDNAm$alnqlet, Vars$aln_qlet),]
Vars_Gadd <- cbind(Var24_TwentyFourDNAm, t(Age24CV))

# Merge GRS data with Gadd data
Vars_All <- merge(Vars_Gadd, GRSdf, by.x = 'aln_qlet', by.y = 'ID')

library(labelled)
Vars_All <- remove_labels(Vars_All)


Age24CV_GRS <- NonGRSEffects('ADA_F24', 'ADA_F24_CV', 'ADAOlink')
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('AXIN1_F24', 'AXIN1_F24_CV', 'AXIN1Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CCL11_F24', 'CCL11_F24_CV', 'CCL11Soma'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CCL19_F24', 'CCL19_F24_CV', 'CCL19Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CCL25_F24', 'CCL25_F24_CV', 'CCL25Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CD6_F24', 'CD6_F24_CV', 'CD6Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CDCP1_F24', 'CDCP1_F24_CV', 'CDCP1Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CSF1_F24', 'CSF1_F24_CV', 'CSF1Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CXCL10_F24', 'CXCL10_F24_CV', 'CXCL10Soma'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CXCL11_F24', 'CXCL11_F24_CV', 'CXCL11Soma'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CXCL9_F24', 'CXCL9_F24_CV', 'CXCL9Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('FGF21_F24', 'FGF21_F24_CV', 'FGF21Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('FGF5_F24', 'FGF5_F24_CV', 'FGF5Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('Flt3L_F24', 'Flt3L_F24_CV', 'Flt3Soma'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL17A_F24', 'IL17A_F24_CV', 'IL17ASoma'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL18R1_F24', 'IL18R1_F24_CV', 'IL18R1Soma1'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL6_F24', 'IL6_F24_CV', 'IL6Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL10_F24', 'IL10_F24_CV', 'IL10Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('OSM_F24', 'OSM_F24_CV', 'OSMSoma'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('SLAMF1_F24', 'SLAMF1_F24_CV', 'SLAMF1Olink'))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('TNFSF14_F24', 'TNFSF14_F24_CV', 'TNFSF14Olink'))



Age24CV_GRS$Protein <- as.factor(Age24CV_GRS$Protein)
test1 <- gather(Age24CV_GRS,  Model, Cor, MethCor, GRSCor, ComCor)
test2 <- gather(Age24CV_GRS,  Model, P, MethP, GRSP, ComP)
resdiff <- test2
resdiff$Diff <- resdiff$ComCor - resdiff$GRSCor
resdiff <- resdiff[resdiff$Model == 'MethP',]

Age24CV_GRS <- NonGRSEffects('ADA_F24', 'ADA_F24_CV', 'ADAOlink',T)
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('AXIN1_F24', 'AXIN1_F24_CV', 'AXIN1Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CCL11_F24', 'CCL11_F24_CV', 'CCL11Soma',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CCL19_F24', 'CCL19_F24_CV', 'CCL19Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CCL25_F24', 'CCL25_F24_CV', 'CCL25Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CD6_F24', 'CD6_F24_CV', 'CD6Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CDCP1_F24', 'CDCP1_F24_CV', 'CDCP1Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CSF1_F24', 'CSF1_F24_CV', 'CSF1Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CXCL10_F24', 'CXCL10_F24_CV', 'CXCL10Soma',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CXCL11_F24', 'CXCL11_F24_CV', 'CXCL11Soma',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('CXCL9_F24', 'CXCL9_F24_CV', 'CXCL9Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('FGF21_F24', 'FGF21_F24_CV', 'FGF21Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('FGF5_F24', 'FGF5_F24_CV', 'FGF5Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('Flt3L_F24', 'Flt3L_F24_CV', 'Flt3Soma',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL17A_F24', 'IL17A_F24_CV', 'IL17ASoma',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL18R1_F24', 'IL18R1_F24_CV', 'IL18R1Soma1',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL6_F24', 'IL6_F24_CV', 'IL6Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('IL10_F24', 'IL10_F24_CV', 'IL10Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('OSM_F24', 'OSM_F24_CV', 'OSMSoma',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('SLAMF1_F24', 'SLAMF1_F24_CV', 'SLAMF1Olink',T))
Age24CV_GRS <- rbind(Age24CV_GRS, NonGRSEffects('TNFSF14_F24', 'TNFSF14_F24_CV', 'TNFSF14Olink',T))
test3 <- gather(Age24CV_GRS, Model, Variance, VE_Episcore, VE_PGS, VE_Combined)

resdiff <- test2
resdiff$Diff <- resdiff$ComCor - resdiff$GRSCor

Moddetails <- cbind(test1, test2, test3)
Moddetails$Variance <- Moddetails$Variance *100
Moddetails$ModelName <- gsub('Cor', '', Moddetails$Model) 
#Moddetails <- Moddetails[c('Protein', 'ModelName', 'Cor', 'P', 'AOVP')]
Moddetails$ModelName <- replace(Moddetails$ModelName, Moddetails$ModelName=='Com', 'Combined')
Moddetails$ModelName <- as.factor(Moddetails$ModelName)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels=c('Meth', 'GRS', 'Combined'))


Moddetails$Protein <- gsub('_F24', '', Moddetails$Protein)


signames <- as.character(unique(Moddetails[Moddetails$AOVP<0.05,]$Protein))
sigtest <- ifelse( (Moddetails$ModelName=='Combined' & Moddetails$AOVP<0.05), '*', '')

Moddetails<- subset(Moddetails, select = !duplicated(names(Moddetails)))
library(ggplot2)
ggplot(Moddetails, aes(x=gsub('_F24', '', Protein), y=Variance, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Variance Explained (%)') + xlab('Protein') + labs(fill='Model') + 
  ggtitle('Age 24') + 
  scale_fill_discrete(breaks = c('Meth', 'GRS', 'Combined'), labels = c('ALSPAC Episcore', 'Xu PGS', 'Combined')) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +
  geom_text(aes(label = sigtest),
            hjust = -0.5,
            color = "black",
            size = 15)


FivePercenters <- Moddetails[Moddetails$ModelName == 'Meth',]
FivePercenters <- FivePercenters[FivePercenters$AOVP<0.05 & FivePercenters$Variance>=5,]
as.character(FivePercenters$Protein)


#------------------------------------------------------------------------------#









