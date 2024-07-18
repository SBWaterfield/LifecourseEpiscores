library(tidyverse)
library(ggplot2)
library(haven)

# Load in DNAm and Olink data
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIESDNAm.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/MothersARIES.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/Data/OlinkVars.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/Data/OlinkMothersVars.Rdata")

# Match variables
Var9_NineDNAm <- Vars[match(NineDNAm$samples$alnqlet, Vars$aln_qlet),]
Var9_NineDNAm[Var9_NineDNAm < 0] <- NA
Var24_TwentyFourDNAm <- Vars[match(TwentyFourDNAm$samples$alnqlet, Vars$aln_qlet),]
Var24_TwentyFourDNAm[Var24_TwentyFourDNAm < 0] <- NA
VarFOM1_FOM1DNAm <- OlinkMothers[match(FOM1DNAm$samples$ALN, OlinkMothers$ALN),]
VarFOM1_FOM1DNAm[VarFOM1_FOM1DNAm < 0] <- NA

FOM1DNAm <- FOM1DNAm$samples
rm(NineDNAm, TwentyFourDNAm)

VarFOM1_FOM1DNAm <- zap_labels(VarFOM1_FOM1DNAm) #Haven labels...

# Load in EPiscore Results
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ALSPACEpsicoresTransferred.RData")

# Load in Training CV results
library(readr)
Age9TrainResults <- read_csv("Age9CVTrainResults.csv")
Age9TrainResults <- Age9TrainResults[,-1]
colnames(Age9TrainResults) <- c('Cor', 'SD', 'Protein')
Age9TrainResults$Protein <- gsub('_F9', '', as.character(Age9TrainResults$Protein))
Age9TrainResults$Transfer <- 'CVTraining'
Age9TrainResults$P <- NA


Age24TrainResults <- read_csv("Age24CVTrainResults.csv")
Age24TrainResults <- Age24TrainResults[,-1]
colnames(Age24TrainResults) <- c('Cor', 'SD', 'Protein')
Age24TrainResults$Protein <- gsub('_F24', '', as.character(Age24TrainResults$Protein))
Age24TrainResults$Transfer <- 'CVTraining'
Age24TrainResults$P <- NA

library(readr)
# Age24TrainResults450k <- read_csv("Age24TrainResults450kALL.csv")
# Age24TrainResults <- read_csv("Age24TrainResults.csv")
# Age24TrainResults <- Age24TrainResults[,-1]
# colnames(Age24TrainResults) <- c('Cor', 'P', 'Protein')
# Age24TrainResults$Protein <- gsub('_F24', '', as.character(Age24TrainResults$Protein))
# Age24TrainResults$Transfer <- 'CVTraining'
# Age24TrainResults$SampleN <- 822

# 
# Age9MaleTrainResults <- read_csv("Age9MaleTrainResults.csv")
# Age9MaleTrainResults <- Age9MaleTrainResults[,-1]
# colnames(Age9MaleTrainResults) <- c('Cor', 'P', 'Protein')
# Age9MaleTrainResults$Protein <- gsub('_F9', '', as.character(Age9MaleTrainResults$Protein))
# Age9MaleTrainResults$Transfer <- 'CVTraining'
# 
# Age24MaleTrainResults <- read_csv("Age24MaleTrainResults.csv")
# Age24MaleTrainResults <- Age24MaleTrainResults[,-1]
# colnames(Age24MaleTrainResults) <- c('Cor', 'P', 'Protein')
# Age24MaleTrainResults$Protein <- gsub('_F24', '', as.character(Age24MaleTrainResults$Protein))
# Age24MaleTrainResults$Transfer <- 'CVTraining'
# 
# Age9FemaleTrainResults <- read_csv("Age9FemaleTrainResults.csv")
# Age9FemaleTrainResults <- Age9FemaleTrainResults[,-1]
# colnames(Age9FemaleTrainResults) <- c('Cor', 'P', 'Protein')
# Age9FemaleTrainResults$Protein <- gsub('_F9', '', as.character(Age9FemaleTrainResults$Protein))
# Age9FemaleTrainResults$Transfer <- 'CVTraining'
# 
# Age24FemaleTrainResults <- read_csv("Age24FemaleTrainResults.csv")
# Age24FemaleTrainResults <- Age24FemaleTrainResults[,-1]
# colnames(Age24FemaleTrainResults) <- c('Cor', 'P', 'Protein')
# Age24FemaleTrainResults$Protein <- gsub('_F24', '', as.character(Age24FemaleTrainResults$Protein))
# Age24FemaleTrainResults$Transfer <- 'CVTraining'

FOM1TrainResults <- read_csv("FOM1CVTrainResults.csv")
FOM1TrainResults <- FOM1TrainResults[,-1]
colnames(FOM1TrainResults) <- c('Cor', 'SD', 'Protein')
FOM1TrainResults$Protein <- gsub('_FOM1', '', as.character(FOM1TrainResults$Protein))
FOM1TrainResults$Transfer <- 'CVTraining'
FOM1TrainResults$P <- NA

#==============================================================================#
# Sort out sample size for each protein

SampleSizes <- function(varmatch, trainresults, suffix){
  
  varmatch <- select(varmatch, contains(suffix))
  df <- data.frame()
  
  for (id in colnames(varmatch)){
    
    vals <- unlist(varmatch[,id])
    N <- sum(!is.na(vals))
    df <- rbind(df, c(gsub(suffix, '', id), N))
    
  }
  
  colnames(df) <- c('Protein', 'SampleN')
  
  df <- left_join(df, trainresults)
  
  return(df)
}


FOM1TrainResults <- SampleSizes(VarFOM1_FOM1DNAm, FOM1TrainResults, '_FOM1')
Age24TrainResults <- SampleSizes(Var24_TwentyFourDNAm, Age24TrainResults, '_F24')
Age9TrainResults <- SampleSizes(Var9_NineDNAm, Age9TrainResults, '_F9')


#==============================================================================#

# Load in ARIES Sample sheets
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Samples.RData")

rm(BirthCell, BirthDNAm, SevenCell, SevenDNAm, NineCell, FifteenCell, FifteenDNAm, TwentyFourCell, Age24_p1, Age24_p2)



EpiscoreCor <- function(Episcores, VarMatch, VarAge, TrainSuffix, TestSuffix, SexSpecific=NULL, Mod=NULL, CV=F){
  
  Resdf <- data.frame()
  
  if(is.null(SexSpecific)){
    for (protein in rownames(Episcores)){
      
      ProteinName <- paste0(gsub(TrainSuffix, "", protein), TestSuffix)
      
      if (ProteinName %in% names(VarMatch)){
      ProtVals <- VarMatch[,ProteinName]
      } else{next}
      ProtVals <- unlist(ProtVals)
      CorrelationResult <- cor.test(t(Episcores[protein,]), ProtVals)
      CorrelationResult <- c(gsub(TrainSuffix, "", protein), CorrelationResult[["estimate"]], CorrelationResult[["p.value"]], CorrelationResult[["parameter"]][["df"]]+2)
      Resdf <- rbind(Resdf, CorrelationResult)
    }
    
    if (!CV){
      Resdf$SD <-NA
    }
    
    colnames(Resdf) <- c('Protein', 'Cor', 'P', 'SampleN', 'SD')
    Resdf$Cor <- as.numeric(Resdf$Cor)
    Resdf$Cor <- signif(Resdf$Cor, 2)
    Resdf$P <- as.numeric(Resdf$P)
    Resdf$P <- signif(Resdf$P, 2)
    Resdf$Transfer <- Mod
    
  }else{
    Sexmatch <- VarAge$Sex == SexSpecific

    for (protein in rownames(Episcores)){
      
      ProteinName <- paste0(gsub(TrainSuffix, "", protein), TestSuffix)
      
      if (!ProteinName %in% names(VarMatch)){next}
      if (sum(!is.na(VarMatch[,ProteinName][Sexmatch])) <10) {next}
      CorrelationResult <- cor.test(t(Episcores[protein,])[Sexmatch], VarMatch[,ProteinName][Sexmatch])
      CorrelationResult <- c(gsub(TrainSuffix, "", protein), CorrelationResult[["estimate"]], CorrelationResult[["p.value"]], CorrelationResult[["parameter"]][["df"]]+2)
      Resdf <- rbind(Resdf, CorrelationResult)
    }
    if (!CV){
      Resdf$SD <-NA
    }
    
    colnames(Resdf) <- c('Protein', 'Cor', 'P', 'SampleN', 'SD')
    Resdf$Cor <- as.numeric(Resdf$Cor)
    Resdf$Cor <- signif(Resdf$Cor, 2)
    Resdf$P <- as.numeric(Resdf$P)
    Resdf$P <- signif(Resdf$P, 2)
    Resdf$Transfer <- Mod
    
    
  }
  
 return(Resdf)  
  
  
}

# Checking transferability of age 24 episcores
Age24EpiscoresResults <-  EpiscoreCor(Age24EpiscoresinAge9, Var9_NineDNAm, NineDNAm,  '_F24', '_F9', Mod = 'Age9')
#Age24EpiscoresResults <- rbind(Age24EpiscoresResults, EpiscoreCor(Age24EpiscoresinAge9, Var9_NineDNAm, NineDNAm,  '_F24', '_F9', 'F', 'Age9Females'))
#Age24EpiscoresResults <- rbind(Age24EpiscoresResults, EpiscoreCor(Age24EpiscoresinAge9, Var9_NineDNAm, NineDNAm,  '_F24', '_F9', 'M', 'Age9Males'))
Age24EpiscoresResults <- rbind(Age24EpiscoresResults, EpiscoreCor(Age24EpiscoresinFOM1, VarFOM1_FOM1DNAm, FOM1DNAm,  '_F24', '_FOM1', Mod = 'Mothers'))
Age24EpiscoresResults <- rbind(Age24EpiscoresResults, Age24TrainResults)
Age24Episcores <- unique(paste0(Age24EpiscoresResults$Protein, '_F24'))
save(Age24EpiscoresResults, file = 'Age24ALSPACEpiscoresFullResults.RData')
protsigs <- unique(Age24EpiscoresResults[Age24EpiscoresResults$P<0.05 & Age24EpiscoresResults$Cor>=0.1 & !Age24EpiscoresResults$Transfer == 'CVTraining',]$Protein)
protCVsigs <- unique(Age24EpiscoresResults[Age24EpiscoresResults$Cor>=0.1 & Age24EpiscoresResults$Transfer == 'CVTraining',]$Protein)
protCVsigs <- intersect(protsigs, protCVsigs)
Age24EpiscoresResults <- Age24EpiscoresResults[Age24EpiscoresResults$Protein %in% protCVsigs,]

Age24EpiscoresResults$Transfer <- ifelse(Age24EpiscoresResults$Transfer == 'CVTraining', 'Age24 (CV)', Age24EpiscoresResults$Transfer)


ggplot(Age24EpiscoresResults, aes(x=Protein, y=Cor, fill = factor(Transfer, levels = c('Age9', 'Age24', 'Mothers', 'Mothers (CV)', 'Age9 (CV)', 'Age24 (CV)')), width=0.6)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin = Cor - SD, ymax = Cor + SD), width = 0.2, position = position_dodge(0.6)) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') + ggtitle('Age 24') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +   labs(fill='Transfer') +
  scale_fill_manual(values = c(Age9 = "#F8766D", Age24 = "#00BA38", Mothers = "#619CFF", CVTraining = 'grey', 
                               'Mothers (CV)' = 'grey', 'Age9 (CV)' = 'grey','Age24 (CV)' = 'grey'))

Age24SigEpiscores <- unique(paste0(Age24EpiscoresResults$Protein, '_F24'))


# Checking transferability of age 9 episcores 
Age9EpiscoresResults <-  EpiscoreCor(Age9EpiscoresinAge24, Var24_TwentyFourDNAm, TwentyFourDNAm,  '_F9', '_F24', Mod = 'Age24')
#Age9EpiscoresResults <- rbind(Age9EpiscoresResults, EpiscoreCor(Age9EpiscoresinAge24, Var24_TwentyFourDNAm, TwentyFourDNAm,  '_F9', '_F24', 'F', 'Age24Females'))
#Age9EpiscoresResults <- rbind(Age9EpiscoresResults, EpiscoreCor(Age9EpiscoresinAge24, Var24_TwentyFourDNAm, TwentyFourDNAm,  '_F9', '_F24', 'M', 'Age24Males'))
Age9EpiscoresResults <- rbind(Age9EpiscoresResults, EpiscoreCor(Age9EpiscoresinFOM1, VarFOM1_FOM1DNAm, FOM1DNAm,  '_F9', '_FOM1', Mod = 'Mothers'))
Age9EpiscoresResults <- rbind(Age9EpiscoresResults, Age9TrainResults)
Age9EpiscoresResults$Cor <- as.numeric(Age9EpiscoresResults$Cor)
Age9EpiscoresResults$P <- as.numeric(Age9EpiscoresResults$P)
Age9Episcores <- unique(paste0(Age9EpiscoresResults$Protein, '_F9'))
save(Age9EpiscoresResults, file = 'Age9ALSPACEpiscoresFullResults.RData')

protsigs <- unique(Age9EpiscoresResults[Age9EpiscoresResults$P<0.05 & Age9EpiscoresResults$Cor>=0.1 & !Age9EpiscoresResults$Transfer == 'CVTraining' ,]$Protein)
protCVsigs <- unique(Age9EpiscoresResults[Age9EpiscoresResults$Cor>=0.1 & Age9EpiscoresResults$Transfer == 'CVTraining',]$Protein)
protCVsigs <- intersect(protsigs, protCVsigs)
Age9EpiscoresResults <- Age9EpiscoresResults[Age9EpiscoresResults$Protein %in% protCVsigs,]

Age9EpiscoresResults$Transfer <- ifelse(Age9EpiscoresResults$Transfer == 'CVTraining', 'Age9 (CV)', Age9EpiscoresResults$Transfer)


ggplot(Age9EpiscoresResults, aes(x=Protein, y=Cor, fill = factor(Transfer, levels = c('Age9', 'Age24', 'Mothers', 'Mothers (CV)', 'Age9 (CV)', 'Age24 (CV)')), width=0.6)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin = Cor - SD, ymax = Cor + SD), width = 0.2, position = position_dodge(0.6)) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') + ggtitle('Age 9') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  +  labs(fill='Transfer') +
  scale_fill_manual(values = c(Age9 = "#F8766D", Age24 = "#00BA38", Mothers = "#619CFF", CVTraining = 'grey', 
                               'Mothers (CV)' = 'grey', 'Age9 (CV)' = 'grey','Age24 (CV)' = 'grey'))
Age9SigEpiscores <- unique(paste0(Age9EpiscoresResults$Protein, '_F9'))


# Checking transferability of mothers episcores
FOM1EpiscoresResults <-  EpiscoreCor(FOM1EpiscoresinAge24, Var24_TwentyFourDNAm, TwentyFourDNAm,  '_FOM1', '_F24', Mod = 'Age24')
FOM1EpiscoresResults <- rbind(FOM1EpiscoresResults, EpiscoreCor(FOM1EpiscoresinAge24, Var24_TwentyFourDNAm, TwentyFourDNAm,  '_FOM1', '_F24', 'F', 'Age24Females'))
FOM1EpiscoresResults <- rbind(FOM1EpiscoresResults, EpiscoreCor(FOM1EpiscoresinAge24, Var24_TwentyFourDNAm, TwentyFourDNAm,  '_FOM1', '_F24', 'M', 'Age24Males'))
FOM1EpiscoresResults <- rbind(FOM1EpiscoresResults, EpiscoreCor(FOM1EpiscoresinAge9, Var9_NineDNAm, NineDNAm,  '_FOM1', '_F9', Mod = 'Age9'))
FOM1EpiscoresResults <- rbind(FOM1EpiscoresResults, EpiscoreCor(FOM1EpiscoresinAge9, Var9_NineDNAm, NineDNAm,  '_FOM1', '_F9', 'F', 'Age9Females'))
FOM1EpiscoresResults <- rbind(FOM1EpiscoresResults, EpiscoreCor(FOM1EpiscoresinAge9, Var9_NineDNAm, NineDNAm,  '_FOM1', '_F9', 'M', 'Age9Males'))
FOM1EpiscoresResults <- rbind(FOM1EpiscoresResults, FOM1TrainResults)
FOM1Episcores <- unique(paste0(FOM1EpiscoresResults$Protein, '_FOM1'))
save(FOM1EpiscoresResults, file = 'AgeFOM1ALSPACEpiscoresFullResults.RData')

#protsigs <- unique(FOM1EpiscoresResults[FOM1EpiscoresResults$P<0.05 & FOM1EpiscoresResults$Cor>=0.1 & !FOM1EpiscoresResults$Transfer == 'CVTraining' ,]$Protein)
#FOM1EpiscoresResults <- FOM1EpiscoresResults[FOM1EpiscoresResults$Protein %in% protsigs,]

FOM1EpiscoresResults <- FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer %in% c('Age9', 'Age24', 'CVTraining'),]
protsigs <- unique(FOM1EpiscoresResults[FOM1EpiscoresResults$P<0.05 & FOM1EpiscoresResults$Cor>=0.1 & !FOM1EpiscoresResults$Transfer == 'CVTraining' ,]$Protein)
protCVsigs <- unique(FOM1EpiscoresResults[FOM1EpiscoresResults$Cor>=0.1 & FOM1EpiscoresResults$Transfer == 'CVTraining',]$Protein)
protCVsigs <- intersect(protsigs, protCVsigs)
FOM1EpiscoresResults <- FOM1EpiscoresResults[FOM1EpiscoresResults$Protein %in% protCVsigs,]

FOM1EpiscoresResults$Transfer <- ifelse(FOM1EpiscoresResults$Transfer == 'CVTraining', 'Mothers (CV)', FOM1EpiscoresResults$Transfer)

ggplot(FOM1EpiscoresResults, aes(x=Protein, y=Cor, fill = factor(Transfer, levels = c('Age9', 'Age24', 'Mothers', 'Mothers (CV)', 'Age9 (CV)', 'Age24 (CV)')), width=0.6)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_errorbar(aes(ymin = Cor - SD, ymax = Cor + SD), width = 0.2, position = position_dodge(0.6)) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') + ggtitle('Middle Age') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=25), axis.title=element_text(size=35),
        legend.text=element_text(size=25), legend.title=element_text(size=35),
        plot.title = element_text(hjust = 0.5, size=25))  + labs(fill='Transfer') +
  scale_fill_manual(values = c(Age9 = "#F8766D", Age24 = "#00BA38", Mothers = "#619CFF", CVTraining = 'grey', 
                               'Mothers (CV)' = 'grey', 'Age9 (CV)' = 'grey','Age24 (CV)' = 'grey'))

FOM1SigEpiscores <- unique(paste0(FOM1EpiscoresResults$Protein, '_FOM1'))


load('AgeFOM1ALSPACEpiscoresFullResults.RData') 
load('Age9ALSPACEpiscoresFullResults.RData') 
load('Age24ALSPACEpiscoresFullResults.RData')
#save(FOM1EpiscoresResults, Age9EpiscoresResults, Age24EpiscoresResults, file = 'TransferResultsPlaceholder.RData')

#===============================================================================#
# N ranges for transfers





#------------------------------------------------------------------------------#
# Is there a difference in transferability of performance by sex?

# All results
load('TransferResultsPlaceholder.RData')
MothersSexTransfer <- FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer %in% c('Age24Males', 'Age24Females'),]
summary(lm(Cor ~ Transfer, data = MothersSexTransfer))
xtest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age24Males',]$Cor
ytest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age24Females',]$Cor
t.test(xtest, ytest)

MothersSexTransfer <- FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer %in% c('Age9Males', 'Age9Females'),]
summary(lm(Cor ~ Transfer, data = MothersSexTransfer))
xtest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age9Males',]$Cor
ytest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age9Females',]$Cor
t.test(xtest, ytest)

#Sig results only
load('TransferResultsPlaceholder.RData')
SigRes <- unique(FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer == 'Age24' & FOM1EpiscoresResults$P<0.05 &
                                        FOM1EpiscoresResults$Cor>0.1 ,]$Protein)
SigRes <- c('CD6', 'CD8A', 'CXCL9','Flt3L', 'IL2B', 'IL18R1', 'TNFB')
FOM1EpiscoresResults <- FOM1EpiscoresResults[FOM1EpiscoresResults$Protein %in% SigRes,]
MothersSexTransfer <- FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer %in% c('Age24Males', 'Age24Females'),]
summary(lm(Cor ~ Transfer, data = MothersSexTransfer))
xtest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age24Males',]$Cor
ytest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age24Females',]$Cor
t.test(xtest, ytest)

load('TransferResultsPlaceholder.RData')

SigRes <- unique(FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer == 'Age9' & FOM1EpiscoresResults$P<0.05 &
                                        FOM1EpiscoresResults$Cor>0.1 ,]$Protein)
SigRes <- c('CD6', 'IL18R1', 'TNFB')
FOM1EpiscoresResults <- FOM1EpiscoresResults[FOM1EpiscoresResults$Protein %in% SigRes,]
MothersSexTransfer <- FOM1EpiscoresResults[FOM1EpiscoresResults$Transfer %in% c('Age9Males', 'Age9Females'),]
summary(lm(Cor ~ Transfer, data = MothersSexTransfer))
xtest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age9Males',]$Cor
ytest <- MothersSexTransfer[MothersSexTransfer$Transfer=='Age9Females',]$Cor
t.test(xtest, ytest)

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Plots to show longitudinal correlations between episcores and proteins.  
load('TransferResultsPlaceholder.RData')


Var9_NineDNAm <- cbind(Var9_NineDNAm, NineDNAm)
Var24_TwentyFourDNAm <- cbind(Var24_TwentyFourDNAm, TwentyFourDNAm)

Results <- Age24EpiscoresResults
VarMatch <- Var24_TwentyFourDNAm
Proteins <- Age24Episcores 
TrainSuffix <- '_F24'
TransferSuffix <- '_F9'

LongEpiscoreProts <- function(Results, VarMatch, Proteins, TrainSuffix, TransferSuffix){
  
  blankdf <- data.frame()
  for (protein in Proteins){
    
    protein <- gsub(TrainSuffix, '', protein)
    Res <- Results[Results$Protein==protein,]   
    Res <- Res[!(Res$Transfer=='Mothers' | Res$Transfer=='CVTraining'),]
    

    ProtCor <- tryCatch({cor.test(VarMatch[, paste0(protein, TrainSuffix)], VarMatch[, paste0(protein, TransferSuffix)]) }, error = function(e){NULL})
    ProtCorAll <- ProtCor[["estimate"]]
   
    
    if(is.null(ProtCor)){next}
    
    Mal <- VarMatch$Sex == 'M'
    
    #ProtCor <- cor.test(VarMatch[, paste0(protein, TrainSuffix)][Mal], VarMatch[, paste0(protein, TransferSuffix)][Mal])
    #ProtCorMal <- ProtCor[["estimate"]]
    
    #ProtCor <- cor.test(VarMatch[, paste0(protein, TrainSuffix)][-Mal], VarMatch[, paste0(protein, TransferSuffix)][-Mal])
    #ProtCorFem <- ProtCor[["estimate"]]
    
    #Res$ProtCor <- c(ProtCorAll, ProtCorFem, ProtCorMal)
    Res$ProtCor <- c(ProtCorAll)
    
    
  blankdf <- rbind(blankdf, Res)
  }
  
  return(blankdf)
 
}





Age24plotRes <- LongEpiscoreProts(Age24EpiscoresResults, Var24_TwentyFourDNAm, Age24Episcores, '_F24', '_F9')
Age9plotRes <- LongEpiscoreProts(Age9EpiscoresResults, Var9_NineDNAm, Age9Episcores, '_F9', '_F24')

# Reduce datasets to only the CV threshold values
Age9values <- Age9TrainResults[Age9TrainResults$Cor>=0.1,]$Protein
Age9values <- Age9values[!is.na(Age9values)]
Age9plotRes <- Age9plotRes[Age9plotRes$Protein %in% Age9values,]


Age24values <- Age24TrainResults[Age24TrainResults$Cor>=0.1,]$Protein
Age24values <- Age24values[!is.na(Age24values)]
Age24plotRes <- Age24plotRes[Age24plotRes$Protein %in% Age24values,]

ggplot(Age24plotRes, aes(x = Cor, y = ProtCor, color = Protein, shape = Transfer)) + 
  geom_point(size=3) + xlab(paste0('Transferred Correlation')) + 
  ylab('Protein @9 - Protein @24 Corerlation') + 
  theme(legend.position="none", axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20)) + 
  geom_abline() + xlim(-0.2, 0.8) + ylim(-0.2, 0.8)

summary(lm(Age24plotRes$Cor ~ Age24plotRes$ProtCor))
cor.test(Age24plotRes$Cor,Age24plotRes$ProtCor)



ggplot(Age9plotRes, aes(x = Cor, y = ProtCor, color = Protein, shape = Transfer)) + 
  geom_point(size=3) + xlab(paste0('Transferred Correlation')) + 
  ylab('Protein @9 - Protein @24 Corerlation') + 
  theme(legend.position="none", axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        axis.text=element_text(size=15), axis.title=element_text(size=20),
        legend.text=element_text(size=15), legend.title=element_text(size=20)) + 
  geom_abline() + xlim(-0.2, 0.8) + ylim(-0.2, 0.8)


summary(lm(Age9plotRes$Cor ~ Age9plotRes$ProtCor))
cor.test(Age9plotRes$Cor,Age9plotRes$ProtCor)

#------------------------------------------------------------------------------#
# Effect of protein variance on model performance
Results <- Age24EpiscoresResults
VarMatch <- Var24_TwentyFourDNAm
Proteins <- Age24Episcores
TrainSuffix <- '_F24'

EpiscoreProtsVars <- function(Results, VarMatch, Proteins, TrainSuffix){
  
  blankdf <- data.frame()
  for (protein in Proteins){
    
    protein <- gsub(TrainSuffix, '', protein)
    Res <- Results[Results$Protein==protein,]   
    Res <- Res[Res$Transfer=='CVTraining',]
    
    ProtCorAll <- var(VarMatch[, paste0(protein, TrainSuffix)], na.rm = T)

    Mal <- VarMatch$Sex == 'M'
    
    #ProtCorMal <- var(VarMatch[, paste0(protein, TrainSuffix)][Mal], na.rm = T)

    #ProtCorFem <- var(VarMatch[, paste0(protein, TrainSuffix)][-Mal], na.rm = T)

    Res$ProtVar <- ProtCorAll
    
    

    
    blankdf <- rbind(blankdf, Res)
  }
  
  return(blankdf)
  
}


Age9VarRes <- EpiscoreProtsVars(Age9EpiscoresResults, Var9_NineDNAm, Age9Episcores, '_F9')
Age24VarRes <- EpiscoreProtsVars(Age24EpiscoresResults, Var24_TwentyFourDNAm, Age24Episcores, '_F24')
FOM1VarRes <- EpiscoreProtsVars(FOM1EpiscoresResults, VarFOM1_FOM1DNAm, FOM1Episcores, '_FOM1')


library(ggplot2)
ggplot(Age9VarRes, aes(x = Cor, y = ProtVar, color = Protein)) + 
  geom_point(size=3) + xlab(paste0('Model trained @ Age 9 CV training performance')) + 
  ylab('Protein @ 9 Variation') + 
  ggtitle(paste0('Relationship between episcore built at age 9 and Protein Variance')) 


summary(lm(Age9VarRes$Cor ~ Age9VarRes$ProtVar))
summary(lm(Age9VarRes[Age9VarRes$Cor>0.10,]$Cor ~ Age9VarRes[Age9VarRes$Cor>0.10,]$ProtVar))


ggplot(Age24VarRes, aes(x = Cor, y = ProtVar, color = Protein)) + 
  geom_point(size=3) + xlab(paste0('Model trained @ Age 24 CV training performance')) + 
  ylab('Protein @ 24 Variation') + 
  ggtitle(paste0('Relationship between episcore built at age 24 and Protein Variance')) 


summary(lm(Age24VarRes$Cor ~ Age24VarRes$ProtVar))
summary(lm(Age24VarRes[Age24VarRes$Cor>0.10,]$Cor ~ Age24VarRes[Age24VarRes$Cor>0.10,]$ProtVar))


ggplot(FOM1VarRes, aes(x = Cor, y = ProtVar, color = Protein)) + 
  geom_point(size=3) + xlab(paste0('Model trained @ FOM1 CV training performance')) + 
  ylab('Protein @ FOM1 Variation') + 
  ggtitle(paste0('Relationship between episcore built at FOM1 and Protein Variance')) 

summary(lm(FOM1VarRes$Cor ~ FOM1VarRes$ProtVar))
summary(lm(FOM1VarRes[FOM1VarRes$Cor>0.10,]$Cor ~ FOM1VarRes[FOM1VarRes$Cor>0.10,]$ProtVar))

#==============================================================================#
# Checking age 24 trained with 450k only compared to EPIC
# CV training results

library(readr)
Age24TrainResults450k <- read_csv("Age24TrainResults450kALL.csv")
Age24TrainResults <- read_csv("Age24TrainResults.csv")

colnames(Age24TrainResults450k) <- c('i', 'r_450k', 'p_450k', 'Protein')
colnames(Age24TrainResults) <- c('i1', 'r_EPIC', 'p_EPIC', 'Protein1')

Age24Compare <- cbind(Age24TrainResults, Age24TrainResults450k)

plot(Age24Compare$r_450k, Age24Compare$r_EPIC)
cor.test(Age24Compare$r_450k, Age24Compare$r_EPIC)



#==============================================================================#
# Results of models that are retrained for multi-CpG models 
library(readr)
Age24Rebuilt <- read_csv("Age24_RebuiltModels.csv")
Age24Rebuilt <- unique(Age24Rebuilt$x)

Age9Rebuilt <- read_csv("Age9_RebuiltModels.csv")
Age9Rebuilt <- unique(Age9Rebuilt$x)

FOM1Rebuilt <- read_csv("FOM1_RebuiltModels.csv")
FOM1Rebuilt <- unique(FOM1Rebuilt$x)
