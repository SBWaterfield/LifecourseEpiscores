library(labelled)

# Load in ALSPAC Olink data
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/OlinkVars.RData")

# Load in ARIES Gadd episcores
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Episcores.Rdata")

# Load in ARIES sample sheets 
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIES_Samples.RData")

# Load in Mothers Data
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/EpiscoresFOM1.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/MothersVars.RData")
load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/antenatalEpiscores.RData")

#Set NA samples
MothersVars[MothersVars < 0] <- NA


# Match samples to episcores
Var9_NineDNAm <- Vars[match(NineDNAm$alnqlet, Vars$aln_qlet),]
Var9_SevenDNAm <- Vars[match(SevenDNAm$alnqlet, Vars$aln_qlet),]
Var9_FifteenDNAm <- Vars[match(FifteenDNAm$alnqlet, Vars$aln_qlet),] #Also used for 24_15 variables
Var24_TwentyFourDNAm <- Vars[match(TwentyFourDNAm$alnqlet, Vars$aln_qlet),]
VarFOM1_AntenatalDNAm <- MothersVars[match(AntenatalDNAm$aln, MothersVars$ALN),]


# transpose episcores
EpiscoreSeven <- data.frame(t(EpiscoreSeven))
EpiscoreNine <- data.frame(t(EpiscoreNine))
EpiscoreFifteen <- data.frame(t(EpiscoreFifteen))
EpiscoreTwentyFour <- data.frame(t(EpiscoreTwentyFour))
EpiscoresFOM1 <- data.frame(t(EpiscoresFOM1))
EpiscoresAntenatal <- data.frame(t(EpiscoreFirst))



# Get matches between episcores and olink data
Age9Olink <- (Var9_NineDNAm[,c(2:90)])
Age24Olink <- (Var24_TwentyFourDNAm[,c(91:179)])
AgeFOM1Olink <- MothersVars[c(1:91)]

Age9OlinkAge7 <- (Var9_SevenDNAm[,c(2:90)])
Age9OlinkAge15 <- (Var9_FifteenDNAm[,c(2:90)])
Age24OlinkAge15 <- (Var9_FifteenDNAm[,c(91:179)])
AgeFOM1OlinkAgeAntenal <- (VarFOM1_AntenatalDNAm[,c(1:91)])


colnames(Age9Olink) <-  gsub("_.*","", colnames(Age9Olink))
colnames(Age24Olink) <-  gsub("_.*","", colnames(Age24Olink))
colnames(AgeFOM1Olink) <-  gsub("_.*","", colnames(AgeFOM1Olink))


colnames(Age9OlinkAge7) <-  gsub("_.*","", colnames(Age9OlinkAge7))
colnames(Age9OlinkAge15) <-  gsub("_.*","", colnames(Age9OlinkAge15))
colnames(Age24OlinkAge15) <-  gsub("_.*","", colnames(Age24OlinkAge15))
colnames(AgeFOM1OlinkAgeAntenal) <-  gsub("_.*","", colnames(AgeFOM1OlinkAgeAntenal))




colnames(Age9Olink) <- paste0(colnames(Age9Olink), "_Olink")
colnames(Age24Olink) <- paste0(colnames(Age24Olink), "_Olink")
colnames(AgeFOM1Olink) <- paste0(colnames(AgeFOM1Olink), "_Olink")


colnames(Age9OlinkAge7) <- paste0(colnames(Age9OlinkAge7), "_Olink")
colnames(Age9OlinkAge15) <- paste0(colnames(Age9OlinkAge15), "_Olink")
colnames(Age24OlinkAge15) <- paste0(colnames(Age24OlinkAge15), "_Olink")
colnames(AgeFOM1OlinkAgeAntenal) <- paste0(colnames(AgeFOM1OlinkAgeAntenal), "_Olink")




matchingnames <- intersect(colnames(Age9Olink), colnames(EpiscoreNine))
matchingnames <- c(matchingnames, 'CXCL11', 'CXCL10', 'MMP1', 'CCL25')



Age9Olink$CXCL11 <- Age9Olink$CXCL11_Olink
Age9Olink$CXCL10 <- Age9Olink$CXCL10_Olink
Age9Olink$MMP1 <- Age9Olink$MMP1_Olink
Age9Olink$CCL25 <- Age9Olink$CCL25_Olink

Age9OlinkAge7$CXCL11 <- Age9OlinkAge7$CXCL11_Olink
Age9OlinkAge7$CXCL10 <- Age9OlinkAge7$CXCL10_Olink
Age9OlinkAge7$MMP1 <- Age9OlinkAge7$MMP1_Olink
Age9OlinkAge7$CCL25 <- Age9OlinkAge7$CCL25_Olink

Age9OlinkAge15$CXCL11 <- Age9OlinkAge15$CXCL11_Olink
Age9OlinkAge15$CXCL10 <- Age9OlinkAge15$CXCL10_Olink
Age9OlinkAge15$MMP1 <- Age9OlinkAge15$MMP1_Olink
Age9OlinkAge15$CCL25 <- Age9OlinkAge15$CCL25_Olink

Age24Olink$CXCL11 <- Age24Olink$CXCL11_Olink
Age24Olink$CXCL10 <- Age24Olink$CXCL10_Olink
Age24Olink$MMP1 <- Age24Olink$MMP1_Olink
Age24Olink$CCL25 <- Age24Olink$CCL25_Olink

Age24OlinkAge15$CXCL11 <- Age24OlinkAge15$CXCL11_Olink
Age24OlinkAge15$CXCL10 <- Age24OlinkAge15$CXCL10_Olink
Age24OlinkAge15$MMP1 <- Age24OlinkAge15$MMP1_Olink
Age24OlinkAge15$CCL25 <- Age24OlinkAge15$CCL25_Olink

AgeFOM1Olink$CXCL11 <- AgeFOM1Olink$CXCL11_Olink
AgeFOM1Olink$CXCL10 <- AgeFOM1Olink$CXCL10_Olink
AgeFOM1Olink$MMP1 <- AgeFOM1Olink$MMP1_Olink
AgeFOM1Olink$CCL25 <- AgeFOM1Olink$CCL25_Olink

AgeFOM1OlinkAgeAntenal$CXCL11 <- AgeFOM1OlinkAgeAntenal$CXCL11_Olink
AgeFOM1OlinkAgeAntenal$CXCL10 <- AgeFOM1OlinkAgeAntenal$CXCL10_Olink
AgeFOM1OlinkAgeAntenal$MMP1 <- AgeFOM1OlinkAgeAntenal$MMP1_Olink
AgeFOM1OlinkAgeAntenal$CCL25 <- AgeFOM1OlinkAgeAntenal$CCL25_Olink

AgeFOM1Olink <- remove_val_labels(AgeFOM1Olink)
AgeFOM1OlinkAgeAntenal <- remove_val_labels(AgeFOM1OlinkAgeAntenal)

# Run correlations between episcore and actual data values
Age9Correlations <- data.frame()
Age7Age9EpiCorrelations <- data.frame()
Age15Age9EpiCorrelations <- data.frame()
Age24Correlations <- data.frame()
Age15Age24EpiCorrelations <- data.frame()
AgeFOM1Correlations <- data.frame()
AgeFOM1AgeAntenalCorrelations <- data.frame()

for (episcore in matchingnames){

  Age9res <- cor.test(EpiscoreNine[,episcore], Age9Olink[,episcore], use = 'complete.obs')
  N <- Age9res[["parameter"]][["df"]]+2
  Age9res <- c(episcore, Age9res$estimate, Age9res$p.value, Age9res[["conf.int"]],N)

  Age24res <- cor.test(EpiscoreTwentyFour[,episcore], Age24Olink[,episcore], use = 'complete.obs')
  N <- Age24res[["parameter"]][["df"]]+2
  Age24res <- c(episcore, Age24res$estimate, Age24res$p.value, Age24res[["conf.int"]],N)

  
  
  AgeFOM1res <- cor.test(EpiscoresFOM1[,episcore], unlist(AgeFOM1Olink[,episcore]), use = 'complete.obs')
  N <- AgeFOM1res[["parameter"]][["df"]]+2
  AgeFOM1res <- c(episcore, AgeFOM1res$estimate, AgeFOM1res$p.value, AgeFOM1res[["conf.int"]],N)
  
  Age7Age9res <- cor.test(EpiscoreSeven[,episcore], Age9OlinkAge7[,episcore], use = 'complete.obs')
  N <- Age7Age9res[["parameter"]][["df"]]+2
  Age7Age9res <- c(episcore, Age7Age9res$estimate, Age7Age9res$p.value, Age7Age9res[["conf.int"]],N)
  
  Age15Age9res <- cor.test(EpiscoreFifteen[,episcore], Age9OlinkAge15[,episcore], use = 'complete.obs')
  N <- Age15Age9res[["parameter"]][["df"]]+2
  Age15Age9res <- c(episcore, Age15Age9res$estimate, Age15Age9res$p.value, Age15Age9res[["conf.int"]],N)
  
  Age15Age24res <- cor.test(EpiscoreFifteen[,episcore], Age24OlinkAge15[,episcore], use = 'complete.obs')
  N <- Age15Age24res[["parameter"]][["df"]]+2
  Age15Age24res <- c(episcore, Age15Age24res$estimate, Age15Age24res$p.value, Age15Age24res[["conf.int"]],N)

  
  Age9Correlations <- rbind(Age9Correlations, Age9res)
  Age24Correlations <- rbind(Age24Correlations, Age24res)
  AgeFOM1Correlations <- rbind(AgeFOM1Correlations, AgeFOM1res)
  
  Age7Age9EpiCorrelations <- rbind(Age7Age9EpiCorrelations, Age7Age9res)
  Age15Age9EpiCorrelations <- rbind(Age15Age9EpiCorrelations, Age15Age9res)
  Age15Age24EpiCorrelations <- rbind(Age15Age24EpiCorrelations, Age15Age24res)
  
  AgeFOM1AgeAntenalres <- cor.test(EpiscoresAntenatal[,episcore], as.numeric(unlist(AgeFOM1OlinkAgeAntenal[,episcore])), use = 'complete.obs')
  N <- AgeFOM1AgeAntenalres[["parameter"]][["df"]]+2
  AgeFOM1AgeAntenalres <- c(episcore, AgeFOM1AgeAntenalres$estimate, AgeFOM1AgeAntenalres$p.value, AgeFOM1AgeAntenalres[["conf.int"]],N)
  AgeFOM1AgeAntenalCorrelations <- rbind(AgeFOM1AgeAntenalCorrelations, AgeFOM1AgeAntenalres)
  

  
}

colnames(Age9Correlations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')
colnames(Age24Correlations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')
colnames(AgeFOM1Correlations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')


colnames(Age7Age9EpiCorrelations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')
colnames(Age15Age9EpiCorrelations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')
colnames(Age15Age24EpiCorrelations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')
colnames(AgeFOM1AgeAntenalCorrelations) <- c('Episcore', 'Cor', 'P', 'LI', 'UI','N')


Age9Correlations$Cor <- as.numeric(Age9Correlations$Cor)
Age24Correlations$Cor <- as.numeric(Age24Correlations$Cor)
AgeFOM1Correlations$Cor <- as.numeric(AgeFOM1Correlations$Cor)

Age7Age9EpiCorrelations$Cor <- as.numeric(Age7Age9EpiCorrelations$Cor)
Age15Age9EpiCorrelations$Cor <- as.numeric(Age15Age9EpiCorrelations$Cor)
Age15Age24EpiCorrelations$Cor <- as.numeric(Age15Age24EpiCorrelations$Cor)
AgeFOM1AgeAntenalCorrelations$Cor <- as.numeric(AgeFOM1AgeAntenalCorrelations$Cor)


Age9Correlations$P <- as.numeric(Age9Correlations$P)
Age24Correlations$P <- as.numeric(Age24Correlations$P)
AgeFOM1Correlations$P <- as.numeric(AgeFOM1Correlations$P)

Age7Age9EpiCorrelations$P <- as.numeric(Age7Age9EpiCorrelations$P)
Age15Age9EpiCorrelations$P <- as.numeric(Age15Age9EpiCorrelations$P)
Age15Age24EpiCorrelations$P <- as.numeric(Age15Age24EpiCorrelations$P)
AgeFOM1AgeAntenalCorrelations$P <- as.numeric(AgeFOM1AgeAntenalCorrelations$P)


Age9Correlations$LI <- as.numeric(Age9Correlations$LI)
Age24Correlations$LI <- as.numeric(Age24Correlations$LI)
AgeFOM1Correlations$LI <- as.numeric(AgeFOM1Correlations$LI)

Age7Age9EpiCorrelations$LI <- as.numeric(Age7Age9EpiCorrelations$LI)
Age15Age9EpiCorrelations$LI <- as.numeric(Age15Age9EpiCorrelations$LI)
Age15Age24EpiCorrelations$LI <- as.numeric(Age15Age24EpiCorrelations$LI)
AgeFOM1AgeAntenalCorrelations$LI <- as.numeric(AgeFOM1AgeAntenalCorrelations$LI)


Age9Correlations$UI <- as.numeric(Age9Correlations$UI)
Age24Correlations$UI <- as.numeric(Age24Correlations$UI)
AgeFOM1Correlations$UI <- as.numeric(AgeFOM1Correlations$UI)

Age7Age9EpiCorrelations$UI <- as.numeric(Age7Age9EpiCorrelations$UI)
Age15Age9EpiCorrelations$UI <- as.numeric(Age15Age9EpiCorrelations$UI)
Age15Age24EpiCorrelations$UI <- as.numeric(Age15Age24EpiCorrelations$UI)
AgeFOM1AgeAntenalCorrelations$UI <- as.numeric(AgeFOM1AgeAntenalCorrelations$UI)



Age9Correlations <- Age9Correlations[order(Age9Correlations$Cor, decreasing = T ),]
Age24Correlations <- Age24Correlations[order(Age24Correlations$Cor, decreasing = T),]
AgeFOM1Correlations <- AgeFOM1Correlations[order(AgeFOM1Correlations$Cor, decreasing = T),]

Age7Age9EpiCorrelations <- Age7Age9EpiCorrelations[order(Age7Age9EpiCorrelations$Cor, decreasing = T),]
Age15Age9EpiCorrelations <- Age15Age9EpiCorrelations[order(Age15Age9EpiCorrelations$Cor, decreasing = T),]
Age15Age24EpiCorrelations <- Age15Age24EpiCorrelations[order(Age15Age24EpiCorrelations$Cor, decreasing = T),]
AgeFOM1AgeAntenalCorrelations <- AgeFOM1AgeAntenalCorrelations[order(AgeFOM1AgeAntenalCorrelations$Cor, decreasing = T),]



Age9Correlations[c('Cor', 'P', 'LI', 'UI')] <- signif(Age9Correlations[c('Cor', 'P', 'LI', 'UI')],2)
Age24Correlations[c('Cor', 'P', 'LI', 'UI')] <- signif(Age24Correlations[c('Cor', 'P', 'LI', 'UI')],2)
AgeFOM1Correlations[c('Cor', 'P', 'LI', 'UI')] <- signif(AgeFOM1Correlations[c('Cor', 'P', 'LI', 'UI')],2)
Age7Age9EpiCorrelations[c('Cor', 'P', 'LI', 'UI')] <- signif(Age7Age9EpiCorrelations[c('Cor', 'P', 'LI', 'UI')],2)
Age15Age9EpiCorrelations[c('Cor', 'P', 'LI', 'UI')] <- signif(Age15Age9EpiCorrelations[c('Cor', 'P', 'LI', 'UI')],2)
Age15Age24EpiCorrelations[c('Cor', 'P', 'LI', 'UI')] <- signif(Age15Age24EpiCorrelations[c('Cor', 'P', 'LI', 'UI')],2)
AgeFOM1AgeAntenalCorrelations[c('Cor', 'P', 'LI', 'UI')] <- signif(AgeFOM1AgeAntenalCorrelations[c('Cor', 'P', 'LI', 'UI')],2)



write.csv(Age9Correlations, file = 'Age9GaddInALSPACCor.csv')
write.csv(Age24Correlations, file = 'Age24GaddInALSPACCor.csv')
write.csv(AgeFOM1Correlations, file = 'AgeFOM1GaddInALSPACCor.csv')
write.csv(Age7Age9EpiCorrelations, file = 'Age7Age9EpiCorrelations.csv')
write.csv(Age15Age9EpiCorrelations, file = 'Age15Age9EpiCorrelations.csv')
write.csv(Age15Age24EpiCorrelations, file = 'Age15Age24EpiCorrelations.csv')
write.csv(AgeFOM1AgeAntenalCorrelations, file = 'AgeFOM1AgeAntenalCorrelations.csv')



# Check how many variables in columns 

GOI9 <- Age9Olink[matchingnames]
test <-sapply(GOI9, function(y) sum(length(which(!is.na(y)))))
max(test) #222
min(test) #222
median(test) #222

GOI24 <- Age24Olink[matchingnames]
test <-sapply(GOI24, function(y) sum(length(which(!is.na(y)))))
max(test) #763
min(test) #763
median(test) #763


GOIFOM1 <- AgeFOM1Olink[matchingnames]
test <-sapply(GOIFOM1, function(y) sum(length(which(!is.na(y)))))
max(test) #622
min(test) #622
median(test) #622

#------------------------------------------------------------------------------#
# Plot data
FOM1Gadd <- read.csv('AgeFOM1GaddInALSPACCor.csv')
Age9Gadd <- read.csv('Age9GaddInALSPACCor.csv')
Age24Gadd <- read.csv('Age24GaddInALSPACCor.csv')
FOM1Gadd$ModelName <- 'Mothers'
Age9Gadd$ModelName <- 'Age9'
Age24Gadd$ModelName <- 'Age24'
Moddetails <- rbind(FOM1Gadd, Age9Gadd, Age24Gadd)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels = c('Age9', 'Age24', 'Mothers'))
Moddetails$Episcore <- factor(Moddetails$Episcore, levels = unique(Moddetails$Episcore))
Moddetails$SameAge <- T
ds5 <- Moddetails[Moddetails$ModelName=='Mothers',]
ModdetailsAll <- Moddetails

library(ggplot2)
ggplot(Moddetails, aes(x=Episcore, y=Cor, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge(width = 0.8), width = 0.8) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') + labs(fill='Age') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  scale_fill_manual(values = c(Age9 = "#F8766D", Age24 = "#00BA38", Mothers = "#619CFF", CVTraining = 'grey', 
                               'Mothers (CV)' = 'grey', 'Age9 (CV)' = 'grey','Age24 (CV)' = 'grey'))

#Plot data of different Episcores in different ages compared to Protein measurement
Age7GaddAge9Olink <-  read.csv('Age7Age9EpiCorrelations.csv')
Age15GaddAge9Olink <-  read.csv('Age15Age9EpiCorrelations.csv')
Age15GaddAge24Olink <-  read.csv('Age15Age24EpiCorrelations.csv')
AgeAntenatalGaddAgeFOM1Olink <- read.csv('AgeFOM1AgeAntenalCorrelations.csv')
Age7GaddAge9Olink$ModelName <- 'Age7Gadd Age9Olink'
Age15GaddAge9Olink$ModelName <- 'Age15Gadd Age9Olink'
Age15GaddAge24Olink$ModelName <- 'Age15Gadd Age24Olink'
AgeAntenatalGaddAgeFOM1Olink$ModelName <- 'AntenatalGadd FOM1Olink'

Moddetails <- rbind(Age7GaddAge9Olink, Age15GaddAge9Olink, Age15GaddAge24Olink)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels = c('Age7Gadd Age9Olink', 'Age15Gadd Age9Olink', 'Age15Gadd Age24Olink'))
Moddetails$SameAge <- F
AgeAntenatalGaddAgeFOM1Olink$SameAge <- F
ModdetailsAll <- rbind(ModdetailsAll, Moddetails, AgeAntenatalGaddAgeFOM1Olink)

library(ggplot2)
ggplot(Moddetails, aes(x=Episcore, y=Cor, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge(width = 0.8), width = 0.8) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 


ds1 <- Moddetails[Moddetails$ModelName=='Age7Gadd Age9Olink',]
ds2 <- Moddetails[Moddetails$ModelName=='Age15Gadd Age9Olink',]
ds3 <- Moddetails[Moddetails$ModelName=='Age15Gadd Age24Olink',]
ds4 <- AgeAntenatalGaddAgeFOM1Olink


ds1$Cor_Age7Gadd_Age9Olink <- ds1$Cor
ds2$Cor_Age15Gadd_Age9Olink <- ds2$Cor
ds3$Cor_Age15Gadd_Age24Olink <- ds3$Cor
ds4$Cor_AntenatalGadd_FOM1Olink <- ds4$Cor
ds5$Cor_FOM1Gadd_FOM1Olink <- ds5$Cor


ds1$LI_Age7Gadd_Age9Olink <- ds1$LI
ds2$LI_Age15Gadd_Age9Olink <- ds2$LI
ds3$LI_Age15Gadd_Age24Olink <- ds3$LI
ds4$LI_AntenatalGadd_FOM1Olink <- ds4$LI
ds5$LI_FOM1Gadd_FOM1Olink <- ds5$LI
ds5$LI_FOM1Gadd_FOM1Olink <- ds5$LI


ds1$UI_Age7Gadd_Age9Olink <- ds1$UI
ds2$UI_Age15Gadd_Age9Olink <- ds2$UI
ds3$UI_Age15Gadd_Age24Olink <- ds3$UI
ds4$UI_AntenatalGadd_FOM1Olink <- ds4$UI
ds5$UI_FOM1Gadd_FOM1Olink <- ds5$UI


Moddetailstrans <- cbind(ds1, ds2, ds3)
Moddetailstrans <- Moddetailstrans[c('Cor_Age7Gadd_Age9Olink', 'Cor_Age15Gadd_Age9Olink', 'Cor_Age15Gadd_Age24Olink',
                                     'LI_Age7Gadd_Age9Olink', 'LI_Age15Gadd_Age9Olink', 'LI_Age15Gadd_Age24Olink',
                                     'UI_Age7Gadd_Age9Olink', 'UI_Age15Gadd_Age9Olink', 'UI_Age15Gadd_Age24Olink',
                                     'Episcore'
                                     )]

axisvals <- c(min(Moddetails$LI), max(Moddetails$UI)) 

# Add linear mmdoel for correlation
library(DEGreport)

ggplot(Moddetailstrans, aes(x=Cor_Age7Gadd_Age9Olink,y=Cor_Age15Gadd_Age9Olink, col=Episcore, fill=Episcore)) +
  geom_pointrange(aes(xmin=LI_Age7Gadd_Age9Olink, xmax=UI_Age7Gadd_Age9Olink)) +
  geom_pointrange(aes(ymin=LI_Age15Gadd_Age9Olink, ymax=UI_Age15Gadd_Age9Olink)) +
  geom_smooth(aes(x=Cor_Age7Gadd_Age9Olink,y=Cor_Age15Gadd_Age9Olink), inherit.aes = F, method = lm) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  geom_vline(xintercept = 0.1,linetype = 2)  +
  geom_vline(xintercept = 0,linetype = 1) +
  xlim(axisvals) + ylim(axisvals) + 
  xlab('Correlation Age7Gadd_Age9Olink (r)') + ylab('Correlation Age15Gadd_Age9Olink (r)') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 
cor.test(Moddetailstrans$Cor_Age7Gadd_Age9Olink, Moddetailstrans$Cor_Age15Gadd_Age9Olink)

ggplot(Moddetailstrans, aes(x=Cor_Age7Gadd_Age9Olink,y=Cor_Age15Gadd_Age24Olink, col=Episcore, fill=Episcore)) +
  geom_pointrange(aes(xmin=LI_Age7Gadd_Age9Olink, xmax=UI_Age7Gadd_Age9Olink)) +
  geom_pointrange(aes(ymin=LI_Age15Gadd_Age24Olink, ymax=UI_Age15Gadd_Age24Olink)) +
  geom_smooth(aes(x=Cor_Age7Gadd_Age9Olink,y=Cor_Age15Gadd_Age24Olink), inherit.aes = F, method = lm) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  geom_vline(xintercept = 0.1,linetype = 2)  +
  geom_vline(xintercept = 0,linetype = 1) +
  xlim(axisvals) + ylim(axisvals) + 
  xlab('Correlation Age7Gadd_Age9Olink (r)') + ylab('Correlation Age15Gadd_Age24Olink (r)') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 
cor.test(Moddetailstrans$Cor_Age7Gadd_Age9Olink, Moddetailstrans$Cor_Age15Gadd_Age24Olink)


ggplot(Moddetailstrans, aes(x=Cor_Age15Gadd_Age9Olink,y=Cor_Age15Gadd_Age24Olink, col=Episcore, fill=Episcore)) +
  geom_pointrange(aes(xmin=LI_Age15Gadd_Age9Olink, xmax=UI_Age15Gadd_Age9Olink)) +
  geom_pointrange(aes(ymin=LI_Age15Gadd_Age24Olink, ymax=UI_Age15Gadd_Age24Olink)) +
  geom_smooth(aes(x=Cor_Age15Gadd_Age9Olink,y=Cor_Age15Gadd_Age24Olink), inherit.aes = F, method = lm) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  geom_vline(xintercept = 0.1,linetype = 2)  +
  geom_vline(xintercept = 0,linetype = 1) +
  xlim(axisvals) + ylim(axisvals) + 
  xlab('Correlation Age15Gadd_Age9Olink (r)') + ylab('Correlation Age15Gadd_Age24Olink (r)') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 
cor.test(Moddetailstrans$Cor_Age15Gadd_Age9Olink, Moddetailstrans$Cor_Age15Gadd_Age24Olink)

# Mothers data
Moddetailstrans <- cbind(ds4, ds5)
Moddetailstrans <- Moddetailstrans[c('Cor_FOM1Gadd_FOM1Olink', 'LI_FOM1Gadd_FOM1Olink', 'UI_FOM1Gadd_FOM1Olink',
                                     'Cor_AntenatalGadd_FOM1Olink', 'LI_AntenatalGadd_FOM1Olink', 'UI_AntenatalGadd_FOM1Olink',
                                     'Episcore'
)]

axisvals <- c(min(min(Moddetailstrans$LI_FOM1Gadd_FOM1Olink), min(Moddetailstrans$LI_AntenatalGadd_FOM1Olink)),
              max(max(Moddetailstrans$UI_FOM1Gadd_FOM1Olink), max(Moddetailstrans$UI_AntenatalGadd_FOM1Olink)))

ggplot(Moddetailstrans, aes(x=Cor_FOM1Gadd_FOM1Olink,y=Cor_AntenatalGadd_FOM1Olink, col=Episcore, fill=Episcore)) +
  geom_pointrange(aes(xmin=LI_FOM1Gadd_FOM1Olink, xmax=UI_FOM1Gadd_FOM1Olink)) +
  geom_pointrange(aes(ymin=LI_AntenatalGadd_FOM1Olink, ymax=UI_AntenatalGadd_FOM1Olink)) +
  geom_smooth(aes(x=Cor_FOM1Gadd_FOM1Olink,y=Cor_AntenatalGadd_FOM1Olink), inherit.aes = F, method = lm) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  geom_vline(xintercept = 0.1,linetype = 2)  +
  geom_vline(xintercept = 0,linetype = 1) +
  xlim(axisvals) + ylim(axisvals) + 
  xlab('Correlation MiddleAgeGadd_MiddleAgeOlink (r)') + ylab('Correlation AntenatalGadd_MiddleAgeOlink (r)') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 
cor.test(Moddetailstrans$Cor_FOM1Gadd_FOM1Olink, Moddetailstrans$Cor_AntenatalGadd_FOM1Olink)


#------------------------------------------------------------------------------#
#All Kids models in one place
Moddetails <- rbind(Age7GaddAge9Olink, Age15GaddAge9Olink, Age15GaddAge24Olink, Age9Gadd, Age24Gadd)
Moddetails$ModelName <- factor(Moddetails$ModelName, levels = c('Age7Gadd Age9Olink', 'Age9', 'Age15Gadd Age9Olink', 'Age15Gadd Age24Olink', 'Age24'))

library(ggplot2)
ggplot(Moddetails, aes(x=Episcore, y=Cor, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge(width = 0.8), width = 0.8) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ggtitle(paste('Gadd Model performance in ALSPAC')) +
  ylab('Correlation (r)') + xlab('Model') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 



#==============================================================================#
# Plotting out how longitudinal projection alters correlation

ModdetailsAll <- ModdetailsAll[!ModdetailsAll$ModelName %in% c('Age15Gadd Age9Olink'),]

# Subset the data for Antenatal, Age9Olink, and Age24Olink model conditions
antenatal_data <- subset(ModdetailsAll, grepl("Antenatal", ModelName))
age9olink_data <- subset(ModdetailsAll, grepl("Age9Olink", ModelName))
age24olink_data <- subset(ModdetailsAll, grepl("Age24Olink", ModelName))
mothers_data <- subset(ModdetailsAll, ModelName == "Mothers")
age9_data <- subset(ModdetailsAll, ModelName == "Age9")
age24_data <- subset(ModdetailsAll, ModelName == "Age24")

# Merge datasets based on the 'Episcore' column
merged_antenatal <- merge(antenatal_data, mothers_data, by = "Episcore", suffixes = c("_Antenatal", "_Mothers"))
merged_age9olink <- merge(age9olink_data, age9_data, by = "Episcore", suffixes = c("_Age9Olink", "_Age9"))
merged_age24olink <- merge(age24olink_data, age24_data, by = "Episcore", suffixes = c("_Age24Olink", "_Age24"))

# Calculate the differences in 'Cor' for each respective 'Episcore' value
cor_diff_antenatal <- merged_antenatal$Cor_Antenatal - merged_antenatal$Cor_Mothers
cor_diff_age9olink <- merged_age9olink$Cor_Age9Olink - merged_age9olink$Cor_Age9
cor_diff_age24olink <- merged_age24olink$Cor_Age24Olink - merged_age24olink$Cor_Age24

# Create a new data frame with the calculated differences
cor_diff_data <- data.frame(
  Episcore = merged_antenatal$Episcore,
  Cor_Diff_Antenatal = cor_diff_antenatal,
  Cor_Diff_Age9Olink = cor_diff_age9olink,
  Cor_Diff_Age24Olink = cor_diff_age24olink
)

colnames(cor_diff_data) <- c('Episcore', 'Antenatal', 'Age7', 'Age15')
# Reshape the dataframe
cor_diff_data_long <- cor_diff_data %>%
  pivot_longer(cols = -Episcore, names_to = "Timepoint", values_to = "Difference")


ggplot(cor_diff_data_long, aes(x=Episcore, y=Difference, fill=Timepoint)) +
  geom_bar(stat='identity', position=position_dodge(width = 0.8), width = 0.8) +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation difference (r)') + xlab('Protein') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 


#==============================================================================#
# Plot the age transfers for age 7 --> 9, age 15 --> 24, Antenatal --> Middle age
AgeTransfers <- rbind(ds1[,1:7], ds3[, 1:7], ds4[,1:7])
AgeTransfers$ModelName <- as.character(AgeTransfers$ModelName)
# Alter model transfer names
AgeTransfers$ModelName <- ifelse(AgeTransfers$ModelName == 'Age7Gadd Age9Olink', 'Age7', AgeTransfers$ModelName)
AgeTransfers$ModelName <- ifelse(AgeTransfers$ModelName == 'Age15Gadd Age24Olink', 'Age15', AgeTransfers$ModelName)
AgeTransfers$ModelName <- ifelse(AgeTransfers$ModelName == 'AntenatalGadd FOM1Olink', 'Antenatal', AgeTransfers$ModelName)
AgeTransfers$ModelName <- factor(AgeTransfers$ModelName, levels = c('Age7', 'Age15', 'Antenatal'))
AgeTransfers <- AgeTransfers[order(-AgeTransfers$Cor), ]
AgeTransfers$Episcore <- factor(AgeTransfers$Episcore, levels = unique(AgeTransfers$Episcore))

library(ggplot2)
ggplot(AgeTransfers, aes(x=Episcore, y=Cor, fill=ModelName)) +
  geom_bar(stat='identity', position=position_dodge(width = 0.8), width = 0.8) +
  geom_hline(yintercept = 0.1,linetype = 2)  +
  geom_hline(yintercept = 0,linetype = 1) +
  ylab('Correlation (r)') + xlab('Protein') + labs(fill='Episcore Age') +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  scale_fill_manual(values = c(Age7 = "#F8766D", Age15 = "#00BA38", Antenatal = "#619CFF", CVTraining = 'grey', 
                               'Mothers (CV)' = 'grey', 'Age9 (CV)' = 'grey','Age24 (CV)' = 'grey'))

#------------------------------------------------------------------------------#
GOI <- c('IFNgamma', 'IL6', 'TNFSF14', 'CD8A', 'OSM')
Vars_place <- remove_val_labels(Vars)
ProtCorResults <- data.frame()
for (episcore in GOI){
  
  episcore9 <- paste0(episcore, '_F9')
  episcore24 <- paste0(episcore, '_F24')
  
  ProtCorrelation <- cor.test(as.numeric(Vars_place[,episcore9]), as.numeric(Vars_place[,episcore24]), use = 'complete.obs')
  ProtCorrelation <- c(episcore, ProtCorrelation$estimate, ProtCorrelation$p.value, ProtCorrelation[["conf.int"]])
  
  
  ProtCorResults <- rbind(ProtCorResults, ProtCorrelation) 
}
colnames(ProtCorResults) <- c('Epicscore', 'Cor', 'P', 'LI', 'UI')






























