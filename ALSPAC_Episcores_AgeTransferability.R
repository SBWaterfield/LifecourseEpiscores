library(tidyverse)
library(ggplot2)
library(aries)
 
# Load in DNAm and Olink data
#load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/ARIESDNAm.RData")
#load("C:/Users/gv20022/OneDrive - University of Bristol/Desktop/MP - PheWAS/OlinkVars.RData")

aries.dir <- aries.dir <- "/projects/MRC-IEU/research/projects/icep2/wp3/006/working/data/aries"

 
load("Data/OlinkVars.RData")

# Load in Data
NineDNAm <- aries.select(aries.dir, time.point="F9")
NineDNAm$meth <- aries.methylation(NineDNAm)

TwentyFourDNAm <- aries.select(aries.dir, time.point="F24")
TwentyFourDNAm$meth <- aries.methylation(TwentyFourDNAm)


#save(NineDNAm, TwentyFourDNAm, file = 'ARIESDNAm.RData')


# Match variables
Var9_NineDNAm <- Vars[match(NineDNAm$samples$alnqlet, Vars$aln_qlet),]
Var24_TwentyFourDNAm <- Vars[match(TwentyFourDNAm$samples$alnqlet, Vars$aln_qlet),]

rm(Vars)

# Episcore function
EpiscoreFun <- function(ARIESDNAm, Coefdata) {
  
  BetaValues <- data.frame(ARIESDNAm[['meth']])
  
  # Reduce Episcore_Coefs to CpGs in beta values
  betanames <- rownames(BetaValues)
  Episcore_Coef <- Coefdata[Coefdata$Variable %in% betanames,] 
  
  
  # Empty table to add episcores (Number of samples +1 for protein name)
  EpiscoreTable <- data.frame()
  ProteinNames <- c()
  
  # Loop to create episcores in siNET data
  
  
  for(Protein in GOI){
    # Load in Protein specific episcore details
    df <- data.frame(Episcore_Coef[Episcore_Coef$`Protein`==Protein,])  
    
    #Grab protein specific CpG sites
    cglist <- unique(df$Variable)
    CGcoef <- df$Coefficient
    CGbeta <- BetaValues[cglist,]
    
    # calculate episcore
    Episcore <- CGcoef*CGbeta
    Episcore <- colSums(Episcore, na.rm = T)

    # Add to Table of data
    newrow <- Episcore
    #newrow <- c(Protein, Episcore)
    EpiscoreTable <- rbind(EpiscoreTable, newrow)
    
    # Get protein name for rows
    ProteinNames <- c(ProteinNames, Protein)
  }
  
  # Column names
  colnames(EpiscoreTable) <- colnames(BetaValues)
  # Episcore_Cols <- c('Episcore')
  # Episcore_Cols <- append(colnames(BetaValues), Episcore_Cols)
  # colnames(EpiscoreTable) <- Episcore_Cols
  
  # Make episcores gene name column the row name
  EpiscoreTable <- data.frame(EpiscoreTable, row.names = ProteinNames)
  # EpiscoreTable <- data.frame(EpiscoreTable[,-1], row.names = EpiscoreTable[,1])
  # i <- c(1:ncol(EpiscoreTable))    
  # EpiscoreTable[ , i] <- apply(EpiscoreTable[ , i], 2, function(x) as.numeric(as.character(x)))
  #  
  return(EpiscoreTable)
}


# calculate episcores of interest
library(readr)
#------------------------------------------------------------------------------#
#Sex specific models into age groups
# 
# #Age 24 male models into other datasets
# Age24MaleIntoFemaleCoef <- read_csv("Age24MaleIntoFemaleCoef.csv")
# Age24MaleIntoFemaleCoef <- Age24MaleIntoFemaleCoef[Age24MaleIntoFemaleCoef$Variable!= '(Intercept)',]
# Age24MaleIntoFemaleCoef <- Age24MaleIntoFemaleCoef[Age24MaleIntoFemaleCoef$Coefficient!= 0,]
# Age24MaleIntoFemaleCoef <- Age24MaleIntoFemaleCoef %>% select(-`...1`)
# Age24MaleIntoFemaleCoef <- Age24MaleIntoFemaleCoef[!duplicated(Age24MaleIntoFemaleCoef),]
# 
# GOI <- unique(Age24MaleIntoFemaleCoef$Protein)
# Age24MaleEpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, Age24MaleIntoFemaleCoef)
# save(Age24MaleEpiscoresinAge24, file='Age24MaleEpiscoresinAge24.RData')
# 
# Age24MaleEpiscoresinAge9 <- EpiscoreFun(NineDNAm, Age24MaleIntoFemaleCoef)
# save(Age24MaleEpiscoresinAge9, file='Age24MaleEpiscoresinAge9.RData')
# 
# #Age 24 female models
# Age24FemaleIntoMaleCoef <- read_csv("Age24FemaleIntoMaleCoef.csv")
# Age24FemaleIntoMaleCoef <- Age24FemaleIntoMaleCoef[Age24FemaleIntoMaleCoef$Variable!= '(Intercept)',]
# Age24FemaleIntoMaleCoef <- Age24FemaleIntoMaleCoef[Age24FemaleIntoMaleCoef$Coefficient!= 0,]
# Age24FemaleIntoMaleCoef <- Age24FemaleIntoMaleCoef %>% select(-`...1`)
# 
# Age24FemaleIntoMaleCoef <- Age24FemaleIntoMaleCoef[!duplicated(Age24FemaleIntoMaleCoef),]
# 
# GOI <- unique(Age24FemaleIntoMaleCoef$Protein)
# Age24FemaleEpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, Age24FemaleIntoMaleCoef)
# save(Age24FemaleEpiscoresinAge24, file='Age24FemaleEpiscoresinAge24.RData')
# 
# Age24FemaleEpiscoresinAge9 <- EpiscoreFun(NineDNAm, Age24FemaleIntoMaleCoef)
# save(Age24FemaleEpiscoresinAge9, file='Age24FemaleEpiscoresinAge9.RData')
# 
# # Age 9 Male models
# Age9MaleIntoFemaleCoef <- read_csv("Age9MaleIntoFemaleCoef.csv")
# Age9MaleIntoFemaleCoef <- Age9MaleIntoFemaleCoef[Age9MaleIntoFemaleCoef$Variable!= '(Intercept)',]
# Age9MaleIntoFemaleCoef <- Age9MaleIntoFemaleCoef[Age9MaleIntoFemaleCoef$Coefficient!= 0,]
# Age9MaleIntoFemaleCoef <- Age9MaleIntoFemaleCoef %>% select(-`...1`)
# 
# Age9MaleIntoFemaleCoef <- Age9MaleIntoFemaleCoef[!duplicated(Age9MaleIntoFemaleCoef),]
# 
# GOI <- unique(Age9MaleIntoFemaleCoef$Protein)
# Age9MaleEpiscoresinAge9 <- EpiscoreFun(NineDNAm, Age9MaleIntoFemaleCoef)
# save(Age9MaleEpiscoresinAge9, file='Age9MaleEpiscoresinAge9.RData')
# 
# Age9MaleEpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, Age9MaleIntoFemaleCoef)
# save(Age9MaleEpiscoresinAge24, file='Age9MaleEpiscoresinAge24.RData')
# 
# # Age 9 female models
# Age9FemaleIntoMaleCoef <- read_csv("Age9FemaleIntoMaleCoef.csv")
# Age9FemaleIntoMaleCoef <- Age9FemaleIntoMaleCoef[Age9FemaleIntoMaleCoef$Variable!= '(Intercept)',]
# Age9FemaleIntoMaleCoef <- Age9FemaleIntoMaleCoef[Age9FemaleIntoMaleCoef$Coefficient!= 0,]
# Age9FemaleIntoMaleCoef <- Age9FemaleIntoMaleCoef %>% select(-`...1`)
# 
# Age9FemaleIntoMaleCoef <- Age9FemaleIntoMaleCoef[!duplicated(Age9FemaleIntoMaleCoef),]
# 
# 
# GOI <- unique(Age9FemaleIntoMaleCoef$Protein)
# Age9FemaleEpiscoresinAge9 <- EpiscoreFun(NineDNAm, Age9FemaleIntoMaleCoef)
# save(Age9FemaleEpiscoresinAge9, file='Age9FemaleEpiscoresinAge9.RData')
# 
# Age9FemaleEpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, Age9FemaleIntoMaleCoef)
# save(Age9FemaleEpiscoresinAge24, file='Age9FemaleEpiscoresinAge24.RData')
#------------------------------------------------------------------------------#
# Age specific Models into other age groups

#Age 9 Models
Age9Coef <- read_csv("Age9ResidProtNoSexChrCoeftabletest.csv")
Age9Coef <- Age9Coef[Age9Coef$Variable!= '(Intercept)',]
Age9Coef <- Age9Coef[Age9Coef$Coefficient!= 0,]
Age9Coef <- Age9Coef %>% select(-`...1`)

Age9Coef <- Age9Coef[!duplicated(Age9Coef),]

GOI <- unique(Age9Coef$Protein)
Age9EpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, Age9Coef)
save(Age9EpiscoresinAge24, file='Age9EpiscoresinAge24.RData')

Age9EpiscoresinAge9 <- EpiscoreFun(NineDNAm, Age9Coef)
save(Age9EpiscoresinAge9, file='Age9EpiscoresinAge9.RData')

# Age 24 models
Age24Coef <- read_csv("Age24ResidProtNoSexChrCoeftabletest.csv")
Age24Coef <- Age24Coef[Age24Coef$Variable!= '(Intercept)',]
Age24Coef <- Age24Coef[Age24Coef$Coefficient!= 0,]
Age24Coef <- Age24Coef %>% select(-`...1`)

Age24Coef <- Age24Coef[!duplicated(Age24Coef),]

Age24Coef450k <- read_csv("Age24ResidProtNoSexChrCoeftabletest450k.csv")
Age24Coef450k <- Age24Coef450k[Age24Coef450k$Variable!= '(Intercept)',]
Age24Coef450k <- Age24Coef450k[Age24Coef450k$Coefficient!= 0,]
Age24Coef450k <- Age24Coef450k %>% select(-`...1`)

Age24Coef450k <- Age24Coef450k[!duplicated(Age24Coef450k),]

Age24Coef <- rbind(Age24Coef, Age24Coef450k)

GOI <- unique(Age24Coef$Protein)
Age24EpiscoresinAge9 <- EpiscoreFun(NineDNAm, Age24Coef)
save(Age24EpiscoresinAge9, file='Age24EpiscoresinAge9.RData')

Age24EpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, Age24Coef)
save(Age24EpiscoresinAge24, file='Age24EpiscoresinAge24.RData')
#------------------------------------------------------------------------------#
# mothers models into age and sex specific groups
FOM1Coef <- read_csv("FOM1ResidProtCoeftable.csv")
FOM1Coef <- FOM1Coef[FOM1Coef$Variable!= '(Intercept)',]
FOM1Coef <- FOM1Coef[FOM1Coef$Coefficient!= 0,]
FOM1Coef <- FOM1Coef %>% select(-`...1`)

FOM1Coef <- FOM1Coef[!duplicated(FOM1Coef),]

GOI <- unique(FOM1Coef$Protein)

FOM1EpiscoresinAge24 <- EpiscoreFun(TwentyFourDNAm, FOM1Coef)
save(FOM1EpiscoresinAge24, file='FOM1EpiscoresinAge24.RData')

FOM1EpiscoresinAge9 <- EpiscoreFun(NineDNAm, FOM1Coef)
save(FOM1EpiscoresinAge9, file='FOM1EpiscoresinAge9.RData')


#------------------------------------------------------------------------------#
rm(NineDNAm, TwentyFourDNAm)

#------------------------------------------------------------------------------#
# All Models into Mothers
load("Data/OlinkMothersVars.RData")
Vars <- OlinkMothers

FOM1DNAm <- aries.select(aries.dir, time.point="FOM")
FOM1DNAm$meth <- aries.methylation(FOM1DNAm)

VarFOM1_FOM1DNAm <- Vars[match(FOM1DNAm$samples$aln, Vars$aln),]

GOI <- unique(FOM1Coef$Protein)
FOM1EpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, FOM1Coef)
save(FOM1EpiscoresinFOM1, file='FOM1EpiscoresinFOM1.RData')


GOI <- unique(Age9Coef$Protein)
Age9EpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, Age9Coef)
save(Age9EpiscoresinAge24, file='Age9EpiscoresinFOM1.RData')

GOI <- unique(Age24Coef$Protein)
Age24EpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, Age24Coef)
save(Age24EpiscoresinFOM1, file='Age24EpiscoresinFOM1.RData')


# GOI <- unique(Age24MaleIntoFemaleCoef$Protein)
# Age24MaleEpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, Age24MaleIntoFemaleCoef)
# save(Age24MaleEpiscoresinFOM1, file='Age24MaleEpiscoresinFOM1.RData')
# 
# 
# GOI <- unique(Age24FemaleIntoMaleCoef$Protein)
# Age24FemaleEpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, Age24FemaleIntoMaleCoef)
# save(Age24FemaleEpiscoresinFOM1, file='Age24FemaleEpiscoresinFOM1.RData')
# 
# 
# GOI <- unique(Age9MaleIntoFemaleCoef$Protein)
# Age9MaleEpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, Age9MaleIntoFemaleCoef)
# save(Age9MaleEpiscoresinFOM1, file='Age9MaleEpiscoresinFOM1.RData')
# 
# 
# GOI <- unique(Age9FemaleIntoMaleCoef$Protein)
# Age9FemaleEpiscoresinFOM1 <- EpiscoreFun(FOM1DNAm, Age9FemaleIntoMaleCoef)
# save(Age9FemaleEpiscoresinFOM1, file='Age9FemaleEpiscoresinFOM1.RData')



rm(FOM1DNAm, Var9_NineDNAm, Var24_TwentyFourDNAm, VarFOM1_FOM1DNAm, Age9Coef, Age24Coef,
   Age24MaleIntoFemaleCoef, Age24FemaleIntoMaleCoef, Age9MaleIntoFemaleCoef, Age9FemaleIntoMaleCoef,
   FOM1Coef)

save.image('ALSPACEpsicoresTransferred.RData')
