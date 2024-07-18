library(aries)
library(haven)
library(limma)
library(ggplot2)

getwd()
#------------------------------------------------------------------------------#
# ARIES data location
aries.dir <- "/projects/MRC-IEU/research/projects/icep2/wp3/006/working/data/aries"

# Load in the child datapoints 
FirstDNAm <- aries.select(aries.dir, time.point="antenatal")
FirstDNAm$meth <- aries.methylation(FirstDNAm)
test <- FirstDNAm$meth

# Load in the ALSPAC variables
mothdata1 <- read_dta('Data/b_4f.dta')
mothdata2 <- read_dta('Data/c_8a.dta')
mothdata3 <- read_dta('Data/d_4b.dta')
mothdata4 <- read_dta('Data/OB_1b.dta')

mothdatajoin <- merge(x=mothdata1,y=mothdata2,by="aln",all=TRUE)
mothdatajoin <- merge(x=mothdatajoin,y=mothdata3,by="aln",all=TRUE)
mothdatajoin <- merge(x=mothdatajoin,y=mothdata4,by="aln",all=TRUE)

mumbloods1 <- read_dta('Data/Mother_samples.dta')
mumbloods2 <- read_dta('Data/mother_metabolomics.dta')
mumbloods <- merge(x=mumbloods1,y=mumbloods2,by="aln",all=TRUE)


# FOM1 is matched to FOM ARIES data
VarFirst <- merge(x=mothdatajoin,y=mumbloods,by="aln",all=TRUE)
FOM1 <- read_dta('Data/FOM1.dta')




#Load in Gadd et al data coefficients
Episcore_Coef <- read.csv("Data/Episcore_Coef.csv")

#------------------------------------------------------------------------------#
# Build Gadd et all episcores in ALSPAC data


# Array specific information
Soma <-  Episcore_Coef[Episcore_Coef$Panel == 'SomaScan',]
Olink <- Episcore_Coef[Episcore_Coef$Panel == 'Olink',]
#ArrayIntersect <- intersect(Soma$`Gene Name`, Olink$`Gene Name`) # 5 shared proteins

# Rename Olink proteins
Olink$'Gene.Name' <- paste(Olink$'Gene.Name', '_Olink', sep = '')

# Combine dataframes
Episcore_Coef <- rbind(Soma, Olink)
Episcore_Genes <- unique(Episcore_Coef$`Gene.Name`) #108, 5 of which are doubled

ARIESDNAm <- FirstDNAm
# Episcore function
EpiscoreFun <- function(ARIESDNAm) {
  
  BetaValues <- data.frame(ARIESDNAm[['meth']])
  
  # Reduce Episcore_Coefs to CpGs in beta values
  betanames <- rownames(BetaValues)
  Episcore_Coef <- Episcore_Coef[Episcore_Coef$CpG.Site %in% betanames,] 
  
  
  # Empty table to add episcores (Number of samples +1 for protein name)
  EpiscoreTable <- data.frame()
  
  
  # Loop to create episcores in siNET data
  
  
  for(Protein in Episcore_Genes){
    
    # Load in Protein specific episcore details
    df <- data.frame(Episcore_Coef[Episcore_Coef$`Gene.Name`==Protein,])  
    
    #Grab protein specific CpG sites
    cglist <- unique(df$CpG.Site)
    CGcoef <- df$CpG.Coeficient
    CGbeta <- BetaValues[cglist,]
    
    # calculate episcore
    Episcore <- CGcoef*CGbeta
    #if (length(Episcore)>977){
      Episcore <- colSums(Episcore, na.rm = T)
    #}
    
    # Add to Table of data
    newrow <- c(Protein, Episcore)
    EpiscoreTable <- rbind(EpiscoreTable, newrow)
  }
  # Column names
  Episcore_Cols <- c('Episcore')
  Episcore_Cols <- append(Episcore_Cols, colnames(BetaValues))
  colnames(EpiscoreTable) <- Episcore_Cols
  
  # Make episcores gene name column the row name
  EpiscoreTable <- data.frame(EpiscoreTable[,-1], row.names = EpiscoreTable[,1])
  i <- c(1:ncol(EpiscoreTable))    
  EpiscoreTable[ , i] <- apply(EpiscoreTable[ , i], 2, function(x) as.numeric(as.character(x)))
  
  return(EpiscoreTable)
}

# Get Episcores for each age group
EpiscoreFirst <- EpiscoreFun(FirstDNAm)
AntenatalDNAm <- FirstDNAm$samples
save(EpiscoreFirst, VarFirst, AntenatalDNAm, file = 'antenatalEpiscores.RData')