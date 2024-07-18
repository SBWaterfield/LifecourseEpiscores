library(aries)
library(haven)
library(limma)
library(glmnet)
library(glmnetUtils)
library(ggplot2)
library(tidyverse)
library(caret)
library(dplyr)
library(zoo)
library(fastDummies)
library(stringr)

homedir <- getwd()
#------------------------------------------------------------------------------#
# # ARIES data location
aries.dir <- aries.dir <- "/projects/MRC-IEU/research/projects/icep2/wp3/006/working/data/aries"

# Load in Data
#NineDNAm <- aries.select(aries.dir, time.point="F9")
#NineDNAm$meth <- aries.methylation(NineDNAm)

FOM1DNAm <- aries.select(aries.dir, time.point="FOM")
FOM1DNAm$meth <- aries.methylation(FOM1DNAm)

#------------------------------------------------------------------------------#
# Load in Olink data 
# ChildSampleNames <- read.csv("Data/ChildSampleNames.csv")

# # Child protein (olink) data (F9 & F24)
# Vars <- read_dta('Data/Child_bloods_6a.dta')
# Vars <- data.frame(Vars)
# 
# ChildVars <- ChildSampleNames[str_detect(ChildSampleNames$`Var.Label`, "NPX"), ] 
# ChildVars <- ChildVars$`Var.name`
# ChildVars <- c('aln_qlet', ChildVars)
# 
# Vars <- Vars %>% select(ChildVars)
# #Vars <- Vars[, ChildVars]
# 
# # Clean ALSPAC data
# Vars[Vars < 0] <- NA # Remove -10 (and other) NA placeholders
# Vars$aln_qlet <-gsub("_","",as.character(Vars$aln_qlet))

#save(Vars, file = 'OlinkVars.RData')
load("Data/OlinkMothersVars.RData")

# age match variables
#Var9_NineDNAm <- Vars[match(NineDNAm$samples$alnqlet, Vars$aln_qlet),]
VarFOM1_FOM1DNAm <- OlinkMothers[match(FOM1DNAm$samples$ALN, OlinkMothers$ALN),]
VarFOM1_FOM1DNAm[VarFOM1_FOM1DNAm < 0] <- NA # Remove -10 (and other) NA placeholders
MothersVars <- VarFOM1_FOM1DNAm
#save(MothersVars, file = 'MothersVars.RData')
#------------------------------------------------------------------------------#
# Function for running protein panel ~ ARIES DNAm penalized regression 

PenRegFun <- function(ARIESDNAm, VarMatch, Age, trainsplit=0.8){

  # Combine data
  samples <- ARIESDNAm[['samples']]
  samples <- cbind(samples, VarMatch)
  
  #Covariates
  #samples$sex <- factor(samples$sex)
  # sex <- samples$sex
  # sex <- ifelse(sex == 'M', 0, 1)
  # sex <- as.factor(sex)
  # samples$batch <- factor(samples$plate)
  # cellcounts <- ARIESDNAm$cell.counts[["blood-gse35069-complete"]]
  # 
  
  #samples <- cbind(samples, cellcounts)
  #samples <- samples[c('sex', 'batch', 'Bcell', 'Mono', 'Neu', 'CD4T', 'CD8T', 'Eos', 'NK')]
  
  # Get DNAm and Impute NA values
  #methdata <- data.frame(t(ARIESDNAm[["meth"]]))
  #methdata <- na.aggregate(methdata)
  #print('Imputation complete')
  #save(methdata, file = "Data/FOM1Methdata.RData")
  load("Data/FOM1Methdata.RData")
  
  # Get proteins of interest 
  ProteinList <- samples %>% select(matches(Age))
  ProteinList <- colnames(ProteinList)
   
  load('EpiscoreModelNames.RData')
  #ProteinList <- EpiscoreNames[['FOM1Episcores']] 
  
  # initialise empty results tables
  Coeftable <- data.frame()
  Testtable <- data.frame()
  Traintable <- data.frame()
  TrainPredsData <- list()
  CVTraintable <- data.frame()
  
  
  counter <- 0
  # Loop over each protein
  for (protein in ProteinList){
  
    #Print progress details
    counter <- counter+1
    #cat('Protein No: ',counter,' modelling', protein)
    
    # Make sure splits are the same for each protein
    set.seed(100)
    
    
    # residualize proteins controlling for sex
    protval <- samples[,protein]
    
    # Make dataframe of methylation and outcome protein
    analysisdf <- cbind(protval, methdata)
    analysisdf <- analysisdf[!is.na(analysisdf[,1]),]

    #remove columns where there is onyl NA values
    analysisdf <- analysisdf[ , colSums(is.na(analysisdf)) < nrow(analysisdf)]
    
    # Split into test and train data
    training.samples <- createDataPartition(analysisdf[,1], p = trainsplit, list = FALSE)
    training.samples <- data.frame(training.samples)
    train.data  <- analysisdf[training.samples$Resample1, ]
    # test.data <- analysisdf[-training.samples$Resample1, ]
  
  
    # Run Model on training data
    x <- subset (train.data, select = -1)
    x <- as.matrix(x)
    y <- matrix(train.data[,1])
    
    elnet.cv <- cva.glmnet(x, y, family="gaussian", nfolds=5, keep=T)
    
    get_model_params <- function(fit) {
      alpha <- fit$alpha
      lambdaMin <- sapply(fit$modlist, `[[`, "lambda.min")
      lambdaSE <- sapply(fit$modlist, `[[`, "lambda.1se")
      error <- sapply(fit$modlist, function(mod) {min(mod$cvm)})
      best <- which.min(error)
      data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
                 lambdaSE = lambdaSE[best], eror = error[best])
    }
     
    modelparams <- get_model_params(elnet.cv)
    

    fit <- glmnet(x, y, family = "gaussian", alpha = modelparams$alpha, lambda = modelparams$lambdaMin, standardize = F) 
    fitcoefficients <- coef(fit)
    fitcoefficients = data.frame(Variable = fitcoefficients@Dimnames[[1]][fitcoefficients@i+1], Coefficient = fitcoefficients@x,
                                 Protein = protein)
    # Get training results
    alphalist <- elnet.cv$alpha
    BestAlpha <- get_model_params(elnet.cv)[[1]]
    BestAlpha <- match(BestAlpha, alphalist)
    
    lambdalist <- elnet.cv[["modlist"]][[BestAlpha]]$lambda
    BestLambda <- match(elnet.cv[["modlist"]][[BestAlpha]]$lambda.min, lambdalist)
    TrainPreds <- elnet.cv[["modlist"]][[BestAlpha]][["fit.preval"]]
    TrainPreds <- TrainPreds[, BestLambda]
    
    df <- TrainPreds
    
    
    Trainingfit <- cor.test(TrainPreds, y)
    TrainRes <- c(Trainingfit[["estimate"]], Trainingfit[["p.value"]], protein)
    
    #CV train results 
    ID <- elnet.cv[["modlist"]][[BestAlpha]][["foldid"]]
    
    fold_correlations <- numeric(length(unique(ID)))
    
    # Loop over each fold
    for (fold_id in unique(ID)) {
      # Subset predictions and actual Y values for the current fold
      fold_preds <- TrainPreds[ID == fold_id]
      fold_actual <- y[ID == fold_id]
      
      # Calculate correlation coefficient for the current fold
      fold_correlations[fold_id] <- cor(fold_actual, fold_preds)
    }
    
    # Calculate the mean correlation coefficient across all folds
    cv_correlation <- mean(fold_correlations)
    cv_SD <- sd(fold_correlations)
    
    CVTrainRes <- c(cv_correlation, cv_SD,  protein)
    
    
    
    fitcoefficients <- coef(fit)
    fitcoefficients = data.frame(Variable = fitcoefficients@Dimnames[[1]][fitcoefficients@i+1], Coefficient = fitcoefficients@x,
                                 Protein = protein)
    
    counter = 1
    # Retrain if model is reduced to intercept
    while (nrow(fitcoefficients)<3){
      
      counter = counter + 1
      lambdalist <- elnet.cv[["modlist"]][[BestAlpha]]$lambda
      lambdalist <- lambdalist[counter]
      fit <- glmnet(x, y, family = "gaussian", alpha = modelparams$alpha, lambda =  lambdalist, standardize = F) 
      
      
      # Force a slightly smaller lambda
      TrainPreds <- elnet.cv[["modlist"]][[BestAlpha]][["fit.preval"]]
      TrainPreds <- TrainPreds[, counter] #Take second largest lambda
      
      df <- TrainPreds
      
      Trainingfit <- cor.test(TrainPreds, y)
      TrainRes <- c(Trainingfit[["estimate"]], Trainingfit[["p.value"]], protein)
      
      #CV train results 
      ID <- elnet.cv[["modlist"]][[BestAlpha]][["foldid"]]
      
      fold_correlations <- numeric(length(unique(ID)))
      
      # Loop over each fold
      for (fold_id in unique(ID)) {
        # Subset predictions and actual Y values for the current fold
        fold_preds <- TrainPreds[ID == fold_id]
        fold_actual <- y[ID == fold_id]
        
        # Calculate correlation coefficient for the current fold
        fold_correlations[fold_id] <- cor(fold_actual, fold_preds)
      }
      
      # Calculate the mean correlation coefficient across all folds
      cv_correlation <- mean(fold_correlations)
      cv_SD <- sd(fold_correlations)
      
      CVTrainRes <- c(cv_correlation, cv_SD,  protein)
      
      fitcoefficients <- coef(fit)
      fitcoefficients = data.frame(Variable = fitcoefficients@Dimnames[[1]][fitcoefficients@i+1], Coefficient = fitcoefficients@x,
                                   Protein = protein)
      
    }
    # Make predictions on test data
    # x.test <- subset (test.data, select = -1)
    # x.test <- as.matrix(x.test)
    # predictions <- predict(fit, x.test)
    
    # Model performance metrics
   # modelmetrics <- data.frame(
   #  RMSE = RMSE(predictions, test.data[, protein]),
   #  Rsquare = R2(predictions, test.data[, protein])
   #  )


    #run model correlation, if it fails return an empty dataframe
    # modelcor <- tryCatch({cor.test(predictions, test.data[,1])}, error=function(e)data.frame())
    # 
    # # Get model data if the model ran 
    # if (is.null(ncol(modelcor))) {
    #   modelcor <- cbind(modelcor$estimate, modelcor$p.value)
    #   modelcor <- data.frame(modelcor)
    #   modelcor$Protein <- protein
    # } else {
    # modelcor <- data.frame()
    # }

    # Add data to results tables
    Coeftable <- rbind(Coeftable, fitcoefficients)
    #Testtable <- rbind(Testtable, modelcor)
    Traintable <- rbind(Traintable, TrainRes)
    TrainPredsData[[protein]] <- df
    CVTraintable <- rbind(CVTraintable, CVTrainRes)
    
  
    
    #Save to csv in case of timeout failure
    write.csv(Coeftable,"FOM1ResidProtCoeftable.csv")
    #write.csv(Testtable,"FOM1ResidProtCOrtable.csv")
    write.csv(Traintable,"FOM1TrainResults.csv")
    write.csv(CVTraintable,"FOM1CVTrainResults.csv")
    
  }
  
  #listdf <- list(Coeftable, Testtable)
  listdf <- list(Coeftable, TrainPredsData)
  return(listdf)
  
}

#------------------------------------------------------------------------------#


FOM1PenReg <- PenRegFun(FOM1DNAm, VarFOM1_FOM1DNAm, '_FOM1', trainsplit=1)
save(FOM1PenReg, file = 'PenRegFOM1ResidProt.Rdata')

