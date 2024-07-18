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
NineDNAm <- aries.select(aries.dir, time.point="F9")
NineDNAm$meth <- aries.methylation(NineDNAm)

#TwentyFourDNAm <- aries.select(aries.dir, time.point="F24")
#TwentyFourDNAm$meth <- aries.methylation(TwentyFourDNAm)

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
load("Data/OlinkVars.RData")

# age match variables
Var9_NineDNAm <- Vars[match(NineDNAm$samples$alnqlet, Vars$aln_qlet),]
#Var24_TwentyFourDNAm <- Vars[match(TwentyFourDNAm$samples$alnqlet, Vars$aln_qlet),]

#------------------------------------------------------------------------------#
# Function for running protein panel ~ ARIES DNAm penalized regression 

PenRegFun <- function(ARIESDNAm, VarMatch, Age, trainsplit=0.8){

  # Combine data
  samples <- ARIESDNAm[['samples']]
  samples <- cbind(samples, VarMatch)
  
  #Covariates
  #samples$sex <- factor(samples$sex)
  sex <- samples$sex
  sex <- ifelse(sex == 'M', 0, 1)
  sex <- as.factor(sex)
  samples$batch <- factor(samples$plate)
  cellcounts <- ARIESDNAm$cell.counts[["blood-gse35069-complete"]]
  
  
  samples <- cbind(samples, cellcounts)
  #samples <- samples[c('sex', 'batch', 'Bcell', 'Mono', 'Neu', 'CD4T', 'CD8T', 'Eos', 'NK')]
  
  # Get DNAm and Impute NA values
  #methdata <- data.frame(t(ARIESDNAm[["meth"]]))
  #methdata <- na.aggregate(methdata)
  #print('Imputation complete')
  #save(methdata, file = "Age9Methdata.RData")
  load("Data/Age9Methdata.RData")
  
  # remove sex Chr
  load('Data/sexcg.RData')
  methdata <- methdata[,!names(methdata) %in% sexanno]
  
  # Get proteins of interest 
  ProteinList <- samples %>% select(matches(Age))
  ProteinList <- colnames(ProteinList)
    
  load('EpiscoreModelNames.RData')
  #ProteinList <- EpiscoreNames[['Age9Episcores']]
  
  # initialise empty results tables
  Coeftable <- data.frame()
  Testtable <- data.frame()
  Traintable <- data.frame()  
  CVTraintable <- data.frame()
  TrainPredsData <- list()

  
  
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
    residprot <- resid(lm(protval ~ sex, na.action = na.exclude))
    residprot <- as.vector(residprot)
    #save(residprot, protval,  file = 'Residtest.RData')
    
    
    # Make dataframe of methylation and outcome protein
    analysisdf <- cbind(residprot, methdata)
    analysisdf <- analysisdf[!is.na(analysisdf[,1]),]
    #analysisdf <- data.frame(analysisdf[complete.cases(analysisdf),])
    if(nrow(analysisdf) >10){
    
    # Split into test and train data
    training.samples <- createDataPartition(analysisdf[,1], p = trainsplit, list = FALSE)
    training.samples <- data.frame(training.samples)
    train.data  <- analysisdf[training.samples$Resample1, ]
    #test.data <- analysisdf[-training.samples$Resample1, ]
  
  
    # Run Model on training data
    x <- subset (train.data, select = -1)
    x <- as.matrix(x)
    y <- matrix(train.data[,1])
    
    elnet.cv <- cva.glmnet(x, y, family="gaussian", nfolds=5, keep=T)
    #save(elnet.cv, file='ELNETFix.RData')
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
    
    # Get training results
    alphalist <- elnet.cv$alpha
    BestAlpha <- get_model_params(elnet.cv)[[1]]
    BestAlpha <- match(BestAlpha, alphalist)
    
    lambdalist <- elnet.cv[["modlist"]][[BestAlpha]]$lambda
    BestLambda <- match(elnet.cv[["modlist"]][[BestAlpha]]$lambda.min, lambdalist)
    
    # Force a non zero lambda
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
    CVTraintable <- rbind(CVTraintable, CVTrainRes)
    
    TrainPredsData[[protein]] <- df
    
    
    #Save to csv in case of timeout failure
    write.csv(Coeftable,"Age9ResidProtNoSexChrCoeftabletest.csv")
    #write.csv(Testtable,"Age9ResidProtNoSexChrCpGtabletest.csv")
    write.csv(Traintable,"Age9TrainResults.csv")
    write.csv(CVTraintable,"Age9CVTrainResults.csv")
    

    
    # If not enough data, go to next iteration of loop 
    }else{next}
  }
  
  listdf <- list(Coeftable, Traintable, TrainPredsData)
  #listdf <- list(Coeftable, Testtable)
  return(listdf)
  
}

#------------------------------------------------------------------------------#

Age9PenReg <- PenRegFun(NineDNAm, Var9_NineDNAm, '_F9', trainsplit = 1.0)
save(Age9PenReg, file = 'PenRegAge9ResidProtNoSexChr_FullSamples.Rdata')
#Age24PenReg <- PenRegFun(TwentyFourDNAm, Var24_TwentyFourDNAm, '_F24', trainsplit=0.5)
#save(Age24PenReg, file = 'PenRegTest2.Rdata')

