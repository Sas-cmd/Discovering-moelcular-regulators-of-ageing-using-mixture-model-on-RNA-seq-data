#@Wrapper to be used to generate clustering at the gene level 
#@dat = Matrix that with gene genes as rows and columns as samples. Please ensure that genes are expressed in at least 90% of samples
#@modes = total number of modes to be fit to the data - the higher the number of modes the longer it will take
#@Partion_name = writes clustering result as a csv file using the partion_name selected 
#@ once clustering is done an exact test can be carried out to determine if clustering is due to sample grouping 


#Script to be run 
FinMM <- MMRun(dat = geneExp, modes = 1:5, Partion_name = Partion_name)

#libraries
library(Biobase)
library(SummarizedExperiment)
library(mclust)
library(RmixmodCombi)
library(MixAll)
library(data.table)
library(tidyverse)

#Functions
BestMM <- function(dat = dat, no_of_mode = modes ){
    
    #Mclust
    BICMclust <- NULL
    dat[is.na(dat)] <- 0
    try({mclustMixed <- Mclust(dat, G = no_of_mode)
    BICMclust <- as.numeric((-2*(mclustMixed$loglik) + mclustMixed$df *log(mclustMixed$n)))})
    
    
    #MixAll Gamma
    BICMixAll <- NULL
    dat[dat == 0] <- NA
    try({MixAllMixed <- clusterGamma(dat, 
                                     nbCluster = no_of_mode, 
                                     models = clusterGammaNames("all", "equal", "free", "free", "equal"),
                                     strategy =  clusterStrategy(), 
                                     criterion = "BIC", nbCore = 0)
    BICMixAll <- as.numeric(MixAllMixed@criterion)})
    
    ##RmixModCombi
    BICRmix <- NULL
    dat[is.na(dat)] <- 0
    try({Rmix <- mixmodCombi(data = dat , nbCluster = no_of_mode, mixmodOutput = NULL,
                             criterion = c("BIC", "ICL"), models = mixmodGaussianModel())
    BICRmix <- as.numeric(-2*(Rmix@mixmodOutput@bestResult@likelihood) + Rmix@mixmodOutput@bestResult@parameters@nbFreeParam *log(Rmix@mixmodOutput@nbSample))})
    
    BestBIC <- min(c(BICMclust, BICMixAll, BICRmix))
    
    if(BestBIC == BICMclust){
        BestOutput <- mclustMixed$classification
        BestModel <- mclustMixed$modelName
    } else if( BestBIC == BICRmix){
        BestOutput <- Rmix@mixmodOutput@bestResult@partition
        BestModel <- Rmix@mixmodOutput@bestResult@model
    } else {
        BestOutput <- MixAllMixed@zi
        BestModel <- MixAllMixed@component@modelName
    }
    
    Mainoutput <- list("Classification" = BestOutput, "BIC" = BestBIC, "Model" = BestModel)
    return(Mainoutput) 
}
MMRun <- function( dat = dat, modes = x, Partion_name = Partion_name){
    AllBIC <- data.table()
    Allpartion <- NULL
    for (j in 1:nrow(dat)) {
        mm <- BestMM(dat = dat[j,], no_of_mode = modes)
        t1 <- as.matrix(mm$Classification)
        colnames(t1) <- rownames(dat)[j]
        t2 <- data.table(t(c(rownames(dat)[j], 
                             mm$BIC, mm$Model)))
        Allpartion <- as.matrix(cbind(Allpartion, t1))
        AllBIC <- rbindlist(list(AllBIC, t2))
        print(j)
    }
    
    
    
    rownames(Allpartion) <- colnames(dat)
    Allpartion <- as.data.frame(Allpartion)
    
    
    write.csv(Allpartion, file = sprintf("BestModel3_%s.csv", Partion_name)) 
    write.csv(AllBIC, file = sprintf("BestModel_BIC3_%s.csv", Partion_name)) 
    
    finlist <- list(Allpartion, AllBIC)
    names(finlist) <- c(paste0("BestModel_", Partion_name), 
                        paste0("BestModel_BIC_", Partion_name))
    
    return(finlist)
}