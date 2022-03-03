## ---------------------------
##
## Script name: Discovering-moelcular-regulators-of-ageing-using-mixture-model-on-RNA-seq-data Figure 4
##
## Author: Sasdekumar Loganthan
##
## Date Created: 2021-10-24
##
## Email: s.loganathan@uq.edu.au
##
## ---------------------------
## Packages required
library(VGAM)
library(tidyverse)
library(data.table)
library(dplyr)
library(mclust)
library(RmixmodCombi)
library(MixAll)

set.seed(2353535)

# Functions ---------------------------

## load up our functions into memory
#Selecting the best MM
BestMM <- function(dat = dat, no_of_mode = modes ){
    #Mclust
    BICMclust <- NULL
    dat[is.na(dat)] <- 0
    try({mclustMixed <- Mclust(dat, G = no_of_mode)
    BICMclust <- as.numeric((-2*(mclustMixed$loglik) + mclustMixed$df *log(mclustMixed$n)))})
    
    
    #MixAll Gamma
    BICMixAll <- NULL
    dat[dat == 0] <- NA
    dat <- as.data.frame(dat)
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
    
    try(if(BestBIC == BICMclust){
        BestOutput <- mclustMixed$classification
        BestModel <- mclustMixed$modelName
    } else if( BestBIC == BICRmix){
        BestOutput <- Rmix@mixmodOutput@bestResult@partition
        BestModel <- Rmix@mixmodOutput@bestResult@model
    } else {
        BestOutput <- MixAllMixed@zi
        BestModel <- MixAllMixed@component@modelName
    })
    
    Mainoutput <- list("Classification" = BestOutput, "BIC" = BestBIC, "Model" = BestModel)
    return(Mainoutput) 
} #dat = dataframe, no_of_modes = total numbe of modes required
MMRun <- function( dat = dat, modes = x, Partion_name = Partion_name){
    AllBIC <- data.table()
    Allpartion <- NULL
    for (j in 1:nrow(dat)) {
        mm <- BestMM(dat = dat[j,], no_of_mode = modes)
        t1 <- as.matrix(mm$Classification)
        colnames(t1) <- rownames(dat)[j]
        t2 <- data.table(t(c(rownames(dat)[j], mm$BIC, mm$Model)))
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
} #dat = dataframe, no_of_modes = total numbe of modes required, Partion_name = Name to save output files

#Data Prep and plotting of Figure 4A-----------
t1 <- read.csv("~/Best_Model_selected_NON_Wrapper.csv")

PackageTypeTable <- as.data.frame(with(t1, table(Tissue, Package)))

PackageTypeTable$Tissue <- gsub("_", " ", PackageTypeTable$Tissue, fixed=TRUE)
PackageTypeTable$Tissue <- tools::toTitleCase(PackageTypeTable$Tissue)

dat.p <- reshape2::melt(PackageTypeTable)

ggplot(dat.p, aes(x = Tissue, y = value , fill= Package)) + geom_bar(stat='identity')   + theme_classic() + ylim(0, 30000) +
    scale_fill_manual(values= wes_palette("GrandBudapest2", n = 2)) + 
    ylab("Numer of Genes")+ theme(axis.text.x = element_text(angle = 90, hjust=1)) 
# Data Prep for Figure 4B ---------------------------  
set.seed(2353535)

obj = readRDS("C:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

fint <- data.frame()
for(i in 1:length(TissueType)){#length(TissueType){
    Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
    pos <- which(colnames(edat) %in% rownames(Tissue)) 
    edat.Tiss <- edat[,pos]
    edat.Tiss[edat.Tiss == 0] <- NA
    #Removing genes that have less and equal to 10% zero values
    edat.Tiss <- edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
    edat.sel <- as.data.frame(edat.Tiss)
    
    
    
    for(j in 1:nrow(edat.sel)){
        dt <- as.data.frame(t(edat.sel[j,]))
        colnames(dt) <- "y2"
        
        fit2 <- vglm(y2 ~ 1, gammaR, data = dt, trace = TRUE, crit = "coef")
        
        t1 <- data.frame(
            rate = Coef(fit2)[1],
            shape = Coef(fit2)[2],
            mean = mean(dt$y2, na.ra = TRUE),
            Gene = rownames(edat.sel[j,]),
            Tissue = TissueType[i])
        rownames(t1) = NULL  
        fint <- rbindlist(list(fint, t1))
    }
}

write.csv(fint, file = "Gamma_parameters for GTEx.csv")

#-----------------------------------------------------#
"Selecting quantile values for each the parameters, shape & rate"
#fint <- read_csv("~/Gamma_parameters for GTEx.csv")

#Simulating for shape parameter 
forshape <- NULL
for(i in 1:4){
    t1 <- fint %>% 
        filter(quartile == i) %>% sample_n(125)
    forshape <- rbind(forshape, t1)        
}

shapedat100 <- as.matrix(t(apply(forshape, 1 , function (t){
    x <- rgamma(100, shape = forshape$shape, rate = forshape$rate)})))
rownames(shapedat100) <- forshape$Gene

shapedat400 <- as.matrix(t(apply(forshape, 1 , function (t){
    x <- rgamma(400, shape = forshape$shape, rate = forshape$rate)})))
rownames(shapedat400) <- forshape$Gene


FinS100 <- MMRun(dat = shapedat100, modes = 1:3, Partion_name = "Shape100")


FinS400 <- MMRun(dat = shapedat400, modes = 1:3, Partion_name = "Shape400")

#Selecting quantile for rate
fint2 <- fint %>% 
    mutate(quartiler = ntile(rate, 4)) %>% 
    ungroup

forrate <- NULL
for(i in 1:4){
    t1 <- fint2 %>% 
        filter(quartiler == i) %>% sample_n(125)
    forrate <- rbind(forrate, t1)        
}

ratedat100 <- as.matrix(t(apply(forrate, 1 , function (t){
    x <- rgamma(100, shape = forrate$shape, rate = forrate$rate)})))
rownames(shapedat100) <- forrate$Gene

ratedat400 <- as.matrix(t(apply(forrate, 1 , function (t){
    x <- rgamma(400, shape = forrate$shape, rate = forrate$rate)})))

FinR100 <- MMRun(dat = ratedat100, modes = 1:3, Partion_name = "Rate100")

FinR400 <- MMRun(dat = ratedat400, modes = 1:3, Partion_name = "Rate400") 

#Plotting of Figure 4B---------------------
temp <- read_csv("D:/Ageing Mixture models project/Gamma distribution simulation/Gamma simulation running MM/BIC&Model/BestModel_BIC3_Rate100.csv")

Gaussian <- length(grep("Gaussian", temp$V2))
Gamma <- length(grep("gamma", temp$V2 ))  
Mclust <- nrow(temp) - Gaussian + Gamma   

Rate100 <- data.frame(Gaussian, Gamma, Mclust)
rownames(Rate100) <- "Rate100"

temp <- read_csv("D:/Ageing Mixture models project/Gamma distribution simulation/Gamma simulation running MM/BIC&Model/BestModel_BIC3_Rate400.csv")

Gaussian <- length(grep("Gaussian", temp$V2))
Gamma <- length(grep("gamma", temp$V2 ))  
Mclust <- nrow(temp) - Gaussian + Gamma   

Rate400 <- data.frame(Gaussian, Gamma, Mclust)
rownames(Rate400) <- "Rate400"

temp <- read_csv("D:/Ageing Mixture models project/Gamma distribution simulation/Gamma simulation running MM/BIC&Model/BestModel_BIC3_Shape100.csv")

Gaussian <- length(grep("Gaussian", temp$V3))
Gamma <- length(grep("gamma", temp$V3 ))  
Mclust <- nrow(temp) - Gaussian + Gamma   

Shape100 <- data.frame(Gaussian, Gamma, Mclust)
rownames(Shape100) <- "Shape100"

temp <- read_csv("D:/Ageing Mixture models project/Gamma distribution simulation/Gamma simulation running MM/BIC&Model/BestModel_BIC3_Shape400.csv")

Gaussian <- length(grep("Gaussian", temp$V3))
Gamma <- length(grep("gamma", temp$V3 ))  
Mclust <- nrow(temp) - Gaussian + Gamma   

Shape400 <- data.frame(Gaussian, Gamma, Mclust)
rownames(Shape400) <- "Shape400"

FinalPackageforSim <- as.matrix(rbind(Shape100, Shape400, Rate100, Rate400))

colnames(FinalPackageforSim) <- c("RmixModCombi", "MixAllGamma", "Mclust")

dat <- reshape2::melt(FinalPackageforSim)


ggplot(dat, aes(x = Var1, y = value, fill= Var2)) + 
    geom_bar(stat='identity')  + 
    scale_fill_manual(values=wes_palette(n = 3, 
    name = "IsleofDogs1")) +
    labs(
        x = "Parameters for the gamma distribution simulation",
        y = "Number of simulated genes", 
        fill = "Type of Package") +
    theme_classic() 

