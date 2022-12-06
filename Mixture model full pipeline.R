#Complete start to finish of mixture models pipeline
##### Main utilitis that needs to be run before running#####
library(Biobase)
library(SummarizedExperiment)
library(mclust)
library(RmixmodCombi)
library(MixAll)
library(data.table)
library(tidyverse)


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
fread_plus <- function(dat.read) {
    fread(dat.read, fill = TRUE, drop = 1) %>% 
        mutate(filename = dat.read)
}
rbindlist_fread <- function(path, pattern = "*.csv") {
    files = list.files(path, pattern, full.names = TRUE)
    rbindlist(lapply(files, function(x) fread_plus(x)))
}
func <- function(i){
    df <- read.csv(paste0(path,l[i]))
    names <- unlist(strsplit(l[i], "_"))
    df["BodyPart"] <- names[1]
    df["Name"] <- names[2]
    return(df)
}
AgeGroupFactors <- function(dat){
    dat$AGE[dat$AGE == 25] <- 1
    dat$AGE[dat$AGE == 35] <- 1
    dat$AGE[dat$AGE == 45] <- 2
    dat$AGE[dat$AGE == 55] <- 2
    dat$AGE[dat$AGE == 65] <- 3
    dat$AGE[dat$AGE == 75] <- 3
    
    return(dat)
}

ForExactTest <- function(data.new){
    
    data.new <- AgeGroupFactors(data.new)
    ft <- as.matrix(with(data.new, table(AGE, Gene)))  
    pt <- fisher.test(ft, workspace = 2e8)
    if(ncol(pt) == 1) next
    pv <- pt$p.value
    output1 <- data.table(t(c(as.numeric(pt$p.value), 
                              x[i]$Gene, x[i]$Package, 
                              as.numeric(x[i]$Best_BIC), 
                              as.numeric(length(unique(t1$Gene))))))
    colnames(output1) <- c("p.value", "Gene", "Package", "Best_BIC", "No_of_Modes")
    return(output1)
}


####Running of mixture models####

#Load data
obj = readRDS("gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

New.list <- list()
for(i in c(22, 23, 29, 34, 35)){
    Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
    pos <- which(colnames(edat) %in% rownames(Tissue)) 
    edat.Tiss <- edat[,pos]
    edat.Tiss[edat.Tiss == 0] <- NA
    #Removing genes that have less and equal to 10% zero values
    edat.Tiss <- edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
    edat.sel <- as.matrix(edat.Tiss)
    
    FinMM <- MMRun(dat = edat.sel, modes = 1:3, Partion_name = TissueType[i])
    New.list[[i]] <- FinMM #if dataset is small, or else save clustering and BIC as csv 
}

t2 = rbindlist_fread("D:/Ageing Mixture models project/Best Model_wrapper/BIC")

colnames(t2) = c("Genes", "BIC", "Model", "filename")



cleanName = function(tissues){
    tt = sub(".*/", "", tissues)
    tt2 = sub("\\..*", "", tt)
    str_remove_all(tt2, "BestModel3_BIC_")
}

t2$tissues = cleanName(t2$filename)

t2$Package = NA
for(i in 1:length(t2$Model)){
    tx = t2$Model[i]   
    if (grepl("Gaussian", tx)) {
        t2$Package[i] = "RmixModCombi"
    } else if(grepl("gamma", t2$Model)) {
        t2$Package[i] = "MixAllGamma"
    }else {t2$Package[i] = "Mclust"}
}

write.csv(t2, "Best_Model_Fin.csv")

#####Exact test #########
"Two exact tests needs to be carried out, one for age and the other for gender. Once the two tests are done we will next use remove genes that are common between the two"

#Ideally the data should have been saved rather than relying on 


temp <- list.files("D:/Ageing Mixture models project/1. New mixture model with wrapper/Best Model_wrapper/Best model", pattern="*.csv", full.names=TRUE)

Classlist <- lapply(temp, fread, header = TRUE)

names(Classlist) = TissueType

Classlist = lapply(Classlist,  function(x) {x <- x %>% select(-V1)})

Finaloutput <- data.table()
for(z in 1:length(TissueType)){
    d1 = Classlist[[z]]
    d1Gene = as.matrix(d1[,-c(1:2)])
    d1age = d1[,1:2]
    col_names <- colnames(d1Gene)
    for(j in 1:ncol(d1Gene)){
        d1Fin = cbind(d1age, d1Gene[,j])
        colnames(d1Fin) <- c("Patients", "AGE", "Gene")
        
        if(length(unique(d1Fin$Gene)) == 1) next
        
        d1Fin <- AgeGroupFactors(d1Fin)
        ft <- as.matrix(with(d1Fin, table(AGE, Gene))) 
        pt <- fisher.test(ft, workspace = 2e8)
        pv <- pt$p.value
        output1 <- data.table(t(c(as.numeric(pt$p.value), 
                                  col_names[j], 
                                  as.numeric(
                                      length(unique(
                                          d1Fin$Gene))))))
        colnames(output1) <- c("p.value", "Gene", 
                               "No_of_Modes")
        
        Finaloutput <- rbindlist(list(Finaloutput, output1))
        print(i)
    
    }
    write.csv(Finaloutput, file = sprintf("BestModel3_Exact_test_%s.csv", name[z]))
}

#Exact test for gender 
Finaloutput <- data.table()
for(z in 1:length(TissueType)){
    d1 = Classlist[[z]]
    d1Gene = as.matrix(d1[,-c(1:2)])
    col_names <- colnames(d1Gene)
    d1Gender = d1[,1]
    colnames(d1Gender) = "SAMPID"
    t6 = pdat %>% filter(our_subtypes == TissueType[z])
    t6 = subset(t6, select = c(SAMPID, GENDER))
    d1Gender = merge(d1Gender, t6, by = "SAMPID")
    for(j in 1:ncol(d1Gene)){
        d1Fin = cbind(d1Gender, d1Gene[,j])
        colnames(d1Fin) <- c("Patients", "GENDER", "Gene")
        
        if(length(unique(d1Fin$Gene)) == 1) next
        
        ft <- as.matrix(with(d1Fin, table(GENDER, Gene))) 
        pt <- fisher.test(ft, workspace = 2e8)
        pv <- pt$p.value
        output1 <- data.table(t(c(as.numeric(pt$p.value), 
                                  col_names[j], 
                                  as.numeric(
                                      length(unique(
                                          d1Fin$Gene))))))
        colnames(output1) <- c("p.value", "Gene", 
                               "No_of_Modes")
        
        Finaloutput <- rbindlist(list(Finaloutput, output1))
        print(i)
        
    }
    write.csv(Finaloutput, file = sprintf("BestModel3_Exact_test_GENDER_%s.csv", name[z]))
}

#####Padjust test #########
#Padjust for age
temp <- list.files("E:/Best Model_NON_wrapper/Exact Test_AGE", pattern="*.csv", full.names=TRUE)

AgeExactList <- lapply(temp, fread, header = TRUE)

obj = readRDS("C:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

AgePadjust <- data.frame()
for(i in 1:length(name)){
    x <- AgeExactList[i][[1]]
    p_adj <- p.adjust(x$p.value, method = "BH")
    xpadj <- cbind(x, p_adj)
    xpadjfin <- xpadj %>% dplyr::filter(p_adj <= 0.05)
    AgePadjust <- rbindlist(list(AgePadjust, xpadjfin))
}

write.csv(AgePadjust, file = "AgePadjust.csv")

#Padjust for Gender

temp <- list.files("E:/Best Model_NON_wrapper/Exact Test_GENDER", pattern="*.csv", full.names=TRUE)
GenderExactList <- lapply(temp, fread, header = TRUE)

#Did this manually because some tissues only have one gender.
#TO DO -- AUTOMATE THIS SECTION
GenderName <- levels(as.factor(c("adipose_subcutaneous", "adipose_visceral_(omentum)", "adrenal_gland", "artery_aorta", "artery_coronary", "artery_tibial", "brain-0", "brain-1", "brain-2", "breast_mammary_tissue", "cells_ebv-transformed_lymphocytes", "cells_transformed_fibroblasts", "colon_sigmoid", "colon_transverse", "esophagus_gastroesophageal_junction", "esophagus_mucosa", "esophagus_muscularis", "heart_atrial_appendage","heart_left_ventricle", "kidney_cortex", "liver", "lung", "minor_salivary_gland", "muscle_skeletal","nerve_tibial","pancreas", "pituitary", "skin", "small_intestine_terminal_ileum", "spleen", "stomach", "thyroid","whole_blood")))      

GenderPadjust <- data.frame()
for(i in 1:length(GenderName)){
    x <- GenderExactList[i][[1]]
    p_adj <- p.adjust(x$p.value, method = "BH")
    xpadj <- cbind(x, p_adj)
    xpadjfin <- xpadj %>% dplyr::filter(p_adj <= 0.05)
    GenderPadjust <- rbindlist(list(GenderPadjust, xpadjfin))
}

write.csv(GenderPadjust, file = "GenderPadjust.csv")

#####Exact test Gender and age comparison #########
FinalAll <- data.frame()
for(i in 1:38){
    EAge <- ExactTestAge %>% filter(Tissue == name[i])
    EGen <- ExactTestGender %>% filter(Tissue == name[i])
    com <- base::intersect(EAge$Gene, EGen$Gene)
    fin <- EAge[!EAge$Gene %in% com,]
    FinalAll <- rbindlist(list(FinalAll, fin))
}

write.csv(FinalAll, "ExactTestAGE_afterGender.csv")
###Main section####
library(data.table)
library(dplyr)
library(grex)
library(ggplot2)
library(naniar)
library(clusterProfiler)
library(org.Hs.eg.db)
library(wesanderson)
maintr <- rbindlist_fread("D:/Ageing Mixture models project/1. New mixture model with wrapper/Best Model_wrapper/Best BIC model")

GTExClean <- function(x){ 
    clean <- cleanid(as.character(x$SelectedGenes))
    SelectedGenes1 <- grex(clean)
    SelectedGenes3 <- cbind(SelectedGenes1, x)
    SelectedGenes2 <- SelectedGenes3[!is.na(SelectedGenes3$entrez_id),]
}

StdDEGenes <- fread("D:/Ageing Mixture models project/1. New mixture model with wrapper/EdgeR/DEAllTissues_gender_All.csv")

StdDEGenesClean <- GTExClean(StdDEGenes)
colnames(StdDEGenesClean)[9] = "SelectedGenes.1"


#write.csv(StdDEGenesClean, "All_EdgeR_Genes.csv")

#Genes detected by edgeR
check2 <- as.data.frame(with(StdDEGenesClean, table(Tissue)))

#write.csv(check2, "AllEdgeRGenes.csv")

#Genes detected by DESeq2
StdDEGenesDESEq <- fread("D:/Ageing Mixture models project/1. New mixture model with wrapper/DESeq2/DESeq2AllTissues_gender.csv")
StdDEGenesDESEq = StdDEGenesDESEq[,-1]
StdDEGenesDESEqClean = GTExClean(StdDEGenesDESEq)

#write.csv(StdDEGenesDESEqClean, "All_DESeq2_Genes.csv")

check3 <- as.data.frame(with(StdDEGenesDESEqClean, 
                             table(Tissue)))

#write.csv(check3, "AllDESeq2Genes.csv")

GTExClean <- function(x){ 
    clean <- cleanid(as.character(x$Gene))
    SelectedGenes1 <- grex(clean)
    SelectedGenes3 <- cbind(SelectedGenes1, x)
    SelectedGenes2 <- SelectedGenes3[!is.na(SelectedGenes3$entrez_id),]
}




MMGenes_GenderFil <- fread("D:/Ageing Mixture models project/1. New mixture model with wrapper/ExactTestAGE_afterGender.csv")

MMGenes_GenderFil <- MMGenes_GenderFil[,-1]

MMGenesClean <- GTExClean(MMGenes_GenderFil)

#Genes detected by mixture models to determine freq of genes
check2 <- as.data.frame(with(MMGenesClean, table(Tissue)))
#write.csv(check2, "AllMMGenes_gender_clean.csv")



#Modes and packages of genes detected by mixture models
new2 <- MMGenes_GenderFil[MMGenes_GenderFil$Gene %in% MMGenesClean$Gene,]

check2 <- as.data.frame(with(MMGenesClean, table(Tissue, Package, No_of_Modes )))

#write.csv(check2, "AllmodesBymixture models.csv")

#Genes that are similar between edgeR and DESeq2

All.all = data.frame()
for(i in 1:38){
    SGenes1 = StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    SGenes2 = StdDEGenesDESEqClean %>% dplyr::filter(Tissue == name[i])
    com = intersect(SGenes1$entrez_id, SGenes2$entrez_id)
    Allgenes = SGenes1[SGenes1$entrez_id %in% com,]
    All.all = rbindlist(list(All.all, Allgenes))
}

#write.csv(All.all, "Common_genes_edgeR_vs_DESeq2.csv")


#Unique genes MM genes in comparison to edgeR
FinalAll <- data.frame()
for(i in 1:38){
    SGenes <- StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    EAge <- MMGenesClean %>% filter(Tissue == name[i])
    com <- base::setdiff(EAge$entrez_id, SGenes$entrez_id)
    Gene <- EAge[EAge$entrez_id %in% com,]
    FinalAll <- rbindlist(list(FinalAll, Gene))
}

write.csv(FinalAll, "UniqueMM_GenderRemoved.csv")

check3 <- as.data.frame(with(FinalAll, table(Tissue)))


#Unique MM genes compared to edgeR and DEseq2

FinalAll1 <- data.frame()
for(i in 1:38){
    SGenes <- StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    EAge <- MMGenesClean %>% filter(Tissue == name[i])
    com <- base::setdiff(EAge$entrez_id, SGenes$entrez_id)
    Gene <- EAge[EAge$entrez_id %in% com,]
    SGenes2 = StdDEGenesDESEqClean %>% dplyr::filter(Tissue == name[i])
    com2 <- base::setdiff(Gene$entrez_id, SGenes2$entrez_id)
    Gene2 <- Gene[Gene$entrez_id %in% com2,]
    FinalAll1 <- rbindlist(list(FinalAll1, Gene2))
}

write.csv(FinalAll1, "UniqueMM_GenderRemoved_edgeR_DESeq2.csv")
#All mixture model genes
GeneCountAllMM <- as.data.frame(with(MMGenesClean, table(Tissue)))
colnames(GeneCountAllMM)[2] <- "MM Genes"
write.csv(GeneCountAllMM, "Genetable_MMALL_GenderRemoved.csv")

#EdgeR Genes
GeneCountEdgeR <- as.data.frame(with(StdDEGenesClean, table(Tissue)))
colnames(GeneCountEdgeR)[2] <- "MM + edgeR Genes"
write.csv(GeneCountEdgeR, "Genetable_EdgeR_GenderRemoved.csv")

#DESeq2 genes
GeneCountDESeq2 = as.data.frame(with(StdDEGenesDESEqClean, table(Tissue)))
colnames(GeneCountDESeq2)[2] <- "MM + DESeq2 Genes"
write.csv(GeneCountDESeq2, 
          "Genetable_DESeq2_GenderRemoved.csv" )

#EdgeR Unique Genes 
GeneCountEdgeRUnique  <- data.frame()
for(i in 1:38){
    SGenes <- StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    EAge <- MMGenesClean %>% filter(Tissue == name[i])
    com <- base::setdiff(SGenes$entrez_id, EAge$entrez_id)
    Gene <- SGenes[SGenes$entrez_id %in% com,]
    SGenes2 <- StdDEGenesDESEqClean %>% 
        dplyr::filter(Tissue == name[i])
    com2 = base::setdiff(Gene$entrez_id, SGenes2$entrez_id)
    Gene2 = Gene[Gene$entrez_id %in% com2,]
    GeneCountEdgeRUnique <- rbindlist(list(GeneCountEdgeRUnique, Gene))
}

GeneCountEdgeRUniquetable <- as.data.frame(with(GeneCountEdgeRUnique, table(Tissue)))

colnames(GeneCountEdgeRUniquetable)[2] <- "edgeR unique Genes"
write.csv(GeneCountEdgeRUniquetable, "Genetable_edgeRUnique_GenderRemoved.csv")

#DESeq2 Unique Genes 
GeneCountDESeq2Unique  <- data.frame()
for(i in 1:38){
    SGenes <- StdDEGenesDESEqClean %>% dplyr::filter(Tissue == name[i])
    EAge <- MMGenesClean %>% filter(Tissue == name[i])
    com <- base::setdiff(SGenes$entrez_id, EAge$entrez_id)
    Gene <- SGenes[SGenes$entrez_id %in% com,]
    SGenes2 <- StdDEGenesClean %>% 
        dplyr::filter(Tissue == name[i])
    com2 = base::setdiff(Gene$entrez_id, SGenes2$entrez_id)
    Gene2 = Gene[Gene$entrez_id %in% com2,]
    GeneCountDESeq2Unique <- rbindlist(
        list(GeneCountDESeq2Unique, Gene))
}

GeneCountDESeq2Uniquetable <- as.data.frame(with(GeneCountDESeq2Unique, table(Tissue)))

colnames(GeneCountDESeq2Uniquetable)[2] <- "DESeq2 unique Genes"
write.csv(GeneCountDESeq2Uniquetable, "Genetable_DESeq2Unique_GenderRemoved.csv")


#Mixture model unique Genes
GeneCount <- as.data.frame(with(FinalAll1, table(Tissue)))
colnames(GeneCount)[2] <- "Mixture model unique Genes"
write.csv(GeneCount, "Genetable_UniqueMM_GenderRemoved.csv")

Fulltable.1 <- full_join(GeneCountAllMM, GeneCountEdgeR, by = "Tissue")

Fulltable.2 = full_join(Fulltable.1, GeneCountDESeq2,
                        by = "Tissue")

Fulltable.3 <- full_join(Fulltable.2, GeneCount, by = "Tissue")

Fulltable.4 <- full_join(Fulltable.3, 
                         GeneCountEdgeRUniquetable, 
                         by = "Tissue")
Fulltable = full_join(Fulltable.4,
                      GeneCountDESeq2Uniquetable,
                      by = "Tissue")

Fulltable$Percent_MMU_genes <- (Fulltable$`Mixture model unique Genes`/Fulltable$`MM Genes`)*100

#Mereging with total number of donors

Number_of_Donors <- as.data.frame(table(pdat$our_subtypes))

colnames(Number_of_Donors) <- c("Tissue", "Number of Donors")

Fulltable <- full_join(Fulltable, Number_of_Donors, by = "Tissue")


Fulltable$Tissue <- gsub("_", " ", Fulltable$Tissue, fixed=TRUE)
Fulltable$Tissue <- tools::toTitleCase(Fulltable$Tissue)
#colnames(Fulltable)[5] <- "edgeR Unique"

write.csv(Fulltable, file = "Mixture model paper FIG 1.csv")

plottable <- Fulltable[,c(1,8)]

plottable[is.na(plottable)] <- 0

plottable$Tissue = with(plottable, reorder(Tissue, Percent_MMU_genes))

"-------Plotting of figure 2A----------"
ggplot(plottable, aes(x = Tissue, y = Percent_MMU_genes)) + 
    geom_bar(stat='identity', fill = wes_palette(n = 1, name = "Moonrise3"))  + 
    labs(x = "Tissue Type", y = "% of MMU Genes") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, 
                                     vjust = 0.45)) 
"Venn diagram for figure 2B to 2D--------------------------"
library(VennDiagram)
#Whole blood 
MMUWB <- MMGenesClean %>% 
    filter(Tissue == "whole_blood") %>% pull(entrez_id) %>% as.numeric

AllWB <- StdDEGenesClean %>% 
    filter(Tissue == "whole_blood") %>% pull(entrez_id) %>% as.numeric

AllDESeq = StdDEGenesDESEqClean %>% filter(Tissue == "whole_blood") %>% pull(entrez_id) %>% as.numeric


venn.diagram(
    x = list(MMUWB, AllWB, AllDESeq),
    category.names = c("MM genes" , 
                       "edgeR genes", "DESeq2 genes"),
    filename = 'Whole blood MMU genes.tiff',
    output=TRUE,
    imagetype="tiff" ,
    height = 6400 , 
    width = 6400 , 
    resolution = 1200,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = wes_palette(n = 3, name = "Moonrise2"),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(5, 150),
    #cat.dist = c(0.02, 0.03),
    cat.fontfamily = "sans",
    cat.col = wes_palette(n = 3, name = "Moonrise2")
)

#Skin
MMuskin <- MMGenesClean %>% 
    filter(Tissue == "skin") %>% pull(entrez_id) %>% as.numeric

Allskin <- StdDEGenesClean %>% 
    filter(Tissue == "skin") %>% pull(entrez_id) %>% as.numeric

AllDESeq = StdDEGenesDESEqClean %>% filter(Tissue == "skin") %>% pull(entrez_id) %>% as.numeric

venn.diagram(
    x = list(MMuskin, Allskin, AllDESeq),
    category.names = c("MM genes" , 
                       "edgeR genes", "DESeq2 genes"),
    filename = 'Skin MMU genes.tiff',
    output=TRUE,
    imagetype="tiff" ,
    height = 6400 , 
    width = 6800 , 
    resolution = 1200,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = wes_palette(n = 3, name = "Chevalier1"),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(5, 1),
    #cat.dist = c(0.02, 0.03),
    cat.fontfamily = "sans",
    cat.col = wes_palette(n = 3, name = "Chevalier1")
)

#Brain-0 
MMUB0 <- MMGenesClean %>% 
    filter(Tissue == "brain-0") %>% pull(entrez_id) %>% as.numeric

AllB0 <- StdDEGenesClean %>% 
    filter(Tissue == "brain-0") %>% pull(entrez_id) %>% as.numeric

AllDESeq = StdDEGenesDESEqClean %>% filter(Tissue == "brain-0") %>% pull(entrez_id) %>% as.numeric

venn.diagram(
    x = list(MMUB0, AllB0, AllDESeq),
    category.names = c("MM genes" , 
                       "edgeR genes", "DESeq2 genes"),
    filename = 'Brain-0 MMU genes.tiff',
    output=TRUE,
    imagetype="tiff" ,
    height = 6400 , 
    width = 6400 , 
    resolution = 1200,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = wes_palette(n = 3, name = "Darjeeling2"),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    #cat.pos = c(5, 150),
    #cat.dist = c(0.02, 0.03),
    cat.fontfamily = "sans",
    cat.col = wes_palette(n = 3, name = "Darjeeling2")
)

####Determining packages from selected genes#####

t1 <- read.csv("D:/Ageing Mixture models project/1. New mixture model with wrapper/Best_Model_Fin.csv")
PackageTypeTable <- as.data.frame(with(t1, table(Tissue, Package)))

PackageTypeTable$Tissue <- gsub("_", " ", PackageTypeTable$Tissue, fixed=TRUE)
PackageTypeTable$Tissue <- tools::toTitleCase(PackageTypeTable$Tissue)

dat.p <- reshape2::melt(PackageTypeTable)

ggplot(dat.p, aes(x = Tissue, y = value , fill= Package)) + geom_bar(stat='identity')   + theme_classic() + ylim(0, 30000) +
    scale_fill_manual(values= wes_palette("GrandBudapest2", n = 2)) + 
    ylab("Numer of Genes")+ theme(axis.text.x = element_text(angle = 90, hjust=1)) 

PackageTypeTable <- as.data.frame(with(MMGenesClean, table(Tissue, Package, No_of_Modes)))

PackageTypeTable$Tissue <- gsub("_", " ", PackageTypeTable$Tissue, fixed=TRUE)
PackageTypeTable$Tissue <- tools::toTitleCase(PackageTypeTable$Tissue)


write.csv(PackageTypeTable, "Selected_MM_Table.csv")
dat.p <- reshape2::melt(PackageTypeTable)
Pubtheme <- function(base_size = 11, base_family = "Helvetica", ticks = TRUE) {
    
    ret <- theme_bw(base_family = , base_size = base_size) + 
        theme(
            axis.line = element_line(color = 'black'),
            axis.title.x = element_text(vjust = -0.3), 
            axis.title.y = element_text(vjust = 0.8),
            legend.background = element_blank(), 
            legend.key = element_blank(), 
            legend.title = element_text(face="plain"),
            panel.background = element_blank(), 
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.background = element_blank(),
            strip.background = element_blank()
        )
    ret              
}

ggplot(dat.p, aes(x = Tissue, y = value , fill= Package)) + geom_bar(stat='identity')   + theme_classic() + 
    scale_fill_manual(values= wes_palette("GrandBudapest1", n = 2)) + 
    ylab("Numer of Genes")+ theme(axis.text.x = element_text(angle = 90)) 


PackageTypeTable <- as.data.frame(with(FinalAll, table(Tissue, Package, No_of_Modes)))

PackageTypeTable$Tissue <- gsub("_", " ", PackageTypeTable$Tissue, fixed=TRUE)
PackageTypeTable$Tissue <- tools::toTitleCase(PackageTypeTable$Tissue)

write.csv(PackageTypeTable, "Selected_MM_Table.csv")
dat.p <- reshape2::melt(PackageTypeTable)
Pubtheme <- function(base_size = 11, base_family = "Helvetica", ticks = TRUE) {
    
    ret <- theme_bw(base_family = , base_size = base_size) + 
        theme(
            axis.line = element_line(color = 'black'),
            axis.title.x = element_text(vjust = -0.3), 
            axis.title.y = element_text(vjust = 0.8),
            legend.background = element_blank(), 
            legend.key = element_blank(), 
            legend.title = element_text(face="plain"),
            panel.background = element_blank(), 
            panel.border = element_blank(),
            panel.grid = element_blank(),
            plot.background = element_blank(),
            strip.background = element_blank()
        )
    ret              
}

ggplot(dat.p, aes(x = Tissue, y = value , fill= Package)) + geom_bar(stat='identity')   + Pubtheme() + 
    scale_fill_manual(values= wes_palette("GrandBudapest1", n = 2)) + 
    ylab("Numer of Genes")+ theme(axis.text.x = element_text(angle = 90)) 



GeneTypeTable <- as.data.frame(with(FinalAll, table(gene_biotype, Tissue)))

dat.p <- reshape2::melt(GeneTypeTable)
ggplot(dat.p, aes(x = Tissue, y = value , fill= gene_biotype)) + geom_bar(stat='identity')  + theme(axis.text.x = element_text(angle = 90))


GeneTypeTable1 <- as.data.frame(table(GeneTypeTable$Tissue))
NewName <- levels(GeneTypeTable1$Var1)

PercentageTable <- data.frame()
for (i in 1:length(NewName)) {
    Tab1 <- GeneTypeTable %>% filter(Tissue == NewName[i])
    Tab1$PCT <- prop.table(Tab1$Freq)
    Tab1$PCT <- Tab1$PCT *100
    PercentageTable <- rbindlist(list(PercentageTable, Tab1))
}

PercentageTable$Freq <- NULL

dat.p <- reshape2::melt(PercentageTable)
ggplot(dat.p, aes(x = Tissue, y = value , fill= gene_biotype)) + geom_bar(stat='identity')  + theme(axis.text.x = element_text(angle = 90))


Test1 <- FinalAll %>% filter(Tissue == NewName[i])
write.csv(Test1, "Test1.csv")


####Comparing to external ageing databases####
#GeneAge database
GeneAgeDB <- read_csv("D:/Ageing Mixture models project/Paper/R/human_genes/genage_human.csv")

#Digitial ageing atals
DAA <- read_csv("D:/Ageing Mixture models project/Paper/R/digital_ageing_atlas_data/digital_ageing_atlas_data1.csv") #4250

DAAHumans <- DAA %>% filter(Species == "Homo sapiens") #3574
DAAHumans[DAAHumans == "#VALUE!"] <- NA
DAAHumans <- DAAHumans[!is.na(DAAHumans$GeneId),]


#Converting Gene symbol to enterzid for DAA humans#
library(org.Hs.eg.db)

DAAGenes <- mapIds(org.Hs.eg.db, DAAHumans$GeneId, 'ENTREZID', 'SYMBOL') 

library(SuperExactTest)

#Mixture model unique genes
x <- list(FinalAll$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes)

total <- sum((length.gene.sets=sapply(x,length)))

res=supertest(x, n=total)

library(RColorBrewer)
cols <- brewer.pal(9, "BuPu")
pal <- colorRampPalette(cols)


plot(res, Layout = "landscape", degree=1:4, 
     color.scale.title = expression(paste(-Log[10],'(',italic(P),')')),
     color.on ='#834C81', color.off ='#fefcfc',
     sort.by = "size", margin = c(0.5, 5, 1, 3), heatmapColor = cols)
dev.print(pdf, 'MMU_CommonGenes_with_databases.pdf')

Databaseintersect <- Reduce(intersect, x)

#52 genes that are common between the databases

DAAGenesnew <- as.data.frame(as.matrix(DAAGenes))
DAAGenesnew$SYMBOL <- rownames(DAAGenesnew)
colnames(DAAGenesnew)[1] <- "ENTREZID"

ComGenes <- DAAGenesnew[DAAGenesnew$ENTREZID %in% Databaseintersect,]

newgenes <- DAAHumans[DAAHumans$GeneId %in% ComGenes$SYMBOL,]

AllNewGenes <- FinalAll[FinalAll$hgnc_symbol %in% newgenes$GeneId,]

base::unique(AllNewGenes$hgnc_symbol)
mainone <- with(AllNewGenes, table(hgnc_symbol, Tissue))
write.csv(mainone, "Common_MMU_genes_with_Databases.csv")

#TweleveCommonGenes <- GeneAgeDB[GeneAgeDB$`entrez gene id` %in% Databaseintersect, ]

#EdgeR unique genes
x <- list(StdDEGenesClean$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes)

total <- sum((length.gene.sets=sapply(x,length)))

res=supertest(x, n=total)

library(RColorBrewer)
cols <- brewer.pal(9, "BuPu")

plot(res, Layout = "landscape", degree=1:4, 
     color.scale.title = expression(paste(-Log[10],'(',italic(P),')')),
     color.on ='#834C81', color.off ='#fefcfc',
     sort.by="size", margin=c(0.5, 5, 1, 3), heatmapColor = cols)

dev.print(pdf, 'EdgeRGene+MMU_genes_common_with_databases.pdf')

Databaseintersect <- Reduce(intersect, x)

#66 genes that are common between the databases

DAAGenesnew <- as.data.frame(as.matrix(DAAGenes))
DAAGenesnew$SYMBOL <- rownames(DAAGenesnew)
colnames(DAAGenesnew)[1] <- "ENTREZID"

ComGenes <- DAAGenesnew[DAAGenesnew$ENTREZID %in% Databaseintersect,]

newgenes <- DAAHumans[DAAHumans$GeneId %in% ComGenes$SYMBOL,]

AllNewGenes <- StdDEGenesClean[StdDEGenesClean$hgnc_symbol %in% newgenes$GeneId,]

base::unique(AllNewGenes$hgnc_symbol)
mainone <- with(AllNewGenes, table(hgnc_symbol, Tissue))
write.csv(mainone, "Common_EdgeR_genes_with_Databases.csv")
#
Databaseintersect <- Reduce(intersect, x)

#66 genes that are common between the databases

#DESeq2 genes superexact test

DAAGenesnew <- as.data.frame(as.matrix(DAAGenes))
DAAGenesnew$SYMBOL <- rownames(DAAGenesnew)
colnames(DAAGenesnew)[1] <- "ENTREZID"

ComGenes <- DAAGenesnew[DAAGenesnew$ENTREZID %in% Databaseintersect,]
dim(ComGenes)
#Finding out number of common genes for MM only and EdgeR that are consitent between datasets

x <- list(FinalAll$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes)
DatabaseintersectMM <- Reduce(intersect, x)

GeneDAAExlcusiveMM <- FinalAll[FinalAll$entrez_id %in% DatabaseintersectMM,]
write.csv(GeneDAAExlcusiveMM, "MMonly_DAAGenes.csv")

uniMMGENES <- setdiff(DatabaseintersectMM, DatabaseintersectEdge)

uniMMGENESread <- FinalAll[FinalAll$entrez_id %in% uniMMGENES,]
write.csv(uniMMGENESread, "Unique_MMonly_DAAGenes.csv")


x <- list(StdDEGenesClean$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes)
DatabaseintersectEdge <- Reduce(intersect, x)

GeneDAAExlcusiveEdage<- StdDEGenesClean[StdDEGenesClean$entrez_id %in% DatabaseintersectMM,]
write.csv(GeneDAAExlcusiveEdage, "EdgeR_DAAGenes.csv")

uniEdgeGENES <- setdiff(DatabaseintersectEdge, DatabaseintersectMM)

uniMMGENESread <- StdDEGenesClean[StdDEGenesClean$entrez_id %in% uniEdgeGENES,]
write.csv(uniMMGENESread, "Unique_EdgeRonly_DAAGenes.csv")




####Pathway over representation####

#Loading the different gene sets

StdDEGenesClean <- fread("D:/Ageing Mixture models project/1. New mixture model with wrapper/EdgeR/All_EdgeR_Genes.csv")
StdDEGenesClean[,1] = NULL


#Genes detected by DESeq2

StdDEGenesDESEqClean = fread("D:/Ageing Mixture models project/All_DESeq2_Genes.csv")
StdDEGenesDESEqClean[,1] = NULL


#GO pathways
GOForGenes <- function(dat){
    grph <- enrichGO(dat$entrez_id, 
                     OrgDb = org.Hs.eg.db , keyType = "ENTREZID", ont = "BP",
                     pvalueCutoff = 0.05, pAdjustMethod = "fdr",
                     qvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 500,
                     readable = FALSE, pool = FALSE)
}


Fint1 <- data.table()
for(i in 1:length(NewName)){
    x <- MMGenesClean %>% filter(Tissue == NewName[i])
    grph <- GOForGenes(x)
    if(nrow(grph) == 0) next
    dotplot(grph, showCategory = 15)
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(NewName[i], nrow(t1))
    Fint1 <- rbindlist(list(Fint1, t1))
    
}

write.csv(Fint1, file = "Go_MM_ALL_Pathways.csv")


Fint1Table <- as.data.frame(table(Fint1$Tissue))

colnames(Fint1Table) <- c("Tissue", "Mixture Model Pathways")

#For edgeR
name = levels(as.factor(unique(StdDEGenesClean$Tissue)))


FintEdge <- data.table()
for(i in 1:length(name)){
    x <- StdDEGenesClean %>% filter(Tissue == name[i])
    if(nrow(x) == 0) next 
    grph <- GOForGenes(x)
    if(nrow(grph) == 0) next
    dotplot(grph, showCategory = 15)
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(name[i], nrow(t1))
    FintEdge <- rbindlist(list(FintEdge, t1))
    
}

write.csv(FintEdge, file = "Go_EdgeR_Pathways.csv")

FintEdgeTable <- as.data.frame(table(FintEdge$Tissue))
colnames(FintEdgeTable) <- c("Tissue", "EdgeR Pathways")


#DESeq2 pathways
FintDESeq2 <- data.table()
for(i in 1:length(name)){
    x <- StdDEGenesDESEqClean %>% filter(Tissue == name[i])
    if(nrow(x) == 0) next 
    grph <- GOForGenes(x)
    if(nrow(grph) == 0) next
    dotplot(grph, showCategory = 15)
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(name[i], nrow(t1))
    FintDESeq2 <- rbindlist(list(FintDESeq2, t1))
    
}

write.csv(FintDESeq2, file = "Go_DESeq2_Pathways.csv")

FintDESeqTable <- as.data.frame(table(FintDESeq2$Tissue))
colnames(FintDESeqTable) <- c("Tissue", "DESeq2 Pathways")

All.all = data.frame()
for(i in 1:38){
    SGenes1 = StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    SGenes2 = StdDEGenesDESEqClean %>% dplyr::filter(Tissue == name[i])
    com = intersect(SGenes1$entrez_id, SGenes2$entrez_id)
    Allgenes = SGenes1[SGenes1$entrez_id %in% com,]
    All.all = rbindlist(list(All.all, Allgenes))
}


AllCommonpathways <- NULL
Commonpathways <- NULL
for (i in 1:length(name)) {
    y <- Fint1 %>% filter(Tissue == name[i])
    if(nrow(y) == 0) next
    x <- FintEdge %>% filter(Tissue == name[i])
    coone <- x[x$ID %in% y$ID,]
    x1 = FintDESeq2 %>% filter(Tissue == name[i])
    coone1 = coone[coone$ID %in% x1$ID]
    AllCommonpathways <- rbind(AllCommonpathways, coone1)
    co <- as.data.frame(length(Reduce(intersect,
                        list(y$Description, 
                        x$Description,
                        x1$Description))))
    co$Tissue <- rep(name[i], nrow(co))
    Commonpathways <- rbind(Commonpathways, co)
    
}

write.csv(AllCommonpathways, "GO_common_pathways_with_ALL_MMGenes.csv")

colnames(Commonpathways) <- c("Mixture model Common Pathways", "Tissue")


FullCommon.1 <- full_join(Fint1Table, FintEdgeTable, by = "Tissue")
Fullcommon.2 = full_join(FullCommon.1, FintDESeqTable, by = "Tissue")


FullCommon1 <- full_join(Fullcommon.2, Commonpathways, by = "Tissue")

write.csv(FullCommon1, file="GO_Number_of_common_pathways.csv")

#Hallmark
library(org.Hs.eg.db)
library(msigdbr)

HallMarkforGenes <- function(dat){
    m_df = msigdbr(species = "Homo sapiens", category = "H")
    m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    grph <- enricher(dat$entrez_id, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                     minGSSize = 20, maxGSSize = 500, 
                     qvalueCutoff = 0.05,TERM2GENE = m_t2g,
                     TERM2NAME = NA)
    return(grph)
}

Fint1 <- data.table()
for(i in 1:length(NewName)){
    x <- MMGenes %>% filter(Tissue == NewName[i])
    grph <- HallMarkforGenes(x)
    if(is.null(grph)) next
    if(nrow(grph) == 0) next
    dotplot(grph, showCategory = 15)
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(NewName[i], nrow(t1))
    Fint1 <- rbindlist(list(Fint1, t1))
    
}

Fint1Table <- as.data.frame(table(Fint1$Tissue))

write.csv(Fint1, file = "Hallmark_MM_ALL_Pathways.csv")

colnames(Fint1Table) <- c("Tissue", "Mixture Model Pathways")

FintEdge <- data.table()
for(i in 1:length(name)){
    x <- StdDEGenesClean %>% filter(Tissue == name[i])
    if(nrow(x) == 0) next 
    grph <- HallMarkforGenes(x)
    if(nrow(grph) == 0) next
    dotplot(grph, showCategory = 15)
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(name[i], nrow(t1))
    FintEdge <- rbindlist(list(FintEdge, t1))
    
}

write.csv(FintEdge, file = "Hallmark_EdgeR_Pathways.csv")

FintEdgeTable <- as.data.frame(table(FintEdge$Tissue))
colnames(FintEdgeTable) <- c("Tissue", "EdgeR Pathways")


AllCommonpathways <- NULL
Commonpathways <- NULL
for (i in 1:length(name)) {
    y <- Fint1 %>% filter(Tissue == name[i])
    if(nrow(y) == 0) next
    x <- FintEdge %>% filter(Tissue == name[i])
    coone <- x[x$ID %in% y$ID,]
    AllCommonpathways <- rbind(AllCommonpathways, coone)
    co <- as.data.frame(length(base::intersect(y$Description , x$Description)))
    co$Tissue <- rep(name[i], nrow(co))
    Commonpathways <- rbind(Commonpathways, co)
    
}

write.csv(AllCommonpathways, "HallMark_common_pathways_MM_ALL.csv")

colnames(Commonpathways) <- c("Mixture model Common Pathways", "Tissue")




FullCommon.1 <- full_join(Fint1Table, FintEdgeTable, by = "Tissue")

FullCommon2 <- full_join(FullCommon.1, Commonpathways, by = "Tissue")

write.csv(FullCommon2, file="HallMark_Number_of_common_pathways.csv")

#WikiPathways

wp2gene <- read.gmt("D:/Ageing Mixture models project/1. New mixture model with wrapper/wikipathways-20200910-gmt-Homo_sapiens.gmt")

#wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


Fint1 <- data.table()
for(i in 1:38){
    x <- MMGenesClean %>% filter(Tissue == name[i])
    grph <- enricher(x$entrez_id, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
    if(is.null(grph)) next
    if(nrow(grph) == 0) next
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(name[i], nrow(t1))
    Fint1 <- rbindlist(list(Fint1, t1))
}

write.csv(Fint1, file = "Wiki_MMonly_Pathways.csv")

Fint1Table <- as.data.frame(table(Fint1$Tissue))

colnames(Fint1Table) <- c("Tissue", "Mixture Model Pathways")

FintEdge <- data.table()
for(i in 1:38){
    x <- StdDEGenesClean %>% filter(Tissue == name[i])
    grph <- enricher(x$entrez_id, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
    if(is.null(grph)) next
    if(nrow(grph) == 0) next
    t1 <- grph@result[grph@result$p.adjust < 0.05,]
    t1[is.na(t1)] <- 0
    t1$Tissue <- rep(name[i], nrow(t1))
    FintEdge <- rbindlist(list(FintEdge, t1))
}

write.csv(FintEdge, file = "Wiki_EdgeR_Pathways.csv")

FintEdgeTable <- as.data.frame(table(FintEdge$Tissue))
colnames(FintEdgeTable) <- c("Tissue", "EdgeR Pathways")


AllCommonpathways <- NULL
Commonpathways <- NULL
for (i in 1:length(name)) {
    y <- Fint1 %>% filter(Tissue == name[i])
    if(nrow(y) == 0) next
    x <- FintEdge %>% filter(Tissue == name[i])
    coone <- x[x$ID %in% y$ID,]
    AllCommonpathways <- rbind(AllCommonpathways, coone)
    co <- as.data.frame(length(base::intersect(y$Description , x$Description)))
    co$Tissue <- rep(name[i], nrow(co))
    Commonpathways <- rbind(Commonpathways, co)
    
}


write.csv(AllCommonpathways, "Wiki_common_pathways.csv")

colnames(Commonpathways) <- c("Mixture model Common Pathways", "Tissue")


FullCommon.1 <- full_join(Fint1Table, FintEdgeTable, by = "Tissue")

FullCommon3 <- full_join(FullCommon.1, Commonpathways, by = "Tissue")

write.csv(FullCommon3, file="WikiPathways_Number_of_common_pathways.csv")

AllCommon.1 <- full_join(FullCommon1, FullCommon2, by = "Tissue")
AllCommon <-  full_join(AllCommon.1, FullCommon3, by = "Tissue")

write.csv(AllCommon, file="Number_of_common_pathways.csv")
library(mclust)
library(RmixmodCombi)
library(MixAll)
library(dplyr)
library(fitdistrplus)

BestMM <- function(dat = dat, no_of_mode = x ){
    
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



size <- levels(as.factor(seq(0.1, 10, by = 0.2)))
means <- levels(as.factor(seq(0.1, 10, by = 0.2)))


m <- matrix(NA, ncol = 50, nrow = 50)
for (i in 1:length(size)) {
    for (j in 1:length(means)) {
        set.seed(1234)
        x <- round(rnbinom(n = 200, size = as.numeric(size[i]), 
                           mu = as.numeric(means[j])), 4)
        be <- BestMM(x, 1:3)
        m[i,j] <-  be$Model
    }
}
write.csv(m, file = "NB_simulation_matrix.csv")

m1 <- read.csv("NB_simulation_matrix.csv")
m1 <- m1[,-1]

mm <- as.data.frame(m1)

mmmm <- matrix(NA, ncol = 50, nrow = 50)
for (i in 1:length(size)) {
    for (j in 1:length(means)){
        x <- mm[i,j]
        if (grepl("gamma", x) == TRUE) {
            tt = 1 
        } else if( grepl("Gaussian", x) == TRUE){
            tt = 2 
        } else {
            tt = 3
        }    
        mmmm[i,j] <- tt
    }
}

mmmm <- as.data.frame(mmmm)
rownames(mmmm) <- seq(0.1, 10, by = 0.2)
colnames(mmmm) <- seq(0.1, 10, by = 0.2)

library(pheatmap)
library(RColorBrewer)
pheatmap(mmmm, cluster_rows = F, cluster_cols = F,
         display_numbers = F,  color = brewer.pal(n = 3, name = "Dark2"),
         breaks = c(0, 1, 2, 3), legend_breaks = c(0, 1, 2, 3))

