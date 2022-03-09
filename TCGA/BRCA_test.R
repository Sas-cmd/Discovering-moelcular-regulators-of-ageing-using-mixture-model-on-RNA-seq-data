library(Biobase)
library(SummarizedExperiment)
library(mclust)
library(RmixmodCombi)
library(MixAll)
library(data.table)
library(dplyr)
library(ggplot2)

setwd("D:/Ageing Mixture models project/1. New mixture model with wrapper/1. Testing against a ref/TCGA_BRCA/TCGA_BRCA")

####Load data and extract count data####
load("D:/Ageing Mixture models project/1. New mixture model with wrapper/1. Testing against a ref/TCGA_BRCA/TCGA_BRCA/TCGA_BRCA.rda")

#saveRDS(data, file = "BRCA.rds")
BRCA.counts <- assays(data)$normalized_count
#saveRDS(BRCA.counts, file = "BRCA_counts.rds")

####Mixture model pipeline####

BRCA.counts[BRCA.counts == 0] <- NA
dim(BRCA.counts)
#19947  1215
BRCA.countsClean <- BRCA.counts[rowSums(is.na(BRCA.counts)) <= floor(length(BRCA.counts[1,])*0.1), ]
dim(BRCA.countsClean)
#15717  1215

#Subsetting by subtypes 

table(data$paper_BRCA_Subtype_PAM50)

# Basal   Her2   LumA   LumB Normal 
# 190     82    562    209     40 


#Sub-setting BRCA-subtype

BRCA_Subtype <- data$paper_BRCA_Subtype_PAM50


BRCA.counts.three <- BRCA.countsClean[,which(BRCA_Subtype %in% c("Basal","Her2","LumB"))]
dim(BRCA.counts.three)
#15717   481
####Running mixture models####
#MM wrapper
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
    
    Mainoutput <- list("Classification" = BestOutput, 
                       "BIC" = BestBIC, "Model" = BestModel)
    return(Mainoutput) 
}

#MM Run wrapper
MMRun <- function( dat = dat, no_of_mode = x, Partion_name = Partion_name, BIC_Name = BIC_Name){
    AllBIC <- data.table()
    Allpartion <- NULL
    for (j in 1:nrow(dat)) {
        mm <- BestMM(dat = dat[j,], no_of_mode = 1:3)
        t1 <- as.matrix(mm$Classification)
        colnames(t1) <- rownames(dat)[j]
        t2 <- data.table(t(c(rownames(dat)[j], mm$BIC, mm$Model)))
        Allpartion <- as.matrix(cbind(Allpartion, t1))
        AllBIC <- rbindlist(list(AllBIC, t2))
        print(j)
    }



    rownames(Allpartion) <- colnames(dat)
    Allpartion <- as.data.frame(Allpartion)

    write.csv(Allpartion, file = Partion_name)
    write.csv(AllBIC, file = BIC_Name)
    
    finlist <- list(Allpartion, AllBIC)
return(finlist)
}

OutputList <- MMRun(BRCA.counts.three, no_of_mode = 1:3, 
          Partion_name = "BRCA_Basal_Her2_LumB.csv", 
          BIC_Name = "BIC_BRCA_Basal_Her2_LumB.csv") 

saveRDS(OutputList, file = "BRCA_BRCA_Basal_Her2_LumB.rds")

####Exact test####
#load data for BIC and partion

load("D:/Ageing Mixture models project/1. New mixture model with wrapper/1. Testing against a ref/TCGA_BRCA/TCGA_BRCA/TCGA_BRCA.rda")

BICall <- read.csv("BIC_BRCA_Basal_Her2_LumB.csv")

maindata <- read.csv("BRCA_Basal_Her2_LumB.csv")
colnames(maindata)[colnames(maindata) == 'X'] <- 'Patient_Id'




#Catagorical values in a dataframe

xx <- as.data.frame(data$paper_BRCA_Subtype_PAM50)
dim(xx)
rownames(xx) <- colnames(data)
xx$Patients_Id <- rownames(xx)
colnames(xx) <- c("Sub_type", "Patient_Id")


ETest <- function(dat = dat, cat.table = cat.table){

    Finaloutput <- data.table()
    for (i in 2:ncol(dat)) {
        oneG <- as.data.frame(dat[,c(1,i)])
        x1 <- merge(oneG, cat.table, 
                    by.x = "Patient_Id", all.x = FALSE)
        colnames(x1) <- c("Patient_Id", "Gene", "Sub_type")
        if(length(unique(x1$Gene)) == 1) next
        t1 <- as.matrix(with(x1, table(Sub_type, Gene)))
        pt <- fisher.test(t1, workspace = 2e8)
        pv <- pt$p.value
        
        output1 <- as.data.frame(t(c(colnames(oneG)[2], 
                                   pv, (ncol(t1)))))
        colnames(output1) <- c("Gene", "p.value", "No_of_Modes")
        
        
        
        Finaloutput <- rbindlist(list(Finaloutput, output1))
        print(paste0(i, "/", ncol(dat)))
    }
 return(Finaloutput)
}

BRCA_etest1 <- ETest(dat = maindata[1:20,], cat.table = xx)
dim(BRCA_etest1)
dim(maindata)
15718 - 15700

write.csv(BRCA_etest1, "BRCA_E.test.csv" )


BRCA_etest1 <- read.csv("BRCA_E.test.csv")
p_adj <- p.adjust(BRCA_etest1$p.value, method = "BH")
BRCA_etestFIN <- cbind(BRCA_etest1, p_adj)

write.csv(BRCA_etestFIN, "BRCA_E_Final.test.csv")

#Mistake with merging, redoing here to make it accurate
BRCA_etest2 = BRCA_etest1[,-c(7:9)]
BICall2 = BICall[,-1]
colnames(BICall2) = c("Gene", "BIC", "Models")
Finout = merge(BRCA_etest2, BICall2, by = "Gene" )
Finout = Finout[,-c(2,3,6)]
write.csv(Finout, file = "BRCA_MM_results.csv")
###############
library(RcmdrMisc)
library(viridis)
library(hrbrthemes)

load("C:/Users/s4504544/Desktop/pam50.rda")
P50 <- as.data.frame(rownames(pam50$centroids))   

xpadjfin1 <- BRCA_etestFIN %>% dplyr::filter(p_adj <= 0.05)

commonGene <- intersect(P50$`rownames(pam50$centroids)`, xpadjfin1$Gene)

NotP50Gene <- P50[!(P50$`rownames(pam50$centroids)` %in% commonGene),]

Genetable <- function(maindata = maindata,
                      cat.table = cat.table,
                      Gene = Gene){
    mainsub <- maindata[,c(1, which(colnames(maindata) %in% Gene))]
    tt1 <- as.data.frame(merge(mainsub, cat.table, by.x = "Patient_Id", all.x = FALSE))
    
    colnames(tt1) <- c("Patient_Id", "Sub_type", Gene)
    tt1$Patient_Id <- NULL
    
    m1 <- as.matrix(table(tt1))
    return(m1)
}

Genetable(maindata, xx, "EGF")


finfin <- data.table()
for (i in 1:length(commonGene)){
    
    ttt <- as.matrix(Genetable(maindata = maindata, cat.table = xx, 
                               Gene = commonGene[i]))
    
    if (nrow(ttt) == 2) {
        ttt1 <- as.data.frame(rowPercents(ttt, digits=2)[1:2,1:3])
    } else {
        ttt1 <- as.data.frame(rowPercents(ttt, digits=2)[1:3,1:3])
    }
    ttt1$Gene <- rep(x = commonGene[i], nrow(ttt1))
    finfin <- rbindlist(list(finfin, ttt1))

}

table(finfin$Gene)

finfin1 <- data.table()
for (i in 1:length(commonGene)) {
    a1 <- finfin %>% filter(Gene == commonGene[i])
    a1$group <- paste0("_C", as.character(c(1:nrow(a1))))
    a1$Gene.group <- paste0(commonGene[i], a1$group )
    a1$group <- NULL
    a1$Gene <- NULL
    a1 <- melt(a1)
    finfin1 <- rbindlist(list(finfin1, a1))
}



ggplot(finfin1, aes(fill=variable, y=value, x= Gene.group)) + 
    geom_bar(position="stack", stat="identity") +
    ggtitle("Clustering of Subtypes by mixture models") +
    theme(axis.text.x = element_text(angle = 90))+
    scale_color_ipsum() +
    scale_fill_ipsum() +
    xlab("")


ggplot(finfin1, aes(x="", y=prop, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none")

finfin$number <- rep(1:3, )

setwd("E:/BRCA_photos")
for (i in 1:nrow(finfin)){
    plpl <- as.data.frame(t(finfin[i, 1:3]))
    plpl$Sub_type <- rownames(plpl)
    ggplot(plpl, aes(x="", y=V1, fill=Sub_type)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) + theme_void() + 
        scale_fill_ipsum() + theme(legend.position = "none")
    ggsave(paste(finfin$Gene[i], i, ".pdf", sep = ""), dpi = 500, scale = 1,
           width = 20, height = 20, units = "cm" )
    print(i)
  
}
    
#######
# for (i in 1:137) {
#     file_name = paste(finfin$Gene[i], i, ".tiff", sep="")
#     tiff(file_name)
#     print(plot_list[[i]])
#     dev.off()
# }
#     pop <- paste0("_Cluster_", rownames(finfin)[i])
#     pop1 <- paste0(finfin$Gene[i], pop)
# 
# 
#     pdf(file= paste0(pop1, ".pdf"))
# 
# 
#     
#         #scale_fill_manual(values=c("#d3d3d3", "#00008b", "#FFB6C1"))
#     
#     dev.off()
#     print(paste0(i, "/", nrow(finfin)))
# }
# 
# #tiff("FileName.tiff", height = 30, width = 20, units= 'cm' ,
#      #compression = "lzw", res = 300)
# #
# ####
# tt1 <- as.data.frame(cbind(maindata$Patient_Id, maindata[,which(colnames(maindata) %in% "ORC6L")]))
# 
# colnames(tt1) <- c("Patient_Id", "ORC6L")
# 
# x1 <- merge(tt1, xx, by.x = "Patient_Id", all.x = FALSE)
# 
# m1 <- as.matrix(with(x1, table(Sub_type, ORC6L)))
# 
# 
# 
# 
# 
# 
# 
# t1 <- as.matrix(with(x1, table(Sub_type, CCNB1)))
# pt <- fisher.test(t1, workspace = 2e8)
# pt$p.value
# 
# brca_subtypes <- TCGAbiolinks::TCGAquery_subtype("brca")
# 
# 
# 
# # opo <- as.data.frame(ppp[, 1, drop = FALSE])
# # opo$Patient_Id <- rownames(opo)
# 
# 
# 
# # ppp <- OutputList[[1]]
# # ppp$Patient_Id <- rownames(ppp)
# 
# x2 <- merge(ppp, xx, by.x = "Patient_Id", all.x = FALSE)
# class(x2$Sub_type)
# 
# opo <- as.data.frame(ppp[, 1, drop = FALSE])
# opo$Patient_Id <- rownames(opo)
# 
# x2 <- merge(opo, xx, by.x = "Patient_Id", all.x = FALSE)
# colnames(x2) <- c("Patient_Id", "Gene", "Sub_type")
# 
# t1 <- as.matrix(with(x2, table(Sub_type, Gene)))
# pt <- fisher.test(t1, workspace = 2e8)
# pv <- pt$p.value
# 
# output1 <- as.data.table(c(colnames(ppp)[1], pv, (ncol(t1)),
#                            BICall[1,]))
# colnames(FinOutput) <- c("Gene", "p.value", "No_of_Modes", "row.no", 
#                          "Gene2", "BIC", "Model")
# 
# Finaloutput <- rbindlist(list(Finaloutput, output1))