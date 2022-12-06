library(DESeq2)
library(Biobase)
library(data.table)


obj = readRDS("D:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType<-unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

setwd("D:/Ageing Mixture models project/1. New mixture model with wrapper/DESeq2")

start.time = Sys.time()
DESeq2ModesAll <- NULL
for(i in c(1:25,27,28,30,31,32,33,35,38)){ #
    edat.sel <- NULL
    Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
    pos <- which(colnames(edat) %in% rownames(Tissue)) 
    edat.Tiss <- edat[,pos]
    edat.Tiss[edat.Tiss == 0] <- NA
    #Removing genes that have less and equal to 10% zero values
    edat.Tiss <- edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
    edat.Tiss[is.na(edat.Tiss)] <- 0
    
    edat.sel <- t(as.matrix(edat.Tiss))
    tn.df <- as.data.frame(cbind(
        Tissue$AGE, Tissue$GENDER, Tissue$DTHHRDY, edat.sel))
    colnames(tn.df)[1:3] <- c("AGE", "GENDER", "DTHHRDY")
    age_map <- c(1,1,2,2,3,3)
    tn.df$AGE <- age_map[tn.df$AGE]
    
    tnn <- t(tn.df[,-c(1:3)])
    
    coldata = tn.df[,1:2]  
    coldata$AGE = as.factor(coldata$AGE)
    coldata$GENDER = as.factor(coldata$GENDER)
    
    dds = DESeqDataSetFromMatrix(countData = tnn,
                                 colData = coldata,
                                 design = ~ GENDER + AGE)
    dds <- DESeq(dds)
    res = results(dds)
    new_DF <- as.data.frame(res[!is.na(res$padj),])
    respadj =  new_DF[new_DF$padj < 0.05,]
    respadj$Tissue <- rep(TissueType[i], nrow(respadj))
    respadj$SelectedGenes <- rownames(respadj)
    
    DESeq2ModesAll <- rbindlist(list(DESeq2ModesAll, respadj))
    print(TissueType[i])
}

for(i in c(26, 29, 34, 36, 37)){ #
    edat.sel <- NULL
    Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
    pos <- which(colnames(edat) %in% rownames(Tissue)) 
    edat.Tiss <- edat[,pos]
    edat.Tiss[edat.Tiss == 0] <- NA
    #Removing genes that have less and equal to 10% zero values
    edat.Tiss <- edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
    edat.Tiss[is.na(edat.Tiss)] <- 0
    
    edat.sel <- t(as.matrix(edat.Tiss))
    tn.df <- as.data.frame(cbind(
        Tissue$AGE, Tissue$GENDER, Tissue$DTHHRDY, edat.sel))
    colnames(tn.df)[1:3] <- c("AGE", "GENDER", "DTHHRDY")
    age_map <- c(1,1,2,2,3,3)
    tn.df$AGE <- age_map[tn.df$AGE]
    
    tnn <- t(tn.df[,-c(1:3)])
    
    coldata = tn.df[,1:2]  
    coldata$AGE = as.factor(coldata$AGE)
    coldata$GENDER = as.factor(coldata$GENDER)
    
    dds = DESeqDataSetFromMatrix(countData = tnn,
                                 colData = coldata,
                                 design = ~ AGE)
    dds <- DESeq(dds)
    res = results(dds)
    new_DF <- as.data.frame(res[!is.na(res$padj),])
    respadj =  new_DF[new_DF$padj < 0.05,]
    respadj$Tissue <- rep(TissueType[i], nrow(respadj))
    respadj$SelectedGenes <- rownames(respadj)
    
    DESeq2ModesAll <- rbindlist(list(DESeq2ModesAll, respadj))
    print(TissueType[i])
}
    end.time <- Sys.time();beepr::beep(3)
    time.taken  <- (end.time - start.time) 
#Time difference of 53.54535 mins

length(unique(DESeq2ModesAll$Tissue))#37
table(DESeq2ModesAll$Tissue)
    
write.csv(DESeq2ModesAll, "DESeq2AllTissues_gender.csv")
    