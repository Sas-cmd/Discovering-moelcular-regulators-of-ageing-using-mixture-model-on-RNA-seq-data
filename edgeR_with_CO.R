## ---------------------------
##
## Script name: edgeR for Discovering-moelcular-regulators-of-ageing-using-mixture-model-on-RNA-seq-data 
##
## Author: Sasdekumar Loganthan
##
## Date Created: 2019-27-10
##
## Email: s.loganathan@uq.edu.au
##
## ---------------------------
## Packages required
library(edgeR)
library(Biobase)

#Data Prep and saving output----------
obj = readRDS("C:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69


TissueType<-unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))


DE3ModesAll <- NULL
for(i in 1:length(TissueType)){ 
    edat.sel <- NULL
    Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
    pos <- which(colnames(edat) %in% rownames(Tissue)) 
    edat.Tiss <- edat[,pos]
    edat.Tiss[edat.Tiss == 0] <- NA
    #Removing genes that have less and equal to 10% zero values
    edat.Tiss <- edat.Tiss[rowSums(
        is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
    
    edat.Tiss[is.na(edat.Tiss)] <- 0
    
    edat.sel <- t(as.matrix(edat.Tiss))
    tn.df <- as.data.frame(cbind(Tissue$AGE, 
                                 Tissue$GENDER, 
                                 Tissue$DTHHRDY, 
                                 edat.sel))
    
    colnames(tn.df)[1:3] <- c("AGE", "GENDER", "DTHHRDY")
    age_map <- c(1,1,2,2,3,3)
    tn.df$AGE <- age_map[tn.df$AGE]
    
    tnn <- t(tn.df[,-c(1:3)])
    
    #y <- DGEList(tnn)  
    age <- factor(tn.df$AGE)
    gender <- factor(tn.df$GENDER)
    dthh <- factor(tn.df$DTHHRDY)
    y <- DGEList(counts = tnn, group = age)
    design <- model.matrix(~gender + age)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef=2:3)
    FDR <- p.adjust(qlf$table$PValue, method = "BH")
    p <- cbind(qlf$table, FDR)
    pp <- p[p$FDR < 0.05,]
    pp$Tissue <- rep(name[i], nrow(pp))
    pp$SelectedGenes <- rownames(pp)
    DE3ModesAll <- rbind(DE3ModesAll, pp)
    print(i)
}

write.csv(DE3ModesAll, file = "DEAllTissues_gender.csv")