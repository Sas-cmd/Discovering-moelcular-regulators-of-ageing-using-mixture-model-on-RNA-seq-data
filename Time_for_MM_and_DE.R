library(data.table)
library(Biobase)
library(tidyverse)
library(mclust)
library(RmixmodCombi)
library(MixAll)

#Functions####
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
AgeGroupFactors <- function(dat){
  dat$AGE[dat$AGE == 25] <- 1
  dat$AGE[dat$AGE == 35] <- 1
  dat$AGE[dat$AGE == 45] <- 2
  dat$AGE[dat$AGE == 55] <- 2
  dat$AGE[dat$AGE == 65] <- 3
  dat$AGE[dat$AGE == 75] <- 3
  
  return(dat)
}
runDESeq2 <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}


runEdgeR <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"], predbeta10=predbeta10[,"pData(e)$conditionB"])
}
algos <- list("DESeq2"=runDESeq2,"edgeR"=runEdgeR)
#Calculating time for 2 largest and 2 smallest GTEx tissues
####Main####
obj = readRDS("D:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

d1 = as.data.frame(table(pdat$our_subtypes))

FinishOut = data.frame()
New.list = list()
for(i in c(7)){
  Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
  pos <- which(colnames(edat) %in% rownames(Tissue)) 
  edat.Tiss <- edat[,pos]
  edat.Tiss[edat.Tiss == 0] <- NA
  #Removing genes that have less and equal to 10% zero values
  edat.Tiss <- edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
  edat.sel <- as.matrix(edat.Tiss)
  edat.sel2 <- t(as.matrix(edat.Tiss))
  tn.df <- as.data.frame(cbind(
    Tissue$AGE, Tissue$GENDER, edat.sel2))
  colnames(tn.df)[1:2] <- c("AGE", "GENDER")
  age_map <- c(1,1,2,2,3,3)
  tn.df$AGE <- age_map[tn.df$AGE]
  tnn <- t(tn.df[,-c(1:2)])
  
   #MMstat = Sys.time()
 #   FinMM <- MMRun(dat = edat.sel, modes = 1:3, Partion_name = TissueType[i])
 #   New.list[[i]] <- FinMM[[1]]
 #  
 #   output1 <- data.table()
 #   d1Gene = New.list[[i]]
 #   d1age = as.data.frame(Tissue$AGE)
 #   colnames(d1age) = "AGE"
 #   age_map <- c(25,35,45,55,65,75)
 #   d1age$AGE <- age_map[d1age$AGE]
 #   col_names <- colnames(d1Gene)
 #   for(j in 1:ncol(d1Gene)){
 #       d1Fin = cbind(rownames(d1Gene), d1age, d1Gene[,j])
 #       colnames(d1Fin) <- c("Patients", "AGE", "Gene")
 #       d1Fin = as.data.frame(d1Fin)
 #       if(length(unique(d1Fin$Gene)) == 1) next
 #       
 #       d1Fin <- AgeGroupFactors(d1Fin)
 #       ft <- as.matrix(with(d1Fin, table(AGE, Gene))) 
 #       pt <- fisher.test(ft, workspace = 2e8)
 #       pv <- pt$p.value
 #       output2 <- data.table(t(c(as.numeric(pt$p.value), 
 #                                 col_names[j], 
 #                                 as.numeric(
 #                                   length(unique(
 #                                     d1Fin$Gene))))))
 #       colnames(output2) <- c("Age_p.value", "Gene", 
 #                              "No_of_Modes")
 #       output1 <- rbindlist(list(output1, output2))
 # }
     # d1Gene = New.list[[i]]
     # d1Gender = as.data.frame(Tissue$GENDER)
     # colnames(d1age) = "GENDER"
       
   #   col_names <- colnames(d1Gene)
   #   outputgender = data.table()
   #   for(j in 1:ncol(d1Gene)){
   #       d1Fin = cbind(rownames(d1Gene), d1Gender, d1Gene[,j])
   #       colnames(d1Fin) <- c("Patients", "GENDER", "Gene")
   #       
   #       if(length(unique(d1Fin$Gene)) == 1) next
   #       
   #       ft <- as.matrix(with(d1Fin, table(GENDER, Gene))) 
   #       pt <- fisher.test(ft, workspace = 2e8)
   #       pv <- pt$p.value
   #       outputgender1 <- data.table(t(c(as.numeric(pt$p.value), 
   #                                 col_names[j], 
   #                                 as.numeric(
   #                                   length(unique(
   #                                     d1Fin$Gene))))))
   #       colnames(outputgender1) <- c("Gender_p.value", "Gene", 
   #                              "No_of_Modes")
   #       
   #       outputgender <- rbindlist(list(outputgender,
   #                                      outputgender1))
   #      
   #   }
   # 
   #   
   # Finaloutput = merge(output1, outputgender, by = "Gene",
   #                     all.x = TRUE)
   # MMend = Sys.time()
  
  #Traditional DE
  TimeDEStart = Sys.time()
  coldata = tn.df[,1:2]  
  coldata$AGE = as.factor(coldata$AGE)
  coldata$GENDER = as.factor(coldata$GENDER)
  tnn[is.na(tnn)] = 0
  dds = DESeqDataSetFromMatrix(countData = tnn,
                               colData = coldata,
                               design = ~ GENDER + AGE)
  dds <- DESeq(dds)
  res = results(dds)
  new_DF <- as.data.frame(res[!is.na(res$padj),])
  respadj =  new_DF[new_DF$padj < 0.05,]
  TimeDEend = Sys.time()
  
  EdgeRStart = Sys.time()
  age <- factor(tn.df$AGE)
  gender <- factor(tn.df$GENDER)
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
  
  EdgeRend = Sys.time()
  
  lastout = data.frame(Tissue = TissueType[i], 
              #MMTotalTime = MMend - MMstat,
              DESeq2Time = TimeDEend - TimeDEStart,
              EdgeRTime = EdgeRend - EdgeRStart)
  FinishOut = rbindlist(list(FinishOut, lastout))
}

