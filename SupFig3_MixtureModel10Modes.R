#To plot this we ran mclust with 10 modes and saved the data separately

MclustMM <- function(dat = dat, no_of_mode = modes ){
  
  #Mclust
  BICMclust <- NULL
  dat[is.na(dat)] <- 0
  try({mclustMixed <- Mclust(dat, G = no_of_mode)
  BICMclust <- as.numeric((-2*(mclustMixed$loglik) + mclustMixed$df *log(mclustMixed$n)))})
  
  BestOutput <- mclustMixed$classification
  BestModel <- mclustMixed$modelName

  Mainoutput <- list(
    "Classification" = mclustMixed$classification, 
    "BIC" = BICMclust, 
    "Model" = mclustMixed$modelName)
  return(Mainoutput) 
}
MMRunMclust <- function( dat = dat, modes = x, Partion_name = Partion_name){
  AllBIC <- data.table()
  Allpartion <- NULL
  for (j in 1:nrow(dat)) {
    mm <- MclustMM(dat = dat[j,], no_of_mode = modes)
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


library(data.table)
library(ggplot2)
library(Biobase)
#Load data
obj = readRDS("C:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

New.list <- list()
for(i in length(TissueType)){
  Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
  pos <- which(colnames(edat) %in% rownames(Tissue)) 
  edat.Tiss <- edat[,pos]
  edat.Tiss[edat.Tiss == 0] <- NA
  #Removing genes that have less and equal to 10% zero values
  edat.Tiss <- edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
  edat.sel <- as.matrix(edat.Tiss)
  
  FinMM <- MMRunMclust(dat = edat.sel, modes = 1:10, Partion_name = TissueType[i])
  New.list[[i]] <- FinMM #if dataset is small, or else save clustering and BIC as csv 
}


df.listf <- list.files("D:/Ageing Mixture models project/1. New mixture model with wrapper/Mclust/Mclust10", pattern="*.csv", full.names=TRUE) 

Mclustlist <- lapply(df.listf, fread, header = TRUE)


toot <- lapply(Mclustlist, apply, 2, function(m){length(unique(m))})


for(i in 1:length(TissueType)){
  if (i == 1){
    df1 <- as.data.frame(toot[i][[1]])
    df1$Gene <- rownames(df1)
    rownames(df1) <- NULL
    df2 <- df1[-c(1:3),]
    df3 <- as.data.frame(table(df2[,1]))
    colnames(df3) <- c("Modes", name[i])
    ModeDis <- df3
    } else {
      df1 <- as.data.frame(toot[i][[1]])
      df1$Gene <- rownames(df1)
      rownames(df1) <- NULL
      df2 <- df1[-c(1:3),]
      df3 <- as.data.frame(table(df2[,1]))
      colnames(df3) <- c("Modes", name[i])
      ModeDis <- merge(ModeDis, df3, by = "Modes", all = TRUE)}
}

dat.m <- reshape2::melt(ModeDis)

ggplot(dat.m, aes(x = variable, y = value, fill= Modes)) + geom_bar(stat='identity')  + theme(axis.text.x = element_text(angle = 90)) + coord_flip()


ModeDis1 <- ModeDis
rownames(ModeDis1) <- ModeDis1$Modes
ModeDis1 <- ModeDis1[,-1]

ModePercThirds <- NULL
for(i in 1:ncol(ModeDis1)){
  One <- ModeDis1[1,i]/sum(ModeDis1[,i], na.rm = T)*100
  Two <- sum(ModeDis1[3,i])/sum(ModeDis1[,i], na.rm = T)*100
  Three <- sum(ModeDis1[4,i])/sum(ModeDis1[,i], na.rm = T)*100
  Four_to_Six <- sum(ModeDis1[5:7,i])/sum(ModeDis1[,i], na.rm = T)*100
  Seven_to_Ten <- sum(ModeDis1[c(2,8:10),i])/sum(ModeDis1[,i], na.rm = T)*100
  Apoe <- rbind(One, Two, Three, Four_to_Six, Seven_to_Ten)
  colnames(Apoe) <- name[i]  
  ModePercThirds<- cbind(ModePercThirds, Apoe)
}

rnames <- rownames(ModePercThirds) 
rnames <- gsub("_", " ", rnames, fixed=TRUE)
rownames(ModePercThirds) <- rnames


cnames <- colnames(ModePercThirds) 
cnames <- gsub("_", " ", cnames, fixed=TRUE)
cnames<- tools::toTitleCase(cnames)
colnames(ModePercThirds) <- cnames


dat.p <- reshape2::melt(ModePercThirds)

library(wesanderson)

ggplot(dat.p, aes(x = Var2, y = value, fill= Var1)) + 
  geom_bar(stat='identity')  + 
  scale_fill_manual(values=wes_palette(n = 5, name= "Moonrise3")) +
  labs(x = "Tissue Type", y = "Percentage of modes", fill = "Number of modes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) 

