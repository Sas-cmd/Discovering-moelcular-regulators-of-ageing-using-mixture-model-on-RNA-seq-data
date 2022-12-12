library(Biobase)
library(data.table)
library(tidyverse)

#Functions####
fread_plus <- function(dat.read) {
  fread(dat.read, fill = TRUE, drop = 1) %>% 
    mutate(filename = dat.read)
}
list_fread <- function(path, pattern = "*.csv") {
  files = list.files(path, pattern, full.names = TRUE)
  lapply(files, function(x) fread_plus(x))
}
rbindlist_fread <- function(path, pattern = "*.csv") {
  files = list.files(path, pattern, full.names = TRUE)
  rbindlist(lapply(files, function(x) fread_plus(x)))
}
cleanNameTissues = function(tissues){
  tt = sub(".*/", "", tissues)
  tt2 = sub("\\..*", "", tt)
  str_remove_all(tt2, "BestModel3_")
}


#####
obj = readRDS("D:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

temp.list = list.files(path = "D:/Ageing Mixture models project/1. New mixture model with wrapper/Best Model_wrapper/Best Model",
                      pattern = ".csv")

read.tissues = cleanNameTissues(temp.list)

t2 = list_fread("D:/Ageing Mixture models project/1. New mixture model with wrapper/Best Model_wrapper/Best Model")

names(t2) = read.tissues
Final.dat = data.frame()
for(i in 1:length(names(t2))){
  Tissue = pdat[which(pdat$our_subtypes == names(t2[i])),]
  pos = which(colnames(edat) %in% rownames(Tissue))
  edat.Tiss = edat[,pos]
  edat.Tiss[edat.Tiss == 0] = NA
  edat.Tiss = edat.Tiss[rowSums(is.na(edat.Tiss)) <= floor(length(edat.Tiss[1,])*0.1), ]
  edat.Tiss[is.na(edat.Tiss)] = 0
  edat.sel <- as.matrix(edat.Tiss)
  
  clus.main = as.data.frame(
    t2[names(t2) == names(t2)[i]][[1]]
  )
  print(names(t2)[i])
  for(j in 3:ncol(clus.main)){
    tt1 = clus.main[,j]
    if(length(unique(tt1)) == 1) next
    x = edat.sel[rownames(edat.sel) == colnames(clus.main)[j],] 
    ttt2 = as.data.frame(cbind(tt1, x))
    colnames(ttt2) = c("Clus", "Exp")
    ttt2$Clus = as.factor(ttt2$Clus)
    aov.res = aov(Exp ~ Clus, data = ttt2)
    x1 = t(as.data.frame(tapply(ttt2$Exp, ttt2$Clus, mean)))
    if(ncol(x1) > 2){
      logFCA_B = log2(x1[,1]) - log2(x1[,2])
      logFCA_C = log2(x1[,1]) - log2(x1[,3])
    }else{
    logFCA_B = log2(x1[,1]) - log2(x1[,2])
    logFCA_C = "NA"
    }
    Cluster.no = as.numeric(length(unique(ttt2$Clus)))
    fin.dat = data.frame(Gene = colnames(clus.main)[j], 
                    TissueType = names(t2)[i],
                    Cluster.no = Cluster.no,
                    logFCA_B = logFCA_B,
                    logFCA_C = logFCA_C,
                    Anova.co_intercept = aov.res$coefficients[1],
                    Anova.co_Clus2 = aov.res$coefficients[2],
                    Anova.co_Clus3 = aov.res$coefficients[3],
                    F.value = summary(aov.res)[[1]][1,4],
                    p.value = summary(aov.res)[[1]][1,5])
    Final.dat = rbindlist(list(Final.dat, fin.dat))
    print(paste0(names(t2)[i], "_", j))
  }
}

write.csv(Final.dat, file = "Fold_change.csv")


Fold.change = fread("Fold_change.csv")
Fold.change$V1 = NULL
colnames(Fold.change)[c(2,10)] = c("Tissue", "p.val.Anova")


Maindata = fread("S2_UniqueMM_Gemes_dectected.csv")
Maindata$V1 = NULL

Fina.data = data.table()
for(i in 1:length(TissueType)){
  x1 = Maindata %>% filter(Tissue == TissueType[i])
  x2 = Fold.change %>% filter(Tissue == TissueType[i])
  Fin.x = merge(x1, x2, by = "Gene", all.x = TRUE)
  Fina.data = rbindlist(list(Fina.data, Fin.x))
}

write.csv(Fina.data, 
          "S2_UniqueMM_Gemes_dectected_Fold_change.csv")

