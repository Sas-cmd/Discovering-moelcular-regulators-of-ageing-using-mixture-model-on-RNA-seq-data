#To produce this image we had to run the mixture models without the wrapper. 

#Each package was run separately results or each of the package were stored separately.

library(data.table)
library(dplyr)
library(tidyverse)

obj = readRDS("C:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))


#Mclust BIC
temp <- list.files("D:/Ageing Mixture models project/Mclust/Mclust3/Mclust BIC", pattern="*.csv", full.names=TRUE)
BICMclustlist <- lapply(temp, fread, header = TRUE)

#MixAll BIC
temp <- list.files("D:/Ageing Mixture models project/MixAll/BIC", pattern="*.csv", full.names=TRUE)
BICMixAlllist <- lapply(temp, fread, header = TRUE)


temp <- list.files("D:/Ageing Mixture models project/RmixCombi/BIC", pattern="*.csv", full.names=TRUE)
BICRmixCombilist <- lapply(temp, fread, header = TRUE)

#Selecting lowest BIC
FinalBIC <- data.table()
for (i in 1:38) {
  
  x <- BICMclustlist[i][[1]]
  y <- BICMixAlllist[i][[1]]
  z <- BICRmixCombilist[i][[1]]
  
  colnames(y) <- c("V1", "Gene", "BIC")
  
  AllBIC <- data.table()
  for (j in 1:nrow(x)) {
    t1 <- x[x$Gene %in% x$Gene[j],]
    t2 <- y[y$Gene %in% x$Gene[j],]
    t3 <- z[z$Gene %in% x$Gene[j],]
    
    RmixMclustDiff <- t1$BIC - t3$BIC
    RmixMixDiff <- t2$BIC - t3$BIC
    AllDiif <- as.data.frame(cbind(RmixMclustDiff, RmixMixDiff, name[i] ))
    if(ncol(AllDiif) == 1) next
    AllBIC <- rbindlist(list(AllBIC, AllDiif))
  }
  FinalBIC <- rbindlist(list(FinalBIC, AllBIC))
}
write.csv(FinalBIC, "BIC_differences_between_packages.csv")

write.csv(FinalBIC, file = "BiC_comparison.csv")

FinalBIC <- read.csv("D:/Ageing Mixture models project/BiC_comparison.csv")

ggplot(FinalBIC, aes(x = V3, y = RmixMixDiff)) + geom_point()

ggplot(FinalBIC, aes(x = V3, y = RmixMclustDiff)) + 
  geom_boxplot() + 
  labs(x = "Tissue Type", y = "BIC score difference between Mclust and RMixModCombi") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90) ) +
  scale_y_continuous(breaks=seq(-500,8500,500))
summary(FinalBIC$RmixMclustDiff)
                          