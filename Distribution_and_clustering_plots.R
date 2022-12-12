library(dplyr)
library(ggplot2)


t1 = read.csv("D:/Ageing Mixture models project/1. New mixture model with wrapper/Best Model_wrapper/Best Model/BestModel3_adipose_subcutaneous.csv")


obj = readRDS("D:\\gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69


TissueType<-unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))
Tissue = pdat[which(pdat$our_subtypes == TissueType[i] ), ] #extracting all the tissuetypes
pos <- which(colnames(edat) %in% rownames(Tissue)) 
edat.Tiss <- edat[,pos]

edat.Tiss1 = as.data.frame(edat.Tiss[rownames(edat.Tiss) == "ENSG00000224081.3",])
colnames(edat.Tiss1) = "expression"
edat.Tiss1$x = seq(1:nrow(edat.Tiss1))
edat.Tiss1$expression = order(edat.Tiss1$expression)
edat.Tiss1$m = rownames(edat.Tiss1)


ggplot(edat.Tiss1, aes(x=x, y= expression) ) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_minimal() + labs(fill="Density")


t1.1 = as.data.frame(t1[,c("ENSG00000224081.3", "Row.names")])
colnames(t1.1) = c('clustering', "m")
t1.1$clustering = as.factor(t1.1$clustering)


data = merge(edat.Tiss1, t1.1, by = "m")

data$expression = order(data$expression )

ggplot(data=data, aes(x=expression, group=clustering, fill=clustering)) +
    geom_density(alpha=.4) + theme_minimal() + 
    scale_fill_manual(values=c("#A7D2CB", "#F2D388", "#874C62"))

    
ggplot(data, aes(x=x, y=expression)) + 
    geom_point(size=3) +
    aes(color = clustering) +
    scale_color_manual(values=c("#A7D2CB", "#F2D388", "#874C62")) + theme_minimal()
    

##
edat.Tiss1 = as.data.frame(edat.Tiss[rownames(edat.Tiss) == "ENSG00000154743.13",])
colnames(edat.Tiss1) = "expression"
edat.Tiss1$x = seq(1:nrow(edat.Tiss1))
edat.Tiss1$expression = order(edat.Tiss1$expression)
edat.Tiss1$m = rownames(edat.Tiss1)


ggplot(edat.Tiss1, aes(x=x, y= expression) ) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_minimal() + labs(fill="Density")


t1.1 = as.data.frame(t1[,c("ENSG00000154743.13", "Row.names")])
colnames(t1.1) = c('clustering', "m")
t1.1$clustering = as.factor(t1.1$clustering)


data = merge(edat.Tiss1, t1.1, by = "m")

data$expression = order(data$expression )

ggplot(data=data, aes(x=expression, group=clustering, fill=clustering)) +
    geom_density(alpha=.4) + theme_minimal() + 
    scale_fill_manual(values=c("#F2D7D9", "#748DA6"))


ggplot(data, aes(x=x, y=expression)) + 
    geom_point(size=3) +
    aes(color = clustering) +
    scale_color_manual(values=c("#F2D7D9", "#748DA6")) + theme_minimal()
