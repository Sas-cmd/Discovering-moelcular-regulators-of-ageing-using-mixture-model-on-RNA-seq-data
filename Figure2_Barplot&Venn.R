## ---------------------------
##
## Script name: Discovering-moelcular-regulators-of-ageing-using-mixture-model-on-RNA-seq-data Figure 2
##
## Author: Sasdekumar Loganthan
##
## Date Created: 2021-06-11
##
## Email: s.loganathan@uq.edu.au
##
## ---------------------------
## Packages required

library(tidyverse)
library(data.table)
library(grex)
library(VennDiagram)
library(wesanderson)

## ---------------------------

## load up our functions into memory
GTExClean <- function(x){ 
    clean <- cleanid(as.character(x$SelectedGenes))
    SelectedGenes1 <- grex(clean)
    SelectedGenes3 <- cbind(SelectedGenes1, x)
    SelectedGenes2 <- SelectedGenes3[!is.na(SelectedGenes3$entrez_id),]
} # To clean Ensemble gene IDs

# Data Prep and Figure 2A ---------------------------  

#Load data
obj = readRDS("~/gtex_portal_normalized.rds") #downloaded and saved into hard-drive
edat <- exprs(obj) # Expression values extracted 30303x8527
pdat <- pData(obj) # Sample level data 8527x69

TissueType <- unique(pdat$our_subtypes)
name <- levels(as.factor(TissueType))

StdDEGenes <- fread("~/DEAllTissues_gender_All.csv")

StdDEGenes <- StdDEGenes[,-1]

StdDEGenesClean <- GTExClean(StdDEGenes)

#write.csv(StdDEGenesClean, "All_EdgeR_Genes.csv")

#Genes detected by edgeR
check2 <- as.data.frame(with(StdDEGenesClean, table(Tissue)))

#write.csv(check2, "AllEdgeRGenes.csv")
MMGenes_GenderFil <- fread("~/ExactTestAGE_afterGender.csv")

MMGenes_GenderFil <- MMGenes_GenderFil[,-1]
colnames(MMGenes_GenderFil)[3] <- "SelectedGenes"

MMGenesClean <- GTExClean(MMGenes_GenderFil)

#Genes detected by mixture models to determine freq of genes
check2 <- as.data.frame(with(MMGenesClean, table(Tissue)))
#write.csv(check2, "AllMMGenes_gender_clean.csv")


#Modes and packages of genes detected by mixture models
new2 <- MMGenes_GenderFil[MMGenes_GenderFil$Gene %in% MMGenesClean$Gene,]

check2 <- as.data.frame(with(MMGenesClean, table(Tissue, Package, No_of_Modes )))

#write.csv(check2, "AllmodesBymixture models.csv")

FinalAll <- data.frame()
for(i in 1:38){
    SGenes <- StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    EAge <- MMGenesClean %>% filter(Tissue == name[i])
    com <- base::setdiff(EAge$entrez_id, SGenes$entrez_id)
    Gene <- EAge[EAge$entrez_id %in% com,]
    FinalAll <- rbindlist(list(FinalAll, Gene))
}

#write.csv(FinalAll, "UniqueMM_GenderRemoved.csv")

check3 <- as.data.frame(with(FinalAll, table(Tissue)))

#All mixture model genes
GeneCountAllMM <- as.data.frame(with(MMGenesClean, table(Tissue)))
colnames(GeneCountAllMM)[2] <- "MM Genes"
#write.csv(GeneCountAllMM, "Genetable_MMALL_GenderRemoved.csv")

#EdgeR Genes
GeneCountEdgeR <- as.data.frame(with(StdDEGenesClean, table(Tissue)))
colnames(GeneCountEdgeR)[2] <- "MM + edgeR Genes"
#write.csv(GeneCountEdgeR, "Genetable_EdgeR_GenderRemoved.csv")

#EdgeR Unique Genes 
GeneCountEdgeRUnique  <- data.frame()
for(i in 1:38){
    SGenes <- StdDEGenesClean %>% dplyr::filter(Tissue == name[i])
    EAge <- MMGenesClean %>% filter(Tissue == name[i])
    com <- base::setdiff(SGenes$entrez_id, EAge$entrez_id)
    Gene <- SGenes[SGenes$entrez_id %in% com,]
    GeneCountEdgeRUnique <- rbindlist(list(GeneCountEdgeRUnique, Gene))
}

GeneCountEdgeRUniquetable <- as.data.frame(with(GeneCountEdgeRUnique, table(Tissue)))

colnames(GeneCountEdgeRUniquetable)[2] <- "edgeR unique Genes"
#write.csv(GeneCount, "Genetable_edgeRUnique_GenderRemoved.csv")


#Mixture model unique Genes
GeneCount <- as.data.frame(with(FinalAll, table(Tissue)))
colnames(GeneCount)[2] <- "Mixture model unique Genes"
#write.csv(GeneCount, "Genetable_UniqueMM_GenderRemoved.csv")

Fulltable.1 <- full_join(GeneCountAllMM, GeneCountEdgeR, by = "Tissue")

Fulltable.2 <- full_join(Fulltable.1, GeneCount, by = "Tissue")

Fulltable <- full_join(Fulltable.2, GeneCountEdgeRUniquetable, by = "Tissue")

Fulltable$Percent_MMU_genes <- (Fulltable$`Mixture model unique Genes`/Fulltable$`MM Genes`)*100

#Mereging with total number of donors

Number_of_Donors <- as.data.frame(table(pdat$our_subtypes))

colnames(Number_of_Donors) <- c("Tissue", "Number of Donors")

Fulltable <- full_join(Fulltable, Number_of_Donors, by = "Tissue")


Fulltable$Tissue <- gsub("_", " ", Fulltable$Tissue, fixed=TRUE)
Fulltable$Tissue <- tools::toTitleCase(Fulltable$Tissue)
colnames(Fulltable)[5] <- "edgeR Unique"

#write.csv(Fulltable, file = "Mixture model paper FIG 1.csv")

plottable <- Fulltable[,c(1,6)]

plottable[is.na(plottable)] <- 0

plottable$Tissue = with(plottable, reorder(Tissue, Percent_MMU_genes))


"-------Plotting of figure 2A----------"
ggplot(plottable, aes(x = Tissue, y = Percent_MMU_genes)) + 
    geom_bar(stat='identity', fill = wes_palette(n = 1, name = "Moonrise3"))  + 
    labs(x = "Tissue Type", y = "% of MMU Genes") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust=1, 
                                     vjust = 0.45)) 



# Venn diagram for figure 2B to 2D--------------------------

#Whole blood 
MMUWB <- MMGenesClean %>% 
    filter(Tissue == "whole_blood") %>% pull(entrez_id) %>% as.numeric

AllWB <- StdDEGenesClean %>% 
    filter(Tissue == "whole_blood") %>% pull(entrez_id) %>% as.numeric


venn.diagram(
    x = list(MMUWB, AllWB),
    category.names = c("Mixture model genes" , "edgeR genes"),
    filename = 'Whole blood MMU genes.tiff',
    output=TRUE,
    imagetype="tiff" ,
    height = 6400 , 
    width = 6400 , 
    resolution = 1200,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = wes_palette(n = 2, name = "Moonrise2"),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(5, 150),
    cat.dist = c(0.02, 0.03),
    cat.fontfamily = "sans",
    cat.col = wes_palette(n = 2, name = "Moonrise2")
)

#Skin
MMuskin <- MMGenesClean %>% 
    filter(Tissue == "skin") %>% pull(entrez_id) %>% as.numeric

Allskin <- StdDEGenesClean %>% 
    filter(Tissue == "skin") %>% pull(entrez_id) %>% as.numeric


venn.diagram(
    x = list(MMuskin, Allskin),
    category.names = c("Mixture model genes" , "edgeR genes"),
    filename = 'Skin MMU genes.tiff',
    output=TRUE,
    imagetype="tiff" ,
    height = 6400 , 
    width = 6800 , 
    resolution = 1200,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = wes_palette(n = 2, name = "Chevalier1"),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(5, 1),
    cat.dist = c(0.02, 0.03),
    cat.fontfamily = "sans",
    cat.col = wes_palette(n = 2, name = "Chevalier1")
)

#Brain-0 
MMUB0 <- MMGenesClean %>% 
    filter(Tissue == "brain-0") %>% pull(entrez_id) %>% as.numeric

AllB0 <- StdDEGenesClean %>% 
    filter(Tissue == "brain-0") %>% pull(entrez_id) %>% as.numeric


venn.diagram(
    x = list(MMUB0, AllB0),
    category.names = c("Mixture model genes" , "edgeR genes"),
    filename = 'Brain-0 MMU genes.tiff',
    output=TRUE,
    imagetype="tiff" ,
    height = 6400 , 
    width = 6400 , 
    resolution = 1200,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = wes_palette(n = 2, name = "Darjeeling2"),
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 0.8,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(5, 150),
    cat.dist = c(0.02, 0.03),
    cat.fontfamily = "sans",
    cat.col = wes_palette(n = 2, name = "Darjeeling2")
)
