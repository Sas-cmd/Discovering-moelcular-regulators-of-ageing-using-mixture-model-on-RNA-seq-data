## ---------------------------
##
## Script name: Discovering-moelcular-regulators-of-ageing-using-mixture-model-on-RNA-seq-data Figure 3
##
## Author: Sasdekumar Loganthan
##
## Date Created: 2021-06-15
##
## Email: s.loganathan@uq.edu.au
##
## ---------------------------
## Packages required
library(data.table)
library(grex)
library(ggplot2)
library(naniar)
library(org.Hs.eg.db)
library(wesanderson)
library(readr)
library(org.Hs.eg.db)
library(tidyverse)
library(SuperExactTest)


#Data Prep and Figure 3A and 3B ---------------------------  
    
#Comparing to external ageing databases
#GeneAge database
GeneAgeDB = read_csv("~/genage_human.csv")
#Digitial ageing atals
DAA <- read_csv("~/digital_ageing_atlas_data1.csv") 

DAAHumans = DAA %>% filter(Species == "Homo sapiens") #3574
DAAHumans[DAAHumans == "#VALUE!"] <- NA
DAAHumans = DAAHumans[!is.na(DAAHumans$GeneId),]

#Converting Gene symbol to enterzid for DAA humans#

DAAGenes = as.data.frame(mapIds(org.Hs.eg.db, DAAHumans$GeneId, 'ENTREZID', 'SYMBOL') )
colnames(DAAGenes)[1] <- "DAA"

#Mixture model unique only genes 
FinalAll = read_csv("~/UniqueMM_GenderRemoved.csv")
#edgeR genes 
StdDEGenesClean <- read.csv("~/All_EdgeR_Genes.csv")

x = list(FinalAll$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes$DAA)
names(x) = c("MMU", "GenAge", "DAA")

total = sum((length.gene.sets=sapply(x,length)))

res=supertest(x, n=total)

"Figure 3A plotting"
cols = wes_palette("Moonrise2", 10, type = "continuous")
plot(res, Layout = "landscape", 
     color.scale.title = expression(paste(-Log[10],
                                          '(',italic(P),')')),
     color.on ='#a1896e', color.off ='#fefcfc', 
     ylim = c(0, 15000),
     sort.by = "size", margin = c(0.5, 5, 1, 3), 
     heatmapColor = cols)

MM.p.val = res$P.value

dev.print(pdf, 'MMU_CommonGenes_with_databases.pdf')

"Figure 3B plotting"

x = list(StdDEGenesClean$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes$DAA)
names(x) = c("edgeR", "GenAge", "DAA")
total = sum((length.gene.sets=sapply(x,length)))

res = supertest(x, n=total)

cols = wes_palette("Moonrise2", 10, type = "continuous")
plot(res, Layout = "landscape",  
     color.scale.title = expression(paste(-Log[10],
                                          '(',italic(P),')')),
     color.on ='#a1896e', color.off ='#fefcfc', 
     ylim = c(0, 15000),
     sort.by="size", margin=c(0.5, 5, 1, 3), 
     heatmapColor = cols)


Std.p.val = res$P.value

dev.print(pdf, 'EdgeRGene+MMU_genes_common_with_databases.pdf')


#Data Prep and Figure 3C -------------------
x = list(StdDEGenesClean$entrez_id, 
          GeneAgeDB$`entrez gene id`, DAAGenes$DAA)
Databaseintersect <- Reduce(intersect, x)
Egenes <- StdDEGenesClean[StdDEGenesClean$entrez_id %in% Databaseintersect,]
Egenes <- as.data.frame(table(Egenes$hgnc_symbol))
colnames(Egenes) <- c("Genes", "Freq")

x = list(FinalAll$entrez_id, GeneAgeDB$`entrez gene id`, DAAGenes$DAA)
Databaseintersect <- Reduce(intersect, x)
Mgenes <- FinalAll[FinalAll$entrez_id %in% Databaseintersect,]
Mgenes <- as.data.frame(table(Mgenes$hgnc_symbol))
colnames(Mgenes) <- c("Genes", "Freq")

FinGenes <- merge(Mgenes, Egenes, by = "Genes", all = TRUE)



colnames(FinGenes) = c("Genes", "MMU", "edgeR")

data <- reshape2::melt(FinGenes)
data[is.na(data)] <- 0

colnames(data) = c("Genes", "Method", "Number_of_Tissues")


data$Genes = with(data, reorder(Genes, Number_of_Tissues))


ggplot(data) + 
    geom_bar(aes(x = Genes, y = Number_of_Tissues, 
                 fill = Method), stat = "identity", 
             alpha = 0.5) + 
    scale_fill_manual(values = wes_palette(n = 2, 
                                name = "Darjeeling2")) +
    labs(x = "Genes", y = "Number of tissues") +
    ylim(0,max(data$Number_of_Tissues)) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) 
