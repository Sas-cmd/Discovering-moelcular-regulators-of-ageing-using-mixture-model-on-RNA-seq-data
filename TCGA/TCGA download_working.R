BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

setwd("C:/New")
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)


GDCdownload(query)
expdat <- GDCprepare(query = query,
                     save = TRUE, 
                     save.filename = "TCGA_BRCA.rda")


table(expdat$definition)


library(Biobase)
library(SummarizedExperiment)

t1 <- expdat@assays@data
t1 <- assays(expdat)$normalized_count

dim(t1)
