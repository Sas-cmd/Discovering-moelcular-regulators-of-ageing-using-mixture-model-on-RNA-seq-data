library("data.table")
library("Biobase")
library("DESeq2")
library("edgeR")
library("mclust")
library("MixAll")
library("RmixmodCombi")

load("DESeq2paper/data/meanDispPairs.RData")


makeSim <- function(n, m, x, beta, meanDispPairs, sf=rep(1,m)) {
  idx <- sample(nrow(meanDispPairs), n, replace=TRUE)
  mu0 <- meanDispPairs[idx,1]
  disp <- meanDispPairs[idx,2]
  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)
  list(mat = matrix(rnbinom(n*m, mu=muMat, size=1/disp), ncol=m),
       disp = disp,
       mu0 = mu0)
}

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
  try({Rmix <- mixmodCombi(data = dat , 
                           nbCluster = no_of_mode, 
                           mixmodOutput = NULL,
                           criterion = c("BIC", "ICL"), 
                           models = mixmodGaussianModel())
  BICRmix <- as.numeric(
    -2*(Rmix@mixmodOutput@bestResult@likelihood) + Rmix@mixmodOutput@bestResult@parameters@nbFreeParam *log(Rmix@mixmodOutput@nbSample))})
  
  if(is.null(BICMclust)){
    BICMclust = 1*10^20}
  if(is.null(BICRmix)){
    BICRmix = 1*10^20} 
  if(is.null(BICMixAll)){
    BICMixAll = 1*10^20}
  
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
    mm <- BestMM(dat = dat[j,], no_of_mode = modes )
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
namesAlgos <- names(algos)

n <- 2000
effSizeLevels <- log2(c(2,3,4))
mLevels <- c(6,8,10)
nreps <- 5

effSizes <- rep(rep(effSizeLevels, each=nreps), times=length(mLevels))
ms <- rep(mLevels, each=nreps * length(effSizeLevels))

library("BiocParallel")
register(SnowParam(workers = 3))

resList <- lapply(seq_along(ms), function(i) {
  set.seed(i)
  m = ms[i]
  es = effSizes[i]
  condition = factor(rep(c("A","B"), each = m/2))
  x = model.matrix(~ condition)
  beta = c(rep(0, n * 8/10), sample(c(-es,es), n * 2/10, TRUE))
  mat = makeSim(n,m,x,beta,meanDispPairs)$mat
  rownames(mat) = paste0("gene", seq(1:nrow(mat)))
  mat[mat == 0] <- NA
  mat = mat[rowSums(is.na(mat)) <= floor(length(mat[1,])*0.1),]
  mat[is.na(mat)] = 0
  e = ExpressionSet(mat, 
                    AnnotatedDataFrame(data.frame(condition)))
  TimeDEStart = Sys.time()
  resTest = lapply(algos, function(f) f(e))
  TimeDEend = Sys.time()
  restest2 = as.matrix(cbind(rownames(mat), 
                             resTest$DESeq2$padj,
                             resTest$edgeR$padj,
                             TimeDEend - TimeDEStart))
  
  colnames(restest2) = c("Gene", "DESeq2", "edgeR", 
                         "Time(sec)")
  
  #For mixture models
  MMstart = Sys.time()
  
  FinMM = MMRun(dat = mat, modes = c(1:2), 
                Partion_name = paste0("sim_",i))
  
  BestBICDF = FinMM[[2]]
  colnames(BestBICDF) = c("Gene", "BIC", "Models")
  
  BestBICDF$Package = lapply(BestBICDF$Models, function(tx){
    if (grepl("Gaussian", tx)) {
      Package = "RmixModCombi"
    } else if(grepl("gamma", tx)) {
      Package = "MixAllGamma"
    }else {Package = "Mclust"}
  })
  
  Clast = FinMM[[1]]
  
  col_names <- colnames(Clast)
  Finaloutput <- data.table()
  for(z in 1:ncol(Clast)){
    d1 = as.data.frame(cbind(Clast[,z], condition))
    colnames(d1) = c("Cluster", "condition")
    
    if(length(unique(d1$Cluster)) == 1) next
      
      ft <- as.matrix(with(d1, table(condition, Cluster))) 
      pt <- fisher.test(ft, workspace = 2e8)
      
      output1 <- data.table(t(c(as.numeric(pt$p.value), 
                                col_names[z], 
                                as.numeric(
                                  length(unique(
                                    d1$Cluster))))))
      colnames(output1) <- c("p.value", "Gene", 
                             "No_of_Modes")
      
      Finaloutput <- rbindlist(list(Finaloutput, output1))
      print(z)
  } 
  MMend = Sys.time()
  FinMMoutput = merge(BestBICDF, Finaloutput, by = "Gene",
                      all.x = TRUE)
  FinMMoutput$Time = MMend - MMstart
  
  Allout = merge(restest2, FinMMoutput, by = "Gene",
                 all.x = TRUE)
  Allout$effSize = es
  Allout$n = m
  
 return(Allout)

    
  
})
res <- do.call(rbind, resList)

save(res, file="results_simulateDE.rds")

