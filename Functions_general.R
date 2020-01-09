#  various functions    


# this function requires log transformed cpm counts as input, using cpm function.
runPCA <- function(NumericMatrix, phenoData){
  
  Num <- t(NumericMatrix)
  
  PCA <- prcomp(Num, center=TRUE, scale.=FALSE) 
  sumPCA <- summary(PCA)
  scores <- PCA$x
  
  # attach sample info to PC scores  
  # make sure sample order matches
  pdatScores <- data.frame(phenoData, scores)
  
  #scree plot, the number of informative PCs = elbow
  #plot(PCA, type="l") 
  
  #separate out the number of PCs that explain most of the variation
  PC1 <- paste("PC1 (", round(100*sumPCA$importance[2, 1], 1),  "%)", sep="")
  PC2 <- paste("PC2 (", round(100*sumPCA$importance[2, 2], 1),  "%)", sep="")
  PC3 <- paste("PC3 (", round(100*sumPCA$importance[2, 3], 1),  "%)", sep="")
  PC4 <- paste("PC4 (", round(100*sumPCA$importance[2, 4], 1),  "%)", sep="")
  PC5 <- paste("PC5 (", round(100*sumPCA$importance[2, 5], 1),  "%)", sep="")
  PC6 <- paste("PC6 (", round(100*sumPCA$importance[2, 6], 1),  "%)", sep="")
  PC7 <- paste("PC7 (", round(100*sumPCA$importance[2, 7], 1),  "%)", sep="")
  PC8 <- paste("PC8 (", round(100*sumPCA$importance[2, 8], 1),  "%)", sep="")
  PC9 <- paste("PC7 (", round(100*sumPCA$importance[2, 9], 1),  "%)", sep="")
  PC10 <- paste("PC8 (", round(100*sumPCA$importance[2, 10], 1),  "%)", sep="")
  
  
  PCimportance <- data.frame(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
  
  PCA.all <- list(PCA=PCA, pdatScores=pdatScores, PCi=PCimportance)
}


# run a PCA on a small number of samples

runPCA.4 <- function(NumericMatrix, phenoData){
  
  Num <- t(NumericMatrix)
  
  PCA <- prcomp(Num, center=TRUE, scale.=FALSE) 
  sumPCA <- summary(PCA)
  scores <- PCA$x
  
  # attach sample info to PC scores  
  # make sure sample order matches
  pdatScores <- data.frame(phenoData, scores)
  
  #scree plot, the number of informative PCs = elbow
  #plot(PCA, type="l") 
  
  #separate out the number of PCs that explain most of the variation
  PC1 <- paste("PC1 (", round(100*sumPCA$importance[2, 1], 1),  "%)", sep="")
  PC2 <- paste("PC2 (", round(100*sumPCA$importance[2, 2], 1),  "%)", sep="")
  PC3 <- paste("PC3 (", round(100*sumPCA$importance[2, 3], 1),  "%)", sep="")
  PC4 <- paste("PC4 (", round(100*sumPCA$importance[2, 4], 1),  "%)", sep="")

  PCimportance <- data.frame(PC1, PC2, PC3, PC4)
  
  PCA.all <- list(PCA=PCA, pdatScores=pdatScores, PCi=PCimportance)
}

