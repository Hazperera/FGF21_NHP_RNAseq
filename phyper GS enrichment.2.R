



##### write a function to filter gene set list based on which genes pass expression filter
  # use list of genes that come out in Module membership table, because those have already passed the filter in WGCNA


filterGeneSets.phyp <- function(MMtable = MMtable, geneSets = geneSets, cutoff=10){
  
  # count how many genes from each gene set exist in MMtable gene list.
  geneSets <- apply(geneSets, 2, as.character)
  NumGenesPer <- apply(geneSets, 2, function(x) sum(x %in% MMtable$HGNCsym))
  keepSets <- NumGenesPer >= cutoff
  setsRemoved <- paste0(paste0(paste0(sum(!keepSets), " sets removed, "), sum(keepSets)), " remaining.")
  geneSets.filt <- as.data.frame(geneSets[,keepSets])
  list(setsRemoved = setsRemoved, geneSets.filtered = geneSets.filt)

}

# # try it out
# geneSets.filt <- filterGeneSets.phyp(MMtable, geneSets)
# 
# geneSets.test <- geneSets.filt$geneSets.filtered



filterGeneSets.ensg <- function(counts, geneSets = geneSets, cutoff=10){
  # counts = count matrix containing ensemble IDs as column "GeneID", that has been through a low-expression filter.
  
  # count how many genes from each gene set exist in low-expression filtered counts matrix.
  NumGenesPer <- sapply(geneSets, function(x) sum(x %in% counts$GeneID))
  keepSets <- NumGenesPer >= cutoff
  setsRemoved <- paste0(paste0(paste0(sum(!keepSets), " sets removed, "), sum(keepSets)), " remaining.")
  geneSets.filt <- geneSets[keepSets]
  list(setsRemoved = setsRemoved, geneSets.filtered = geneSets.filt)
}


##### write a function that tests for gene set enrichment with hypergeometric test

  # Have a table of genes and which module they belong to,
  # and a table of gene sets (columns) 
  # that made the obesity and muscle related gene list for WGNCA. 

annotateModules.table <- function(MMtable=geneInfo, geneSets = geneSets, nTopMod=10){

  ## make a list of modules + genes in each module
  MMtable <- select(MMtable, one_of("moduleColor", "HGNCsymbol"))
  names(MMtable) <- c("moduleColor", "HGNCsym")
  MMtable$HGNCsym <- as.character(MMtable$HGNCsym)

  # write a function to put inside lapply to grab gene names from each module
  mod.names <- levels(as.factor(MMtable$moduleColor))

  getGenes <- function(modName){
    MMtable[MMtable$moduleColor == modName, "HGNCsym"]
  }

  mod.list <- lapply(mod.names, getGenes)
  names(mod.list) <- mod.names



  ## run a loop over the columns of the module table to get a table of gene sets and respective p values for each module

  ## function for one module and one gene set, that will go into apply function over columns of geneSet
  geneSets <- apply(geneSets, 2, as.character)

  # vector of genes that will serve as full list of possible genes  ** Consider what this list should contain... input or all genes?
  fullGeneList <- MMtable$HGNCsym
  fullGeneList <- unique(na.omit(fullGeneList))
 
  oneGeneSet <- geneSets[,1]
  #curModule <- mod.list[[1]]  # this is just for testing while building function...
    
  phyp.mods <- function(oneGeneSet) {
      subset <- curModule
      geneSet <- na.omit(as.character(oneGeneSet))
      fullList <- as.character(fullGeneList)  
      numBlack <- sum(subset %in% geneSet)
      totBlack <- length(geneSet)
      totWhite <- length(fullList)-totBlack
      handful <- length(subset)
      phyper(numBlack, totBlack, totWhite, handful, lower.tail = F) 
  }

  ## testing... testing... apply to whole geneSet table for one module
  # forOneModule <- apply(GeneSets, 2, phyp.mods)

  # build a table to populate with a loop that will generate vector (apply over GeneSets columns) for each module
  phypResults <- matrix(nrow=ncol(geneSets), ncol=length(mod.list))
  colnames(phypResults) <- names(mod.list)
  rownames(phypResults) <- colnames(geneSets)

  # make ze loop  !! the unlist is important !!
  for (i in 1:ncol(phypResults)){
      curModule <- unlist(mod.list[i])
      phypResults[,i] <- apply(geneSets, 2, phyp.mods)
  }

  phypResults2 <- apply(phypResults, c(1,2), function(x) signif(x,2))

  ## Now make a list with top ranked gene sets for each module
  # make an empty matrix to populate
  topRanked <- list()
  for (i in 1:ncol(phypResults2)){
    topRankOneMod <- phypResults2[order(phypResults2[,i]), i, drop=FALSE][1:nTopMod, ,drop=FALSE]
    topRankOneMod <- data.frame(geneSet=rownames(topRankOneMod), pval = topRankOneMod[,1])
    rownames(topRankOneMod) <- NULL
    topRanked[[i]] <- topRankOneMod
    names(topRanked)[[i]] <- colnames(phypResults2)[i]
  }
  
  #return from function
  list(scoreTable = phypResults2, topRanked = topRanked)
}






#### annotate modules using gene sets in list form #####
annotateModules <- function(MMtable=geneInfo, geneSets, nTopMod=10){
  
  ## make a list of modules + genes in each module using ensgID
  MMtable <- select(MMtable, one_of("moduleColor", "GeneID"))
  names(MMtable) <- c("moduleColor", "GeneID")
  MMtable$GeneID <- as.character(MMtable$GeneID)
  
  # write a function to put inside lapply to grab gene names from each module
  mod.names <- levels(as.factor(MMtable$moduleColor))
  
  getGenes <- function(modName){
    MMtable[MMtable$moduleColor == modName, "GeneID"]
  }
  
  mod.list <- lapply(mod.names, getGenes)
  names(mod.list) <- mod.names
  
  
  
  ## run a loop over the columns of the module table to get a table of gene sets and respective p values for each module
  
  ## function for one module and one gene set, that will go into apply function over columns of geneSet
  geneSets <- lapply(geneSets, as.character)
  
  # vector of genes that will serve as full list of possible genes  ** Consider what this list should contain... input or all genes?
  fullGeneList <- MMtable$GeneID
  fullGeneList <- unique(na.omit(fullGeneList))
  
  oneGeneSet <- geneSets[[1]]
  #curModule <- mod.list[[1]]  # this is just for testing while building function...
  
  phyp.mods <- function(oneGeneSet) {
    subset <- curModule
    geneSet <- na.omit(oneGeneSet)
    fullList <- as.character(fullGeneList)  
    numBlack <- sum(subset %in% geneSet)
    totBlack <- length(geneSet)
    totWhite <- length(fullList)-totBlack
    handful <- length(subset)
    phyper(numBlack, totBlack, totWhite, handful, lower.tail = F) 
  }
  
  ## testing... testing... apply to whole geneSet list for one module
  # forOneModule <- lapply(geneSets, phyp.mods)
  
  # build a table to populate with a loop that will generate vector (apply over GeneSets columns) for each module
  phypResults <- matrix(nrow=length(geneSets), ncol=length(mod.list))
  colnames(phypResults) <- names(mod.list)
  rownames(phypResults) <- names(geneSets)
  
  # make ze loop  !! the unlist is important !!
  for (i in 1:ncol(phypResults)){
    curModule <- mod.list[[i]]
    phypResults[,i] <- sapply(geneSets, phyp.mods)
  }
  
  phypResults2 <- apply(phypResults, c(1,2), function(x) signif(x,2))
  
  # calc adjusted p values
  phypResultsFDR <- apply(phypResults2, 2, function(x) p.adjust(x, method="BH", n=length(geneSets)))
  
  ## Now make a list with top ranked gene sets for each module
  # make an empty matrix to populate
  topRanked <- list()
  for (i in 1:ncol(phypResults2)){
    topRankOneMod <- phypResults2[order(phypResults2[,i]), i, drop=FALSE][1:nTopMod, ,drop=FALSE]
    topRankOneMod <- data.frame(geneSet=rownames(topRankOneMod), pval = topRankOneMod[,1])
    rownames(topRankOneMod) <- NULL
    topRanked[[i]] <- topRankOneMod
    names(topRanked)[[i]] <- colnames(phypResults2)[i]
  }

  
  #return from function
  list(pvalTable = phypResults2, FDRtable = phypResultsFDR, topRanked = topRanked)
}



#### Make a cleaned up annotation table for gene modules ####


# # first annotate modules with a gene set collection
# modAnno.OBcan <- annotateModules(AT.MM, geneSets=OBcan.ATexpr, nTopMod=20)


# then feed the result into the clean up function:
trimModuleAnnotation <- function(modAnno, numberNames, FDRcutoff, numTop ){
  
  # make table for FDR and filter down to cutoff
  modAnno.FDR <- as.data.frame(modAnno$FDRtable)
  modAnno.FDR <- modAnno.FDR %>% mutate(GeneSet = rownames(modAnno.FDR))
  
  numberOfMods <- ncol(modAnno.FDR)-1
  
  modAnno.FDR <- gather(modAnno.FDR, "Module", "FDR", 1:numberOfMods) %>%
                  filter(FDR < FDRcutoff) %>%
                  arrange(Module, FDR) %>%
                  mutate(GSmod = paste(GeneSet, Module, sep="_"))
  
  # make table for p value and merge with FDR table
  modAnno.pval <- as.data.frame(modAnno$pvalTable)
  modAnno.pval <- modAnno.pval %>% mutate(GeneSet = rownames(modAnno.pval))
  
  modAnno.pval <- gather(modAnno.pval, "Module", "Pval", 1:numberOfMods) %>%
                  mutate(GSmod = paste(GeneSet, Module, sep="_"))
  
  modAnno.trim <- modAnno.FDR %>% left_join(modAnno.pval[,3:4], by="GSmod") 
    
  table(modAnno.trim$Module)
  length(unique(modAnno.trim$Module))

  # put on module number names
  modAnno.trim <- merge(numberNames[,c("moduleColor", "modNumber")], modAnno.trim, by.x="moduleColor", by.y="Module", all.x=TRUE)
  # could get rid of "grey" non-module, but could also leave in to see if annotated...
  modAnno.trim <- modAnno.trim %>% rename(Module = modNumber) %>%
                  select(Module:FDR, Pval)

  # replace NAs with "no significant annotations"
  modAnno.trim$GeneSet[is.na(modAnno.trim$GeneSet)] <- "No significant overlap"
  noOverlap <- modAnno.trim %>% filter(is.na(FDR))

  # trim further to top # per module, ranked by p value, since there can be quite a few ties in FDR
  modAnno.trim <- modAnno.trim %>% group_by(Module) %>%
                  top_n(n=-numTop, wt=Pval) # n is negative because we want the lowest Pvals for each Module
  # table(modAnno.trim$Module)

  # stick no overlaps back on and arrange by module number
  modAnno.trim <- modAnno.trim %>%
                      bind_rows(noOverlap) %>%
                      arrange(as.numeric(str_split_fixed(Module, "_", 2)[,2]), Pval)
  
  # table(modAnno.trim$Module)
  
  modAnno.trim
}









##### hypergeometric test for gene set enrichment from top DE genes ####

# ATTENTION: rank list by p value and fold change - try gene set enrichment with up and down genes, up genes only, down genes only.

# rankedList = AT.PostVsPre
# geneSets = HallmarkGS[-1,]
# setSizes = c(10,25,50,100,250,500)

GSenrich.TopDEG <- function(rankedList, fullGeneList, geneSets = geneSets, setSizes = c(10,25,50,100,250,500), nTopGS=10){
  
  ## run a loop for each set size from top DE gene list to get table of gene sets and respective p values for each set size

  ## function for one module and one gene set, that will go into apply function over columns of geneSet
  geneSets <- lapply(geneSets, as.character)

  # vector of genes that will serve as full list of possible genes  
  fullGeneList <- unique(na.omit(fullGeneList))

  oneGeneSet <- geneSets[[1]]
  # setSize=setSizes[1]  # this is just for testing while building function...

  phyp.sets <- function(oneGeneSet) {
    subset <- rankedList$GeneID[1:setSizes[i]]
    geneSet <- na.omit(as.character(oneGeneSet))
    fullList <- as.character(fullGeneList)  
    numBlack <- sum(subset %in% geneSet)
    totBlack <- length(geneSet)
    totWhite <- length(fullList)-totBlack
    handful <- length(subset)
    phyper(numBlack, totBlack, totWhite, handful, lower.tail = F) 
  }

  ## testing... testing... apply to whole geneSet table for one module
  # forOneSetSize <- lapply(geneSets, phyp.sets)

  # build a table to populate with a loop that will generate vector of p values (apply over GeneSets columns) for each set size
  phypResults <- matrix(nrow=length(geneSets), ncol=length(setSizes))
  colnames(phypResults) <- setSizes
  rownames(phypResults) <- names(geneSets)

  # make ze loop  !! the unlist is important !!
  for (i in 1:ncol(phypResults)){
    phypResults[,i] <- unlist(lapply(geneSets, phyp.sets))
  }

  phypResults2 <- apply(phypResults, c(1,2), function(x) signif(x,2))
  
  # make a table with adjusted p values - this can be used to make set size vs -log10(pval) line graph
  adj.phypResults <- apply(phypResults, 2, function(x) p.adjust(x, method="BH", n=nrow(phypResults)))
  adj.phypResults <- apply(adj.phypResults, c(1,2), function(x) signif(x,2))
  
  

  ## Now make a list with top ranked gene sets for each set size, using unadjusted p value b/c useful for ranking
  # make an empty matrix to populate
  topRanked <- list()
  for (i in 1:ncol(phypResults2)){
    topRankOneSet <- phypResults2[order(phypResults2[,i]), i, drop=FALSE][1:nTopGS, ,drop=FALSE]
    topRankOneSet <- data.frame(geneSet=rownames(topRankOneSet), pval = topRankOneSet[,1])
    rownames(topRankOneSet) <- NULL
    topRankOneSet$setSize <- setSizes[i]
    topRanked[[i]] <- topRankOneSet
    names(topRanked)[[i]] <- colnames(phypResults2)[i]
  }
  
  #turn into data frame for easy writing
  topRanked.df <- bind_rows(topRanked)
  topRanked.df <- topRanked.df[,c(3,1,2)]
  
  # grab adjusted pvalues and add to table
  topRanked.adj <- list()
  for (i in 1:ncol(adj.phypResults)){
    topRankOneSet <- adj.phypResults[order(adj.phypResults[,i]), i, drop=FALSE][1:nTopGS, ,drop=FALSE]
    topRankOneSet <- data.frame(geneSet=rownames(topRankOneSet), adjPval = topRankOneSet[,1])
    rownames(topRankOneSet) <- NULL
    topRankOneSet$setSize <- setSizes[i]
    topRanked.adj[[i]] <- topRankOneSet
    names(topRanked.adj)[[i]] <- colnames(adj.phypResults)[i]
  }
  #turn into data frame for easy writing
  topRanked.adj.df <- bind_rows(topRanked.adj)
  
  # stick together
  topRanked.df$unique <- paste(topRanked.df$setSize, topRanked.df$geneSet, sep="_")
  topRanked.adj.df$unique <- paste(topRanked.adj.df$setSize, topRanked.adj.df$geneSet, sep="_")
  topRanked.df2 <- merge(topRanked.df, topRanked.adj.df[,c("unique", "adjPval")], by="unique")
  topRanked.df2 <- topRanked.df2 %>% arrange(setSize, pval) %>% select(setSize, geneSet, pval, adjPval)

#return from function
list(phypResults = phypResults2, adj.phypResults=adj.phypResults, topRanked = topRanked.df2)
}


# test!
# AT.PostVsPre.GSE <- GSenrich.TopDEG(AT.PostVsPre, geneSets = HallmarkGS, fullGeneList= AT.PostVsPre$GeneID, setSizes = c(10,25,50,100,250,500), nTopGS=10)


##### function for plotting ####

# phypResults <- AT.PostVsPre.GSE  # for testing

plotPhyper.adj <- function(phypResults, cutoff=.05, title="Enrichment of gene sets within top ranked DE genes"){
  
  adjResults <- as.data.frame(phypResults$adj.phypResults)
  results <- as.data.frame(phypResults$phypResults)
  adjResults$geneSets <- rownames(adjResults)
  L.adjResults <- gather(adjResults, "setSize", "adjPval", 1:(ncol(adjResults)-1))
  # add column for color 
  L.adjResults$sig <- ifelse(L.adjResults$adjPval < cutoff, "sig", "ns")
  sigGS <- filter(L.adjResults, sig=="sig")$geneSets
  L.adjResults$sigColor <- NA
  L.adjResults[L.adjResults$geneSets %in% sigGS, "sigColor"] <- L.adjResults[L.adjResults$geneSets %in% sigGS, "geneSets"]
  
  L.adjResults$setSize <- as.numeric(L.adjResults$setSize)
  
  ggplot(L.adjResults, aes(x=setSize, y=-log10(adjPval), group=geneSets, color=sigColor))+
    geom_line(size=1)+
    scale_color_discrete(name="Gene Set")+
    theme(legend.text = element_text(size=8), plot.title = element_text(size=14))+
    geom_hline(yintercept=-log10(cutoff), linetype="dotted")+
    annotate("text", label=paste0("FDR = ",cutoff), fontface="bold", x=max(L.adjResults$setSize)*.95, y=1.4)+
    ylab("-log10 FDR (BH adjusted p value)")+
    ggtitle(title)
  
}


plotPhyper.pval <- function(phypResults, cutoff=.05, title="Enrichment of gene sets within top ranked DE genes"){
  
  results <- as.data.frame(phypResults$phypResults)
  results$geneSets <- rownames(results)
  L.results <- gather(results, "setSize", "adjPval", 1:(ncol(results)-1))
  # add column for color 
  L.results$sig <- ifelse(L.results$adjPval < cutoff, "sig", "ns")
  sigGS <- filter(L.results, sig=="sig")$geneSets
  L.results$sigColor <- NA
  L.results[L.results$geneSets %in% sigGS, "sigColor"] <- L.results[L.results$geneSets %in% sigGS, "geneSets"]
  
  L.results$setSize <- as.numeric(L.results$setSize)
  
  ggplot(L.results, aes(x=setSize, y=-log10(adjPval), group=geneSets, color=sigColor))+
    geom_line(size=1)+
    scale_color_discrete(name="Gene Set")+
    theme(legend.text = element_text(size=8), plot.title = element_text(size=14))+
    geom_hline(yintercept=-log10(cutoff), linetype="dotted")+
    annotate("text", label=paste0("p = ",cutoff), fontface="bold", x=max(L.results$setSize)*.95, y=1.4)+
    ggtitle(title)+
    ylab("-log10 unadjusted p value")
  
}


# test
# plotPhyper.pval(AT.PostVsPre.GSE, title="Pathways in AT post vs pre FGF21 DE genes")
# plotPhyper.adj(AT.PostVsPre.GSE, title="Pathways in AT post vs pre FGF21 DE genes")


#### phyper for a any list of genes of interest, not ranked/no set sizes #####


GSenrich.GOIs <- function(GOIs, fullGeneList, geneSets){
  
  ### function for one gene set, that will go into apply function over columns of geneSet
  
  geneSets <- lapply(geneSets, function(x) as.character(x))
  oneGeneSet <- geneSets[[1]]
  
  # vector of genes that will serve as full list of possible genes  
  fullGeneList <- unique(na.omit(fullGeneList))

  phyp.GOIs <- function(oneGeneSet) {
    geneSet <- na.omit(as.character(oneGeneSet))
    fullList <- as.character(fullGeneList)  
    numBlack <- sum(GOIs %in% geneSet)
    totBlack <- length(geneSet)
    totWhite <- length(fullList)-totBlack
    handful <- length(GOIs)
    phyper(numBlack, totBlack, totWhite, handful, lower.tail = F) 
  }
  
  ## testing... testing... apply to whole geneSet table for one module
   pvals <- unlist(lapply(geneSets, phyp.GOIs))
   
  # make a table with pvals and adjusted pvals
  results <- data.frame(GeneSet = names(geneSets), pVal=pvals)
  results$adjPval <- p.adjust(results$pVal)
  results[,2:3] <- apply(results[,2:3], c(1,2), function(x) signif(x,2))
  results <- results %>% arrange(pVal)
  
  #return from function
  results
}

# #testing
# test.results <- GSenrich.GOIs(GOIs=GOIs, geneSets=Hallmark.GS, fullGeneList=as.character(genes$GeneID))

# now list genes in each significant or high-ranked gene set...

# GOIs=GOIs
# GSresults<- test.results
# geneSets <- Hallmark.GS
# nTopGS = 10

# gene set enrichment results must be ranked  for this function

genesInSigGS <- function(GOIs, GSresults, geneSets, nTopGS=10){
  
  # grab matching gene sets
  topGS <- geneSets[names(geneSets) %in% as.character(GSresults$GeneSet[1:nTopGS])]
  # put top genesets in same order as results list
  topGS <- topGS[as.character(GSresults$GeneSet[1:nTopGS])]
  
  genesPres.ensg <- list()
  
  for (i in (1:nTopGS)){
  genesPres.ensg[[i]] <- GOIs[GOIs %in% topGS[[i]]]
  }
  
  names(genesPres.ensg) <- names(topGS)
  
  ## ** feature to add: spit out second list of genes translated to HGNCsymbol
  genesPres.hgnc <- list()
  for (i in 1:nTopGS){
    ensg <- data.frame(GeneID=genesPres.ensg[[i]])
    hgnc <- merge(ensg, genes[,c("GeneID", "HGNCsymbol")], by="GeneID")
    hgnc <- hgnc[ensg$GeneID,]
    genesPres.hgnc[[i]] <- as.character(hgnc$HGNCsymbol)
  }
  names(genesPres.hgnc) <- names(topGS)
  
  list(genesPres.ensg=genesPres.ensg, genesPres.hgnc= genesPres.hgnc)
}

#test
 #genesPresent <- genesInSigGS(GOIs, test.results, Hallmark.GS, nTopGS=10)







