
##### plotting themes #####

smallFacetLabels <- theme_bw(14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.key = element_blank(), plot.title = element_text(size=14))+
  theme(strip.text.x = element_text(size=10, margin = margin(.1, 0, .1, 0, "cm")), 
        axis.text = element_text(size=10))

theme.14 <- theme_bw(14) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.key = element_blank(), plot.title = element_text(size=14))

theme.10 <- theme_bw(10) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.key = element_blank(), plot.title = element_text(size=10))

# share a legend between two plots
share_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#use like:
# commlegend <- share_legend(p1)
# p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
#                                p2 + theme(legend.position="none"),
#                                nrow=1), commlegend, nrow=2, heights=c(10, 1))




##### plot gene expression #####

plotTop10.paired <- function(top10=top10, counts=lcpm.normCounts.all, by="Treatment", anno, title="Expression of top 10 DE genes", ncol=3){
  top10counts <- counts[rownames(counts) %in% top10$MMG,]
  top10counts <- merge(MMGtoHGNC, top10counts, by.x="MMG", by.y="row.names")
  top10counts <- top10counts[match(top10$MMG, top10counts$MMG),]
  sampleID <- colnames(top10counts)[-c(1:4)]
  top10counts <- transpose(top10counts)
  colnames(top10counts) <- top10counts[4,]
  top10counts <- top10counts[-c(1:4),]
  top10counts$SampleID <- sampleID
  top10counts <- merge(top10counts, anno[,1:13], by="SampleID")
  L.top10 <- gather(top10counts, "Gene", "Counts", 2:(nrow(top10)+1))
  L.top10$Counts <- as.numeric(L.top10$Counts)
  L.top10$Treatment <- factor(L.top10$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  L.top10$Gene <- factor(L.top10$Gene, levels=c(unique(L.top10$Gene)))
  
  ggplot(L.top10, aes_string(x=by, y="Counts"))+
    geom_line(aes(group=AnimalID), color="grey30")+
    #geom_boxplot(outlier.color=NA, aes(color=Treatment))+
    geom_point(size=3, shape=21, aes(fill=Treatment))+
    stat_summary(fun.y="mean", geom="point", shape=23, size=2, fill="white")+
    facet_wrap(~Gene, ncol=ncol, scales="free_y")+
    ylab("log2 counts per million")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "none", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}


plotTop10.paired.CR <- function(top10=top10, counts=voomCounts.CR, by="Treatment", anno=annoCR, title="Expression of top 10 DE genes"){
  
  top10counts <- counts[rownames(counts) %in% top10$MMG,]
  top10counts <- merge(MMGtoHGNC, top10counts, by.x="MMG", by.y="row.names")
  sampleID <- colnames(top10counts)[-c(1:4)]
  top10counts <- transpose(top10counts)
  colnames(top10counts) <- top10counts[4,]
  top10counts <- top10counts[-c(1:4),]
  top10counts$SampleID <- sampleID
  top10counts <- merge(top10counts, annoCR[,c(1,7:9)], by="SampleID")
  L.top10 <- gather(top10counts, "Gene", "Counts", 2:(nrow(top10)+1))
  L.top10$Counts <- as.numeric(L.top10$Counts)
  L.top10$Treatment <- factor(L.top10$Treatment, levels =c("CalR.Pre", "CalR.Post"))
  L.top10$Gene <- factor(L.top10$Gene, levels=c(unique(L.top10$Gene)))

  ggplot(L.top10, aes(x=Treatment, y=Counts))+
    geom_line(aes(group=AnimalID), color="grey30")+
    #geom_boxplot(outlier.color=NA, aes(color=Treatment))+
    geom_point(size=3, shape=21, aes(fill=Treatment))+
    stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white")+
    facet_wrap(~Gene, ncol=3, scales="free_y")+
    ylab("log2 counts per million")+
    ggtitle("Top DE genes in calorie-restricted soleus")+
    theme(legend.position = "none", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
  
}



plotGOI.paired.HGNC <- function(GOI=GOI, counts=voomCounts.A, 
                                by="Treatment", anno, group="AnimalID", ncol=5,
                                title="Expression of top 10 DE genes", scale="free_7"){
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$HGNCsymbol),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(GOIcounts, anno, by="SampleID")
  L.GOI <- gather(GOIcounts, "Gene", "Counts", 2:(nGOI+1))
  L.GOI$Counts <- as.numeric(L.GOI$Counts)
  L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  L.GOI$Gene <- factor(L.GOI$Gene, levels=c(unique(L.GOI$Gene)))
  
  ggplot(L.GOI, aes_string(x=by, y="Counts"))+
    geom_line(aes_string(group=group), color="grey30")+
    #geom_boxplot(outlier.color=NA, aes(color=Treatment))+
    geom_point(size=2, shape=21, aes(fill=Treatment))+
    scale_fill_manual(values=c("orangered", "deepskyblue4", "olivedrab"))+
    stat_summary(fun.y="mean", geom="point", shape=23, size=2, fill="white")+
    facet_wrap(~Gene, ncol=ncol, scales=scale)+
    ylab("Gene Expression log2(cpm)")+
    xlab("Treatment Time point")+
    scale_x_discrete(labels=c("Pre", "Post", "Washout"))+
    labs(fill="Treatment")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "none", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}


plotGOI.paired.HGNC.1 <- function(GOI=GOI, counts=voomCounts.A, 
                                by="Treatment", anno, group="AnimalID", title=""){
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$HGNCsymbol),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(GOIcounts, anno, by="SampleID")
  L.GOI <- gather(GOIcounts, "Gene", "Counts", 2:(nGOI+1))
  L.GOI$Counts <- as.numeric(L.GOI$Counts)
  L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  L.GOI$Gene <- factor(L.GOI$Gene, levels=c(unique(L.GOI$Gene)))
  
  ggplot(L.GOI, aes_string(x=by, y="Counts"))+
    geom_line(aes_string(group=group), color="grey30")+
    #geom_boxplot(outlier.color=NA, aes(color=Treatment))+
    geom_point(size=2, shape=21, aes(fill=Treatment))+
    scale_fill_manual(values=c("orangered", "deepskyblue4", "olivedrab"))+
    stat_summary(fun.y="mean", geom="point", shape=23, size=2, fill="white")+
    ylab("Expression \nlog2(cpm)")+
    xlab("Time point")+
    scale_x_discrete(labels=c("Pre", "Post", "Washout"))+
    labs(fill="Treatment")+
    ggtitle(title)+
    theme.10+
    theme(legend.position = "none", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}


# GOI=OBtrim.GS$combo_fatty_acid_metabolism[1:30]
# anno=annoA
# title="Fatty acid metab gene set"

plotGOI.paired.ensg <- function(GOI=GOI, counts=voomCounts.A, by="Treatment", anno, ncol=5, title="Expression of top 10 DE genes"){
  counts <- merge(MMGtoHGNC[,2:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$GeneID),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-c(1:2)]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[2,]
  GOIcounts <- GOIcounts[-c(1:2),,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(GOIcounts, anno, by="SampleID")
  L.GOI <- gather(GOIcounts, "Gene", "Counts", 2:(nGOI+1))
  L.GOI$Counts <- as.numeric(L.GOI$Counts)
  L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  L.GOI$Gene <- factor(L.GOI$Gene, levels=c(unique(L.GOI$Gene)))
  
  ggplot(L.GOI, aes_string(x=by, y="Counts"))+
    geom_line(aes(group=AnimalID), color="grey30")+
    #geom_boxplot(outlier.color=NA, aes(color=Treatment))+
    geom_point(size=3, shape=21, aes(fill=Treatment))+
    stat_summary(fun.y="mean", geom="point", shape=23, size=2, fill="white")+
    facet_wrap(~Gene, ncol=ncol, scales="free_y")+
    ylab("log2 counts per million")+
    ggtitle(title)+
    theme(legend.position = "none", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}





# counts=lcpm.normCounts.AT
# by="Bodyweight"
# anno=forcor.AT.BWlag
# title="DEG associated with % change in BW: adipose tissue"


plotGOI.phenCor.HGNC <- function(GOI=GOI$HGNCsymbol, counts=voomCounts.A, pheno="BWpctCh.lag.wash2", anno=AT.PD.T, 
                                 color="Treatment", title="Genes associated with phenotypical data", ncol=4){
  
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$HGNCsymbol),-1]
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(GOIcounts, anno, by="SampleID")
  L.GOI <- gather(GOIcounts, "Gene", "Counts", 2:(length(GOI)+1))
  L.GOI$Counts <- as.numeric(L.GOI$Counts)
  L.GOI$Gene <- factor(L.GOI$Gene, levels=c(unique(L.GOI$Gene)))
  L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  
  
  ggplot(L.GOI, aes_string(x=pheno, y="Counts", color=color))+
    geom_point(size=3)+
    stat_smooth(method="lm", se=FALSE)+
    facet_wrap(~Gene, ncol=ncol,scales="free")+
    ylab("Gene Expression log2(cpm)")+
    labs(color="Treatment")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "bottom")
}



# GOI=HandC2.GS$REACTOME_RAP1_SIGNALLING
# pheno="TAGpctCh_maxCh"

plotGOI.phenCor.ensg <- function(GOI, counts=voomCounts.A, pheno="BWpctCh.lag.wash2", anno=annoAT.PD, 
                                 color="Sex", title="Genes associated with phenotypical data", ncol=4){
  
  counts <- merge(MMGtoHGNC[,c(2:4)], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$GeneID),]
  GOIcounts <- GOIcounts[,-c(1:2)]
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(GOIcounts, anno, by="SampleID")
  L.GOI <- gather(GOIcounts, "Gene", "Counts", 2:(length(GOI)+1))
  L.GOI$Counts <- as.numeric(L.GOI$Counts)
  L.GOI$Gene <- factor(L.GOI$Gene, levels=c(unique(L.GOI$Gene)))
  #remove genes with NA name
  L.GOI <- na.omit(L.GOI)
  
  ggplot(L.GOI, aes_string(x=pheno, y="Counts", color=color))+
    geom_point(size=3)+
    stat_smooth(method="lm", se=FALSE)+
    facet_wrap(~Gene, ncol=ncol,scales="free")+
    ylab("Gene Expression")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "bottom", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}


# GOI=top10$HGNCsymbol
# IDtype="HGNCsymbol"
# pheno="BW.con"
# anno=AT.PD.T
# pointColor="AnimalID"
# color="AnimalID"
# ncol=4

plotGOI.pairedxTime <- function(GOI=GOI, IDtype="HGNCsymbol", counts=voomCounts.A, pheno="BWpctCh.lag.wash2", 
                                group="AnimalID", anno=AT.PD.T, 
                                 pointColor="Treatment", xlab=pheno, title=NULL, ncol=4){
  
  counts <- merge(MMGtoHGNC[,c("MMG",IDtype)], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts[,IDtype]),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(GOIcounts, anno, by="SampleID")
  L.GOI <- gather(GOIcounts, "Gene", "Counts", 2:(nGOI+1))
  L.GOI$Counts <- as.numeric(L.GOI$Counts)
  L.GOI$Gene <- factor(L.GOI$Gene, levels=c(unique(L.GOI$Gene)))
  L.GOI$Treatment <- factor(L.GOI$Treatment, levels=c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  
  ggplot(L.GOI, aes_string(x=pheno, y="Counts"))+
    geom_line(aes_string(group=group), color="grey30")+
    geom_point(size=3, shape=21, aes(fill=Treatment))+
    #scale_fill_manual(values=c("orangered2", "steelblue"))+
    #stat_smooth(method="lm", se=FALSE)+
    facet_wrap(~Gene, ncol=ncol,scales="free")+
    ylab("Gene Expression log2(cpm)")+
    xlab(pheno)+
    labs(fill="Treatment")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "bottom" )
}



# GOI <- c("ANGPTL8", "ELOVL3", "SCD")
# counts = voomCounts.A
# pheno="TAGpctCh_maxCh"
# anno=annoAcor

plotGOIchVphenoCh <-  function(counts, GOI, pheno, anno, ncol=5, scale="free_y",
                               title="Change in gene expression and pheno variable"){
  
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$HGNCsymbol),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(anno[c("SampleID", "AnimalID", "Treatment")], GOIcounts, by="SampleID")
  GOIcounts <- merge(GOIcounts, phenoAnimal, by.x="AnimalID", by.y="AnimalID")
  GOIchange <- gather(GOIcounts, "Gene", "Counts", 4:(nGOI+3) )
  GOIchange <- GOIchange %>% select(-SampleID) %>% spread( "Treatment", "Counts") 
  GOIchange <- mutate(GOIchange, Post.Pre=as.numeric(FGF21.Post)-as.numeric(FGF21.Pre), Wash.Post = as.numeric(FGF21.Wash)-as.numeric(FGF21.Post))
  PostvPre <- GOIchange %>% select(-c(FGF21.Post:FGF21.Wash, Wash.Post)) 

  
  # L.GOI$Counts <- as.numeric(L.GOI$Counts)
  # L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  PostvPre$Gene <- factor(PostvPre$Gene, levels=c(unique(PostvPre$Gene)))
  # 
 
  
  ggplot(PostvPre, aes_string(x=pheno, y="Post.Pre"))+
    geom_point(size=2, color="navyblue")+
    stat_smooth(method="lm", se=FALSE, color="grey20")+
    facet_wrap(~Gene, ncol=ncol,scales=scale)+
    ylab("Change in gene expression")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "none")
}




plotGOIpctChVphenoCh <-  function(counts, GOI, pheno, anno, ncol=5, scale="free_y",
                               title="Change in gene expression and pheno variable"){
  
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$HGNCsymbol),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(anno[c("SampleID", "AnimalID", "Treatment")], GOIcounts, by="SampleID")
  GOIcounts <- merge(GOIcounts, phenoAnimal, by.x="AnimalID", by.y="AnimalID")
  GOIchange <- gather(GOIcounts, "Gene", "Counts", 4:(nGOI+3) )
  GOIchange <- GOIchange %>% select(-SampleID) %>% spread( "Treatment", "Counts") 
  GOIchange <- mutate(GOIchange, Post.Pre=((as.numeric(FGF21.Post)-as.numeric(FGF21.Pre))/as.numeric(FGF21.Post))*100, Wash.Post = as.numeric(FGF21.Wash)-as.numeric(FGF21.Post))
  PostvPre <- GOIchange %>% select(-c(FGF21.Post:FGF21.Wash, Wash.Post)) 
  
  
  # L.GOI$Counts <- as.numeric(L.GOI$Counts)
  # L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  PostvPre$Gene <- factor(PostvPre$Gene, levels=c(unique(PostvPre$Gene)))
  # 
  
  
  ggplot(PostvPre, aes_string(x=pheno, y="Post.Pre"))+
    geom_point(size=2, color="navyblue")+
    stat_smooth(method="lm", se=FALSE, color="grey20")+
    facet_wrap(~Gene, ncol=ncol,scales=scale)+
    ylab("Percent change in gene expression")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "none")
}




plotGOIchVphenoCh.1 <-  function(counts, GOI, pheno, anno){
  
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(GOI, counts$HGNCsymbol),-1]
  GOIcounts <- na.omit(GOIcounts)
  nGOI <- nrow(GOIcounts)
  sampleID <- colnames(GOIcounts)[-1]
  GOIcounts <- transpose(GOIcounts)
  colnames(GOIcounts) <- GOIcounts[1,]
  GOIcounts <- GOIcounts[-1,,drop=FALSE]
  GOIcounts$SampleID <- sampleID
  GOIcounts <- merge(anno[c("SampleID", "AnimalID", "Treatment")], GOIcounts, by="SampleID")
  GOIcounts <- merge(GOIcounts, phenoAnimal, by.x="AnimalID", by.y="AnimalID")
  GOIchange <- gather(GOIcounts, "Gene", "Counts", 4:(nGOI+3) )
  GOIchange <- GOIchange %>% select(-SampleID) %>% spread( "Treatment", "Counts") 
  GOIchange <- mutate(GOIchange, Post.Pre=as.numeric(FGF21.Post)-as.numeric(FGF21.Pre), Wash.Post = as.numeric(FGF21.Wash)-as.numeric(FGF21.Post))
  PostvPre <- GOIchange %>% select(-c(FGF21.Post:FGF21.Wash, Wash.Post)) 
  
  
  # L.GOI$Counts <- as.numeric(L.GOI$Counts)
  # L.GOI$Treatment <- factor(L.GOI$Treatment, levels =c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  PostvPre$Gene <- factor(PostvPre$Gene, levels=c(unique(PostvPre$Gene)))
  # 
  
  
  ggplot(PostvPre, aes_string(x=pheno, y="Post.Pre"))+
    geom_point(size=2, color="navyblue")+
    stat_smooth(method="lm", se=FALSE, color="grey20")+
    ylab("Change in \nexpression")+
    smallFacetLabels+
    theme.10+
    theme(legend.position = "none")
}




#### plot modules #####


# # testing
# Musc.MEs <- read.csv("Results/WGCNA/All Genes/M.allgenes.MEexp.csv")
# names(Musc.MEs)[1] <- "SampleID"
# 
# MEs <- Musc.MEs
# anno <- annoM[,-c(14:46)]

plotModExpr <- function(MEs, anno, title="Module expression", ncol=ncol){
  
  ## merge module eigengenes and sample data
  MEanno <- merge(MEs, anno, by="SampleID")
  
  L.MEanno <- gather(MEanno, "Module", "ME.expr", 2:(ncol(MEs)) )
  L.MEanno$TissTreat <- paste(L.MEanno$Tissue, L.MEanno$Treatment, sep=".")
  # put groups in logical order
  L.MEanno$TissTreat <- factor(L.MEanno$TissTreat, levels = c("Adipose.FGF21.Pre", "Adipose.FGF21.Post", "Adipose.FGF21.Wash", "Gastrocnemius.FGF21.Pre", "Gastrocnemius.FGF21.Post", "Gastrocnemius.FGF21.Wash", "Soleus.CalR.Pre", "Soleus.CalR.Post"))
  L.MEanno$Treatment <- factor(L.MEanno$Treatment, levels = c("FGF21.Pre", "FGF21.Post", "FGF21.Wash", "CalR.Pre", "CalR.Post"))    
  # put modules in order of original list
  L.MEanno$Module <- factor(L.MEanno$Module, levels = unique(L.MEanno$Module))
  
  
  # for treatment group
  ggplot(L.MEanno, aes(x=Treatment, y=ME.expr))+
    #geom_boxplot(outlier.colour = NA)+
    geom_line(aes(group=AnimalID), color="grey30")+
    geom_point(size=2, shape=21, aes(fill=Treatment))+
    scale_fill_manual(values=c("orangered", "deepskyblue4", "olivedrab"), labels=c("Pre", "Post", "Washout"))+
    #geom_point(size=2, position = position_jitter(width=.2))+
    stat_summary(fun.y="mean", geom="point", shape=23, size=2, fill="white")+
    facet_wrap(~Module, ncol=ncol, scales="free_y")+
    xlab("Treatment Time point")+
    scale_x_discrete(labels=c("Pre", "Post", "Washout"))+
    ylab("Module eigengene expression")+
    labs(fill="Treatment")+
    smallFacetLabels+
    theme(legend.position = "top",  axis.text.x = element_text(angle=20, hjust=1, vjust=1))+
    ggtitle(title)
  
}

# #testing
# plotModules(Musc.MEs, annoM[,-c(14:46)])



# testing
# MEs <- AT.MEs[,-ncol(AT.MEs)]
# Mods <- BWxTimeMods
# Phenos <- c("BWpctCh.con", "BWpctCh.lag.maxCh", "BWpctCh.lag.w6", "BWpctCh.lag.wash2")
# anno <- AT.PD.T

plotModPhenoCor <- function(MEs, anno, Mods, Phenos, color=NULL, ncol=5, title="Module Pheno Correlation"){
  
  ## merge module eigengenes, sample data, and phenodata by sample ID to get all time points
  nMods <- length(Mods)
  nPheno <- length(Phenos)
  
  anno <- anno %>% select(one_of(Phenos, "SampleID", "Treatment", "Sex", "AnimalID"))
  MEs <- MEs %>% select(one_of(Mods, "SampleID"))
  ME.PD.T <- merge(anno, MEs, by="SampleID")
  
  L.MEanno <- gather(ME.PD.T, "Measurement", "Value", 2:(nPheno+1) )
  L.MEanno <- gather(L.MEanno, "Module", "ME.expr", 5:(nMods+4)) # have to change this if # of pheno variables changes
  # L.MEanno$TissTreat <- paste(L.MEanno$Tissue, L.MEanno$Treatment, sep=".")
  # # put groups in logical order
  # L.MEanno$TissTreat <- factor(L.MEanno$TissTreat, levels = c("Adipose.FGF21.Pre", "Adipose.FGF21.Post", "Adipose.FGF21.Wash", "Gastrocnemius.FGF21.Pre", "Gastrocnemius.FGF21.Post", "Gastrocnemius.FGF21.Wash", "Soleus.CalR.Pre", "Soleus.CalR.Post"))
  # L.MEanno$Treatment <- factor(L.MEanno$Treatment, levels = c("FGF21.Pre", "FGF21.Post", "FGF21.Wash", "CalR.Pre", "CalR.Post"))    
  # put modules in order of original list
  L.MEanno$Module <- factor(L.MEanno$Module, levels = unique(L.MEanno$Module))
  
 
  ggplot(L.MEanno, aes_string(x="Value", y="ME.expr", color=color))+
  geom_point(size=3)+
  stat_smooth(method="lm", se=FALSE)+
  facet_wrap(~Module, ncol=ncol,scales="free")+
  ylab("Module expression")+
  ggtitle(title)+
  theme(legend.position = "none", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}


# MEs <- AT.MEs
# anno <- AT.PD.T
# Mods <- topBW$Module
# pheno="BW.con"
# pointColor="AnimalID"


plotModPhenoCor.pairedxTime <- function(MEs, anno, Mods, pheno, pointColor="AnimalID", xlab=pheno, title=NULL, ncol=5){
  
  ## merge module eigengenes, sample data, and phenodata by sample ID to get all time points
  nMods <- length(Mods)
  nAnno <- ncol(anno)
  endMods <- nMods+nAnno
  
  MEs <- MEs %>% select(one_of(Mods, "SampleID"))
  ME.PD.T <- merge(anno, MEs, by="SampleID")
  L.MEanno <- gather(ME.PD.T, "Module", "ME.expr", (nAnno+1):endMods) # have to change this if # of pheno variables changes
  # put modules in order of original list
  L.MEanno$Module <- factor(L.MEanno$Module, levels = unique(L.MEanno$Module))
  L.MEanno$Treatment <- factor(L.MEanno$Treatment, levels=c("FGF21.Pre", "FGF21.Post", "FGF21.Wash"))
  
  
  ggplot(L.MEanno, aes_string(x=pheno, y="ME.expr"))+
    geom_line(aes(group=AnimalID), color="grey30")+
    geom_point(size=3, shape=21, aes(fill=Treatment))+
    #scale_fill_manual(values=c("orangered2", "steelblue"))+
    #stat_smooth(method="lm", se=FALSE)+
    facet_wrap(~Module, ncol=ncol,scales="free")+
    ylab("Module eigengene expression")+
    xlab(pheno)+
    labs(fill="Treatment")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "bottom", axis.text.x = element_text(angle=20, hjust=1, vjust=1))
}


# ME.change = AT.ME.PostvPre
# Mods = mods <- c("MEmediumorchid", "MEpink", "MEgrey60", "MEturquoise", "MEplum3", "MEdarkred", "MEcoral3")
# pheno = "BWpctChange_Wash.w2"


plotChangeModPhenoCor <- function(ME.change, Mods, pheno,  ncol=5, scale="free_y",
                                  title="Change in module expression and pheno variable"){
  
  ## merge module eigengenes, sample data, and phenodata by sample ID to get all time points
  nMods <- length(Mods)
  ME.change <- ME.change %>% select(one_of("AnimalID", "Sex", pheno, Mods))
  L.MEchange <- gather(ME.change, "Module", "ME.expr", 4:(nMods+3)) # have to change this if # of pheno variables changes
  # L.MEanno$TissTreat <- paste(L.MEanno$Tissue, L.MEanno$Treatment, sep=".")
  # # put groups in logical order
  # L.MEanno$TissTreat <- factor(L.MEanno$TissTreat, levels = c("Adipose.FGF21.Pre", "Adipose.FGF21.Post", "Adipose.FGF21.Wash", "Gastrocnemius.FGF21.Pre", "Gastrocnemius.FGF21.Post", "Gastrocnemius.FGF21.Wash", "Soleus.CalR.Pre", "Soleus.CalR.Post"))
  # L.MEanno$Treatment <- factor(L.MEanno$Treatment, levels = c("FGF21.Pre", "FGF21.Post", "FGF21.Wash", "CalR.Pre", "CalR.Post"))    
  # put modules in order of original list
  L.MEchange$Module <- factor(L.MEchange$Module, levels = unique(L.MEchange$Module))
  
  
  ggplot(L.MEchange, aes_string(x=pheno, y="ME.expr"))+
    geom_point(size=2, color="navyblue")+
    stat_smooth(method="lm", se=FALSE, color="grey20")+
    facet_wrap(~Module, ncol=ncol, scales=scale)+
    ylab("Module expression")+
    ggtitle(title)+
    smallFacetLabels+
    theme(legend.position = "none")
}


###### for heatmaps #####

# counts <- voomCounts.A
# GOI=Imm1

grabGOIcounts <- function(counts, GOI, anno){
  counts <- merge(MMGtoHGNC[,3:4], counts, by.x="MMG", by.y="row.names")
  GOIcounts <- counts[match(unique(GOI), counts$HGNCsymbol),-1]
  GOIcounts <- unique(na.omit(GOIcounts))
  rownames(GOIcounts) <- GOIcounts$HGNCsymbol
  GOIcounts <- GOIcounts[,-1]
  GOIcounts <- GOIcounts[,match(anno$SampleID,colnames(GOIcounts))]
  GOIcounts
}



