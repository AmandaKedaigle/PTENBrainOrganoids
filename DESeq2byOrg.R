#Takes a Seurat Object and performs differential expression analysis 
#By splitting the object into seperate samples (i.e organoids) and using DESeq2,
#which accounts for noise between samples and then looks for DEGs that overcome that noise
#based on method by Sean Simmons

library(Seurat)
library(DESeq2)
library(writexl)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(viridis)
library(RRHO2)
library(patchwork)


setwd("~/Documents/DorsalKadoshima/PTEN_FinalObjects")


#Function that performs the DESeq calculation by summing counts in each sample
combineDE<-function(seur,id,condition="treat",base="wt",combineOn="orig.ident",
                    minCells=20, #Minimum # of cells that must be in this cluster per sample to keep that sample
                    minBatches=2, #Minimum # of samples you can have per condition. Cannot be lower than 2.
                    minReads=10, #Mininum # of reads to have total per gene to calculate DE for that gene
                    genes=c(),  #Genes to consider for DE analysis, if you don't want to use all expressed genes.
                    form="" #design formula for DESeq2 if you want to include more variables in addition to "condition"
                    )
{
  print("Subsample")
  seur<-subset(seur,idents=c(id))
  
  Idents(seur)=combineOn
  genes.use=rownames(GetAssayData(seur,slot="counts"))
  if(length(genes)>0){genes.use=genes}
  
  print("Combine data per sample")
  data.all=data.frame(row.names = genes.use)
  for(i in levels(Idents(seur))) {
    temp.cells=WhichCells(seur,ident=i)
    if (length(temp.cells)==1) data.temp=(GetAssayData(seur,slot="counts")[genes.use,temp.cells])
    if (length(temp.cells)>1) data.temp=apply(GetAssayData(seur,slot="counts")[genes.use,temp.cells],1,sum)
    data.all=cbind(data.all,data.temp)
    colnames(data.all)[ncol(data.all)]=i
  }
  
  print("Filter samples for minimum cells")
  keepOrgs=names(summary(Idents(seur)))[summary(Idents(seur))>minCells]
  numOrg=length(keepOrgs)
  print(paste("Keeping", numOrg, "samples"))
  data.all=data.all[,keepOrgs]

  extraColumns<-strsplit(form,"+",fixed=T)[[1]]
  val=seur@meta.data[,c(condition,combineOn,extraColumns)]
  val=val[!duplicated(val[,2]),]
  rownames(val)=val[,2]
  keepBatch=as.character(val[keepOrgs,1])
  levels = levels(factor(keepBatch))
  if(length(levels)<2) {
    print("Not enough batches per treatment group with minimum # of cells!")
    return(NULL)
  }
  for (level in levels) {
    if(sum(keepBatch==level)<2) {
      print("Not enough batches per treatment group with minimum # of cells!")
      return(NULL)
    }
  }
  
  print("Save meta data")
  colDat=factor(keepBatch)
  if (base != "") { colDat = relevel(colDat, ref=base)}
  colDat=data.frame(colDat)
  colnames(colDat)="condition"
  rownames(colDat)=colnames(data.all)
  colDat[keepOrgs,extraColumns]=val[keepOrgs,extraColumns]
  
  print("Run DESeq2")
  design= ~ condition
  if(nchar(form)>0){
    design=as.formula(paste("~",form,"+ condition",sep=""))
  }
  print(design)
    
  dds <- DESeqDataSetFromMatrix(countData = data.all,colData = colDat,design = design)
  dds <- dds[ rowSums(counts(dds) > minReads)>=2, ]
  dds <- DESeq(dds)
  out=data.frame(results(dds))
  out=out[order(out$pvalue),]
  return(out)
}


  
#Load Seurat Object
seur = readRDS("PGP1_1m_rep1_celltypes.rds")
seur = subset(seur, subset=treat %in% c("mut", "wt"))
  
#Set Ident to cluster or celltypes, whatever groups you want to seperate before detecting DEGs in that group
Idents(seur) = "CellType"
  
#Set condition to the metadata column you want DEGs between, and base to the "wildtype" or base level of that column
condition = "treat"
base = "wt"
  
#set combineOn to the metadata column that contains the samples (i.e. different organoids)
combineOn = "orig.ident"
  
dir="DEGs_PGP1_1m_wtvMut/"
dir.create(dir)
  
#This loop will run DE analysis for each cluster and save a .xlsx file for each!
ids=levels(seur)
alldegs = list()
for (i in ids) {
  degs <- combineDE(seur, id=i, condition=condition, base=base, combineOn=combineOn)
  degs$gene = rownames(degs)
  #if (length(degs)>0) {
  #  write_xlsx(degs, path=paste0(dir,gsub("/","-",i),".DEGs.xlsx"), format_headers = T)
  #}
  alldegs[[i]] = degs
}
saveRDS(alldegs, paste0(dir,"DEGs.rds"))

  
#Compare to atlas DEGs over time
atlas = readRDS("../DEGs_Celltypes_Over_Time_23d-1.5m.RDS")
names(atlas)[names(atlas)=="Newborn DL PN"] = "Immature DL PN"
names(atlas)[names(atlas)=="Newborn PN"] = "Newborn DL PN"
for (i in names(alldegs)) {
  atlasi = atlas[[i]]
  atlasi$gene = rownames(atlasi)
  degs = alldegs[[i]]
  genes = degs[!is.na(degs$padj) & degs$padj<0.01,"gene"]
  
  #RRHO
  deg.list = data.frame(gene=degs$gene, value = (-1*log10(degs$padj)*ifelse(degs$log2FoldChange>0,1,-1)))
  deg.list = deg.list[!is.na(deg.list$value),]
  
  atlas.list = data.frame(gene=atlasi$gene, value = (-1*log10(atlasi$padj)*ifelse(atlasi$log2FoldChange>0,1,-1)))
  atlas.list = atlas.list[!is.na(atlas.list$value),]
  
  deg.list = deg.list[deg.list$gene %in% atlas.list$gene,]
  atlas.list = atlas.list[atlas.list$gene %in% deg.list$gene,]
  
  RRHO2_obj <-  RRHO2_initialize(deg.list, atlas.list, labels = c(paste0(i, " MT v WT"), paste0(i, " Atlas Over Time")),  boundary=0.03)
  png(paste0(dir, "RRHO-", i, ".png"))
  print(RRHO2_heatmap(RRHO2_obj))
  dev.off()
  
  uu = RRHO2_obj$genelist_uu$gene_list_overlap_uu
  dd = RRHO2_obj$genelist_dd$gene_list_overlap_dd
  
  #diverging bars
  geneOrder = degs[order(degs$log2FoldChange),"gene"]
  degs$gene = factor(degs$gene, levels=geneOrder)
  atlasi$gene = factor(atlasi$gene, levels=geneOrder)
  p1=ggplot(degs[genes,], aes(x=gene, y=log2FoldChange, fill=log2FoldChange)) + 
    geom_bar(stat="identity") +
    #scale_fill_gradientn(values=c(0,0.125,1), colors=c("gray","pink","red"),limits=c(0,8),name="-log(adj. p value)") +
    scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,name="Fold Change") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "bottom")+
    xlab(NULL) + ylab(paste0(i, " MT v WT DEGs"))
  p2= ggplot(atlasi[genes,], aes(x=gene, y=log2FoldChange,fill=-1*log10(padj))) +
    geom_bar(stat="identity") +
    scale_fill_gradientn(values=c(0,0.125,1), colors=c("gray","pink","red"),limits=c(0,-1*log10(min(atlasi[genes,"padj"]))),name="-log(adj. p value)") +
    #scale_fill_gradient2(low="blue",mid="white",high="red",midpoint=0,name="Fold Change",limits=c(-2,2.5)) +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none",axis.text.y = element_blank(), 
          axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
    xlab(NULL) + ylab(paste0(i, " Atlas Over Time DEGs"))
  
  pdf(paste0(dir, "DEGs-",i,"-BarPlot.pdf"), width=5.5, height=3+0.05*length(genes))
  print(p1 + p2 + plot_layout(guides="collect") & theme(legend.position = "bottom"))
  dev.off()
  
}
  

#Read in DEGs for each and find GO terms -----
datasets =  c("DEGs/DEGs_Mito_1m_wtvMut/DEGs.rds","DEGs/DEGs_PGP1_1m_wtvMut/DEGs.rds","DEGs/DEGs_PGP1_1m_wtvKO/DEGs.rds")
names = c("Mito_1m","PGP1_1m","PGP1_1m_wtvKO")

datasets = c("DEGs/DEGs_Mito_3m_rep1_wtvMut/DEGs.rds","DEGs/DEGs_Mito_3m_rep2_wtvMut/DEGs.rds","DEGs/DEGs_PGP1_3m_wtvMut/DEGs.rds","DEGs/DEGs_PGP1_3m_wtvKO/DEGs.rds")
names = c("Mito_3m_r1","Mito_3m_r2","PGP1_3m","PGP1_3m_wtvKO")

datasets = c("DEGs_1m_Celltypes_Over_WT.RDS","DEGs_1m_Celltypes_Over_Mut.RDS", "DEGs_3m_Celltypes_Over_WT.RDS", "DEGs_3m_Celltypes_Over_Mut.RDS")
names = c("WTvWT_1m","MutvMut_1m", "WTvWT_3m","MutvMut_3m")

for (d in 1:length(datasets)) {
allGenesUp = list()
allGenesDown = list()
allGenes = list()
dat = data.frame()
for (ct in c("aRG", "oRG","IP","Newborn DL PN","Immature DL PN",
     "Subcortical progenitors","Subcortical neurons","Cortical hem",
     "Cajal Retzius","PN","CFuPN","CPN","Immature IN")) {
  #for (d in 1:length(datasets)) {
    dat = readRDS(datasets[[d]])
    if (ct %in% names(dat)) {
      res = dat[[ct]]
      genesUp = rownames(res[!is.na(res$padj) & res$padj<0.05 & res$log2FoldChange>0,])
      genesDown = rownames(res[!is.na(res$padj) & res$padj<0.05 & res$log2FoldChange<0,])
      genes = rownames(res[!is.na(res$padj) & res$padj<0.05,])
      allGenesUp[[ct]] = genesUp
      allGenesDown[[ct]] = genesDown
      allGenes[[ct]] = genes
    }
  }
  
 # if(length(allGenesUp) >0 | length(allGenesDown)>0) {
    ckU = compareCluster(allGenesUp, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")
    ckUs = simplify(ckU, cutoff=0.5)
    ckD = compareCluster(allGenesDown, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")
    ckDs = simplify(ckD, cutoff=0.5)
    #ck = compareCluster(allGenes, fun="enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP")
    #ck = simplify(ck, cutoff=0.5)

    #allGenesUni = lapply(allGenes, function(x) {bitr(x, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")$UNIPROT})
    #ck = compareCluster(allGenesUni, fun="enrichKEGG", organism = "hsa", keyType = "uniprot")
    
    p1 = dotplot(ckUs, showCategory=3) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
      ggtitle("Upregulated in PGP1") + scale_color_viridis(option="D", end=0.8)

    p2 = dotplot(ckDs, showCategory=3) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
      ggtitle("Upregulated in Mito210") + scale_color_viridis(end=0.8)
    
    #p3 = dotplot(ck, showCategory=5) + theme(axis.text.x=element_text(angle=60, hjust=1)) + 
    #  ggtitle(paste0("Changed between CellLines ", ct)) + scale_color_viridis(end=0.8)

    pdf(paste0(names[[d]],"-GO.pdf"), width=20, height=9)
    print(p1 + p2)
    dev.off()
}
  #}


