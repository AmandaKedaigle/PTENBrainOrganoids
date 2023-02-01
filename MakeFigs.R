library(Seurat)
library(lme4)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)
library(patchwork)
library(ggpattern)

setwd("~/Documents/DorsalKadoshima/PTEN_FinalObjects/")

objs = c("Mito_1m_rep1_celltypes.rds","Mito_3m_rep1_celltypes.rds","Mito_3m_rep2_celltypes.rds",
         "PGP1_1m_rep1_celltypes.rds","PGP1_3m_celltypes.rds")
names = c("Mito_1m_rep1","Mito_3m_rep1","Mito_3m_rep2","PGP1_1m_rep1","PGP1_3m_rep1")

for (i in 1:length(objs)) {
  seur = readRDS(objs[[i]])
  name = names[[i]]
  Idents(seur) = "CellType"
  levels =  c("aRG","Cajal Retzius","Cortical hem","Cycling Progenitors","IP",
              "Immature DL PN","Newborn DL PN","Subcortical neurons", "Subcortical progenitors", 
              "Unknown", "Choroid plexus/Cortical hem",
              "CFuPN","CPN","oRG","PN","Immature IN",
              "oRG/Astroglia","Astroglia")
  cols = c("#41ae76","#ee8866","#bebada","#bbcc33","#fdb462",
           "#f768a1","#fa9fb5","#77aadd", "#77d3dd",
           "darkgray", "powderblue",
           "#cc6677","#882255","#225522","#aa4499","#332288",
           "#009988","#0077bb")
  cols1 = cols[match(levels(seur),levels)]
  # tiff(filename = paste0("TSNEs/",name,".colors.celltype.tsne.tiff"),res=200,width=600,height=600)
  # print(DimPlot(seur,pt.size=0.1, cols=cols1, reduction = "tsne") + NoAxes() + NoLegend())
  # dev.off()

  numOrgs = length(unique(seur$orig.ident))
  seur$org = paste0(seur$treat, seur$orig.ident)
  seur$org = factor(seur$org, levels = c("wt1", "wt2","wt3","wt4","wt5","wt7","wt8","wt9","mut1","mut2","mut3","mut4","mut5","mut6","ko4","ko5","ko6","ko7","ko8","ko9"))
  
  #TSNEs of indiv orgs
  # tiff(filename = paste0("TSNEs/",name,"eachOrg.colors.celltype.tsne.tiff"),res=200,width=1500,height=600*(numOrgs/3))
  # print(DimPlot(seur,pt.size=0.1, cols=cols1, reduction = "tsne", split.by = "org", ncol = 3) + NoAxes() + NoLegend())
  # dev.off()
  # seur$oRG = ifelse(seur$CellType=="oRG","oRG","Other")
  # seur$oRG[seur$CellType=="aRG"] = "aRG"
  # seur$oRG[seur$CellType=="IP"] = "IP"
  # levels=c("aRG","Other","oRG","IP")
  # seur$oRG = factor(seur$oRG, levels = levels)
  # seur$neu = ifelse(seur$CellType=="Immature DL PN", "Immature DL PN", NA)
  # seur$neu[seur$CellType=="Newborn DL PN"] = "Newborn DL PN"
  # seur$neu[seur$CellType=="CFuPN"] = "CFuPN"
  # seur$neu[seur$CellType=="CPN"] = "CPN"
  # seur$neu[seur$CellType=="PN"] = "PN"
  # levels=c("Immature DL PN","Newborn DL PN","CFuPN","CPN","PN")
  # seur$neu = factor(seur$neu,  levels=levels)
  # tiff(filename = paste0("TSNEs/",name,"eachOrg.neurons.wKO.tsne.tiff"),res=200,width=1500,height=600*(numOrgs/3))
  # print(DimPlot(seur,pt.size=0.1, cols=c("#f768a1","#fa9fb5","#cc6677","#882255","#aa4499")[levels %in% seur$neu], reduction = "tsne", split.by = "org", ncol = 3, group.by = "neu") +
  #         NoAxes() + NoLegend())
  # dev.off()
  # seur = subset(seur, subset=treat!="ko")
  # tiff(filename = paste0("TSNEs/",name,".progen.noKO.tsne.tiff"),res=200,width=600,height=600)
  # print(DimPlot(seur,group.by="oRG",pt.size=0.1, cols=c("#41ae76","lightgray","#225522","#fdb462")[levels %in% seur$oRG], reduction = "tsne") + NoAxes() + NoLegend())
  # dev.off()
  
  #Barcharts of cell type proportions
  # counts = as.matrix(table(seur$CellType,seur$org))
  # counts = counts[,colSums(counts)!=0]
  # counts = t(t(counts)/colSums(counts))
  # counts = melt(counts,varnames = c('CellType','Org'))
  # treat = table(seur$org, seur$treat)
  # treat = treat[rowSums(treat)!=0,]
  # assigns = max.col(treat)
  # assigns = colnames(treat)[assigns]
  # names(assigns) = rownames(treat)
  # counts$treat = NA
  # for (o in 1:length(assigns)) { counts$treat[counts$Org==names(assigns)[[o]]] = assigns[[o]]}
  # counts$treat = factor(counts$treat, levels = c("wt","mut","ko"))
  # 
  # order = c("aRG","IP","Newborn DL PN","Immature DL PN",
  #           "Subcortical progenitors","Subcortical neurons","Cortical hem",
  #           "Cajal Retzius","oRG","PN","CFuPN","CPN","Immature IN", "Unknown")
  # counts$CellType = factor(counts$CellType, levels = order)
  # cols2 = cols[match(order[order %in% counts$CellType],levels)]
  # pdf(paste0("BarPlots/", name,'percent.Cells.in.Prog.wKO.pdf'),height=2, width=3*1.5)
  # print(ggplot(counts[counts$CellType %in% c("oRG","IP","aRG"),],aes(y=value, x=Org,fill=CellType,color=treat)) +
  #         geom_bar(stat="identity",show.legend=T,size=0.8)+
  #         scale_fill_manual(values=cols)+
  #         geom_col_pattern(aes(pattern_alpha=treat),pattern="stripe",pattern_color="black",pattern_size=0.3)+
  #         scale_pattern_alpha_manual(values=c(0,0,1))+
  #         scale_color_manual(values=c("darkgray","black", "black"))+
  #         facet_grid(. ~ CellType) + theme_classic() +
  #         labs(y="Percentage of Cells", x = "Organoid") +
  #         theme(axis.text.x = element_text(angle=70, hjust=1, size=14), legend.position="none"))
  # dev.off()
  
  #Distribution of apoptosis genes
  apop = read.table("../../CellAssignments/txtLists/apoptosis.txt")$V1
  seur = AddModuleScore(seur, list(apop))
  pdf(paste0("Apoptosis/", name,'ApoptosisModuleScore.pdf'),height=5, width=5)
  print(VlnPlot(seur, "Cluster1",split.by="treat", pt.size=0.1) + ggtitle("Apoptosis Module Score"))
  dev.off()
  
  #Vln plots of marker genes
  #Stacked Violin Plot from Min Tang
  # modify_vlnplot<- function(obj, 
  #                           feature, 
  #                           pt.size = 0, 
  #                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
  #                           ...) {
  #   VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
  #     xlab("") + ylab(feature) + ggtitle("") + 
  #     theme(legend.position = "none", 
  #           axis.text.x = element_blank(), 
  #           axis.ticks.x = element_blank(), 
  #           axis.title.y = element_text(size = rel(1), angle = 0), 
  #           axis.text.y = element_text(size = rel(1)), 
  #           plot.margin = plot.margin ) 
  # }
  # ## extract the max value of the y axis
  # extract_max<- function(p){
  #   ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  #   return(ceiling(ymax))
  # }
  # ## main function
  # StackedVlnPlot<- function(obj, features,
  #                           pt.size = 0, 
  #                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
  #                           ...) {
  #   plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  #   
  #   # Add back x-axis title to bottom plot. patchwork is going to support this?
  #   plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
  #     theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.ticks.x = element_line())
  #   
  #   # change the y-axis tick to only max value 
  #   ymaxs<- purrr::map_dbl(plot_list, extract_max)
  #   plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
  #                             scale_y_continuous(breaks = c(y)) + 
  #                             expand_limits(y = y))
  #   
  #   patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  # }
  # 
  # genes1 = c("SOX2","HES1","NEUROG1","NEUROD2","TBR1","PAX3","NEFL","WLS","LHX5")
  # genes3 = c("SOX2","HOPX","EOMES","NEUROD2","BCL11B","LDB2","BHLHE22","SATB2","DLX5")
  # sub = subset(seur, subset=CellType != "Unknown")
  # order = c("aRG","oRG","IP","Newborn DL PN","Immature DL PN",
  #           "Subcortical progenitors","Subcortical neurons","Cortical hem",
  #           "Cajal Retzius","PN","CFuPN","CPN","Immature IN", "Unknown")
  # sub$CellType = factor(sub$CellType, levels = order)
  # cols2 = cols[match(order[order %in% sub$CellType],levels)]
  # 
  # pdf(paste0("VlnPlots/", name,'markerGenes.pdf'),height=10, width=0.75*length(unique(sub$CellType)))
  # print(StackedVlnPlot(sub, features = if(i %in% c(1,4)){genes1} else {genes3},
  #               cols=cols2, group.by = "CellType", ) + NoAxes())
  # dev.off()
  
  # #cell type proportion function from ASD github
  # meta=seur@meta.data[,c("CellType","treat","orig.ident")] #treat contains wtVMut info, and orig.ident contains organoid identity
  # meta$orig.ident = factor(paste0(meta$treat,"_",meta$orig.ident))
  # meta = meta[meta$CellType != "Unknown",] #remove unknown cells
  # meta$CellType = factor(gsub(" ",".",meta$CellType))
  # ko = F
  # if (length(unique(meta$treat)) == 3) {
  #   metaAll = meta
  #   meta = meta[meta$treat %in% c("wt","mut"),]
  #   ko = T
  # }
  # meta$treat = factor(meta$treat, levels=c("wt","mut")) #set wt as base
  # 
  # ret = lmm_celltype(meta, celltype = "CellType")
  # ret$adj.pval = p.adjust(ret$pval, method="BH") #get adjusted p value amount of regressions we ran
  # ret$dataset = name
  # ret$CellType = rownames(ret)
  # ret$sig = ret$adj.pval<0.05
  # 
  # 
  # saveRDS(ret,file=paste0("CellTypeComp/",name,".CellComposition.LMM.NoUnk.rds"))
  # write.table(ret,file=paste0("CellTypeComp/",name,".CellComposition.LMM.NoUnk.txt"),sep="\t",quote=F,row.names=T)
  # 
  # if (ko) {
  #   meta = metaAll[metaAll$treat %in% c("wt","ko"),]
  #   meta$treat = factor(meta$treat, levels=c("wt","ko")) #set wt as base
  #   ret = lmm_celltype(meta, celltype = "CellType")
  #   ret$adj.pval = p.adjust(ret$pval, method="BH") #get adjusted p value amount of regressions we ran
  #   ret$dataset = paste0(name, ".KOvWT")
  #   ret$CellType = rownames(ret)
  #   ret$sig = ret$adj.pval<0.05
  # 
  #   saveRDS(ret,file=paste0("CellTypeComp/",name,"KOvWT.CellComposition.LMM.NoUnk.rds"))
  #   write.table(ret,file=paste0("CellTypeComp/",name,"KOvWT.CellComposition.LMM.NoUnk.txt"),sep="\t",quote=F,row.names=T)
  # }
  
  # #pathways
  # pathfiles = list.files(path="Pathways/", pattern="*.txt")
  # paths = lapply(pathfiles, function(x) {read.delim(paste0("Pathways/",x), header=F, sep="\t")[[2]]})
  # paths = lapply(paths, function(x){sapply(x, function(j) {sapply(strsplit(as.character(j), ";"),'[',1)})})
  # seur = AddModuleScore(seur, paths)
  # colnames(seur@meta.data)[(ncol(seur@meta.data)-length(paths)+1):ncol(seur@meta.data)] = sapply(pathfiles,tools::file_path_sans_ext)
  # 
  # #Plot distributions (inspired by Super Plots Lord et. al. 2020)
  # df = seur@meta.data[,c("org", "treat","CellType", sapply(pathfiles,tools::file_path_sans_ext))]
  # dfm = melt(df)
  # ReplicateAverages <- dfm %>% group_by(variable,treat,CellType,org) %>% summarize(mean(value))
  # colnames(ReplicateAverages)[[5]] = 'value'
  # pdf(paste0("Pathways/",name,".PTEN-Kegg-Pathways.pdf"), width=15, height=15)
  # print(ggplot(dfm, aes(x=treat,y=value)) +
  #   geom_violin(aes(fill=treat)) +
  #   stat_summary(fun=median, geom="crossbar")+
  #   #scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  #   geom_beeswarm(data=ReplicateAverages, size=3, cex = 4)  +
  #   theme_classic() +
  #   facet_grid(CellType ~ variable) +
  #   labs(y="Module Expression Score", x=NULL, color="Organoid Mean Value"))
  # dev.off()
  # 
  # #Calculate significance with mixed model
  # df$treat = factor(df$treat, levels = c("wt","mut","ko"))
  # ret=c()
  # for (c in unique(df$CellType)) {
  #   dfc = df[df$CellType==c,]
  #   for (k in 4:ncol(df)) {
  #     colnames(dfc)[k] = "col"
  #     res=lmer(col ~ treat + (1|org),data=dfc)
  #     res2=lmer(col ~ (1|org),data=dfc)
  #     colnames(dfc)[k] = colnames(df)[k]
  #     coef=fixef(res)["treatmut"]
  #     OR=exp(coef)
  #     CI=confint(res,parm="treatmut",method="Wald")
  #     CI_OR=exp(CI)
  #     anov=anova(res,res2)
  #     pv=anov$"Pr(>Chisq)"[2]
  #     resf=as.numeric(c(coef,OR,CI,CI_OR,pv))
  #     resf = c(resf, c)
  #     resf = c(resf, colnames(df)[k])
  #     ret=rbind(ret,resf)
  #   }
  # }
  # ret=data.frame(ret)
  # colnames(ret)=c("coef","OR","CI_coef_low","CI_coef_high",
  #                 "CI_OR_low","CI_OR_high","pval", "CellType", "Pathway")
  # ret$adj.pval = p.adjust(ret$pval, method="BH")
  # if (sum(ret$adj.pval<0.2)>0) {print(paste0("SIGNIFICANT PATHWAY:",name))}
  # write.table(ret, paste0("Pathways/", name, ".PTEN.Kegg.Pathways.pvals.csv"), sep=",")
  
  #Cacoa
  # orgs = names(table(seur$org)[table(seur$org)>0])
  # sample.groups = sapply(orgs, function(x) {substr(x,1, nchar(x)-1)})
  # seur$orgV = as.character(seur$org)
  # ko = "ko" %in% sample.groups
  # 
  # if (ko) {
  #   seur1 = subset(seur, subset=treat %in% c("wt","mut"))
  #   sample.groups1 = sample.groups[sample.groups %in% c("wt","mut")]
  # } else {
  #   seur1 = seur
  #   sample.groups1 = sample.groups
  # }
  # 
  # cao <- Cacoa$new(seur1, sample.groups=sample.groups1, cell.groups=seur1$CellType, sample.per.cell=seur1$orgV, 
  #                  ref.level="wt", target.level="mut", graph.name="RNA_snn", embedding=seur1@reductions$tsne@cell.embeddings)
  # 
  # cao$estimateCellLoadings()
  # cao$estimateExpressionShiftMagnitudes()
  # cao$plotCellLoadings(show.pvals=TRUE)
  # ggsave(paste0("cacoa/",name,".cacoa.cellLoadings.pdf"), width=8, height=8)
  # cao$plotExpressionShiftMagnitudes()
  # ggsave(paste0("cacoa/",name,".cacoa.expShifts.pdf"), width=8, height=6)
  
}

#Cell composition figures
retAll = data.frame()
for (i in c(2,3,5)) {  #1 and 4 1m, 2 3 5 3m
  ret = readRDS(paste0("CellTypeComp/",names[[i]],".CellComposition.LMM.NoUnk.rds"))
  retAll = rbind(retAll, ret)
  #if (file.exists(paste0("CellTypeComp/",names[[i]],"KOvWT.CellComposition.LMM.NoUnk.rds"))) {
  #    ret = readRDS(paste0("CellTypeComp/",names[[i]],"KOvWT.CellComposition.LMM.NoUnk.rds"))
  #    retAll = rbind(retAll, ret)
  #}
}
retAll$CellType = factor(retAll$CellType, levels = c("Cycling.Progenitors","aRG","IP","Newborn.DL.PN",
                                                   "Immature.DL.PN","Subcortical.progenitors","Subcortical.neurons","Cortical.hem","Cajal.Retzius",
                                                   "oRG","PN","CFuPN","CPN","Immature.IN")) #order of y axis
cortical = retAll[retAll$CellType %in% c("Cycling.Progenitors","aRG","IP","Newborn.DL.PN",
                                         "Immature.DL.PN","oRG","PN","CFuPN","CPN"),]
ggplot(cortical,aes(x=dataset,y=CellType,fill=coef,size=-log(adj.pval,10),color=sig)) +
  geom_point(shape=21) + scale_fill_gradient2(name="Regression Coeff.") +
  scale_color_manual(values=c("gray","black"), guide=FALSE) +
  xlab("") + ylab("") + coord_flip() + theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.box="horizontal")
ggsave("CellTypeComp/3m.Cortical.CellComposition.pdf", width = 7.8, height=2.7)


#DEG correlation with atlas genes
#Compare to atlas DEGs over time via spearman corr of signed logFC
atlas = readRDS("DEGs/DEGs_Celltypes_Over_Time_2-4m.RDS")
names(atlas)[names(atlas)=="Newborn DL PN"] = "Immature DL PN"
names(atlas)[names(atlas)=="Newborn PN"] = "Newborn DL PN"

rhos = data.frame()
for (d in c(2,3,5)) { #1 and 4 1m, 2 3 5 3m
  alldegs = readRDS(paste0("DEGs/DEGs_",names[[d]],"/DEGs.rds"))
  for (i in names(alldegs)) {
    atlasi = atlas[[i]]
    atlasi$gene = rownames(atlasi)
    degs = alldegs[[i]]
    
    deg.list = -1*log10(degs$padj)*ifelse(degs$log2FoldChange>0,1,-1)
    names(deg.list) = degs$gene
    
    atlas.list = -1*log10(atlasi$padj)*ifelse(atlasi$log2FoldChange>0,1,-1)
    names(atlas.list) = atlasi$gene
    
    atlas.list = atlas.list[names(deg.list)]
    good = !is.na(atlas.list) & !is.na(deg.list)
    atlas.list = atlas.list[good]
    deg.list = deg.list[good]
    
    rho = cor.test(deg.list, atlas.list, method = "spearman")$estimate
    rhos[names[[d]], i] = rho
    
    dat = data.frame("PTEN" = as.character(deg.list), "Atlas" = as.character(atlas.list), "Gene" = names(deg.list))
    dat[,1:2] = lapply(dat[,1:2], rank)
    #dat = dat[dat$Gene %in% c("BCL11B","FEZF2","CRYM","TLE4","LDB2","SOX5","GADD45B","LMO3","ID2"),]
    dat = dat[dat$Gene %in% c("SATB2","CUX1","CUX2","DKK3","INHBA","BTG1","EPHA3","TMTC4","NNMT","CAV1","CHN2"),]
    png(paste0("DEGs/rankcorr_CPNTargets",names[[d]],i,".png"), width=600, height=500)
    print(ggplot(dat, aes(PTEN, Atlas)) + geom_point(col="darkgray") +
     theme_classic() +
     #geom_text(y = 2450, x = 200, label=paste0("rho=", signif(rho, digits=3))))
     geom_label_repel(aes(label=Gene)))
    dev.off()
  }
}

rhos$dataset = rownames(rhos)
rhosm = melt(rhos)
ggplot(rhosm, aes(x=variable, y=dataset, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2( na.value = "lightgray") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("DEGs/DEGsvAtlas_spearmanCorr_3m.pdf", width=6, height=3)
