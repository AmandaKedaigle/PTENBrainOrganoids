#Code for uploading data to single cell portal

library(Seurat)
library(Matrix)

setwd("~/Documents/DorsalKadoshima/PTEN_FinalObjects/")

objs = c("Mito_1m_rep1_celltypes.rds","Mito_3m_rep1_celltypes.rds","Mito_3m_rep2_celltypes.rds",
         "PGP1_1m_rep1_celltypes.rds","PGP1_3m_celltypes.rds")
names = c("Mito_1m_rep1","Mito_3m_rep1","Mito_3m_rep2","PGP1_1m_rep1","PGP1_3m_rep1")


for (i in 1:length(names)) {
seur = readRDS(objs[[i]])
sampleName = names[[i]]
savedir = "SCP"
#dir.create(savedir)

#Add metadata by SCP conventions
seur$biosample_id = paste0(seur$treat, "_", seur$orig.ident)

seur$donor_id = ifelse(startsWith(sampleName, "Mito"), "Mito 210", "PGP1")

seur$sex = "male"
#Hues is female

seur$species = "NCBITaxon:9606"
seur$species__ontology_label = "Homo sapiens"
seur$disease = "MONDO:0005260"
seur$disease__ontology_label = "autism"
seur$organ = "UBERON:0000955"
seur$organ__ontology_label = "brain"
seur$library_preparation_protocol = "EFO:0008995"
seur$library_preparation_protocol__ontology_label = "10X sequencing"
seur$biosample_type = "DerivedType_Organoid"

rownames(seur@meta.data) = paste0(rownames(seur@meta.data), sampleName)

#meta - study wide vs cluster based ----
swkeep=c("biosample_id", "donor_id", "nFeature_RNA","nCount_RNA",
       "sex", "species", "species__ontology_label", "disease", "disease__ontology_label", "organ",
       "organ__ontology_label", "library_preparation_protocol", "library_preparation_protocol__ontology_label", "biosample_type")
line1="NAME\tbiosample_id\tdonor_id\tnGene\tnUMI\tsex\tspecies\tspecies__ontology_label\tdisease\tdisease__ontology_label\torgan\torgan__ontology_label\tlibrary_preparation_protocol\tlibrary_preparation_protocol__ontology_label\tbiosample_type"
line2="TYPE\tgroup\tgroup\tnumeric\tnumeric\tgroup\tgroup\tgroup\tgroup\tgroup\tgroup\tgroup\tgroup\tgroup\tgroup"

meta=seur@meta.data[,swkeep]
saveMeta=paste(savedir,"/meta_", sampleName, ".txt",sep="")
write(c(line1,line2),saveMeta,sep="\n")
write.table(meta,saveMeta,sep="\t",col.names=F,quote=F,append=T)

saveMeta=paste(savedir,"/meta_all.txt",sep="")
#write(c(line1,line2),saveMeta,sep="\n")
write.table(meta,saveMeta,sep="\t",col.names=F,quote=F,append=T)

cbkeep = c("CellType","treat","percent.mito","percent.ribo","seurat_clusters")
cbLabels = c("CellType","treatment","percent.mito", "percent.ribo", "Cluster")
cbTypes = c("group","group","numeric", "numeric", "group")

umap=as.data.frame(seur@reductions$tsne@cell.embeddings)
umap[,cbLabels] = seur@meta.data[rownames(umap),cbkeep]
rownames(umap) = paste0(rownames(umap), sampleName)
saveUMAP=paste(savedir,"/umap_",sampleName,".txt",sep="")
line1=paste(c("NAME\tX\tY", cbLabels), collapse="\t")
line2=paste(c("TYPE\tnumeric\tnumeric",cbTypes), collapse="\t")
write(c(line1,line2),saveUMAP,sep="\n")
write.table(umap,saveUMAP,sep="\t",col.names=F,quote=F,append=T)

saveData=paste(savedir,"/expression_counts_",sampleName,".txt",sep="")
writeMM(GetAssayData(seur, slot = "counts"),saveData)
genes = rownames(GetAssayData(seur, slot = "counts"))
barcodes = colnames(GetAssayData(seur, slot = "counts"))
barcodes = paste0(barcodes, sampleName)
write.table(genes,paste0(savedir,"/expression_",sampleName,"_genes.txt"),col.names=F,row.names=F,quote=F)
write.table(barcodes, paste0(savedir,"/expression_",sampleName,"_barcodes.txt"),col.names=F,row.names=F,quote=F)

saveData=paste(savedir,"/expression_data_",sampleName,".txt",sep="")
writeMM(GetAssayData(seur, slot = "data"),saveData)
genes = rownames(GetAssayData(seur, slot = "data"))
barcodes = colnames(GetAssayData(seur, slot = "data"))
barcodes = paste0(barcodes, sampleName)
write.table(genes,paste0(savedir,"/expression_data_",sampleName,"_genes.txt"),col.names=F,row.names=F,quote=F)
write.table(barcodes, paste0(savedir,"/expression_data_",sampleName,"_barcodes.txt"),col.names=F,row.names=F,quote=F)

}
