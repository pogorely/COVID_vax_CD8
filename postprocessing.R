source("functions.R")

##part 1: seurat analysis

#Load aggregated 10X data and create Seurat object
cov1.data <- Read10X(data.dir = "libs_aggregate/filtered_feature_bc_matrix/")
cov1 <- CreateSeuratObject(counts = cov1.data$`Gene Expression`)
#Remove low quality cells
cov1[["percent.mt"]] <- PercentageFeatureSet(cov1, pattern = "^MT-")
plot1 <- FeatureScatter(cov1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cov1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cov1 <- subset(cov1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
cov1<- NormalizeData(cov1, normalization.method = "LogNormalize", scale.factor = 10000)

#Find 2000 variable genes excluding TCR and BCR genes
markers.remove<-fread("markers.remove.csv", stringsAsFactors = F)
cov.m_filt<-cov1
cov.m_filt <- FindVariableFeatures(cov.m_filt, selection.method = "vst", nfeatures = 2000+135)
VariableFeatures(cov.m_filt) <- VariableFeatures(cov.m_filt)[!(VariableFeatures(cov.m_filt) %in% markers.remove$V2)]

#Regress cell cycle variability
all.genes <- rownames(cov.m_filt)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cov.m_filt <- CellCycleScoring(cov.m_filt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cov.m_filt<- ScaleData(cov.m_filt, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
cov.m_filt <- RunPCA(cov.m_filt, features = VariableFeatures(object = cov.m_filt) , verbose = FALSE)
cov.m_filt <- FindNeighbors(cov.m_filt, dims = 1:15)
cov.m_filt <- FindClusters(cov.m_filt, resolution = 0.5, verbose = FALSE)
cov.m_filt <- RunUMAP(cov.m_filt, dims = 1:15)
cov.m_filt <- RunTSNE(cov.m_filt, dims = 1:15)
cov.m_filt$replicate<-(sapply(X = strsplit(colnames(cov.m_filt), split = "-"), FUN = "[", 2))
saveRDS(cov.m_filt, file = "10x_aggr_step1.rds")

umapCoord <- as.data.frame(Embeddings(object = cov.m_filt[["umap"]]))
tsneCoord<-as.data.frame(Embeddings(object = cov.m_filt[["tsne"]]))
tmp<-cbind(umapCoord, batch=cov.m_filt$replicate, cluster=cov.m_filt$seurat_clusters, tsneCoord, nCount=cov.m_filt$nCount_RNA, nFeature=cov.m_filt$nFeature_RNA)
write.csv(tmp, file="10x_aggr_step1.csv")

#Remove B cells and monocytes
aggr<-readRDS("10x_aggr_step1.rds")

Idents(aggr) <- "seurat_clusters"
cl.markers <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

aggr.s<-subset(aggr, idents = c(8,10, 11), invert=T)

#find carrier cell barcodes
inf_analysis<-add_citeseq_tcrseq_metadata("10x_aggr_step1.csv")
writeLines(inf_analysis[donor=="1814",barcode,],con = "barcodes_1814.tsv")

#Remove donor #1814 (carrier cells) and all remaining non cd8 cells 
Idents(aggr.s) <- "orig.ident"
d1814<-fread("barcodes_1814.tsv", stringsAsFactors = F, header = F)

aggr.f<-aggr.s[,!colnames(aggr.s) %in% d1814$V1]
Idents(aggr.f) <- "orig.ident"
cd8<-subset(x = aggr.f, subset = (CD8A>0.5|CD8B>0.5), invert=F)
cd8 <- FindNeighbors(cd8, dims = 1:15)
cd8 <- FindClusters(cd8, resolution = 0.5, verbose = FALSE)
cd8 <- RunUMAP(cd8, dims = 1:15)
cd8 <- RunTSNE(cd8, dims = 1:15)
DimPlot(cd8, label=T)
cd8.markers <- FindAllMarkers(cd8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

umapCoord <- as.data.frame(Embeddings(object = cd8[["umap"]]))
tsneCoord<-as.data.frame(Embeddings(object = cd8[["tsne"]]))

cd8$CD8A<-FetchData(cd8, vars = "CD8A")
cd8$CD8B<-FetchData(cd8, vars = "CD8B")
cd8$IL7R<-FetchData(cd8, vars = "IL7R")
cd8$s100a4<-FetchData(cd8, vars = "S100A4")
cd8$cd27<-FetchData(cd8, vars = "CD27")
cd8$cd4<-FetchData(cd8, vars = "CD4")
cd8$CCR7<-FetchData(cd8, vars = "CCR7")
cd8$TCF7<-FetchData(cd8, vars = "TCF7")
cd8$GZMB<-FetchData(cd8, vars = "GZMB")
cd8$GZMH<-FetchData(cd8, vars = "GZMH")
cd8$SELL<-FetchData(cd8, vars = "SELL")
cd8$GZMK<-FetchData(cd8, vars = "GZMK")
tmp<-cbind(umapCoord, batch=cd8$replicate, cluster=cd8$seurat_clusters, tsneCoord, nCount=cd8$nCount_RNA, nFeature=cd8$nFeature_RNA,
           CD8A=cd8$CD8A, CD8B=cd8$CD8B, IL7R=cd8$IL7R, s100a4=cd8$s100a4, cd27=cd8$cd27, cd4=cd8$cd4,
           ccr7=cd8$CCR7, tcf7=cd8$TCF7, gzmb=cd8$GZMB, gzmh=cd8$GZMH, sell=cd8$SELL, gzmk=cd8$GZMK)
write.csv(tmp, file="cd8_only.csv")
saveRDS(cd8, file = "cd8_only.rds")

##Part 2: add Citeseq, TCR, dextramer assignments for cd8 cells identified on previous step. 

inf_analysis8<-add_citeseq_tcrseq_metadata("cd8_only.csv")
legend<-fread("metadata/dextramer_legend_infvax.txt")
setkey(legend,batch,hash,dex)
#fix for double dextramer assignments
tmp<-character()
for (i in 1:nrow(inf_analysis8[grepl("_",bestdex16),]))
{
  ind<-inf_analysis8[grepl("_",bestdex16),,][i,.(batch,hash,dex=as.integer(unlist(strsplit(bestdex16,split = "_")))),]
  tmp[i]<-(legend[ind,paste(dex_id,collapse="|"),])
  inf_analysis8[grepl("_",bestdex16),,][i,dextr:=legend[ind,paste(dex_id,collapse="|"),],]
} 
inf_analysis8[grepl("_",bestdex16),dextr:=tmp,]

#filter low confidence dextramer assignments
inf_analysis8[bestdex16_max<4,dextr:=NA,]
inf_analysis8[bestdex16_max<4,spike:=NA,]

inf_analysis8[!is.na(cdr3b_nt)&!is.na(cdr3a_nt),clonotype_id:=.GRP,.(batch,donor,cdr3b_nt,cdr3a_nt)]
inf_analysis8[!is.na(cdr3b_nt),clonotype_id_beta:=.GRP,.(batch,donor,cdr3b_nt)]

#dextr_clone: imputation of dextramer assignment for a clonotype.
inf_analysis8[,dextr_clone:=dextr,.(clonotype_id,donor)]
inf_analysis8[!is.na(clonotype_id),dextr_clone:=names(sort(-table(dextr)))[1],.(clonotype_id,donor)]
#spike_clone: imputation of spike/non-spike assignment for a clonotype.
inf_analysis8$spike_clone<-FALSE
inf_analysis8[dextr_clone%in%c("A24_NYN","B15_NAF","B15_NQF","B44_AEA","B44_AEV","A24_QYI","A02_YLQ","A01_LTD","B15_NAF|B15_NQF","B15_NQF|B15_NAF"),spike_clone:=TRUE,]
inf_analysis8[is.na(dextr_clone),spike_clone:=NA,]

#epitope column: fixed NQF/NAF double assignments, filtered other umbiguous assignments
inf_analysis8[dextr=="B15_NAF|B15_NQF",dextr:="B15_NQF|B15_NAF",]
inf_analysis8[dextr%in%c("B15_NAF|B15_NQF","B15_NQF|B15_NAF"),spike:=TRUE,]
inf_analysis8[,epitope:=dextr_clone,]
inf_analysis8[epitope%in%c("B15_NAF","B15_NQF","B15_NAF|B15_NQF","B15_NQF|B15_NAF"),epitope:="B15_NQF_NAF",]
inf_analysis8[grepl("|",epitope,fixed=T),epitope:=NA,]
inf_analysis8[is.na(epitope),spike_clone:=NA,]

#write out final table
write.tsv(inf_analysis8,fname="cd8_only_dextr.tsv")