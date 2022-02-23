source("functions.R")

##part 1: seurat analysis

#Load aggregated 10X data and create Seurat object
cov1.data <- Read10X(data.dir = "libs_aggregate_rev/filtered_feature_bc_matrix/")
cov1 <- CreateSeuratObject(counts = cov1.data$`Gene Expression`)
#Remove low quality cells
cov1[["percent.mt"]] <- PercentageFeatureSet(cov1, pattern = "^MT-")
plot1 <- FeatureScatter(cov1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cov1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
cov1 <- subset(cov1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
cov1<- NormalizeData(cov1, normalization.method = "LogNormalize", scale.factor = 10000)

#Find 2000 variable genes excluding TCR and BCR genes
markers.remove<-fread("markers.remove.csv", stringsAsFactors = F)
cov.m_filt<-cov1
cov.m_filt <- FindVariableFeatures(cov.m_filt, selection.method = "vst", nfeatures = 2000+210)
VariableFeatures(cov.m_filt) <- VariableFeatures(cov.m_filt)[!(VariableFeatures(cov.m_filt) %in% markers.remove$V2)]

#Regress cell cycle variability
all.genes <- rownames(cov.m_filt)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cov.m_filt <- CellCycleScoring(cov.m_filt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cov.m_filt<- ScaleData(cov.m_filt, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
cov.m_filt <- RunPCA(cov.m_filt, features = VariableFeatures(object = cov.m_filt) , verbose = FALSE)
cov.m_filt <- FindNeighbors(cov.m_filt, dims = 1:20)
cov.m_filt <- FindClusters(cov.m_filt, resolution = 0.5, verbose = FALSE)
cov.m_filt <- RunUMAP(cov.m_filt, dims = 1:20)
cov.m_filt <- RunTSNE(cov.m_filt, dims = 1:20)
cov.m_filt$replicate<-(sapply(X = strsplit(colnames(cov.m_filt), split = "-"), FUN = "[", 2))
saveRDS(cov.m_filt, file = "10x_aggr_step1.rds")

umapCoord <- as.data.frame(Embeddings(object = cov.m_filt[["umap"]]))
tsneCoord<-as.data.frame(Embeddings(object = cov.m_filt[["tsne"]]))
tmp<-cbind(umapCoord, batch=cov.m_filt$replicate, cluster=cov.m_filt$seurat_clusters, tsneCoord, nCount=cov.m_filt$nCount_RNA, nFeature=cov.m_filt$nFeature_RNA)
write.csv(tmp, file="10x_aggr_step1.csv")

#Remove B cells and monocytes
aggr<-readRDS("10x_aggr_step1.rds")

Idents(aggr) <- "seurat_clusters"
DimPlot(aggr, label=T)
cl.markers <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

aggr.s<-subset(aggr, idents = c(11,12,13,14), invert=T)

#find carrier cell barcodes, or load them from file to save time.
#inf_analysis<-add_citeseq_tcrseq_metadata_revision("10x_aggr_step1.csv")
#writeLines(inf_analysis[donor=="1814",barcode,],con = "barcodes_1814.tsv") this is empty now

#Remove donor #1814 (carrier cells) and all remaining non cd8 cells 
Idents(aggr.s) <- "orig.ident"
d1814<-fread("barcodes_1814.tsv", stringsAsFactors = F, header = F)

aggr.f<-aggr.s[,!colnames(aggr.s) %in% d1814$V1]

Idents(aggr.f) <- "orig.ident"
cd8b<-subset(x = aggr.f, subset = (CD8A>0.5|CD8B>0.5), invert=F)
cd8b <- FindNeighbors(cd8b, dims = 1:20)
cd8b <- FindClusters(cd8b, resolution = 0.5, verbose = FALSE)
cd8b <- RunUMAP(cd8b, dims = 1:20)
cd8b <- RunTSNE(cd8b, dims = 1:20)
DimPlot(cd8b, label=T)
cd8b.markers <- FindAllMarkers(cd8b, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cd8<-subset(cd8b, idents = c(11), invert=T)
cd8 <- FindNeighbors(cd8, dims = 1:20)
cd8 <- FindClusters(cd8, resolution = 0.5, verbose = FALSE)
cd8 <- RunUMAP(cd8, dims = 1:20)
cd8 <- RunTSNE(cd8, dims = 1:20)
DimPlot(cd8, label=T)
cd8.markers <- FindAllMarkers(cd8, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cd8.markers, file="all.cd8.markers.csv")

umapCoord <- as.data.frame(Embeddings(object = cd8[["umap"]]))
tsneCoord<-as.data.frame(Embeddings(object = cd8[["tsne"]]))

DefaultAssay(cd8) <- "RNA"
cd8$CD8A<-FetchData(cd8, vars = "CD8A")
cd8$CD8B<-FetchData(cd8, vars = "CD8B")
cd8$IL7R<-FetchData(cd8, vars = "IL7R")
cd8$GNLY<-FetchData(cd8, vars = "GNLY")
cd8$cd27<-FetchData(cd8, vars = "CD27")
cd8$cd4<-FetchData(cd8, vars = "CD4")
cd8$CCR7<-FetchData(cd8, vars = "CCR7")
cd8$TCF7<-FetchData(cd8, vars = "TCF7")
cd8$GZMB<-FetchData(cd8, vars = "GZMB")
cd8$GZMH<-FetchData(cd8, vars = "GZMH")
cd8$SELL<-FetchData(cd8, vars = "SELL")
cd8$GZMK<-FetchData(cd8, vars = "GZMK")
cd8$NKG7<-FetchData(cd8, vars = "NKG7")
cd8$TIGIT<-FetchData(cd8, vars = "TIGIT")
cd8$PDCD1<-FetchData(cd8, vars = "PDCD1")
cd8$CTLA4<-FetchData(cd8, vars = "CTLA4")
cd8$LAG3<-FetchData(cd8, vars = "LAG3")
cd8$TOX<-FetchData(cd8, vars = "TOX")
cd8$EOMES<-FetchData(cd8, vars = "EOMES")
cd8$TBX21<-FetchData(cd8, vars = "TBX21")


tmp<-cbind(umapCoord, batch=cd8$replicate, cluster=cd8$seurat_clusters, tsneCoord, nCount=cd8$nCount_RNA, nFeature=cd8$nFeature_RNA,
           CD8A=cd8$CD8A, CD8B=cd8$CD8B, IL7R=cd8$IL7R, GNLY=cd8$GNLY, cd27=cd8$cd27, cd4=cd8$cd4,
           ccr7=cd8$CCR7, tcf7=cd8$TCF7, gzmb=cd8$GZMB, gzmh=cd8$GZMH, sell=cd8$SELL, gzmk=cd8$GZMK,
           nkg7=cd8$NKG7, TIGIT=cd8$TIGIT, PDCD1=cd8$PDCD1, ctla4=cd8$CTLA4, lag3=cd8$LAG3,
           TOX=cd8$TOX, EOMES=cd8$EOMES, TBX21=cd8$TBX21)
write.csv(tmp, file="cd8_only_genes.csv")
saveRDS(cd8, file = "cd8_only.rds")
#write.csv(cd8,file="cd8_only.csv")

##part two: add CITEseq, TCRseq, dextramer info

inf_analysis_rev<-add_citeseq_tcrseq_metadata_revision("cd8_only_genes.csv")
legend<-fread("metadata_rev/meta_revision.tsv")
legend$batch<-paste0("batch_",legend$batch)
setkey(legend,batch,hash,dex)

#write separate dextramers separated by pipe.

tmp<-character()
tmpd<-character()
tmps<-character()
for (i in 1:nrow(inf_analysis_rev[grepl("FB.*FB",bestdex16),])) 
{   
  ind<-inf_analysis_rev[grepl("FB.*FB",bestdex16),,][i,.(batch,hash,dex=as.integer(unlist(strsplit(gsub("FB_dex","",bestdex16),split = "_")))),]
  tmp[i]<-(legend[ind,paste(dex_id,collapse="|"),])
  tmpd[i]<-(legend[ind,paste(unique(donor),collapse="|"),])
  tmps[i]<-(legend[ind,paste(unique(sample_category),collapse="|"),])
} 
inf_analysis_rev[grepl("FB.*FB",bestdex16),dextr:=tmp,]
inf_analysis_rev[grepl("FB.*FB",bestdex16),donor:=tmpd,]
inf_analysis_rev[grepl("FB.*FB",bestdex16),sample_category:=tmps,]

#filter low confidence dextramer assignments
inf_analysis_rev[bestdex16_max<4,dextr:=NA,]
inf_analysis_rev[bestdex16_max<4,spike:=NA,]
inf_analysis_rev[bestdex16_max<4,donor:=NA,]
inf_analysis_rev[bestdex16_max<4,sample_category:=NA,]

inf_analysis_rev[!is.na(cdr3b_nt)&!is.na(cdr3a_nt),clonotype_id:=.GRP,.(batch,cdr3b_nt,cdr3a_nt)] 
inf_analysis_rev[!is.na(cdr3b_nt),clonotype_id_beta:=.GRP,.(batch,cdr3b_nt)] 


#Alternative imputation code:
inf_analysis_rev[,dextr_clone:=dextr,] #worst case: dextr clone is same as dextr. 
inf_analysis_rev[,donor_clone:=donor,]
inf_analysis_rev[,sample_category_clone:=sample_category,]
for (clone in unique(inf_analysis_rev$clonotype_id))
{
  if (inf_analysis_rev[(clonotype_id==clone)&!is.na(clonotype_id),sum(table(dextr)==max(table(dextr)))==1,])#notie condition
  {
    inf_analysis_rev[(clonotype_id==clone)&!is.na(clonotype_id),dextr_clone:=names(sort(-table(dextr)))[1],]
    inf_analysis_rev[(clonotype_id==clone)&!is.na(clonotype_id),donor_clone:=names(sort(-table(donor)))[1],]
    inf_analysis_rev[(clonotype_id==clone)&!is.na(clonotype_id),sample_category_clone:=names(sort(-table(sample_category)))[1],]
  }
}

#spike_clone: imputation of spike/non-spike assignment for a clonotype.
inf_analysis_rev$spike_clone<-FALSE
inf_analysis_rev[dextr_clone%in%c("A24_NYN","B15_NQK_A","B15_NQK_Q","B44_AEA","B44_AEV","A24_QYI","A02_YLQ","A01_LTD","B15_NQK_A|B15_NQK_Q","B15_NQK_Q|B15_NQK_A"),spike_clone:=TRUE,]
inf_analysis_rev[is.na(dextr_clone),spike_clone:=NA,]

#epitope column: fixed NQF/NAF double assignments, filtered other ambiguous assignments
inf_analysis_rev[dextr=="B15_NQK_A|B15_NQK_Q",dextr:="B15_NQK_Q|B15_NQK_A",]
inf_analysis_rev[dextr%in%c("B15_NQK_A|B15_NQK_Q","B15_NQK_Q|B15_NQK_A"),spike:=TRUE,]
inf_analysis_rev[,epitope:=dextr_clone,]
inf_analysis_rev[epitope%in%c("B15_NQK_A","B15_NQK_Q","B15_NQK_A|B15_NQK_Q","B15_NQK_Q|B15_NQK_A"),epitope:="B15_NQK",]
inf_analysis_rev[grepl("|",epitope,fixed=T),epitope:=NA,] 
inf_analysis_rev[is.na(epitope),spike_clone:=NA,]
inf_analysis_rev[grepl("|",donor,fixed=T),donor:=NA,]
inf_analysis_rev[grepl("|",donor_clone,fixed=T),donor_clone:=NA,]

legend_donors<-fread("metadata_rev/donor_barcode_legend.txt")
legend_donors[,donor_category:=paste0(new_name,"_",new_category),]
setkey(legend_donors,donor_category)

inf_analysis_rev$internal_sample_id<-legend_donors[paste0(inf_analysis_rev$donor_clone,"_",inf_analysis_rev$sample_category_clone),ID,]
inf_analysis_rev$internal_HLA_group<-legend_donors[paste0(inf_analysis_rev$donor_clone,"_",inf_analysis_rev$sample_category_clone),HLA,]
inf_analysis_rev$donor_category<-legend_donors[paste0(inf_analysis_rev$donor_clone,"_",inf_analysis_rev$sample_category_clone),Category,]
inf_analysis_rev<-inf_analysis_rev[!is.na(donor_clone)&!is.na(dextr_clone)&donor_clone!="NA"&!is.na(internal_HLA_group),,]

write.tsv(inf_analysis_rev,fname="cd8_only_dextr_rev_clean.tsv")


