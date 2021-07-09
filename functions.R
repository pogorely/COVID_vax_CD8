library(data.table)
library(igraph)
library(stringr)
library(Seurat)
library(viridis)

add_citeseq_tcrseq_metadata<-function(path_to_aggr, postfix="aggr"){
  ctg<-fread(path_to_aggr)  
  bcds<-fread("libs_aggregate/filtered_feature_bc_matrix/barcodes.tsv.gz",header=F)
  ftrs<-fread("libs_aggregate/filtered_feature_bc_matrix/features.tsv.gz",header=F)
  mtr<-fread("libs_aggregate/filtered_feature_bc_matrix/matrix.mtx.gz",skip=3,fill=F)
  TCRs<-fread("libs_aggregate/filtered_contig_annotations.csv")
  
  ctg$cdr3b<-NA
  ctg$cdr3b_nt<-NA
  ctg$vb<-NA
  ctg$jb<-NA
  ctg$cdr3a<-NA
  ctg$cdr3a_nt<-NA
  ctg$va<-NA
  ctg$ja<-NA
  ctg$cdr3a2<-NA
  ctg$cdr3a2_nt<-NA
  ctg$va2<-NA
  ctg$ja2<-NA
  
  #get cell UMI and add TCRbeta and TCRalpha with largest UMI for this cell barcode
  for (i in 1:nrow(ctg))
  {
    if ((i%%100)==0)print(i)
    ctg$cdr3b[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$cdr3
    ctg$cdr3b_nt[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$cdr3_nt
    ctg$vb[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$v_gene
    ctg$jb[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$j_gene
    ctg$cdr3a[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$cdr3
    ctg$cdr3a_nt[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$cdr3_nt
    ctg$va[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$v_gene
    ctg$ja[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$j_gene
    
    ctg$cdr3a2[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$cdr3
    ctg$cdr3a2_nt[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$cdr3_nt
    ctg$va2[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$v_gene
    ctg$ja2[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$j_gene
  }
  setnames(ctg,"V1","barcode")
  #add dextramer barcodes and cell hashing statistics
  ftrs_ind<-ftrs[,((.N-29):.N),]#all ftrs_location:
  mtrdex<-mtr[V1%in%ftrs_ind,,]
  mtrdex$dex<-ftrs[mtrdex$V1,V2,]
  mtr_casted<-dcast(mtrdex, V2~ dex,value.var = "V3")
  mtr_casted$barcode<-bcds[mtr_casted$V2]$V1
  mtr_casted<-na_to0(mtr_casted[,c("barcode",c(paste0("FB_dex",1:17),paste0("FB_hash",1:10),"FB_CD3","FB_CCR7","FB_CD45RA"))])
  
  dexb<-merge(ctg,mtr_casted,by="barcode",all=T) # now contain all relevant counts. 
  #lets add batch info. 
  dexb$batch<-paste0("batch_",sapply(strsplit(dexb$barcode,split = "-",fixed=T),"[[",2))
  #fix for batch1 donor hashing issues
  dexb[batch=="batch_1",FB_hash7:=FB_CD3,]#in batch 1 7 is absent, but CD3 is there instead, so we substitute.  
  #get best_hashtag/batch columns
  hashes<-grep(pattern = "hash",names(dexb),value = T)
  dexs2<-grep(pattern = "dex",names(dexb),value = T)[1:16]
  whichx<-function(x)
  {
    if (!is.na(sum(x))){
      if (sum(x)>0){
        base::which(x)
      }
      else{return(NA)}
    }
    else{return(NA)}
  }
  dexb$hash<-(apply(prop.table(as.matrix(dexb[,..hashes]),margin = 1),MARGIN = 1,function(x)whichx(x>0.5)))
  dexb$bestdex16<-(apply(prop.table(as.matrix(dexb[,..dexs2]),margin = 1),MARGIN = 1,function(x)paste0(whichx(x>0.3),collapse="_")))
  #add donors/dex metadata.
  legend<-fread("metadata/dextramer_legend_infvax.txt")
  legend_donor<-fread("metadata/donor_legend_infvax.txt")
  
  setkey(legend_donor,batch, hash)
  setkey(x = legend,batch,hash,dex)
  
  dexb$spike<-legend[dexb[,.(batch,hash,dex=as.integer(bestdex16)),],spike,]
  dexb$dextr<-legend[dexb[,.(batch,hash,dex=as.integer(bestdex16)),],dex_id,]
  dexb$category=legend_donor[dexb[,.(batch,hash),],category,]
  dexb$donor=legend_donor[dexb[,.(batch,hash),],donor,]
  dexb<-dexb[!is.na(UMAP_1),]
  
  dexb$bestdex16_max<-(apply((as.matrix(dexb[,..dexs2])),MARGIN = 1,max))
  dexb<-dexb[!is.na(hash),,]
  dexb
}

write.tsv<-function(x,fname,rows=F){
  write.table(x,sep="\t",quote = F,row.names = rows,file = fname)
}

na_to0 <- function (x) {
  x[is.na(x)]<-0
  x
}
