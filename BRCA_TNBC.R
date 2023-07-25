rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
set.seed(999)
options(stringsAsFactors = F)
setwd("~/")
#####read raw data from 10X in GSE169246
# Breast cancer TNBC treated
 data1 <- Matrix::readMM("./GSE169246/matrix.mtx.gz")
  head(data1)
  c1 <- fread("./GSE169246/barcodes.tsv.gz",header=F,data.table=F)
  r1 <- fread("./GSE169246/features.tsv.gz",header=F,data.table=F)
  head(c1)
  head(r1)
  colnames(data1) <- c1$V1
  rownames(data1) <- r1$V1
  meta_data <- fread("./GSE169246/GSE169246.metadata.csv.gz",header = T)
  head(meta_data)
  rownames(meta_data) <- meta_data$`Cell barcode`
  data_treated <- CreateSeuratObject(counts = data1,meta.data = meta_data)
  data_treated1 <- subset(x = data_treated, subset = Origin == "t")
  table(data_treated1$Sample)
  ####### tumor tisse
  name_tumor <- c("Post_P002_t","Post_P003_t","Post_P005_t","Post_P012_t","Post_P013_t","Post_P017_t","Post_P018_t","Post_P019_t","Post_P022_t","Post_P023_t","Post_P025_t","Pre_P002_t","Pre_P004_t","Pre_P005_t","Pre_P007_t","Pre_P0>
  for (i in 1:length(name_tumor)) {
    print(name_tumor[i])
    assign(name_tumor[i],subset(x = data_treated, subset = Sample == name_tumor[i]))
  }
  data_treated1 <-  merge(Post_P002_t,y=c(Post_P003_t,Post_P005_t,Post_P012_t,Post_P013_t,Post_P017_t,Post_P018_t,Post_P019_t,Post_P022_t,Post_P023_t,Post_P025_t,Pre_P002_t,Pre_P004_t,Pre_P005_t,Pre_P007_t,Pre_P012_t,Pre_P013_t,Pre>
  temp1 <- subset(data_treated1,subset=Tissue=="breast")
  data_treated1<- temp1
  table(data_treated1@meta.data$orig.ident)

  table(data_treated1@meta.data$Sample)

  data_treated1@meta.data$orig.ident <- data_treated1$Sample

  table(data_treated1@meta.data$orig.ident)
  table(data_treated1@meta.data$Sample)

  data_treated1[["percent.Hb"]] <- PercentageFeatureSet(data_treated1, pattern = "^HB")
  data_treated1[["percent.mt"]] <- PercentageFeatureSet(data_treated1, pattern = "^MT-")
  VlnPlot(data_treated1, features = c( "percent.mt","percent.Hb","nCount_RNA", "nFeature_RNA"),ncol =4,pt.size =0.01)
  VlnPlot(data_treated1, features = c( "percent.mt","percent.Hb","nCount_RNA", "nFeature_RNA"),ncol =4,pt.size =0.01,group.by = "Sample")
  data_treated_tumor <- data_treated1
  save(data_treated_tumor,file = "./rep/temp_data/TNBC.tumor.raw.Rdata")


  data_treated1 <- subset(x = data_treated, subset = Origin == "b")
  table(data_treated1$Sample)


