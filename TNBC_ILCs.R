table(sce.all$celltype)
ILCs_cell <- subset(x = sce.all, subset = celltype == "IlCs")
table(ILCs_cell$celltype)
table(ILCs_cell$Tissue)
ILCs_cell <- subset(x = ILCs_cell, subset = Tissue == "blood")
table(ILCs_cell$celltype)

as.vector(as.data.frame(table(ILCs_cell$orig.ident))[,1])[-1]



{
  a.list <- SplitObject(ILCs_cell, split.by = "orig.ident")
  file_name <- ls(a.list)
  for (i in 1:length(ls(a.list))) {
    print(file_name[i])
    gene.exp <- as.matrix(eval(parse(text =paste('a.list','$',file_name[i],"@assays[['RNA']]@counts",sep = ""))))
    dataan <- as(as.matrix(gene.exp), "dgCMatrix")
    assign(file_name[i],CreateSeuratObject(
      counts = dataan, project =names(a.list[i]),
      meta.data =eval(parse(text =paste('a.list',"$",file_name[i],"@meta.data",sep = "")))))
  }
  #ILCs_cell <-merge(Post_P008_b,y=c(Post_P013_b,Post_P018_b,Post_P018_t,Post_P020_b,Post_P022_b,Post_P022_t,Post_P023_b,Post_P023_t,Post_P024_b,Post_P025_b,Post_P025_t,Pre_P008_b,Pre_P013_b,Pre_P018_b,Pre_P018_t,Pre_P020_b,Pre_P02>

  ILCs_cell <-merge(eval(parse(text =as.vector(as.data.frame(table(ILCs_cell$orig.ident))[,1])[1] )),
                    y=c(Post_P013_b,Post_P018_b,Post_P020_b,Post_P022_b,Post_P023_b,Post_P024_b,Post_P025_b,Pre_P008_b,Pre_P013_b,Pre_P018_b,
                        Pre_P020_b,Pre_P022_b,Pre_P023_b,Pre_P024_b,Pre_P025_b),project = "CCA")

  table(ILCs_cell$celltype)
  table(ILCs_cell$orig.ident)

}

ILCs_cell <- NormalizeData(ILCs_cell, normalization.method = "LogNormalize", scale.factor = 10000)
ILCs_cell <- FindVariableFeatures(ILCs_cell, selection.method = "vst", nfeatures = 2000)
ILCs_cell <- ScaleData(ILCs_cell, verbose = T)
ILCs_cell <- RunPCA(ILCs_cell, features = VariableFeatures(object = ILCs_cell))
DimPlot(ILCs_cell, reduction = "pca")
ILCs_cell <- RunTSNE(ILCs_cell, reduction = "pca", dims = 1:30)
ILCs_cell <- RunUMAP(ILCs_cell, reduction = "pca", dims = 1:30)
ILCs_cell <- FindNeighbors(ILCs_cell, reduction = "pca", dims = 1:30)
#ILCs_cell <- FindClusters(ILCs_cell, resolution =0.6)###yuan 0.05
ILCs_cell <- FindClusters(ILCs_cell, resolution =0.12)###yuan 0.05
#ILCs_cell <- FindClusters(ILCs_cell, resolution =5)###yuan 0.05

DimPlot(ILCs_cell, reduction = "tsne",label = T)

