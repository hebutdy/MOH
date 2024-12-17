library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)

pbmc <- Read10X('/home/u202232805049/repro/CCST-main/transform/mouse/C3/temp') #这里填写解压后的目录
pbmc <- CreateSeuratObject(pbmc) #转化为Seurat对象
Pattern = '^mt-' #线粒体基因的名字，根据实际去匹配，人是^MT-开头，小鼠^mt-，大鼠^Mt-
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = Pattern)
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt <  15) #进行细胞过滤，这里是要求每个细胞里面至少有200个的基因，但不大于8000个基因，线粒体基因含量不大于15%。这里的过滤就根据实际来，没有绝对的标准。
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) #进行log归一化
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #寻找高变基因
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) #根据基因数进行数据的缩放
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) #进行PCA降纬，特征是我们前面找到的高变基因
pbmc <- FindNeighbors(pbmc, dims = 1:30) #寻找临近的特征点
pbmc <- FindClusters(object = pbmc, verbose = T, resolution = 0.7) #分群，关键在于resolution调整分辨率，在0～1之间，约大就分越多分群

# 查看
ls("package:celldex")
ref <- celldex::MouseRNAseqData()

norm_count = GetAssayData(pbmc, slot="data")
pbmc_anno <- SingleR(test = norm_count, ref = ref, labels = ref$label.main, cluster = pbmc$seurat_clusters)  #这里我们选择按照分群去注释

celltype = data.frame(ClusterID = rownames(pbmc_anno),
                      celltype = pbmc_anno$labels,
                      stringsAsFactors = F)