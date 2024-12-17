
library(MAESTRO)
library(SingleR)

library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(scater)
library(bluster)
library(GenomeInfoDb)
library(IRanges)
library(rtracklayer)
library(devtools)
library(reticulate)
library(Seurat)
library(Signac)
library(cluster)
library(igraph)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tinytex)
library(tidyverse)
library(parallel)
library(data.table)
library(Hmisc)
library(patchwork)
library(GenomicRanges)
library(factoextra)

## Setup path and environment
source("/root/autodl-tmp/MOH/scRNA_scATAC_mouse.r")
# set python environment
Sys.setenv(RETICULATE_PYTHON = "/root/miniconda3/envs/R/bin/python")
use_python("/root/miniconda3/envs/R/bin/python")
py_config()
lisa_path <- "/root/autodl-tmp/MOHr/lisa_output/"
jaspar_path <- "/root/autodl-tmp/MOH/lisa_path/"

lymph_obj <- ReadData(h5Path = "/root/autodl-tmp/MOH/dataset/mouse3/L2", data_type = "scRNA_scATAC", min_cell = 0.05, dataFormat = "h5")



ATAC_gene_peak <- CalGenePeakScore(peak_count_matrix = lymph_obj@assays$ATAC@counts,organism = "GRCm38")

GAS_obj <- calculate_GAS_v1(ATAC_gene_peak = ATAC_gene_peak, obj = lymph_obj, method = "wnn")
GAS <- GAS_obj[[1]]
lymph_obj <- GAS_obj[[2]]


cell_gene_matrix <- as.matrix(GAS)


cell_similarity_matrix <- dist(t(cell_gene_matrix))
cell_similarity_matrix <- 1 - cell_similarity_matrix
cell_embedding <- cmdscale(cell_similarity_matrix, k = 128)
HGT_result <- run_HGT(GAS = as.matrix(GAS),result_dir='/root/autodl-tmp/MOH/result/mouse3', data_type='scRNA_scATAC', envPath=NULL, lr=0.2, epoch=30, n_hid=128, n_heads=16)
cell_hgt_matrix <- HGT_result[['cell_hgt_matrix']]
merged_embedding <- (cell_embedding + cell_hgt_matrix) / 2
cell_hgt_matrix <- merged_embedding



rownames(cell_hgt_matrix) <- colnames(GAS)

lymph_obj <- lymph_obj[, colnames(GAS)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]

## ST
CCST_embed_matrix <- read.csv("/root/autodl-tmp/MOH/chapter4/L2/Combined_representations.csv", header = TRUE)

CCST_embed_matrix <- as.matrix(CCST_embed_matrix)


cosine_similarity_matrix <- matrix(0, nrow = nrow(CCST_embed_matrix), ncol = nrow(cell_hgt_matrix))

max_cols <- max(ncol(cell_hgt_matrix), ncol(CCST_embed_matrix))
# [1] 3563  128
cell_hgt_matrix <- cbind(cell_hgt_matrix, matrix(0, nrow = nrow(cell_hgt_matrix), ncol = max_cols - ncol(cell_hgt_matrix)))
# [1] 2432  128
CCST_embed_matrix <- cbind(CCST_embed_matrix, matrix(0, nrow = nrow(CCST_embed_matrix), ncol = max_cols - ncol(CCST_embed_matrix)))


new_CCST_embed_matrix <- CCST_embed_matrix
new_deep_embed_matrix <- cell_hgt_matrix

dim(new_CCST_embed_matrix)
dim(new_deep_embed_matrix)

# 计算余弦相似度矩阵
cosine_similarity_matrix <- proxy::simil(x = cell_hgt_matrix, y = CCST_embed_matrix, method = "cosine")


for (i in 1:nrow(cosine_similarity_matrix)) {
    max_similarity_index <- which.max(cosine_similarity_matrix[i, ])

    if (cosine_similarity_matrix[i, max_similarity_index] > 0.3) {
      # 归一化并加权平均
      normalized_embed1 <- cell_hgt_matrix[i, ] / sum(cell_hgt_matrix[i, ])
      normalized_embed2 <- CCST_embed_matrix[max_similarity_index, ] / sum(CCST_embed_matrix[max_similarity_index, ])
      weighted_average <- (normalized_embed1 + normalized_embed2) / 2
#       weighted_average <- (cell_hgt_matrix[i, ] + CCST_embed_matrix[max_similarity_index, ]) / 2
      # 反归一化并更新新生成的嵌入矩阵
      new_deep_embed_matrix[i, ] <- weighted_average * sum(cell_hgt_matrix[i, ])
      new_CCST_embed_matrix[max_similarity_index, ] <- weighted_average * sum(CCST_embed_matrix[max_similarity_index, ])
#       new_deep_embed_matrix[i, ] <- weighted_average
#         new_CCST_embed_matrix[max_similarity_index, ] <- weighted_average
    }
}

for (j in 1:ncol(cosine_similarity_matrix)) {
    max_similarity_index <- which.max(cosine_similarity_matrix[, j])

    if (cosine_similarity_matrix[max_similarity_index, j] > 0.3) {
      # 归一化并加权平均
      normalized_embed1 <- cell_hgt_matrix[max_similarity_index, ] / sum(cell_hgt_matrix[max_similarity_index, ])
      normalized_embed2 <- CCST_embed_matrix[j, ] / sum(CCST_embed_matrix[j, ])
      weighted_average <- (normalized_embed1 + normalized_embed2) / 2
#       weighted_average <- (cell_hgt_matrix[max_similarity_index, ] + CCST_embed_matrix[j, ]) / 2

      # 反归一化并更新新生成的嵌入矩阵
      new_deep_embed_matrix[max_similarity_index, ] <- weighted_average * sum(cell_hgt_matrix[max_similarity_index, ])
      new_CCST_embed_matrix[j, ] <- weighted_average * sum(CCST_embed_matrix[j, ])
#       new_deep_embed_matrix[max_similarity_index, ] <- weighted_average
#         new_CCST_embed_matrix[j, ] <- weighted_average
    }
 }



HGT_embedding <-
  CreateDimReducObject(embeddings = new_deep_embed_matrix,
                       key = "HGT_",
                       assay = "RNA")


lymph_obj@reductions[['HGT']] <- HGT_embedding
lymph_obj <-
  FindVariableFeatures(lymph_obj, selection.method = "vst", nfeatures = 2000)



lymph_obj <-
  RunUMAP(
    lymph_obj,
    reduction = 'HGT',
    dims = 1:ncol(new_deep_embed_matrix),
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_"
  )
lymph_obj <-
  FindNeighbors(lymph_obj,
                reduction = "HGT",
                dims = 1:ncol(new_deep_embed_matrix))

lymph_obj <- FindClusters(lymph_obj,
                          reduction.type="HGT",
                          algorithm= 1,
                          resolution = 0.15,
                          verbose = T,
                          graph.name = "RNA_snn")

current_clusters <- as.numeric(as.character(lymph_obj$seurat_clusters))


graph.out <- as.factor(lymph_obj$seurat_clusters)

DefaultAssay(lymph_obj) <- "RNA"


##后续正常的下游实验流程
