library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)
library(dplyr)
library(patchwork)

seurat <- Read10X('/home/u202232805049/repro/CCST-main/transform/mouse/C3/temp')
seurat <- CreateSeuratObject(seurat)

refdata <- celldex::MouseRNAseqData()

testdata <- GetAssayData(seurat, slot="counts")

cellpred <- SingleR(test = testdata,
                    ref = refdata,
                    labels = refdata$label.main)

celltype = data.frame(ClusterID = rownames(cellpred),
                      celltype = cellpred$labels,
                      stringsAsFactors = F)

write.csv(celltype, "metadata.csv")