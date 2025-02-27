#--------------load packages-----------------------------
library(Seurat)
library(tidyverse)
library(Cairo)
library(sc.utils)
library(clusterProfiler)
library(cowplot)
library(dittoSeq)
library(ggpubr)
library(SeuratWrappers)
library(ggsci)
library(velocyto.R)
library(scCustomize)

#--------------0 Quality Control--------------------------------------------
#' read loom and GEO expression mtx file, combine and merge into one seurat object
#' @param loom_path loom file from RNA velocyto
#' @param exprSet_path exprSet data from GEO
#' @param project_name project name in seurat
#' @author Yi Xiong
generate_seurat <- function(loom_path = NULL, exprSet_path = NULL, project_name = NULL,min.cells = 0, min.features = 0, tag = "Med2-Peri:"){
  ldat <- ReadVelocity(file = loom_path)
  ##modify the cell name
  # tmp = ldat$spliced
  # cells = colnames(tmp)
  # make sure the consistency of barcode between velocyto and cell ranger
  ldat <- lapply(ldat,function(x) {
    colnames(x) <-  gsub(tag,"",colnames(x))
    colnames(x) <- gsub("x","-1",colnames(x))
    x
  })
  
  Seurat_obj <- as.Seurat(x = ldat)#
  Seurat.data <- Read10X(data.dir = exprSet_path)
  
  Seurat_GEO <- CreateSeuratObject(Seurat.data, min.cells = min.cells, min.features = min.features, project = project_name)
  
  Seurat_obj[['RNA']] <- CreateAssayObject(counts = Seurat_GEO@assays$RNA@counts)
  orig.ident <- Idents(object = Seurat_GEO)
  Idents(object = Seurat_obj) <- orig.ident
  orig.ident_chr = as.character(orig.ident);names(orig.ident_chr) = colnames(Seurat_obj)
  Seurat_obj <- AddMetaData(
    object = Seurat_obj,
    metadata = orig.ident_chr,
    col.name = 'sample'
  )
  # s1 = colnames(E10)
  # s2 = colnames(E10_GEO)
  # length(intersect(s1,s2))
  return(Seurat_obj)
}

library(data.table)
metadata <- fread("./metadata.csv",data.table = F)

seu_list <- map(c("WSY_S_HSCR_X",
                  "WSY_S_HSCR_K",
                  "LMH_S_HSCR_X",
                  "LMH_S_HSCR_K",
                  "JBX_S_HSCR_X",
                  "JBX_S_HSCR_K",
                  "GZY_S_HSCR_X",
                  "GZY_S_HSCR_K",
                  "LZL_K",
                  "LZL_X",
                  "LZY_K",
                  "LZY_X",
                  "YRZ_HD_K",
                  "YRZ_HD_X"),function(i){
                    seu <- generate_seurat(loom_path = paste0("./cellranger/",i,"/03.cellranger/",i,"/velocyto/",i,".loom"), 
                                           exprSet_path = paste0("./cellranger/",i,"/03.cellranger/",i,"/outs/filtered_feature_bc_matrix/"), 
                                           project_name = i,
                                           tag = paste0(i,":"))
                  })

grep( "^MT-", rownames(seu_list[[1]]), value = T)
#seu_list[[1]]
#Quality Control for each matrix
for (i in 1:length(seu_list)){
  DefaultAssay(seu_list[[i]]) <- "RNA"
  seu_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seu_list[[i]], pattern = "^MT-")
  meta <- seu_list[[i]]@meta.data
  
  meta <- meta %>% mutate(cell.id = rownames(.)) %>% left_join(metadata, by = "sample")
  rownames(meta) = meta$cell.id
  seu_list[[i]]@meta.data <-  meta
}

for (i in 1:length(seu_list)){
  print(head(seu_list[[i]][[]]))
}

#Visualise QC metrics as violin plot
violinplot <- list()
for (i in 1:length(seu_list)){
  violinplot[[i]] <- VlnPlot(seu_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.05)
}

cowplot::plot_grid(plotlist = violinplot)

#batch QC
for (i in 1:length(seu_list)){
  seu_list[[i]] <- subset(seu_list[[i]],subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 50)
}

#print filter result
for (i in 1:length(seu_list)) {
  print(paste0(names(seu_list)[i]," ",seu_list[[i]]@assays$RNA@counts@Dim))
}


for (i in 1:length(seu_list)) {
  print(i)
  seu_list[[i]] <- NormalizeData(seu_list[[i]], verbose = FALSE)
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

###Find the anchors
anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:30, 
                                  reduction = "cca")
# Run the standard workflow for visualization and clustering
#Integrate Datasets into a new matrix
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)

ElbowPlot(integrated, ndims = 30, reduction = "pca")


# remove doublet-------------------------------
scrublet = reticulate::import("scrublet")
scr <- scrublet$Scrublet(counts_matrix = t(as.matrix(integrated@assays$RNA@counts)),
                         expected_doublet_rate = 0.12)
scr_out = scr$scrub_doublets()

integrated$doulet_scores = scr_out[[1]]
integrated$predicted_doulets = scr_out[[2]]

####1 cluster #############
ndim = 20;resolution = 1
integrated <- FindNeighbors(integrated, dims = 1:ndim,
                            features = NULL, k.param = 20,
                            compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                            annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                            force.recalc = FALSE)
integrated <- FindClusters(integrated, resolution = resolution, algorithm = 1)

#Run UMAP
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:ndim, 
                      n.neighbors = 30L, n.components = 2L,
                      min.dist = 0.5, spread = 1)
head(integrated[[]])
DimPlot(integrated, reduction = "umap", label = TRUE, label.size = 3,
        repel = TRUE) + NoLegend()

DimPlot(integrated, reduction = "umap",group.by = "Individual_Name", label = TRUE, label.size = 3,
        repel = TRUE) + NoLegend()



########### 2 cell type prediction--------------------------------
library(SeuratDisk)

Convert("./raw_data/Full_obj_log_counts_soupx_v2.h5ad", ".h5seurat")
# This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory

# This .d5seurat object can then be read in manually
seuratObject <- LoadH5Seurat("example_dir/example_ad.h5Seurat")

options(reticulate.conda_binary = "/home/byron/miniconda3/envs/scrna/")


library(reticulate)
# reticulate::use_miniconda(condaenv = "/home/byron/miniconda3/envs/scrna/", conda = "auto", required = NULL)
options(reticulate.conda_binary = "/home/byron/miniconda3/bin/conda", SCP_env_name = "SCP_env")
library(SCP)
library(Seurat)
sc <- import("scanpy")
adata <- sc$read_h5ad("./raw_data/Full_obj_log_counts_soupx_v2.h5ad")
srt <- adata_to_srt(adata)
srt

DefaultAssay(integrated) <- "RNA"
integrated <- RunKNNPredict(
  srt_query = integrated, srt_ref = srt,
  ref_group = "category", filter_lowfreq = 20
)
head(integrated[[]])
integrated$pred.cell.type.major = integrated$KNNPredict_classification
CellDimPlot(srt = integrated, group.by = "pred.cell.type.major", reduction = "umap", label = TRUE)

integrated <- RunKNNPredict(
  srt_query = integrated, srt_ref = srt,
  ref_group = "Integrated_05", filter_lowfreq = 20
)

integrated$pred.cell.type.minor = integrated$KNNPredict_classification

integrated$pred.cell.type.minor_v2 = dplyr::recode(integrated$pred.cell.type.minor, 
                                                   'CLC+ Mast cell' = "Mast cells",
                                                   'Mast cell' = "Mast cells",
                                                   'BEST2+ Goblet cell' = "Goblet cells", 
                                                   'BEST4+ epithelial' = "BEST4+ epithelial",
                                                   Colonocyte = "Enterocytes",
                                                   'Contractile pericyte (PLN+)' = "Pericytes", 
                                                   'Pericyte' = "Pericytes",
                                                   'Distal progenitor' = "Intestinal Stem cells",
                                                   'Goblet cell' = "Goblet cells", 
                                                   'Microfold cell' = "Intestinal Stem cells",
                                                   'Progenitor (NEUROG3+)' = "Enteroendocrine cells", 
                                                   'EC cells (NPW+)' = "EC", 
                                                   'Proximal progenitor' = "Intestinal Stem cells",
                                                   'Stem cells' = "Intestinal Stem cells",
                                                   'Stem cells' = "Intestinal Stem cells",
                                                   'Stromal 2 (NPY+)' = "Stromal cells",
                                                   'arterial capillary' = "EC (capillary)",
                                                   'Branch A4 (IN)' = "Neurons",
                                                   'angiogenic pericyte' = "Pericytes",
                                                   cDC1 = "DCs",
                                                   cDC2 = "DCs",
                                                   pDC = "DCs",
                                                   CLP = "B cells",
                                                   'Cycling B cell' = "B cells (cycling)",
                                                   'cycling EC' = "EC (cycling)",
                                                   "Cycling plasma cell" = "Plasma cell (cycling)",
                                                   "Differentiating glia" = "glia",
                                                   "Adult Glia" = "glia",
                                                   "CLP" = "B cells",
                                                   'DZ GC cell' = "B cells",
                                                   'LEC3 (ADGRG3+)' = "Lymphatic endothelium",
                                                   'FCRL4+ Memory B' = "B cells",
                                                   'Fetal arterial EC' = "EC (arterial)",
                                                   'Fetal venous EC' = "EC (venous)",
                                                   'GC B cell' = "B cells",
                                                   'Germ' = "B cells",
                                                   'L cells (PYY+)' = "L cells",
                                                   'IgM plasma cell' = "Plasma cell (IgM)",
                                                   'Immature pericyte' = "Pericytes",
                                                   'LYVE1+ Macrophage' = "Macrophages",
                                                   'LZ GC cell' = "B cells",
                                                   'ILC3' = "Lymphoid cells",
                                                   'ILCP' = "Lymphoid cells",
                                                   'LTi-like NCR- ILC3' = "Lymphoid cells",
                                                   'LTi-like NCR+ ILC3' = "Lymphoid cells",
                                                   'Lymphoid DC' = "Lymphoid cells",
                                                   'Mature arterial EC' = "EC (arterial)",
                                                   'Mature venous EC' = "EC (venous)",
                                                   'Memory B' = "B cells (memory)",
                                                   'Mesothelium' = "Mesenchymal cells",
                                                   'MMP9+ Inflammatory macrophage' = "Macrophages",
                                                   'Monocytes' = "Monocytes",
                                                   'myofibroblast' = "Stromal cells",
                                                   'Naive B' = "B cells",
                                                   'Pro-B' = "B cells",
                                                   'gdT' = "T cells (gd)",
                                                   'TRGV4 gdT' = "T cells (gd)",
                                                   'TRGV2 gdT' = "T cells (gd)",
                                                   'TRDV2/TRGV9 gdT' = "T cells (gd)",
                                                   'Treg' = "T cells",
                                                   'Glia 3 (BCAN+)' = "glia",
                                                   'Glia 1 (DHH+)' = "glia",
                                                   'SELL+ CD8 T' = "T cells (CD4)",
                                                   'SELL+ CD4 T' = "T cells (CD4)",
                                                   'Activated CD4 T' = "T cells (CD4)",
                                                   'Activated CD8 T' = "T cells (CD8)",
                                                   'CD8 Tmem' = "T cells (CD8)",
                                                   'CX3CR1+ CD8 Tmem' = "T cells (CD8)",
                                                   'NK cell' = "NK cells",
                                                   'NK T cell' = "NK cells",
                                                   'MAIT cell' = "T cells",
                                                   'T reticular' = "T cells",
                                                   'Tfh' = "T cells",
                                                   'Th17' = "T cells",
                                                   'STAT1+ Naive B' = "B cells",
                                                   'Branch B3 (IPAN)' = "Neurons",
                                                   'Mesoderm 2 (ZEB2+)' = "Mesenchymal cells",
                                                   'Mesothelium' = "Mesenchymal cells",
                                                   'Stromal 1 (ADAMDEC1+)' = "Stromal cells",
                                                   'Stromal 1 (CCL11+)' = "Stromal cells",
                                                   'Stromal 2 (NPY+)' = "Stromal cells",
                                                   'Stromal 3 (C7+)' = "Stromal cells",
                                                   'Stromal 3 (KCNN3+)' = "Stromal cells",
                                                   'mLN Stroma (FMO2+)' = "Stromal cells",
                                                   'Transitional Stromal 3 (C3+)' = "Stromal cells",
                                                   'SMC (PLPP2+)' = "Stromal cells",
                                                   'venous capillary' = "EC (capillary)")

integrated$pred.cell.type.major_v2 = dplyr::recode(integrated$pred.cell.type.minor, 
                                                   'CLC+ Mast cell' = "T cells",
                                                   'Mast cell' = "T cells",
                                                   'BEST2+ Goblet cell' = "Epithelial", 
                                                   'BEST4+ epithelial' = "Epithelial",
                                                   Colonocyte = "Epithelial",
                                                   'Contractile pericyte (PLN+)' = "Mesenchymal", 
                                                   'Pericyte' = "Mesenchymal",
                                                   'Distal progenitor' = "Epithelial",
                                                   'Goblet cell' = "Epithelial", 
                                                   'Microfold cell' = "Epithelial",
                                                   'Progenitor (NEUROG3+)' = "Epithelial", 
                                                   'EC cells (NPW+)' = "Endothelial", 
                                                   'Proximal progenitor' = "Epithelial",
                                                   'Stem cells' = "Epithelial",
                                                   'Stromal 2 (NPY+)' = "Mesenchymal",
                                                   'arterial capillary' = "Endothelial",
                                                   'Branch A4 (IN)' = "Neuronal",
                                                   'angiogenic pericyte' = "Mesenchymal",
                                                   cDC1 = "Myeloid",
                                                   cDC2 = "Myeloid",
                                                   pDC = "Myeloid",
                                                   CLP = "B cells",
                                                   'Cycling B cell' = "B cells",
                                                   'cycling EC' = "Endothelial",
                                                   "Cycling plasma cell" = "B cells",
                                                   "Differentiating glia" = "Neuronal",
                                                   "Adult Glia" = "Neuronal",
                                                   'DZ GC cell' = "B cells",
                                                   'LEC3 (ADGRG3+)' = "Endothelial",
                                                   'FCRL4+ Memory B' = "B cells",
                                                   'Fetal arterial EC' = "Endothelial",
                                                   'Fetal venous EC' = "Endothelial",
                                                   'GC B cell' = "B cells",
                                                   'Germ' = "B cells",
                                                   'L cells (PYY+)' = "Epithelial",
                                                   'IgM plasma cell' = "B cells",
                                                   'Immature pericyte' = "Mesenchymal",
                                                   'LYVE1+ Macrophage' = "Myeloid",
                                                   'LZ GC cell' = "B cells",
                                                   'ILC3' = "T cells",
                                                   'ILCP' = "T cells",
                                                   'LTi-like NCR- ILC3' = "T cells",
                                                   'LTi-like NCR+ ILC3' = "T cells",
                                                   'Lymphoid DC' = "Myeloid",
                                                   'Mature arterial EC' = "Endothelial",
                                                   'Mature venous EC' = "Endothelial",
                                                   'Memory B' = "B cells",
                                                   'Mesothelium' = "Mesenchymal",
                                                   'MMP9+ Inflammatory macrophage' = "Myeloid",
                                                   'Monocytes' = "Myeloid",
                                                   'myofibroblast' = "Mesenchymal",
                                                   'Naive B' = "B cells",
                                                   'Pro-B' = "B cells",
                                                   'gdT' = "T cells",
                                                   'TRGV4 gdT' = "T cells",
                                                   'TRGV2 gdT' = "T cells",
                                                   'TRDV2/TRGV9 gdT' = "T cells",
                                                   'Treg' = "T cells",
                                                   'Glia 3 (BCAN+)' = "Neuronal",
                                                   'Glia 1 (DHH+)' = "Neuronal",
                                                   'SELL+ CD8 T' = "T cells",
                                                   'SELL+ CD4 T' = "T cells",
                                                   'Activated CD4 T' = "T cells",
                                                   'Activated CD8 T' = "T cells",
                                                   'CD8 Tmem' = "T cells",
                                                   'CX3CR1+ CD8 Tmem' = "T cells",
                                                   'NK cell' = "T cells",
                                                   'NK T cell' = "T cells",
                                                   'MAIT cell' = "T cells",
                                                   'T reticular' = "T cells",
                                                   'Tfh' = "T cells",
                                                   'Th17' = "T cells",
                                                   'STAT1+ Naive B' = "B cells",
                                                   'Branch B3 (IPAN)' = "Neuronal",
                                                   'Mesoderm 2 (ZEB2+)' = "Mesenchymal",
                                                   'Mesothelium' = "Mesenchymal",
                                                   'Stromal 1 (ADAMDEC1+)' = "Mesenchymal",
                                                   'Stromal 1 (CCL11+)' = "Mesenchymal",
                                                   'Stromal 2 (NPY+)' = "Mesenchymal",
                                                   'Stromal 3 (C7+)' = "Mesenchymal",
                                                   'Stromal 3 (KCNN3+)' = "Mesenchymal",
                                                   'mLN Stroma (FMO2+)' = "Mesenchymal",
                                                   'Transitional Stromal 3 (C3+)' = "Mesenchymal",
                                                   'SMC (PLPP2+)' = "Mesenchymal",
                                                   'venous capillary' = "Endothelial",
                                                   'Activated T' = "T cells",
                                                   'Macrophages' = "Myeloid",
                                                   'Enteroendocrine cells' = "Epithelial",
                                                   TA = "Epithelial",
                                                   Tuft = "Epithelial")


integrated$Sample_Condition = paste0(integrated$Condition,"_",integrated$Region)

DimPlot(integrated, reduction = "umap",group.by = "pred.cell.type.minor", label = TRUE, label.size = 3,
        repel = TRUE) + NoLegend()

DimPlot(integrated, reduction = "umap",group.by = "pred.cell.type.minor_v2", label = TRUE, label.size = 3,
        repel = TRUE) + NoLegend()

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
integrated <- CellCycleScoring(integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


#### 3 Cell type refinement---------------------------------------------------------
head(integrated[[]])

## refine major cell types
table(integrated$seurat_clusters,integrated$pred.cell.type.major)
Idents(integrated) = "integrated_snn_res.1"
integrated <- RenameIdents(object = integrated, 
                           "0" = "B/plasma cell, red blood cells",
                           "1" = "Endothelial",
                           "2" = "Mesenchymal",
                           "3" = "Mesenchymal",
                           "4" = "Endothelial",
                           "5" = "Mesenchymal",
                           "6" = "T cells",
                           "7" = "B/plasma cell, red blood cells",
                           "8" = "Epithelial",
                           "9" = "T cells",
                           "10" = "Epithelial",
                           "11" = "Epithelial",
                           "12" = "Mesenchymal",
                           "13" = "Epithelial",
                           "14" = "Mesenchymal",
                           "15" = "B/plasma cell, red blood cells",
                           "16" = "Myeloid",
                           "17" = "Mesenchymal",
                           "18" = "B/plasma cell, red blood cells",
                           "19" = "Epithelial",
                           "20" = "Epithelial",
                           "21" = "Epithelial",
                           "22" = "Neuronal",
                           "23" = "B/plasma cell, red blood cells",
                           "24" = "Epithelial",
                           "25" = "Epithelial",
                           "26" = "Mesenchymal",
                           "27" = "Mesenchymal",
                           "28" = "Epithelial",
                           "29" = "Epithelial", 
                           "30" = "Myeloid",
                           "31" = "T cells",
                           "32" = "Epithelial",
                           "33" = "Epithelial",
                           "34" = "B/plasma cell, red blood cells",
                           "35" = "Epithelial/Neuronal",
                           "36" = "Mesenchymal",
                           "37" = "Endothelial",
                           "38" = "Mesenchymal",
                           "39" = "B/plasma cell, red blood cells",
                           "40" = "T cells",
                           "41" = "Mesenchymal",
                           "42" = "Epithelial")

integrated$cell.type.major_v3 <- Idents(integrated)

CellDimPlot(integrated,
            group.by = "cell.type.major_v3", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

## 3.1 refine Mesenchymal ----------------------
mes <- subset(integrated,subset = cell.type.major_v3 == "Mesenchymal")
DefaultAssay(mes) <- "RNA"
ndim = 20;resolution = 0.8

mes <- ScaleData(mes)
mes <- RunPCA(mes, npcs = 30, verbose = FALSE)

mes <- FindNeighbors(mes, dims = 1:ndim,
                     features = NULL, k.param = 20,
                     compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                     annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                     force.recalc = FALSE)
mes <- FindClusters(mes, resolution = 1.2, algorithm = 1)

#Run UMAP
mes <- RunUMAP(mes, reduction = "pca", dims = 1:ndim, 
               n.neighbors = 30L, n.components = 2L,
               min.dist = 1.2, spread = 1)

CellDimPlot(mes,
            group.by = "seurat_clusters", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(mes,
            group.by = "RNA_snn_res.1.2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(mes,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(mes$RNA_snn_res.1.2,mes$pred.cell.type.minor))

Idents(mes) = "RNA_snn_res.1.2"
mes <- RenameIdents(object = mes, 
                    "0" = "Contractile pericyte (PLN+)",
                    "1" = "Stromal 1 (CCL11+/ADAMDEC1+)",
                    "2" = "Stromal 3 (C7+)",
                    "3" = "Stromal 1 (CCL11+/ADAMDEC1+)",
                    "4" = "Stromal 2 (NPY+)",
                    "5" = "Stromal 3 (KCNN3+)",
                    "6" = "myofibroblast",
                    "7" = "Pericyte",
                    "8" = "Stromal 1 (CCL11+)",
                    "9" = "Transitional Stromal 3 (C3+)",
                    "10" = "Transitional Stromal 3 (C3+)",
                    "11" = "angiogenic pericyte",
                    "12" = "Stromal 2 (NPY+)",
                    "13" = "Contractile pericyte (PLN+)",
                    "14" = "Stromal 2 (NPY+)",
                    "15" = "Transitional Stromal 3 (C3+)",
                    "16" = "Stromal 1 (CCL11+)",
                    "17" = "SMC (PLPP2+)",
                    "18" = "Contractile pericyte (PLN+)",
                    "19" = "Transitional Stromal 3 (C3+)",
                    "20" = "Stromal 4 (MMP1+)",
                    "21" = "T reticular",
                    "22" = "Pericytes/Stromal Transition",
                    "23" = "Pericytes/Stromal Transition",
                    "24" = "Cycling stroma",
                    "25" = "Mesenchymal-like EC",
                    "26" = "Contractile pericyte (PLN+)")

mes$cell.type.minor_v1 <- Idents(mes)
CellDimPlot(mes,
            group.by = "cell.type.minor_v1", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)


## 3.2 refine B plasmal RBC ----------------------
b_plasm <- subset(integrated,subset = cell.type.major_v3 == "B/plasma cell, red blood cells")
DefaultAssay(b_plasm) <- "RNA"

ndim = 20;resolution = 1.2

b_plasm <- ScaleData(b_plasm)
b_plasm <- RunPCA(b_plasm, npcs = 30, verbose = FALSE)

b_plasm <- FindNeighbors(b_plasm, dims = 1:ndim,
                         features = NULL, k.param = 20,
                         compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                         annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                         force.recalc = FALSE)
b_plasm <- FindClusters(b_plasm, resolution = resolution, algorithm = 1)

#Run UMAP
b_plasm <- RunUMAP(b_plasm, reduction = "pca", dims = 1:ndim, 
                   n.neighbors = 30L, n.components = 2L,
                   min.dist = 1.2, spread = 1)

CellDimPlot(b_plasm,
            group.by = "RNA_snn_res.1.2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(b_plasm,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(b_plasm$RNA_snn_res.1.2,b_plasm$pred.cell.type.minor))

Idents(b_plasm) = "RNA_snn_res.1.2"
b_plasm <- RenameIdents(object = b_plasm, 
                        "0" = "Naive B",
                        "1" = "Naive B",
                        "2" = "Naive B",
                        "3" = "Naive B",
                        "4" = "LZ GC cell",
                        "5" = "Cycling B cell",
                        "6" = "Naive B",
                        "7" = "Naive B",
                        "8" = "IgA plasma cell",
                        "9" = "IgA plasma cell",
                        "10" = "DZ GC cell",
                        "11" = "Naive B",
                        "12" = "Cycling B cell",
                        "13" = "Cycling B cell",
                        "14" = "Memory B",
                        "15" = "Naive B",
                        "16" = "Naive B",
                        "17" = "IgG plasma cell",
                        "18" = "Cycling B cell",
                        "19" = "Undefined",
                        "20" = "IgA plasma cell")

## 3.3 refine Endothelial ----------------------
endo <- subset(integrated,subset = cell.type.major_v3 == "Endothelial")
DefaultAssay(endo) <- "RNA"

ndim = 20;resolution = 1.2

endo <- ScaleData(endo)
endo <- RunPCA(endo, npcs = 30, verbose = FALSE)

endo <- FindNeighbors(endo, dims = 1:ndim,
                      features = NULL, k.param = 20,
                      compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                      annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                      force.recalc = FALSE)
endo <- FindClusters(endo, resolution = resolution, algorithm = 1)

#Run UMAP
endo <- RunUMAP(endo, reduction = "pca", dims = 1:ndim, 
                n.neighbors = 30L, n.components = 2L,
                min.dist = 1.2, spread = 1)

CellDimPlot(endo,
            group.by = "RNA_snn_res.1.2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(endo,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(endo$RNA_snn_res.1.2,endo$pred.cell.type.minor))

Idents(endo) = "RNA_snn_res.1.2"
endo <- RenameIdents(object = endo, 
                     "0" = "arterial capillary EC",
                     "1" = "Mature venous EC",
                     "2" = "Mature venous EC",
                     "3" = "Mature arterial EC",
                     "4" = "Mature venous EC",
                     "5" = "arterial capillary EC",
                     "6" = "Mature venous EC",
                     "7" = "arterial capillary EC",
                     "8" = "arterial capillary EC",
                     "9" = "Mature venous EC",
                     "10" = "venous capillary EC",
                     "11" = "venous capillary EC",
                     "12" = "Mature venous EC",
                     "13" = "arterial capillary EC",
                     "14" = "arterial capillary EC",
                     "15" = "venous capillary EC",
                     "16" = "Fetal venous EC",
                     "17" = "Mature venous EC",
                     "18" = "arterial capillary EC",
                     "19" = "LEC",
                     "20" = "Fetal arterial EC",
                     "21" = "arterial capillary EC")

endo$cell.type.minor_v1 <- Idents(endo)

CellDimPlot(endo,
            group.by = "cell.type.minor_v1", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

## 3.4 refine T cells ----------------------
t_cells <- subset(integrated,subset = cell.type.major_v3 == "T cells")
DefaultAssay(t_cells) <- "RNA"

ndim = 20;resolution = 1.2

t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells, npcs = 30, verbose = FALSE)

t_cells <- FindNeighbors(t_cells, dims = 1:ndim,
                         features = NULL, k.param = 20,
                         compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                         annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                         force.recalc = FALSE)
t_cells <- FindClusters(t_cells, resolution = resolution, algorithm = 1)

#Run UMAP
t_cells <- RunUMAP(t_cells, reduction = "pca", dims = 1:ndim, 
                   n.neighbors = 30L, n.components = 2L,
                   min.dist = 1.2, spread = 1)


a <- as.data.frame(table(t_cells$RNA_snn_res.1.2,t_cells$pred.cell.type.minor))

Idents(t_cells) = "RNA_snn_res.1.2"
t_cells <- RenameIdents(object = t_cells, 
                        "0" = "SELL+ CD4 T",
                        "1" = "gdT",
                        "2" = "NK cell",
                        "3" = "Th17",
                        "4" = "gdT",
                        "5" = "Tfh",
                        "6" = "CD8 T",
                        "7" = "SELL+ CD4 T",
                        "8" = "ILC3",
                        "9" = "Treg",
                        "10" = "CD4 T",
                        "11" = "Th17",
                        "12" = "CD8 Tmem",
                        "13" = "gdT",
                        "14" = "Undefined",
                        "15" = "Undefined",
                        "16" = "Undefined")

t_cells$cell.type.minor_v1 <- Idents(t_cells)

## 3.5 refine Myeloids ----------------------
Myeloid <- subset(integrated,subset = cell.type.major_v3 == "Myeloid")
DefaultAssay(Myeloid) <- "RNA"

ndim = 20;resolution = 1.2

Myeloid <- ScaleData(Myeloid)
Myeloid <- RunPCA(Myeloid, npcs = 30, verbose = FALSE)

Myeloid <- FindNeighbors(Myeloid, dims = 1:ndim,
                         features = NULL, k.param = 20,
                         compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                         annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                         force.recalc = FALSE)
Myeloid <- FindClusters(Myeloid, resolution = resolution, algorithm = 1)

#Run UMAP
Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:ndim, 
                   n.neighbors = 30L, n.components = 2L,
                   min.dist = 1.2, spread = 1)

CellDimPlot(Myeloid,
            group.by = "RNA_snn_res.1.2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(Myeloid,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(Myeloid$RNA_snn_res.1.2,Myeloid$pred.cell.type.minor))

Idents(Myeloid) = "RNA_snn_res.1.2"
Myeloid <- RenameIdents(object = Myeloid, 
                        "0" = "Macrophages",
                        "1" = "Mast cell",
                        "2" = "Mast cell",
                        "3" = "LYVE1+ Macrophages",
                        "4" = "cDC2",
                        "5" = "Macrophages",
                        "6" = "Monocytes",
                        "7" = "Monocytes",
                        "8" = "cDC1",
                        "9" = "cDC2",
                        "10" = "Mast cell",
                        "11" = "Macrophages",
                        "12" = "pDC",
                        "13" = "MMP9+ Inflammatory macrophage",
                        "14" = "Undefined",
                        "15" = "Mast cell",
                        "16" = "Macrophages",
                        "17" = "Monocytes/LYVE1+ Macrophages",
                        "18" = "Lymphoid DC",
                        "19" = "Monocytes",
                        "20" = "Undefined")

Myeloid$cell.type.minor_v1 <- Idents(Myeloid)


## 3.6 refine Neuronal ----------------------
Neuronal <- subset(integrated,subset = cell.type.major_v3 == "Neuronal")
DefaultAssay(Neuronal) <- "RNA"

ndim = 20;resolution = 1.2

Neuronal <- ScaleData(Neuronal)
Neuronal <- RunPCA(Neuronal, npcs = 30, verbose = FALSE)

Neuronal <- FindNeighbors(Neuronal, dims = 1:ndim,
                          features = NULL, k.param = 20,
                          compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                          annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                          force.recalc = FALSE)
Neuronal <- FindClusters(Neuronal, resolution = resolution, algorithm = 1)

#Run UMAP
Neuronal <- RunUMAP(Neuronal, reduction = "pca", dims = 1:ndim, 
                    n.neighbors = 30L, n.components = 2L,
                    min.dist = 1.2, spread = 1)

CellDimPlot(Neuronal,
            group.by = "RNA_snn_res.1.2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(Neuronal,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(Neuronal$RNA_snn_res.1.2,Neuronal$pred.cell.type.minor))

Idents(Neuronal) = "RNA_snn_res.1.2"
Neuronal <- RenameIdents(object = Neuronal, 
                         "0" = "Glia",
                         "1" = "Glia",
                         "2" = "Glia 1 (DHH+)",
                         "3" = "Glia",
                         "4" = "Glia 1 (DHH+)",
                         "5" = "Glia 3 (BCAN+)",
                         "6" = "Glia",
                         "7" = "Glia",
                         "8" = "Glia 3 (BCAN+)",
                         "9" = "Glia 1 (DHH+)",
                         "10" = "Glia 2 (ELN+)",
                         "11" = "Glia 2 (ELN+)",
                         "12" = "Glia 2 (ELN+)",
                         "13" = "Neuroblast",
                         "14" = "Neuroblast")

Neuronal$cell.type.minor_v1 <- Idents(Neuronal)


## 3.7 refine Epi Neuronal ----------------------
Epi_neu <- subset(integrated,subset = cell.type.major_v3 == "Epithelial/Neuronal")
DefaultAssay(Epi_neu) <- "RNA"

ndim = 20;resolution = 1.2

Epi_neu <- ScaleData(Epi_neu)
Epi_neu <- RunPCA(Epi_neu, npcs = 30, verbose = FALSE)

Epi_neu <- FindNeighbors(Epi_neu, dims = 1:ndim,
                         features = NULL, k.param = 20,
                         compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                         annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                         force.recalc = FALSE)
Epi_neu <- FindClusters(Epi_neu, resolution = resolution, algorithm = 1)

#Run UMAP
Epi_neu <- RunUMAP(Epi_neu, reduction = "pca", dims = 1:ndim, 
                   n.neighbors = 30L, n.components = 2L,
                   min.dist = 1.2, spread = 1)

CellDimPlot(Epi_neu,
            group.by = "RNA_snn_res.1.2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(Epi_neu,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(Epi_neu$RNA_snn_res.1.2,Epi_neu$pred.cell.type.minor))

Idents(Epi_neu) = "RNA_snn_res.1.2"
Epi_neu <- RenameIdents(object = Epi_neu, 
                        "0" = "Progenitor (NEUROG3+)",
                        "1" = "Progenitor (NEUROG3+)",
                        "2" = "EC cells (NPW+)",
                        "3" = "EC cells (NPW+)",
                        "4" = "L cells (PYY+)",
                        "5" = "Branch A4 (IN)",
                        "6" = "Colonocyte",
                        "7" = "Branch B3 (IPAN)",
                        "8" = "Tuft",
                        "9" = "Branch A4 (IN)",
                        "10" = "Paneth",
                        "11" = "Glia",
                        "12" = "Undefined")

Epi_neu$cell.type.minor_v1 <- Idents(Epi_neu)


## 3.8 refine Epithelial ----------------------
Epithelial <- subset(integrated,subset = cell.type.major_v3 == "Epithelial")
DefaultAssay(Epithelial) <- "RNA"

ndim = 20;resolution = 1

Epithelial <- ScaleData(Epithelial)
Epithelial <- RunPCA(Epithelial, npcs = 30, verbose = FALSE)

Epithelial <- FindNeighbors(Epithelial, dims = 1:ndim,
                            features = NULL, k.param = 20,
                            compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                            annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                            force.recalc = FALSE)
Epithelial <- FindClusters(Epithelial, resolution = resolution, algorithm = 1)

#Run UMAP
Epithelial <- RunUMAP(Epithelial, reduction = "pca", dims = 1:ndim, 
                      n.neighbors = 30L, n.components = 2L,
                      min.dist = 1.2, spread = 1)

CellDimPlot(Epithelial,
            group.by = "RNA_snn_res.1", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

CellDimPlot(Epithelial,
            group.by = "pred.cell.type.minor_v2", 
            label = T,
            label_insitu = TRUE,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right"
)

a <- as.data.frame(table(Epithelial$RNA_snn_res.1,Epithelial$pred.cell.type.minor))

Idents(Epithelial) = "RNA_snn_res.1"
Epithelial <- RenameIdents(object = Epithelial, 
                           "0" = "TA",
                           "1" = "TA",
                           "2" = "TA",
                           "3" = "BEST2+ Goblet cell",
                           "4" = "TA",
                           "5" = "Stem cells",
                           "6" = "BEST2+ Goblet cell",
                           "7" = "Stem cells",
                           "8" = "Colonocyte",
                           "9" = "TA",
                           "10" = "TA",
                           "11" = "BEST4+ epithelial",
                           "12" = "TA",
                           "13" = "TA",
                           "14" = "TA",
                           "15" = "Cycling B cell",
                           "16" = "TA",
                           "17" = "TA",
                           "18" = "BEST2+ Goblet cell",
                           "19" = "Tuft",
                           "20" = "TA",
                           "21" = "TA",
                           "22" = "TA",
                           "23" = "TA",
                           "24" = "Microfold cell",
                           "25" = "TA")

Epithelial$cell.type.minor_v1 <- Idents(Epithelial)


## 3.9 merge cell types ----------------------------
metadata_celltypes <- rbind(b_plasm@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            endo@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            Epi_neu@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            Epithelial@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            mes@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            Myeloid@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            Neuronal@meta.data[,c("orig.ident","cell.type.minor_v1")],
                            t_cells@meta.data[,c("orig.ident","cell.type.minor_v1")])

metadata_celltypes <- metadata_celltypes[colnames(integrated),]

integrated$cell.type.minor_v1 = metadata_celltypes$cell.type.minor_v1

saveRDS(integrated,file = "./output/Data/integrated_annot_v4.Rds")
