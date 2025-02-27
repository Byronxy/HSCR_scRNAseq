options(reticulate.conda_binary = "/home/bio/miniconda3/condabin/conda", SCP_env_name = "SCP_env")
library(SCP)
library(cowplot)

library(Seurat)
library(reticulate)
library(CytoTRACE)
library(Cairo)
library(scCustomize)
'%ni%' <- Negate('%in%')
integrated <- readRDS(file = "./output/Data/integrated_annot_v4.Rds")
DefaultAssay(integrated) = 'RNA'

## figure 1b----------------
integrated$cell.type.major_v5 = factor(integrated$cell.type.major_v5,levels = c("Epithelial",
                                                                                "Stromal",
                                                                                "Fibroblasts/Pericytes",
                                                                                "Neuronal",
                                                                                "Endothelial",
                                                                                "B/plasma cells",
                                                                                "Myeloid",
                                                                                "T cell"))
genes <- c("EPCAM","FABP1","MT1G","APOA1",
           "PDGFRA","BMP4","F3","VIM",
           "HHIP","KCNJ8","ABCC9","BGN","THY1",
           "NPM1","S100B","PLP1",
           "CLDN5","CDH5","PECAM1","RAMP2",
           "CD79B","JCHAIN","IGHA1","IGHG1",
           "S100A4","LYZ","HLA-DRB1",
           "CD3D","CD69","CD8B","CD8A","CD4","FOS")
DotPlot_Custom(integrated, 
               # group.by = "cell.type.major_v5",
               features = genes)


## figure 1c------------------------
meta_data = integrated@meta.data %>%
  dplyr::select(Individual_Name,
                Sample_Month,
                Sample_Condition,
                cell.type.major_v5)


library(ggplot2)
library(ggh4x)

meta_data$type <- factor(meta_data$type, c("pre", "post"))

p1 <- ggplot(meta_data, aes(
  y = Individual_Name
)) +
  geom_bar(aes(fill = cell.type.major_v5), position = "fill") +
  scale_fill_manual(values = cluster_colors) +
  ggnewscale::new_scale_fill() +
  geom_col(data = dplyr::distinct(meta_data, Individual_Name, Sample_Month, Sample_Condition), aes(x = -.3, fill = Sample_Month), width = .9) +
  geom_col(data = dplyr::distinct(meta_data, Individual_Name, Sample_Month, Sample_Condition), aes(x = -.15), fill = "white", color = "white") +
  geom_text(data = dplyr::distinct(meta_data, Individual_Name, Sample_Month, Sample_Condition), aes(x = -.075, label = Individual_Name)) +
  scale_fill_manual(values = month_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_nested(Sample_Condition + Sample_Month ~ .,
               scale = "free_y", space = "free_y",
               strip = strip_nested(
                 by_layer_y = TRUE,
                 background_y = list(element_rect(), element_blank()),
                 text_y = list(element_text(), element_blank())
               )
  )

p1

p2 <- ggplot() +
  geom_text(data = data.frame(
    y = 21,
    x = c(-.225, -0.075, .5),
    label = c("Cohort", "ID", "Cell type proporation")
  ), aes(x = x, y = y, label = label), fontface = "bold") +
  scale_x_continuous(expand = c(0, 0), limits = c(-.3, 1)) +
  coord_cartesian(clip = "off")

library(patchwork)


p2 + p1 & plot_layout(heights = c(1, 20)) &
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank()
  ) &
  labs(x = NULL, y = NULL)

## figure 1d----------------------
circ_data <- prepare_circlize_data(integrated, scale = 0.7)
set.seed(1234)

plot_circlize(circ_data,do.label = F, pt.size = 0.05, 
              col.use = cluster_colors ,
              bg.color = 'white', 
              kde2d.n = 1000, 
              repel = T, 
              label.cex = 0.1)

add_track(circ_data, 
          group = "Sample_Month", 
          colors = month_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, 
          group = "Condition",
          colors = condition_colors, track_num = 3)
add_track(circ_data, 
          group = "Sample_Condition",
          colors = sample_condition_colors, track_num = 4)

## figure 2a------------------
epi <- readRDS( "./output/Data/epi_annot_v4_monocle.Rds")
CellDimPlot(
  srt = epi, group.by = "cell.type.minor_v2",
  reduction = "UMAP", theme_use = "theme_blank",
  palcolor = epi_colors
)


## figure 2b--------------------
cellcycle <- c("TOP2A","MKI67","CENPF","HMGB2","SMC4")
inflammation <- c("TNF","CD44","IL6","IL10","IL23A")
genes <- c(cellcycle,inflammation)

plot_df = clust.genes_sub %>% filter(feature %in% c(genes))
plot_df$feature = factor(plot_df$feature, levels = genes)
library(ggpubr)
p = ggplot(data = plot_df, aes(x = feature, y = group)) +
  geom_point(aes(size = -log10(padj), fill = logFC), shape = 21, color = 'black', stroke = 1) +
  labs(x = '', y = '') +
  scale_size_area(name = '-log10(q-value)') +
  scale_fill_gradient2(low = '#08519c', 
                       mid = 'white', 
                       high = 'red',
                       limits = c(-0.3, 0),
                       name = 'logFC\n (HD vs. HAD)',
                       na.value = '#08519c') +
  theme_pubr() +
  RotatedAxis()

## figure 2c-----------------------
df = meta %>%
  select(Condition,cell.type.cycle,cell.type.minor_v2) %>%
  group_by(Condition,cell.type.minor_v2,cell.type.cycle) %>%
  summarise (n = n()) %>%
  mutate(freq = n / sum(n))
df = df %>% filter(cell.type.cycle %in% c("Cycling"))
# df$freq = df$freq*100

library(RColorBrewer)
p = df %>%
  filter(cell.type.minor_v2 %in% c("Stem cells","TA","Goblet cells","Colonocyte")) %>%
  mutate(cell.type.minor_v2 = factor(cell.type.minor_v2,levels = c("Stem cells","TA","Goblet cells","Colonocyte"))) %>%
  ggplot(aes(x = cell.type.minor_v2,y = Condition,fill = freq)) +
  geom_jjpie(aes(piefill = freq),
             color = "black",
             # strock = 1,
             width = 1) +
  coord_fixed() +
  # scale_fill_gradient2(low = '#d53d4e',mid = '#fefcbf',high = '#3589bd',
  #                      midpoint = 0.2,
  #                      name = 'Cycling cell \nProportion') +
  scale_fill_distiller(palette = "Spectral") +
  xlab('') + ylab('') +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', size = 15),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'gray', linewidth = 2),
        axis.title.x = element_text(family = 'sans', face = 'bold', size = 15),
        axis.text = element_text(family = 'sans', colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.ticks = element_line(linewidth = 1))


# figure 2h-----------------------------

## monocle
Mat_epi <- epi_had@assays$RNA@counts
#Generate barcode file.
Barcode_epi <- epi_had@meta.data
#Generate feature file.
Feature_epi <- rownames(epi_had)
Feature_epi <- as.data.frame(Feature_epi)
colnames(Feature_epi) <- c("gene_short_name")
rownames(Feature_epi) <- Feature_epi$gene_short_name
#Establish the monocle2 object
pd<-new("AnnotatedDataFrame",data = Barcode_epi)
fd<-new("AnnotatedDataFrame",data = Feature_epi)
Monocle_epi <- newCellDataSet(Mat_epi,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size()) 
#Estimate Size Factors and Dispersions.
Monocle_epi <- estimateSizeFactors(Monocle_epi)
Monocle_epi <- estimateDispersions(Monocle_epi)
#Generate intestine differentiation for trajectory reconstruction
Gene_list <- c("SPINK4", "AGR2", "REG4", "MUC1", "TFF3", "SPDEF", "MUC2", "AQP8", "MS4A12", "SLC26A3", "CA2", "CA1", "KRT20", "GUCA2B", "CEACAM1", "RPS20", "RPS12", "RPS19", "RPS29",
               "CFTR", "CD44", "PROM1","MLLT10", "PTPLAD1", "OLFM4", "CDCA7", "RGMB", "PTPRO", "CDK6", "RNF43", "DNMT3A", "EZH2", "STMN1", "METTL3", "LGR5", "ASCL2")
HSMM_data <- Monocle_epi

# HSMM_data <- setOrderingFilter(HSMM_data,Gene_list)

epi_had <- FindVariableFeatures(epi_had)
var.genes <- VariableFeatures(epi_had)
HSMM_data <- setOrderingFilter(HSMM_data, var.genes)

HSMM_data <- reduceDimension(HSMM_data, max_components = 2,
                             norm_method = "log", pseudo_expr = 1,
                             reduction_method = 'DDRTree', verbose = T,
                             reducedModelFormulaStr = "~nCount_RNA + sample + percent.mt + S.Score + G2M.Score",
                             # auto_param_selection = F,
                             cores = 10)

HSMM_data <- orderCells(HSMM_data,root_state = 11)

Mat_epi <- epi_hd@assays$RNA@counts
#Generate barcode file.
Barcode_epi <- epi_hd@meta.data
#Generate feature file.
Feature_epi <- rownames(epi_hd)
Feature_epi <- as.data.frame(Feature_epi)
colnames(Feature_epi) <- c("gene_short_name")
rownames(Feature_epi) <- Feature_epi$gene_short_name
#Establish the monocle2 object
pd<-new("AnnotatedDataFrame",data = Barcode_epi)
fd<-new("AnnotatedDataFrame",data = Feature_epi)
Monocle_epi <- newCellDataSet(Mat_epi,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size()) 
#Estimate Size Factors and Dispersions.
Monocle_epi <- estimateSizeFactors(Monocle_epi)
Monocle_epi <- estimateDispersions(Monocle_epi)
#Generate intestine differentiation for trajectory reconstruction
Gene_list <- c("SPINK4", "AGR2", "REG4", "MUC1", "TFF3", "SPDEF", "MUC2", "AQP8", "MS4A12", "SLC26A3", "CA2", "CA1", "KRT20", "GUCA2B", "CEACAM1", "RPS20", "RPS12", "RPS19", "RPS29",
               "CFTR", "CD44", "PROM1","MLLT10", "PTPLAD1", "OLFM4", "CDCA7", "RGMB", "PTPRO", "CDK6", "RNF43", "DNMT3A", "EZH2", "STMN1", "METTL3", "LGR5", "ASCL2")
HSMM_data <- Monocle_epi

HSMM_data <- setOrderingFilter(HSMM_data,Gene_list)

HSMM_data <- reduceDimension(HSMM_data, max_components = 2,
                             norm_method = "log", pseudo_expr = 1,
                             reduction_method = 'DDRTree', verbose = T,
                             reducedModelFormulaStr = "~nCount_RNA + sample + percent.mt + S.Score + G2M.Score",
                             # auto_param_selection = F,
                             cores = 10)

HSMM_data <- orderCells(HSMM_data)

HSMM_data_had <- readRDS(file = "./output/Data/monocle2_epi_had.Rds")
HSMM_data_hd <- readRDS(file = "./output/Data/monocle2_epi_hd.Rds")

cds.list <- list()
cds.list[['HAD']] = HSMM_data_had
cds.list[['HD']] = HSMM_data_hd

# DTP analysis follow the paper: A pooled single-cell genetic screen identifies regulatory checkpoints in the
# continuum of the epithelial-to-mesenchymal transition
source("./code/DTP_function.R")

# Identify genes that are expressed in at least 50 of cells 
expressed_genes.list <- list()
expressed_genes.list[["HAD"]] <- row.names(fData(HSMM_data_had)[Matrix::rowSums(Biobase::exprs(HSMM_data_had) > 0) > 50 ,])
length(expressed_genes.list[["HAD"]])
expressed_genes.list[["HD"]] <- row.names(fData(HSMM_data_hd)[Matrix::rowSums(Biobase::exprs(HSMM_data_hd) > 0) > 50 ,])
length(expressed_genes.list[["HD"]])

expressed_genes <- unique(union(expressed_genes.list[["HAD"]], expressed_genes.list[["HD"]]))
length(expressed_genes)

# Use dynamic time warping to align HAD and HD pseudospatial trajectories and create a cds object of aligned trajectories
HD.to.HAD.aligned.cds <- getDTWcds(cds.list[['HD']],cds.list[['HAD']], 
                                   ref = "HAD", query = "HD", 
                                   expressed_genes = expressed_genes, cores = 5)
cds.aligned.list <- list()
cds.aligned.list[["HD to HAD"]] <- HD.to.HAD.aligned.cds

head(pData(cds.aligned.list[["HD to HAD"]]))

for(alignment in names(cds.aligned.list)){
  
  cds.aligned.list[[alignment]] <- preprocess_cds(cds.aligned.list[[alignment]])
  
}


# Identify genes that are differentially expressed across pseudospace as a function of treatment
aligned.pseudospace.DEG.test.list <- list()

for(alignment in names(cds.aligned.list)) {
  
  aligned.pseudospace.DEG.test.list[[alignment]] <- differentialGeneTest(cds.aligned.list[[alignment]][expressed_genes], 
                                                                         fullModelFormulaStr="~sm.ns(Pseudotime, df=3)*Cell.Type", 
                                                                         reducedModelFormulaStr="~sm.ns(Pseudotime, df=3) + Cell.Type", cores = 5)
  
}

Pseudospatial.aligned.sig.genes.list <- list()

df = aligned.pseudospace.DEG.test.list[[1]]

for(sample in names(aligned.pseudospace.DEG.test.list)){
  
  Pseudospatial.aligned.sig.genes.list[[sample]] <- row.names(subset(aligned.pseudospace.DEG.test.list[[sample]], 
                                                                     qval <= 1e-10))
  print(sample)
  print(length(Pseudospatial.aligned.sig.genes.list[[sample]]))
}

df = data.frame(gene = pseudospatial_aligned_sig_maxRank_genes)

## Visualization trajectory pseudotime
pseudospatial_aligned_sig_maxRank_genes = Pseudospatial.aligned.sig.genes.list[["HD to HAD"]]

aligned.pseudospace.smoothed.combined.exprs.list <- list()

for(alignment in names(cds.aligned.list)){
  
  aligned.pseudospace.smoothed.combined.exprs.list[[alignment]] <- genSmoothCurves(cds.aligned.list[[alignment]][pseudospatial_aligned_sig_maxRank_genes], 
                                                                                   new_data=data.frame(Pseudotime=rep(seq(0,99, by=1), 2), 
                                                                                                       Cell.Type=c(rep(unique(pData(cds.aligned.list[[alignment]])$Cell.Type)[1], 100), 
                                                                                                                   rep(unique(pData(cds.aligned.list[[alignment]])$Cell.Type)[2], 100))), 
                                                                                   trend_formula="~sm.ns(Pseudotime, df=3)*Cell.Type",
                                                                                   cores=1)
  
}

for(alignment in names(aligned.pseudospace.smoothed.combined.exprs.list)){
  
  aligned.pseudospace.smoothed.combined.exprs.list[[alignment]] <- aligned.pseudospace.smoothed.combined.exprs.list[[alignment]][rowSums(is.na(aligned.pseudospace.smoothed.combined.exprs.list[[alignment]])) == 0,]
  
}

aligned.pseudospace.smoothed.combined.exprs.scaled.list <- list()

for(alignment in names(aligned.pseudospace.smoothed.combined.exprs.list)){
  
  aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]] <- scale(t(scale(t(aligned.pseudospace.smoothed.combined.exprs.list[[alignment]]))))
  aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]][aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]] > 3] <- 3
  aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]][aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]] < -3] <- -3
  aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]][is.na(aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]])] <- 0
  
}

row_dist.list <- list()

for(alignment in names(aligned.pseudospace.smoothed.combined.exprs.scaled.list)){
  
  row_dist.list[[alignment]] <- as.dist((1 - cor(Matrix::t(aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]])))/2)
  row_dist.list[[alignment]][is.na(row_dist.list[[alignment]])] <- 1
  
}

dim(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['HD to HAD']])

saveRDS(aligned.pseudospace.smoothed.combined.exprs.scaled.list,file = "./output/Data/aligned.pseudospace.smoothed.combined.exprs.scaled.list.Rds")

aligned.pseudospace.smoothed.combined.exprs.scaled.list <- readRDS(file = "./output/Data/aligned.pseudospace.smoothed.combined.exprs.scaled.list.Rds")

library(ComplexHeatmap)
Gene_list[Gene_list %in% pseudospatial_aligned_sig_maxRank_genes]
genes_to_show <- c("SPINK4", "AGR2",   "REG4",   "TFF3",   "SPDEF",  "MUC2",   "CA2",    "CA1",    "KRT20",  "GUCA2B", "RPS20" , "RPS12" , "RPS19",  "RPS29" , "CFTR",
                   "CD44",   "OLFM4",  "CDCA7",  "RGMB",   "STMN1" , "ASCL2")
mark_at = which(rownames(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['HD to HAD']]) %in% genes_to_show)
ha = rowAnnotation(foo = anno_mark(at = mark_at, labels = genes_to_show))

heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
annCol = data.frame(Cell.Type=c(rep(unique(pData(cds.aligned.list[[alignment]])$Cell.Type)[1], 100), 
                                rep(unique(pData(cds.aligned.list[[alignment]])$Cell.Type)[2], 100)))

annColors <- list()
rownames(annCol) = colnames(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['HD to HAD']])
annColors[["Cell.Type"]] <- c("HD" = "#ee0301" ,
                              "HAD" =  "#00198a")
p <- ComplexHeatmap::pheatmap(mat = aligned.pseudospace.smoothed.combined.exprs.scaled.list[['HD to HAD']],
                              cluster_cols=FALSE, 
                              cluster_rows=TRUE,
                              show_rownames=F, 
                              show_colnames=F, 
                              #color = NMF:::ccRamp(x = heatmap.BlBkRd,n = 64),
                              color = circlize::colorRamp2(c(-2, 0, 2), c("cyan", "black", "yellow")),
                              clustering_distance_rows=row_dist.list[["HD to HAD"]],
                              clustering_method = "ward.D2",
                              cutree_rows=num_clusters,
                              right_annotation = ha,
                              annotation_col = annCol[c("Cell.Type")],
                              annotation_colors = annColors[c("Cell.Type")],
                              # cellheight = 0.1, # 热图高度固定
                              # cellwidth = 0.1
                              # row_labels = genes_to_show
)




# figure 3a-----------
epi <- subset(integrated,subset = cell.type.major_v5 == "Epithelial")
stem <- subset(epi,subset = cell.type.minor_v2 == "Stem cells")

stem <- ScaleData(stem)
stem <- RunPCA(stem, npcs = 30, verbose = FALSE)

stem <- FindNeighbors(stem, dims = 1:20,
                      features = NULL, k.param = 20,
                      compute.SNN = TRUE, prune.SNN = 1/15, nn.method = "rann",
                      annoy.metric = "euclidean", nn.eps = 0, verbose = TRUE,
                      force.recalc = FALSE)
stem <- FindClusters(stem, resolution = 0.8, algorithm = 1)
# Run UMAP
stem <- RunUMAP(stem, reduction = "pca", dims = 1:20, 
                n.neighbors = 30L, n.components = 2L,
                min.dist = 1.2, spread = 1)

stem <- RenameIdents(object = stem, 
                     "0" = "C5_OLFM4lo/FBAP1+",
                     "1" = "C4_OLFM4+",
                     "2" = "C3_OLFM4+/LGR5lo",
                     "3" = "C3_OLFM4+/LGR5lo",
                     "4" = "C2_LGR5hi",
                     "5" = "C1_CLU+",
                     "6" = "C4_OLFM4+",
                     "7" = "C2_LGR5hi",
                     "8" = "C5_OLFM4lo/FBAP1+",
                     "9" = "C3_OLFM4+/LGR5lo",
                     "10" = "C4_OLFM4+",
                     "11" = "C4_OLFM4+")

CellDimPlot(stem,
            group.by = c("cell.type.minor_v3"), 
            label = F,
            label_insitu = F,
            reduction = "UMAP", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right",
            palcolor = stem_colors
)

# figure 3c------
# remotes::install_github("jr-leary7/SCISSORS")
stem_phate <- SCISSORS::RunPHATE(
  object = stem,
  n.components = 2,
  n.PC = 20,
  mds.method = "smacof",
  dist.metric = "cosine",
  random.seed = 312
)

CellDimPlot(stem_phate,
            group.by = c("cell.type.minor_v3"), 
            label = F,
            label_insitu = F,
            reduction = "phate", show_stat = FALSE,
            theme_use = "theme_blank", 
            legend.position = "right",
            palcolor = stem_colors
)

# figure 5a -----------------------
integrated.list <- SplitObject(integrated,"Sample_Condition")
liana_integrated.list <- lapply(integrated.list,liana_wrap)
liana_integrated.list_aggr <-  lapply(liana_integrated.list,function(x){x %>% liana_aggregate()})

pval_list <- map(c(1:4),function(x){
  liana_trunc <- liana_integrated.list_aggr[[x]] %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.05) # note that these pvals are already corrected
  
  heat_freq(liana_trunc)
  freqs <- liana_trunc %>% group_by(source, target) %>% 
    summarise(COUNT = n(),.groups = "keep") 
  freqs$group = names(liana_integrated.list_aggr)[x]
  return(freqs)
})

pval_df <- do.call(rbind,pval_list)
pval_df = pval_df %>%
  filter(source != target) 
pval_df$merge = paste(pval_df$source, pval_df$target, sep = " | ")

p = pval_df %>%  
  mutate(group = factor(group ,levels = c("Control_D","Control_R","HD_D","HD_R"))) %>% 
  filter(target %in% c("Epithelial")) %>%
  ggplot(aes(x = group, y = merge, color = group,size = COUNT,fill = COUNT)) +
  geom_point(shape = 21, color = 'black', stroke = 1) +
  scale_fill_distiller(palette = "RdYlBu") +
  scale_size_continuous(range = c(3, 15)) + # Adjusts the range of point sizes
  labs(x = "", y = "") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', size = 15),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'black', linewidth = 1),
        axis.title.x = element_text(family = 'sans', face = 'bold', size = 15),
        axis.text = element_text(family = 'sans', colour = 'black'),
        axis.ticks = element_line(linewidth = 1)) +
  RotatedAxis()

# figure 5h ---------------------
pdata %>%
  filter(source %in% c("Stromal1","Stromal3","Stromal2","Stromal4","Stromal3_Transitional")) %>%
  ggplot(aes(x = target, y = interaction)) +
  geom_point(aes(size = specificity, fill = magnitude), shape = 21, color = 'black', stroke = 0.5) +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(.~source,nrow = 1) +
  xlab('') + ylab('') +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(family = 'sans', size = 15),
        panel.background = element_rect(fill = 'white'),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'gray', linewidth = 2),
        axis.title.x = element_text(family = 'sans', face = 'bold', size = 15),
        axis.text = element_text(family = 'sans', colour = 'black'),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.ticks = element_line(linewidth = 1))

# figure 6b-------
library(UCell)
library(msigdbr)
library(fgsea)
m_df <-  msigdbr(species = "Homo sapiens")
m_df_hallmark <- m_df %>% filter(str_detect(gs_name, "HALLMARK"))
unique(m_df_hallmark$gs_name)
m_df_sub <- m_df %>% filter(gs_name %in% unique(m_df_hallmark$gs_name))
m_df_sub_list <- m_df_sub %>% 
  split(x = .$gene_symbol, f= .$gs_name)

stroma <- AddModuleScore_UCell(stroma, features = m_df_sub_list)
FeatureStatPlot(stroma, 
                stat.by = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING_UCell",
                            "HALLMARK_INFLAMMATORY_RESPONSE_UCell"), 
                group.by = "cell.type.minor_v3",
                palcolor = stroma_colors,
                # split.by = "Condition",
                add_box = TRUE
                # pairwise_method = "t.test",
                #comparisons = T
)

# figure 7d-------
ggplot(markers_s4, aes(rank, -log10(p_val_adj))) + 
  geom_point(size=3, color=markers_s4$pt.col) + 
  ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(rank, -log10(p_val_adj), label=gene), size=3) +
  ylab("-log10(p-val)") + 
  ggpubr::theme_pubr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = .5)
  )



