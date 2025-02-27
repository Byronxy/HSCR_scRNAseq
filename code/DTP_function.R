source("./script/Pseudospace_support_functions.R")

getDTWcds <- function (query_cds, ref_cds, ref, query, expressed_genes, cores) 
{
  alignment_genes <- intersect(row.names(subset(fData(ref_cds), 
                                                use_for_ordering)), row.names(subset(fData(query_cds), 
                                                                                     use_for_ordering)))
  ref_align_cds <- ref_cds[alignment_genes]
  query_align_cds <- query_cds[alignment_genes]
  pData(ref_align_cds)$cell_id <- row.names(pData(ref_align_cds))
  pData(ref_align_cds)$Pseudotime <- 100 * pData(ref_align_cds)$Pseudotime/max(pData(ref_align_cds)$Pseudotime)
  ref_align_cds <- ref_align_cds[alignment_genes, as.character(arrange(pData(ref_align_cds), 
                                                                       Pseudotime)$cell_id)]
  pData(query_align_cds)$cell_id <- row.names(pData(query_align_cds))
  pData(query_align_cds)$Pseudotime <- 100 * pData(query_align_cds)$Pseudotime/max(pData(query_align_cds)$Pseudotime)
  query_align_cds <- query_align_cds[alignment_genes, as.character(arrange(pData(query_align_cds), 
                                                                           Pseudotime)$cell_id)]
  smoothed_ref_exprs <- genSmoothCurves(ref_align_cds[alignment_genes], 
                                        data.frame(Pseudotime = seq(0, 100, by = 1)), cores = cores)
  smoothed_ref_exprs <- smoothed_ref_exprs[rowSums(is.na(smoothed_ref_exprs)) == 
                                             0, ]
  vst_smoothed_ref_exprs <- vstExprs(ref_cds, expr_matrix = smoothed_ref_exprs)
  smoothed_query_exprs <- genSmoothCurves(query_align_cds[alignment_genes], 
                                          data.frame(Pseudotime = seq(0, 100, by = 1)), cores = cores)
  smoothed_query_exprs <- smoothed_query_exprs[rowSums(is.na(smoothed_query_exprs)) == 
                                                 0, ]
  vst_smoothed_query_exprs <- vstExprs(query_cds, expr_matrix = smoothed_query_exprs)
  alignment_genes <- intersect(row.names(vst_smoothed_ref_exprs), 
                               row.names(vst_smoothed_query_exprs))
  ref_matrix <- t(scale(t(vst_smoothed_ref_exprs[alignment_genes, 
  ])))
  query_matrix <- t(scale(t(vst_smoothed_query_exprs[alignment_genes, 
  ])))
  
  ref_query_dtw <- align_cells(ref_matrix, query_matrix, step_pattern = rabinerJuangStepPattern(3, 
                                                                                                "c"), open.begin = F, open.end = F)
  align_res <- warp_pseudotime(ref_align_cds, query_align_cds, 
                               ref_query_dtw)
  query_ref_aligned <- align_res$query_cds
  pData(query_ref_aligned)$Pseudotime <- pData(query_ref_aligned)$Alignment_Pseudotime
  ref_aligned_cell_ids <- setdiff(row.names(pData(ref_align_cds)), 
                                  "duplicate_root")
  query_aligned_cell_ids <- setdiff(row.names(pData(query_align_cds)), 
                                    "duplicate_root")
  combined_exprs <- cbind(exprs(query_cds[expressed_genes, 
                                          query_aligned_cell_ids]), exprs(ref_cds[expressed_genes, 
                                                                                  ref_aligned_cell_ids]))
  pData_ref <- pData(ref_align_cds)[, c("Condition", "cell.type.minor_v2", 
                                        "Pseudotime")]
  pData_ref$Cell.Type <- ref
  pData_query_aligned <- pData(query_ref_aligned)[, c("Condition", 
                                                      "cell.type.minor_v2", "Pseudotime")]
  pData_query_aligned$Cell.Type <- query
  combined_pData <- rbind(pData_query_aligned, pData_ref)
  combined_pData <- combined_pData[colnames(combined_exprs), 
  ]
  combined_pd <- new("AnnotatedDataFrame", data = combined_pData)
  
  fd <- new("AnnotatedDataFrame", data = fData(ref_cds)[row.names(combined_exprs), 
                                                        1:2])
  
  ref_queryToRef_combined_cds <- newCellDataSet(combined_exprs, 
                                                phenoData = combined_pd, 
                                                featureData = fd, 
                                                expressionFamily = negbinomial.size(), 
                                                lowerDetectionLimit = 1)
  pData(ref_queryToRef_combined_cds)$cell_id <- row.names(pData(ref_queryToRef_combined_cds))
  return(ref_queryToRef_combined_cds)
}

compare_cell_types_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)*Cell.Type", 
                                               min_expr = NULL, cell_size = 0.75, cell_alpha = 1, line_size = 1, 
                                               nrow = NULL, ncol = 1, panel_order = NULL, color_by = "Cell.Type", 
                                               shade_by = NULL, df = 3, maxit = 300, relative_expr = TRUE, 
                                               pseudocount = 0) 
{
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- Biobase::exprs(cds_subset)
    cds_exprs <- cds_exprs + pseudocount
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- t(t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- melt(ceiling(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- Biobase::exprs(cds_subset)
    cds_exprs <- cds_exprs + pseudocount
    cds_exprs <- melt(as.matrix(cds_exprs))
  }
  colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  if (integer_expression) {
    cds_exprs$expression <- cds_exprs$expression
  }
  else {
    cds_exprs$expression <- log10(cds_exprs$expression)
  }
  if (is.null(cds_exprs$gene_short_name) == FALSE) {
    cds_exprs$gene_label <- as.character(cds_exprs$gene_short_name)
    cds_exprs$gene_label[is.na(cds_exprs$gene_label)] <- cds_exprs$gene_id
  }
  else {
    cds_exprs$gene_label <- cds_exprs$gene_id
  }
  cds_exprs$gene_label <- factor(cds_exprs$gene_label)
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, 
                         Cell.Type = pData(cds_subset)$Cell.Type)
  model_expectation <- genSmoothCurves(cds_subset, cores = 1, 
                                       trend_formula = trend_formula, relative_expr = T, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  cds_exprs$expectation <- apply(cds_exprs, 1, function(x) model_expectation[x[2], 
                                                                             x[1]])
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (is.null(panel_order) == FALSE) {
    cds_exprs$gene_label <- factor(cds_exprs$gene_label, 
                                   levels = panel_order)
  }
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), 
                        alpha = I(cell_alpha))
  }
  else {
    q <- q + geom_point(size = I(cell_size))
  }
  q <- q + geom_line(aes(Pseudotime, expectation, color = Cell.Type), 
                     size = line_size)
  q <- q + scale_y_log10() + facet_wrap(~gene_label, nrow = nrow, 
                                        ncol = ncol, scales = "free_y")
  q <- q + ylab("Expression") + xlab("Pseudospace")
  q <- q + monocle:::monocle_theme_opts()
  q
}

delta_auc_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)*Cell.Type", 
                                      min_expr = NULL, df = 3, maxit = 300, relative_expr = TRUE, 
                                      pseudocount = 0, cores = 1) 
{
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- Biobase::exprs(cds_subset)
    cds_exprs <- cds_exprs + pseudocount
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- t(t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- melt(ceiling(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- Biobase::exprs(cds_subset)
    cds_exprs <- cds_exprs + pseudocount
    cds_exprs <- melt(as.matrix(cds_exprs))
  }
  colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
                     by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  if (integer_expression) {
    cds_exprs$expression <- cds_exprs$expression
  }
  else {
    cds_exprs$expression <- log10(cds_exprs$expression)
  }
  if (is.null(cds_exprs$gene_short_name) == FALSE) {
    cds_exprs$gene_label <- as.character(cds_exprs$gene_short_name)
    cds_exprs$gene_label[is.na(cds_exprs$gene_label)] <- cds_exprs$gene_id
  }
  else {
    cds_exprs$gene_label <- cds_exprs$gene_id
  }
  cds_exprs$gene_label <- factor(cds_exprs$gene_label)
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, 
                         Cell.Type = pData(cds_subset)$Cell.Type)
  model_expectation <- genSmoothCurves(cds_subset, cores = cores, 
                                       trend_formula = trend_formula, relative_expr = T, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  cds_exprs$expectation <- apply(cds_exprs, 1, function(x) model_expectation[x[2], 
                                                                             x[1]])
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  pseudotime_max = max(cds_exprs$Pseudotime)
  cds_exprs.auc = cds_exprs %>% dplyr::select(id, Cell.Type, 
                                              Pseudotime, expectation) %>% dplyr::mutate(pseudotime_quantile = cut(Pseudotime, 
                                                                                                                   breaks = seq(0, pseudotime_max, pseudotime_max/5))) %>% 
    unique() %>% dplyr::group_by(id, Cell.Type, pseudotime_quantile) %>% 
    dplyr::summarize(AUC = MESS::auc(Pseudotime, expectation, 
                                     type = "spline")) %>% dplyr::group_by(id, pseudotime_quantile) %>% 
    dplyr::summarize(auc_difference = (AUC[1] - AUC[2])/(AUC[1] + 
                                                           AUC[2])) %>% dplyr::arrange(desc(auc_difference))
  cds_exprs.auc$pseudotime_quantile = as.numeric(factor(cds_exprs.auc$pseudotime_quantile, 
                                                        levels = sort(unique(cds_exprs.auc$pseudotime_quantile))))
  cds_exprs.auc_rank = tidyr::spread(cds_exprs.auc, key = pseudotime_quantile, 
                                     value = auc_difference)
  return(as.data.frame(cds_exprs.auc_rank))
}

get_gsea_sig_results <- function(GSAhyper_list, qval_cutoff){
  
  GSA_hyper_results.list <- list()
  
  for(cluster in names(GSAhyper_list)){
    
    GSAhyper_df <- as.data.frame(GSAhyper_list[[cluster]]$p.adj)
    GSAhyper_df$gene_set <- row.names(GSAhyper_df)
    colnames(GSAhyper_df) <- c("qval","gene_set")
    
    GSA_hyper_results.list[[cluster]] <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>% 
      mutate(gene_set = factor(gene_set, levels = gene_set))
    
  }
  
  return(GSA_hyper_results.list)
  
}
