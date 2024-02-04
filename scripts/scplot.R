
#######################
# N cells per cluster #
#######################


get_cell_number<- function(obj, space, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  write_tsv(df, here(glue::glue("results/tsv/{space}/{space}_number_of_cells_per_cluster.tsv")))
  
  number_of_cells_per_cluster_wider<- df  %>% 
    pivot_wider(names_from = seurat_clusters, values_from = n)
  
  write_tsv(df, here(glue::glue("results/tsv/{space}/{space}_number_of_cells_per_cluster_wide.tsv")))
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y")
  
  ggsave(filename = here(glue::glue("results/figures/{space}/final_annotation/{space}_number_of_cells_per_cluster.pdf")),
         plot = p1, width = 8, height = 50,
         limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =percent)) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y")
  
  ggsave(filename = here(glue::glue("results/figures/{space}/final_annotation/{space}_percent_of_cells_per_cluster.pdf")), plot = p2, width = 8, height = 50,
         limitsize = FALSE)
}


#######################
# Get Remain Clusters #
#######################

#>41 cells with >= 3 samples

cluster_remain_41_exclusive<- function(obj, space, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>% count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  cluster_remain<- c()
  for (i in 0:(length(number_of_cells_per_cluster$seurat_clusters  %>% unique) -1)){
    n_sample<- number_of_cells_per_cluster  %>% filter(seurat_clusters == i, n > 41)  %>% pull(sample_id)  %>% length()
    if (n_sample >= 3){
      cluster_remain<- c(cluster_remain, i)
    }
  }
  print(space)
  print(cluster_remain)
  
  number_of_cells_per_cluster_wider<- number_of_cells_per_cluster  %>%
    pivot_wider(names_from = seurat_clusters, values_from = n)
  
  number_of_cells_per_cluster_wider_remain<- number_of_cells_per_cluster  %>%
    filter(seurat_clusters %in% cluster_remain)  %>%
    pivot_wider(names_from = seurat_clusters, values_from = n)
  write_tsv(number_of_cells_per_cluster_wider, here("./results/tsv/", paste0(space, "/", space, "_number_of_cells_all_clusters_wide.tsv")))
  write_tsv(number_of_cells_per_cluster_wider_remain, here("./results/tsv/", paste0(space, "/", space, "_number_of_cells_remain_clusters_wide.tsv")))
  
  return(cluster_remain)
}



###################
# Remove Clusters #
###################

remove_clusters<- function(obj, obj_clusters, space, space_name, seurat_clusters, 
                           outdir = "./", fsuffix = "clusters_removed"){
  
  print(paste0("Export data: ", outdir, "data/objects/", space, "/", space, "_",  fsuffix ,".rds"))
  print(paste0("Export file: ", outdir, "results/figures/", space, "/", space, "_",  fsuffix ,".pdf"))
  # Remove clusters
  obj$filter <- ifelse(obj$seurat_clusters %in% obj_clusters, "Keep", "Remove")
  
  obj_removed<- subset(x = obj, subset = filter == "Keep")
  
  obj_removed$seurat_clusters<- factor(x = obj_removed$seurat_clusters,
                                       levels = obj_clusters)
  
  P_umap<- DimPlot(obj_removed, reduction = "umap",
                   label = TRUE, pt.size = 0.2 ) + 
    labs(title = NULL, y = NULL, x = NULL) +
    theme(text = element_text(size = 20))
  
  #return(P_umap)
  #Save rds
  
  ggsave(paste0(outdir, "results/figures/", space, "/", space, "_",  fsuffix ,"_umap.pdf"), P_umap, width = 12, height = 8)
  
  saveRDS(obj_removed, paste0(outdir, "data/objects/", space, "/", space, "_",  fsuffix ,".rds"))
  
  ExportSeurat(obj_removed, paste0(outdir, "data/objects/", space, "/", space, "_", fsuffix ,".bcs"), overwrite=TRUE)
  
}



########################
# Bi-Clustered Heatmap #
########################


plot_bi_clustered_heatmap <- function(obj, ident = "seurat_clusters", slot = "data",
                                      space, outdir = "./", fsuffix = "clusters_removed"){
  set.seed(123)
  clusters <- obj@meta.data %>% pull(ident) %>% unique
  markers <- presto::wilcoxauc(obj, ident , assay = slot,
                               groups_use = clusters)
  write_tsv(markers,
            paste0(outdir, "results/tsv/", space, "/", space, "_",  fsuffix ,"_all_markers.tsv"))
  
  
  top_markers<- markers  %>% filter(padj <=0.05, abs(logFC) >=0.8)
  
  write_tsv(top_markers,
            paste0(outdir, "results/tsv/", space, "/", space, "_",  fsuffix ,"_top_markers(padj<=0.05, abs(logFC) >= 0.8.tsv"))
  
  top_markers_unique <- top_markers %>% pull(feature) %>% unique()
  
  # top_10_markers <- top_markers(markers, n = 10, auc_min = .6,
  #                               pct_in_min = 20, pct_out_max = 30)
  # top10 <- top_10_markers[,-1]  %>% unlist()  %>% unique() %>% na.omit()
  
  
  # top_5_markers <- top_markers(markers, n = 5, auc_min = .6,
  #                              pct_in_min = 20, pct_out_max = 30)
  # 
  # top5 <- top_5_markers[,-1]  %>% unlist()  %>% unique() %>% na.omit()
  
  
  top_10_markers<- top_markers %>% group_by(group) %>% arrange(padj, -abs(logFC), .by_group = TRUE)  %>% slice_head(n = 10)
  
  top10<- top_10_markers  %>% pull(feature)  %>% unique()
  
  write_tsv(top_10_markers,
            paste0(outdir, "results/tsv/", space, "/", space, "_",  fsuffix ,"_top10_markers(padj<=0.05, abs(logFC) >= 0.8.tsv"))
  
  top_5_markers<- top_markers %>% group_by(group) %>% arrange(padj, -abs(logFC), .by_group = TRUE)  %>% slice_head(n = 5)
  
  top5<- top_5_markers  %>% pull(feature)  %>% unique()
  
  write_tsv(top_5_markers,
            paste0(outdir, "results/tsv/", space, "/", space, "_",  fsuffix ,"_top5_markers(padj<=0.05, abs(logFC) >= 0.8.tsv"))
  
  
  ## Built data for plotting
  obj_sub <- obj[, sample.int(ncol(obj), size = 10000)]
  
  all_cells <-  obj_sub@meta.data$seurat_clusters %in% clusters
  
  all_mat<- obj_sub[["RNA"]]@data[top_markers_unique, all_cells] %>% as.matrix()
  top_5_mat<- obj_sub[["RNA"]]@data[top5, all_cells] %>% as.matrix()
  top_10_mat<- obj_sub[["RNA"]]@data[top10, all_cells] %>% as.matrix()
  
  ## scale the rows
  all_mat<- t(scale(t(all_mat)))
  top_5_mat<- t(scale(t(top_5_mat)))
  top_10_mat<- t(scale(t(top_10_mat)))
  
  ## Heatmap annotations
  #cluster_anno<- obj_sub@meta.data$seurat_clusters[all_cells]
  cluster_id <- obj_sub$seurat_clusters[all_cells]
  genotype <-  obj_sub$genotype[all_cells]
  age <- obj_sub$age[all_cells]
  len <- cluster_id  %>% sort()  %>% unique()  %>% length()
  
  genotype <- factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+"))
  age <- factor(age, levels= c("6mos", "14mos", "18mos", "sick"))
  
  reorder_df<- data.frame(cluster_id = cluster_id,
                          genotype = genotype,
                          age = age)
  
  reorder_df<- reorder_df %>% mutate(ind = row_number()) %>%
    arrange(cluster_id, genotype, age)
  
  P_cluster = createPalette(len,  c("#ff0000", "#00ff00", "#0000ff"))
  P4 <- scales::hue_pal()(4)
  P4.2 = brewer.pal(4, "Dark2")
  
  column_ha<- HeatmapAnnotation(cluster_id = cluster_id[reorder_df$ind],
                                genotype = genotype[reorder_df$ind],
                                age = age[reorder_df$ind],
                                col = list(cluster_id = setNames(P_cluster, levels(cluster_id)),
                                           genotype = setNames(P4, levels(genotype)),
                                           age = setNames(P4.2, levels(age))),
                                na_col = "grey")
  column_split_df = data.frame(cluster_id = cluster_id[reorder_df$ind],
                               genotype = genotype[reorder_df$ind])
  
  column_ha_cluster<- HeatmapAnnotation(
    cluster_id = cluster_id,
    col = list(cluster_id = setNames(P_cluster, levels(cluster_id))),
    na_col = "grey",
    show_annotation_name = TRUE)
  
  
  ## Colors
  color_range_all_low = quantile(all_mat, 0.1)
  color_range_top5_low = quantile(top_5_mat, 0.1)
  color_range_top10_low = quantile(top_10_mat, 0.1)
  
  color_range_all_high = quantile(all_mat, 0.95)
  color_range_top5_high = quantile(top_5_mat, 0.95)
  color_range_top10_high = quantile(top_10_mat, 0.95)
  
  # col_fun_all = circlize::colorRamp2(c(color_range_all_low,0,1), viridis(3))
  # col_fun_top5 = circlize::colorRamp2(c(color_range_top5_low, 0, 1), viridis(3))
  # col_fun_top10 = circlize::colorRamp2(c(color_range_top10_low, 0, 1), viridis(3))
  
  col_fun_all = circlize::colorRamp2(c(-1,0,2), c("#440154FF", "#238A8DFF", "#FDE725FF"))
  col_fun_top5 = circlize::colorRamp2(c(-1, 0, 2), c("#440154FF", "#238A8DFF", "#FDE725FF"))
  col_fun_top10 = circlize::colorRamp2(c(-1, 0, 2), c("#440154FF", "#238A8DFF", "#FDE725FF"))
  
  
  p1 = Heatmap(all_mat, name = "Expression",
               column_split = factor(cluster_id),
               top_annotation = column_ha_cluster,
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               show_column_names = FALSE,
               cluster_column_slices = TRUE,
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               column_title_gp = gpar(fontsize = 6),
               column_title_rot = 90,
               row_km = 20,
               column_gap = unit(0.5, "mm"),
               border = NA,
               row_names_gp = gpar(fontsize = 5),
               col = col_fun_all,
               use_raster = TRUE,
               raster_by_magick = FALSE,
               raster_quality = 8)
  
  p2 = Heatmap(top_5_mat, name = "Expression",
               column_split = factor(cluster_id),
               top_annotation = column_ha_cluster,
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               show_column_names = FALSE,
               cluster_column_slices = TRUE,
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               column_title_gp = gpar(fontsize = 5),
               column_gap = unit(0.3, "mm"),
               border = NA,
               row_names_gp = gpar(fontsize = 5),
               col = col_fun_top5,
               use_raster = TRUE,
               raster_by_magick = FALSE,
               raster_quality = 8)
  
  p3 = Heatmap(top_10_mat, name = "Expression",
               column_split = factor(cluster_id),
               top_annotation = column_ha_cluster,
               cluster_columns = TRUE,
               show_column_dend = FALSE,
               show_column_names = FALSE,
               cluster_column_slices = TRUE,
               cluster_rows = TRUE,
               show_row_dend = FALSE,
               column_title_gp = gpar(fontsize = 5),
               column_gap = unit(0.3, "mm"),
               border = NA,
               row_names_gp = gpar(fontsize = 5),
               col = col_fun_top10,
               use_raster = TRUE,
               raster_by_magick = FALSE,
               raster_quality = 8)
  return(list(p1 = p1, p2 = p2, p3 = p3))
}


#########
# Anova #
#########

get_cell_number_between_age<- function(obj, space, vjust = 3, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_age_anova.pdf")),
         plot = p1, width = 10, height = 30, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =percent)) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_age_anova.pdf")),
         plot = p2, width = 10, height = 30, limitsize = FALSE)
}


get_cell_number_between_genotype<- function(obj, space, vjust= 3, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y =n)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_genotype_anova.pdf")),
         plot = p1, width = 10, height = 30, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y =percent)) +
    scale_y_continuous(labels = scales::percent) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)), vjust = vjust) +
    theme_bw()
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_genotype_anova.pdf")),
         plot = p2, width = 10, height = 30, limitsize = FALSE)
}


##############################
# Custom Comparison with pva #
##############################
###############
# By Genotype #
###############

genotype_comparisons <- list(c("Bcl6tg/+", "CD70-/-;Bcl6tg/+"),
                             c("WT", "Bcl6tg/+"),
                             c("CD70-/-", "CD70-/-;Bcl6tg/+"),
                             c("WT", "CD70-/-;Bcl6tg/+"))

get_cell_number_between_genotype_custom<- function(obj, space, vjust = 0, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = n)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Genotypes", y = paste0("Number of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_genotype_wilcox.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = percent)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Genotypes", y = paste0("Percent of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_genotype_wilcox.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}


##########
# By Age #
##########
age_comparisons <- list( c("6mos", "14mos"),
                         c("14mos", "18mos"),
                         c("6mos", "18mos"))

get_cell_number_between_age_custom<- function(obj, space, vjust = 0, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Age", y = paste0("Number of cells per cluster"))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_age_wilcox.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y = percent)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", vjust = vjust) +
    labs(x = "Age", y = paste0("Percent of cells per cluster"))
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_age_wilcox.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}

##############################
# Custom Comparison with sig #
##############################

###############
# By Genotype #
###############

genotype_comparisons <- list(c("Bcl6tg/+", "CD70-/-;Bcl6tg/+"),
                             c("WT", "Bcl6tg/+"),
                             c("CD70-/-", "CD70-/-;Bcl6tg/+"),
                             c("WT", "CD70-/-;Bcl6tg/+"))

get_cell_number_between_genotype_custom_sig<- function(obj, space, vjust = 0.5, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = n)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Genotypes", y = paste0("Number of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_genotype_wilcox_symbol.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = genotype, y = percent)) +
    geom_boxplot(aes(fill = genotype), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ age, scale= "free_y") +
    stat_compare_means(comparisons = genotype_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Genotypes", y = paste0("Percent of cells per cluster "))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_genotype_wilcox_symbol.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}


##########
# By Age #
##########
age_comparisons <- list( c("6mos", "14mos"),
                         c("14mos", "18mos"),
                         c("6mos", "18mos"))

get_cell_number_between_age_custom_sig<- function(obj, space, vjust = 0.5, seurat_clusters){
  number_of_cells_per_cluster<- obj@meta.data %>%
    count(sample_id, age, genotype, seurat_clusters, .drop = FALSE)
  
  df<- left_join(number_of_cells_per_cluster, total_number_of_cells_per_sample) %>%
    mutate(percent = n/total_number)
  
  p1<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y =n)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Age", y = paste0("Number of cells per cluster"))
  
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_number_of_cells_per_cluster_between_age_wilcox_symbol.pdf")),
         plot = p1, width = 10, height = 35, limitsize = FALSE)
  
  p2<- df %>%
    mutate(age = factor(age, levels = c("6mos", "14mos", "18mos", "sick"))) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "CD70-/-", "Bcl6tg/+", "CD70-/-;Bcl6tg/+")))  %>%
    ggplot(aes(x = age, y = percent)) +
    geom_boxplot(aes(fill = age), outlier.colour = "NA") +
    geom_jitter(width = 0.2) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(seurat_clusters ~ genotype, scale= "free_y") +
    stat_compare_means(comparisons = age_comparisons, method= "wilcox.test", label = "p.signif",
                       vjust = vjust,
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.5),
                                          symbols = c("****", "***", "**", "*", "nsig"))) +
    labs(x = "Age", y = paste0("Percent of cells per cluster"))
  
  ggsave(filename = here(glue::glue("results/figures/{space}/{space}_cluster_removed_percent_of_cells_per_cluster_between_age_wilcox_symbol.pdf")),
         plot = p2, width = 10, height = 35, limitsize = FALSE)
}

# Extract scaled matrix
GetMatrixFromSeurat<- function(obj, features, group, ...) {
    p<- DotPlot(object = obj, features = features, group.by = group, ...)
    df<- p$data
  
    exp_mat<-df %>% 
        select(-pct.exp, -avg.exp) %>%  
        arrange(id) %>%
        pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
        as.data.frame() 

    row.names(exp_mat) <- exp_mat$features.plot  
    exp_mat <- exp_mat[,-1] %>% as.matrix()
    
    percent_mat<-df %>% 
        select(-avg.exp, -avg.exp.scaled) %>%  
        arrange(id) %>%
        pivot_wider(names_from = id, values_from = pct.exp) %>% 
        as.data.frame()
  
   row.names(percent_mat) <- percent_mat$features.plot  
   percent_mat <- percent_mat[,-1] %>% as.matrix()
   
   percent_mat <- percent_mat[complete.cases(exp_mat), ]
   exp_mat <- exp_mat[complete.cases(exp_mat), ] 


   if (!identical(dim(exp_mat), dim(percent_mat))) {
     stop("the dimension of the two matrice should be the same!")
   }
  
   if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
     stop("column names of the two matrice should be the same!")
   }
  
   if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
     stop("column names of the two matrice should be the same!")
   }
  
   return(list(exp_mat = exp_mat, percent_mat = percent_mat))
  
}

# Extract unscaled matrix
GetMatrixFromSeurat2<- function (obj, features, group, ...) 
{
  p <- DotPlot(object = obj, features = features, group.by = group, 
               ...)
  df <- p$data
  exp_mat <- df %>% select(-pct.exp, -avg.exp.scaled) %>% arrange(id) %>% 
    pivot_wider(names_from = id, values_from = avg.exp) %>% 
    as.data.frame()
  row.names(exp_mat) <- exp_mat$features.plot
  exp_mat <- exp_mat[, -1] %>% as.matrix()
  percent_mat <- df %>% select(-avg.exp, -avg.exp.scaled) %>% 
    arrange(id) %>% pivot_wider(names_from = id, values_from = pct.exp) %>% 
    as.data.frame()
  row.names(percent_mat) <- percent_mat$features.plot
  percent_mat <- percent_mat[, -1] %>% as.matrix()
  percent_mat <- percent_mat[complete.cases(exp_mat), ]
  exp_mat <- exp_mat[complete.cases(exp_mat), ]
  print(ncol(exp_mat))
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  if (!all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  if (!all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}

MakeClusterDotPlot<- function(exp_mat, percent_mat, col_fun, legend_title = "expression",
                              column_title = "Clustered Dotplot", row_km = 4, 
                              row_names_font_size = 3, column_ha, cluster_id = NULL,
                              legend_dot_size = unit(2,"mm"), cluster_rows=TRUE, cluster_columns = TRUE, ...){
  
  if (!is.null(cluster_id)){
    exp_mat <- exp_mat[,cluster_id]
    colnames(exp_mat) <- gsub(".*_", "", cluster_id)
    
    percent_mat <- percent_mat[,cluster_id]
    colnames(percent_mat) <- gsub(".*_", "", cluster_id)
  }
    
    layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = NA, fill = NA))
    grid.circle(x=x,y=y,r= sqrt(pindex(percent_mat, i, j)/100) * legend_dot_size,
                gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}
  
  hp<- Heatmap(exp_mat,
               heatmap_legend_param=list(title= legend_title),
               column_title = column_title, 
               col=col_fun,
               cluster_rows = cluster_rows,
               cluster_columns = cluster_columns,
               rect_gp = gpar(type = "none"),
               layer_fun = layer_fun,
               row_names_gp = gpar(fontsize = row_names_font_size),
               row_km = row_km ,
               border = "black",
               top_annotation = column_ha,
               ...)
  # see https://github.com/jokergoo/ComplexHeatmap/issues/737 for retain k-means
  
  return(hp)
}


