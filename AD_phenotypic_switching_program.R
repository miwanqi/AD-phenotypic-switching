
rm(list = ls())
gc()


dir <- "D:\\3_subcluster"
dir_data <- "D:\\3_亚聚类\\20231024_亚群rds"
cells <- c('SMC','Fib','EC','Macro','Mo','Neut','DC','Mast','B','T_NK')



i=1
cell <- cells[i]
dir_cell <- paste0(dir,"\\",cell)
dir.create(dir_cell)
setwd(dir_cell)
scRNA_sub <- readRDS(paste0(dir_data,"\\",cells[i],"_subcluster.rds"))

###____修改细胞亚群名称####
cluster2celltype <- c('SMC_1' = "ModSMC",
                      'SMC_2' = "Pro_SynSMC_1",
                      'SMC_3' = 'ConSMC',
                      'SMC_4' = "Pro_SynSMC_2",
                      'SMC_5' = 'SynSMC1',
                      'SMC_6' = 'SynSMC2')
scRNA_sub[['cluster']] = unname(cluster2celltype[scRNA_sub@meta.data$cluster])
scRNA_sub <- RenameIdents(scRNA_sub,cluster2celltype)
scRNA_sub$cluster <- factor(scRNA_sub$cluster,levels = c("ConSMC","ModSMC","Pro_SynSMC_1","Pro_SynSMC_2","SynSMC1","SynSMC2"))
scRNA_sub@active.ident <- factor(scRNA_sub@active.ident,levels = c("ConSMC","ModSMC","Pro_SynSMC_1","Pro_SynSMC_2","SynSMC1","SynSMC2"))
saveRDS(scRNA_sub,file = paste0(dir_data,"\\",cells[i],"_subcluster.rds"))


####____分化####

library(monocle)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(ggsignif)
library(patchwork)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
library(ClusterGVis)
##########
cds <- readRDS('7_SMC_monocle2.pseudotime.cds')

monocle_time <- as.data.frame(cds@phenoData@data)
monocle_time$cell <- rownames(monocle_time)
monocle_time <- monocle_time[,c(23,21)]
colnames(monocle_time)[2] <- "Monocle2_Pseudotime"
write.table(monocle_time,"7_monocle_time.txt",sep = "\t",row.names = F,col.names = T)


##########

library(destiny) # 加载 destiny...
library(Biobase)
library(grid)
library(rgl)

raw_data <- as.matrix(GetAssayData(scRNA_sub))
anno_mat <- scRNA_sub@meta.data
set <- ExpressionSet(raw_data, phenoData=AnnotatedDataFrame(anno_mat))
dmap <- DiffusionMap(set)
cor<-eigenvectors(dmap)[,1:3]
scRNA_sub_start <- subset(scRNA_sub,idents ="ConSMC" )
cor_V<-cor[colnames(scRNA_sub_start),]
cor_V_s <- cor_V[which(cor_V[,1] > -0.015 & cor_V[,1] < -0.01 & cor_V[,2] < -0.02  &  cor_V[,3] > 0.02),]

dpt <- DPT(dmap,tips = which(colnames(set) == "A1C_CACTTCAGACCTACAAGGGGCGGATAA"))

save(dmap,file = "7_DC_dmap.RData")
save(dpt,file = "7_DC_dpt.RData")

dc_data <- plot(dpt,dcs = c(1,-2), pch = 20)
dc_time <- dc_data$data
dc_time$cell <- rownames(dc_time)
dc_time <- dc_time[,c(6,3)]
colnames(dc_time)[2] <- "DC_Pseudotime"
write.table(dc_time,"7_DC_time.txt",sep = "\t",row.names = F,col.names = T)


#########
library(SCORPIUS)

exprData <- t(as.matrix(scRNA_sub@assays$RNA@data))
cluster <- scRNA_sub$cluster


space <- reduce_dimensionality(exprData, "spearman")
draw_trajectory_plot(space, cluster, contour = TRUE)
traj <- infer_trajectory(space)

gimp <- gene_importances(
  exprData, 
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 
gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(exprData[,gene_sel])

time <- traj$time

save(list = c('space','traj','gimp','time'),file = "7_SCORPIUS.RData")

SCORPIUS_time <- as.data.frame(time)
SCORPIUS_time$cell <- rownames(SCORPIUS_time)
SCORPIUS_time <- SCORPIUS_time[,c(2,1)]
colnames(SCORPIUS_time)[2] <- "SCORPIUS_Pseudotime"
write.table(SCORPIUS_time,"7_SCORPIUS_time.txt",sep = "\t",row.names = F,col.names = T)


###########
install.packages("qs")
scRNA_sub <- Seurat::NormalizeData(scRNA_sub,verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000) %>%
  ScaleData(verbose = FALSE)
scale.data <- scRNA_sub@assays$RNA@scale.data
scale.gene <- rownames(scale.data)
counts <- scRNA_sub@assays$RNA@counts
counts <- counts[scale.gene,]

sim <- SingleCellExperiment(assays = List(counts = counts)) 

# umap reduction
umap = scRNA_sub@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
reducedDims(sim) = SimpleList(UMAP = umap)

# metadata
meta = scRNA_sub@meta.data
# colData(sim)相当于meta.data，但他不是data.frame格式
# 所以需要分步赋予
colData(sim)$sampleId = meta$orig.ident
colData(sim)$cluster = meta$cluster

slingshot <- slingshot(sim, 
                 clusterLabels = 'cluster',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 start.clus= "ConSMC",  # 选择起点,"ConSMC"
                 end.clus = NULL     # 这里我们不指定终点
)     
colnames(colData(slingshot))
saveRDS(slingshot,"7_slingshot.rds")

lin1 <- getLineages(slingshot, 
                    clusterLabels = "cluster", 
                    start.clus = 'ConSMC',#可指定起始细胞簇
                    #end.clus=c("Comp","Col15a1","Ccl19", "Coch", "Cxcl12", "Fbln1", "Bmp4", "Npnt", "Hhip"),#可指定终点细胞簇
                    reducedDim = "UMAP")
plot(reducedDims(slingshot)$UMAP,col = brewer.pal(10,'Paired')[scRNA_sub$cluster],pch=16,asp=1)
lines(SlingshotDataSet(lin1), lwd=2,col = 'black',type = 'lineages')

slingshot_time <- as.data.frame(slingshot@colData@listData[["slingshot"]]@assays@data@listData[["pseudotime"]])

slingshot_time$cell <- rownames(slingshot_time)
slingshot_time <- slingshot_time[,c(3,1,2)]
colnames(slingshot_time) <- c("cell","Slingshot_Lineage1_Pseudotime","Slingshot_Lineage2_Pseudotime")
write.table(slingshot_time,"7_slingshot_time.txt",sep = "\t",row.names = F,col.names = T)


###############
monocle_time <- read.table("7_monocle_time.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
dc_time <- read.table("7_DC_time.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
SCORPIUS_time <- read.table("7_SCORPIUS_time.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)

slingshot_time <- read.table("7_slingshot_time.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)


phenotype <- as.data.frame(scRNA_sub$cluster)
phenotype$cell <- rownames(phenotype)
colnames(phenotype)[1] <- "phenotype"

times <- list(monocle_time,dc_time, SCORPIUS_time, slingshot_time,phenotype) %>%
  purrr::reduce(full_join, by = "cell")
rownames(times) <- times$cell
times <- times[,-1]



####________轨迹分析工作流封装函数####
trajectory_analysis_workflow <- function(dt, trajectories, 
                                         trajectory_threshold = trajectory_threshold,
                                         pseudotime_trim_quantile = 0.1,
                                         color_sub = color_sub,
                                         seurat_obj = seurat_obj,
                                         cell_order = cell_order,
                                         start_node = start_node,  # 新增起始节点参数
                                         plot_result = TRUE) {
  
  # 加载必要的程序包
  library(igraph)    # 图论分析
  library(dplyr)     # 数据操作
  library(tidyr)     # 数据整理
  library(reshape2)  # 数据重塑
  library(ggplot2)
  library(bezier)
  
  
  
  # 1. 伪时间值整合模块
  
  ## 使所有方法获取的伪时间分布在[0,1]之间
  normalize_pseudotime <- function(x) {
    (rank(x, na.last = "keep") - 1) / (sum(!is.na(x)) - 1)
  }
  
  ##  跨方法归一化
  dt_norm <- dt %>% mutate(across(1:(ncol(dt)-1), normalize_pseudotime))
  
  # 鲁棒性整合伪时间
  ## 抗离群值：排除1个最小和1个最大观测;容错处理：自动忽略NA值，当有效值≥3时仍可计算;保留信息：使用中间50%数据的均值
  dt_integrated <- dt_norm %>%
    rowwise() %>%  ## 逐行计算模式
    mutate(IntegratedTime = {
      # 提取方法的值（列序号1-5）
      method_values <- unlist(cur_data()[1:(ncol(dt)-1)])
      valid_values <- method_values[!is.na(method_values)]
      n_valid <- length(valid_values)
      
      if (n_valid >= 3) {
        sorted <- sort(valid_values)
        # 截取中间50%数据（排除首尾各1个）
        trimmed <- sorted[-c(1, n_valid)]  
        mean(trimmed)
      } else if (n_valid >= 1) {
        # 有效值不足时直接取均值
        mean(valid_values)
      } else {
        NA_real_
      }
    }) %>%  ##生成每个细胞（行）的共识伪时间,对5个方法的归一化值排序,截取排序后的第2-4个元素（排除最小/最大值）,计算截尾均值
    ungroup()
  rownames(dt_integrated) <- rownames(dt_norm)
  
  
  # 2. 轨迹结构整合模块
  build_transition_matrix <- function(trajs) {
    # 解包所有轨迹中的边
    ## 遍历每个方法的轨迹数据
    all_edges <- unlist(lapply(trajs, function(single_method_trajs){
      unlist(lapply(single_method_trajs, function(x){## 处理单条轨迹字符串
        paths <- unlist(strsplit(x, ">")) ## 分割节点
        if(length(paths) < 2) return(character(0)) ## 过滤无效轨迹
        sapply(1:(length(paths)-1), function(i) paste0(paths[i], ">", paths[i+1])) ## 生成连续节点对作为有向边
      }))
    }))
    
    
    # 构建频率统计表
    freq_table <- if(length(all_edges) > 0) {
      data.frame(edge = all_edges) %>% 
        separate(edge, c("from", "to"), sep = ">",extra = "drop",fill = "right")
    } else {
      data.frame(from = character(), to = character())
    }
    
    ## 计算有向线段出现的次数
    freq_table_count <- as.data.frame(table(freq_table$from,freq_table$to))
    colnames(freq_table_count) <- c(colnames(freq_table),"n")
    
    # 补全所有可能的节点组合
    full_grid <- expand.grid(from = cell_order, to =cell_order)
    freq_table_count <- full_grid %>%
      left_join(freq_table_count, by = c("from", "to")) %>%
      mutate(n = coalesce(n, 0L))
    
    
    
    # 转换为邻接矩阵
    freq_matrix <- reshape2::acast(freq_table_count, from ~ to, drop = FALSE)
    return(freq_matrix)
  }
  
  
  #留一法
  n_methods <- length(trajectories)
  # 生成留一法组合（每次排除一个方法）
  loo_combinations <- map(1:n_methods, ~trajectories[-.x])
  
  # 生成每个LOO组合的转移矩阵
  matrix_list <- map(loo_combinations, function(comb) {
    # 计算当前组合的转移矩阵
    comb_matrix <- comb %>% 
      map(~ build_transition_matrix(.x)) %>%  # 每个方法生成矩阵
      reduce(`+`)                             # 累加当前组合的矩阵
    
    # 在循环内过滤低频边（阈值设为方法数的1/2）
    local_threshold <- ceiling(n_methods * 0.5)
    filtered_matrix <- ifelse(comb_matrix >= local_threshold, 1, 0)
    
    filtered_matrix
  })
  
  
  # 汇总所有LOO矩阵
  combined_matrix <- reduce(matrix_list, `+`)
  
  # 整合所有方法的转移邻接矩阵
  #combined_matrix <- lapply(trajectories, build_transition_matrix) %>% Reduce(`+`, .)
  
  # 过滤低频转移边
  filtered_matrix <- ifelse(combined_matrix >= trajectory_threshold, 1, 0)
  
  # 构建共识轨迹图
  consensus_graph <- graph_from_adjacency_matrix(
    filtered_matrix, 
    mode = "directed",
    diag = FALSE
  )
  
  
  # 3. 轨迹排序优化模块
  # 环路检测与处理
  ## 输入各个Method的轨迹（如Method1的B->E->C->D->A）
  ## cell_order应为生物学意义的正确顺序（如干细胞到终末分化细胞的路径）
  handle_cycles <- function(graph, cell_order) {
    if (!is_dag(graph)) { ## 验证图是否为有向无环图
      # 方案A：移除形成环路的边
      ## 删除低权重边
      fas <- feedback_arc_set(graph, weights = 1/E(graph)$weight)
      acyclic_graph <- delete_edges(graph, fas)
      
      # 若仍存在环路，启用方案B（生物学约束）
      ## 过滤不在cell_order中的节点
      if (!is_dag(acyclic_graph)) {
        edge_ends <- ends(consensus_graph, E(consensus_graph), names = TRUE)
        from_nodes <- edge_ends[,1]
        to_nodes <- edge_ends[,2]
        
        valid_edges <- (from_nodes %in% cell_order) & (to_nodes %in% cell_order)
        from_idx <- match(from_nodes[valid_edges], cell_order)
        to_idx <- match(to_nodes[valid_edges], cell_order)
        
        good_edges <- (to_idx > from_idx) & !is.na(from_idx) & !is.na(to_idx)
        edge_ids <- which(valid_edges)[good_edges]
        
        acyclic_graph <- subgraph.edges(graph, edge_ids, delete.vertices = FALSE)
        
        if (!is_dag(acyclic_graph)) stop("无法移除环路")
      }
      return(acyclic_graph)
    } else {
      return(graph)
    }
  }
  
  acyclic_graph <- handle_cycles(consensus_graph, cell_order)
  
 
  
  
  # 4. 拓扑排序
  
  # 拓扑排序获取全局顺序
  topological_order <- if (start_node %in% V(acyclic_graph)$name) {
    c(start_node, setdiff(topo_sort(acyclic_graph)$name, start_node))
  } else {
    topo_sort(acyclic_graph)$name
  }
  
  
  
  # 按细胞类型整合排序
  dt_ordered <- data.frame()
  for (ct in topological_order) {
    subset_cells <- dt_integrated %>%
      filter(phenotype == ct) %>%
      arrange(IntegratedTime) ##根据整合后伪时间对数据帧的行进行排序
    dt_ordered <- bind_rows(dt_ordered, subset_cells)
  }
  
  # 4. 结果验证模块
  result <- list(
    consensus_graph = acyclic_graph,
    ordered_data = dt_ordered,
    combined_matrix = combined_matrix,
    topological_order = topological_order,
    final_trajectory = paste(unique(dt_ordered$phenotype), collapse = ">")
  )
  seurat_obj$integrated_time <- dt_integrated$IntegratedTime[match(rownames(seurat_obj@meta.data), rownames(dt_integrated))]
  
  # 5. 结果可视化
  if (plot_result && !is.null(seurat_obj)) {
    # 获取轨迹图的边信息
    edge_list <- as_edgelist(acyclic_graph)
    
    
    # 生成可视化数据
    umap_data <- Embeddings(seurat_obj, "umap") %>%
      as.data.frame() %>%
      mutate(
        cell_id = rownames(seurat_obj@meta.data),  # 关键修正：添加细胞ID
        phenotype = seurat_obj$cluster,
        pseudotime = seurat_obj$integrated_time
      )
    
    
    # 计算表型中心坐标（使用加权中位数）
    phenotype_centroids <- umap_data %>%
      group_by(phenotype) %>%
      summarise(
        UMAP_1_centroid = median(UMAP_1[phenotype %in% V(acyclic_graph)$name]),  # 仅包含轨迹细胞
        UMAP_2_centroid = median(UMAP_2[phenotype %in% V(acyclic_graph)$name]),
        .groups = 'drop'
      )
    
    
    
    # 生成轨迹曲线数据
    trajectory_curves <- map(1:nrow(edge_list), function(i) {
      from_node <- edge_list[i, 1]
      to_node <- edge_list[i, 2]
      
      # 获取端点对应表型（假设节点名为细胞ID）
      from_phenotype <- as.character(umap_data$phenotype[umap_data$phenotype == from_node]%>% na.omit() %>% first())
      to_phenotype <- as.character(umap_data$phenotype[umap_data$phenotype == to_node]%>% na.omit() %>% first())
      
      # 获取标准化坐标
      from_pt <- phenotype_centroids %>% 
        filter(phenotype == from_phenotype) %>%
        select(UMAP_1_centroid, UMAP_2_centroid) %>%
        unlist()
      
      to_pt <- phenotype_centroids %>% 
        filter(phenotype == to_phenotype) %>%
        select(UMAP_1_centroid, UMAP_2_centroid) %>%
        unlist()
      
      # 生成分段贝塞尔曲线（含多个箭头）
      bezier_df <- data.frame(
        bezier::bezier(
          t = seq(0, 1, length.out = 30),  # 增加路径点密度
          p = rbind(
            from_pt, 
            from_pt + (to_pt - from_pt)*0.3,  # 动态调整控制点
            to_pt
          )
        )
      ) %>% setNames(c("x", "y"))
      
      
      # 生成箭头标记位置
      arrow_positions <- seq(0.3, 0.9, by = 0.2)  # 沿路径添加多个箭头
      arrow_df <- map_df(arrow_positions, function(pos) {
        idx <- floor(pos * nrow(bezier_df))
        data.frame(
          x = bezier_df$x[idx],
          y = bezier_df$y[idx],
          xend = bezier_df$x[idx + 2],  # 方向向量
          yend = bezier_df$y[idx + 2]
        )
      })
      
      list(
        curve_data = tibble(
          path_id = i,
          x = bezier_df$x,
          y = bezier_df$y
        ),
        arrow_data = tibble(
          path_id = rep(i, nrow(arrow_df)),  # 确保与箭头数匹配
          x = arrow_df$x,
          y = arrow_df$y,
          xend = arrow_df$xend,
          yend = arrow_df$yend
        )
      )
      
    }) 
    
    # 分别提取曲线和箭头数据
    curve_df <- map_dfr(trajectory_curves, ~ .x$curve_data)
    arrow_df <- map_dfr(trajectory_curves, ~ .x$arrow_data)
    
    
    #生成权重标签数据
    weight_labels <- map_dfr(1:nrow(edge_list), function(i) {
      from <- edge_list[i, 1]
      to <- edge_list[i, 2]
      weight <- combined_matrix[from, to]
      
      # 获取曲线中段点作为标签位置
      curve_sub <- curve_df %>% filter(path_id == i)
      # 安全索引计算
      n_points <- nrow(curve_sub)
      mid_idx <- max(ceiling(n_points * 0.4), 2)  # 确保至少第二个点
      prev_idx <- max(mid_idx - 1, 1)             # 防止负索引
      next_idx <- min(mid_idx + 1, n_points)      # 防止越界
      
      # 生成标签位置数据
      curve_sub %>% 
        slice(mid_idx) %>% 
        mutate(
          x_prev = curve_sub$x[prev_idx],
          y_prev = curve_sub$y[prev_idx],
          x_next = curve_sub$x[next_idx],
          y_next = curve_sub$y[next_idx],
          angle = atan2(y_next - y_prev, x_next - x_prev) * 180/pi,
          weight_label = sprintf("%d", weight)
        ) %>% 
        select(x, y, angle, weight_label)
      
    })
    
    
    # 分图层绘制
    p <- ggplot(umap_data, aes(UMAP_1, UMAP_2)) +
      geom_point(aes(color = pseudotime), size = 0.5, alpha = 0.8) +
      
      # 绘制主轨迹曲线
      geom_path(
        data = curve_df,
        aes(x, y, group = path_id),
        color = "grey40",
        linewidth = 0.8,
        alpha = 0.9
      ) +
      
      # 绘制动态箭头
      geom_segment(
        data = arrow_df,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "#2B4050",
        linewidth = 0.5,
        arrow = arrow(
          angle = 30,
          length = unit(0.2, "cm"),
          type = "closed"
        )
      ) +
      
      # 增强表型标签
      geom_label(
        data = phenotype_centroids,
        aes(UMAP_1_centroid, UMAP_2_centroid, 
            label = str_replace(phenotype, "_", "\n"),  # 自动换行
            fill = phenotype),
        color = "white",
        fontface = "bold",
        size = 4.5,
        label.padding = unit(0.3, "lines"),
        label.r = unit(0.15, "lines"),
        alpha = 0.9
      ) +
      # 权重标签
      geom_text(
        data = weight_labels,
        aes(x = x, y = y, label = weight_label),
        color = "red",
        size = 4.5,
        hjust = 0.5,
        vjust = 0.5
      ) +
      scale_color_viridis()+
      scale_fill_manual(values = color_sub) +
      coord_fixed(ratio = 1,expand = FALSE) +
      theme_classic(base_size = 14) +
      labs(
        title = "Integrated Trajectory in UMAP Space",
        color = "Pseudotime",
        x = "UMAP-1", 
        y = "UMAP-2"
      ) +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right",
        legend.box = "horizontal",
        panel.border = element_rect(fill = NA, color = "grey70")
      )
    
    print(p)
    ggsave("7_integrated_trajectory_umap.pdf", p, width = 10, height = 8)
    
  }
  
  return(result)
}



########
# 示例用法
# 加载必要的程序包
library(igraph)
library(dplyr)


trajectories <- list(
  Method1 = c("ConSMC>ModSMC>SynSMC2", "ConSMC>SynSMC1>Pro_SynSMC_2>Pro_SynSMC_1"),
  Method2 = c("ConSMC>ModSMC","ConSMC>SynSMC1","ConSMC>SynSMC2","ConSMC>Pro_SynSMC_1","ConSMC>Pro_SynSMC_2","ConSMC>ModSMC"),
  Method3 = c("ConSMC>ModSMC","ConSMC>SynSMC1>SynSMC2","ConSMC>Pro_SynSMC_1>Pro_SynSMC_2"),
  Method4 = c("SynSMC2>SynSMC1>ConSMC>ModSMC>Pro_SynSMC_2>Pro_SynSMC_1"),
  Method5 = c("ConSMC>Pro_SynSMC_1>Pro_SynSMC_2", "ConSMC>SynSMC1>SynSMC2","ConSMC>ModSMC>SynSMC2"),
  Method6 = c("ConSMC>SynSMC1","SynSMC2>ModSMC","SynSMC2>SynSMC1","SynSMC2>Pro_SynSMC_2>Pro_SynSMC_1")
)

# 调用函数，指定起始节点为"B"
result <- trajectory_analysis_workflow(times, trajectories, 
                                       start_node = "ConSMC",
                                       cell_order = c("ConSMC","ModSMC","Pro_SynSMC_1","Pro_SynSMC_2","SynSMC1", "SynSMC2"),
                                       trajectory_threshold = 3,
                                       color_sub = color_SMC,
                                       seurat_obj = scRNA_sub
                                       )
print(result$final_trajectory)

saveRDS(result,"7_intergrated_result.rds")
