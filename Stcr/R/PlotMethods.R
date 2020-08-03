#' Scatter plot of a binary meta.data on a reduced dimension 
#'
#' @param sset StcrSet object which contains cellName
#' @param group_by uses meta.data of cell from sset
#' @param cell_ind A vector of cells to be highlighted
#' @param reduction Which dimensionality reduction to use.
#' @param color1 primary color, default navy
#' @param color2 backgound cell color, default grey
#'
#' @export PlotTCR
#'
#' @examples
#' \dontrun{
#' plt = PlotTCR(sset, group_by = )
#' print(plt)
#' plt = PlotTCR(sset, group_by = )
#' print(plt)
#' }
PlotTCR <- function(sset,group_by = NULL,cell_ind = NULL,reduction = "umap",color1 ="navy",color2= "grey"){
  
  TCR = sset@cellName
  names(TCR) = sset@cellName
  if(!is.vector(cell_ind)){
    if(is.null(group_by)){
      stop("group_by parameter is missing. Use binary data for from sset@meta.data")
    } else {
      eval(parse(text = paste0("sset@meta.data$",group_by,"=as.factor(sset@meta.data$",group_by,")")))
      if(nlevels(eval(parse(text = paste0("sset@meta.data$",group_by)))) > 2){
        stop("data is not binary")
      } else {
        if(eval(parse(text = paste0("\"FALSE\" == (levels(sset@meta.data$",group_by,")[1])"))) || eval(parse(text = paste0("\"TRUE\" == (levels(sset@meta.data$",group_by,")[1])"))) ){
          cell_ind = rownames(eval(parse(text = paste0("subset(sset@meta.data, ",group_by,"==\"", levels(eval(parse(text = paste0("as.factor(sset@meta.data$",group_by,")"))))[2],"\")"))))
        } else {
          cell_ind = rownames(eval(parse(text = paste0("subset(sset@meta.data, ",group_by,"==\"", levels(eval(parse(text = paste0("as.factor(sset@meta.data$",group_by,")"))))[1],"\")"))))
        }
        TCR[sset@cellName %in% cell_ind] = "1"
        TCR[!(sset@cellName %in% cell_ind)] = "2"
      }
      
    }
  } else {
    TCR[sset@cellName %in% cell_ind] = "1"
    TCR[!(sset@cellName %in% cell_ind)] = "2"
  }

  sortedind=rev(order(TCR))
  if(reduction == "umap"){
    df = data.frame(x=sset@reductions$umap@cell.embeddings[sortedind,1],
                    y=sset@reductions$umap@cell.embeddings[sortedind,2],
                    col = TCR[sortedind])
  } else if(reduction == "tsne"){
    df = data.frame(x=sset@reductions$tsne@cell.embeddings[sortedind,1],
                    y=sset@reductions$tsne@cell.embeddings[sortedind,2],
                    col = TCR[sortedind])
  } else if(reduction == "pca"){
    df = data.frame(x=sset@reductions$pca@cell.embeddings[sortedind,1],
                    y=sset@reductions$pca@cell.embeddings[sortedind,2],
                    col = TCR[sortedind])
  } else {
    stop("define your reduced dimension")
  }
  
  plt = ggplot(df) + 
    geom_point(aes(x=x,y=y,colour=col),size=3) +
    scale_colour_manual(values=c(color1,color2)) +
    theme(
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      plot.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(), 
      axis.line = element_blank(),
      axis.ticks=element_blank(),
      legend.position = "none"
    )
  print(plt)
}


#' Plot of (non)canonical frequencies by cluster
#'
#' PlotFrequency takes either canonical_ind or gene regex string from chains and
#' returns a ggplot object of (non)canonical frequencies in each of the clusters.
#' All gene regex string are combinded with "AND" operation to assign canonical_ind.
#'
#' @param CDR3PairMat 
#' @param ident It can be cluster or celltype
#' @param reverse switch to cannonical / non-cannonical
#' @param selected_ind a vector of selected cell index
#' @param chain1_v_gene ex) c("TRAV11","TRAJ18"...) used to assign canonical and noncanonical
#' @param chain1_j_gene regex string for chain1_j_gene. default value is NULL. ex "TRAV11" 
#' @param chain1_d_gene regex string for chain1_d_gene. default value is NULL. ex "TRAJ18"
#' @param chain2_v_gene regex string for chain2_v_gene. default value is NULL. ex "(TRBV13|TRBV29|TRBV1\\*)"
#' @param chain2_j_gene regex string for chain2_j_gene. default value is NULL.
#' @param chain2_d_gene regex string for chain2_d_gene. default value is NULL.
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' 
#' test_cdr3_pair_mat = SetCDR3PairMat(sset_test, test_cdr3_mat, chain1 = c("TRA"), chain2 = c("TRB"))
#' PlotFrequency <- function(test_cdr3_pair_mat, ident = sset_test@ident , selected_ind= NULL,
#'  chain1_v_gene = NULL, chain1_d_gene = NULL, chain1_j_gene = NULL, 
#'  chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL)
#' If you have a list of canonical cell index
#' PlotFrequency <- function(test_cdr3_pair_mat, ident = sset_test@ident , canonical_ind= canonical_ind)
#'
#' }
PlotFrequency <- function(CDR3PairMat, ident = NULL, reverse = TRUE, selected_ind = NULL,
  chain1_v_gene = NULL, chain1_d_gene = NULL, chain1_j_gene = NULL, 
  chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL){
  mat_cdr3_pair = CDR3PairMat[[1]]
  df_clonotype = CDR3PairMat[[2]]
  ident = as.factor(ident)
  df_mat_cdr3_pair = data.frame(mat_cdr3_pair)
  df_mat_cdr3_pair$ident = ""
  for(i in levels(ident)){
    cell_ind = names(ident[ident %in% c(i)])
    df_mat_cdr3_pair$ident[rownames(df_mat_cdr3_pair) %in% cell_ind] = i
  }

  for(i in 1:nrow(df_clonotype)){
    df_clonotype$ident[i] =  paste( unique(df_mat_cdr3_pair[mat_cdr3_pair[,i] > 0,'ident']), collapse = ",")
  }


  df_clonotype2 = subset(df_clonotype, n > 0)
  df_clonotype2$selected = FALSE
  if(is.null(selected_ind)){
    if(is.null(chain1_v_gene) & is.null(chain1_d_gene) & is.null(chain1_j_gene) & 
        is.null(chain2_v_gene) & is.null(chain2_d_gene) & is.null(chain2_j_gene)){
      stop("neither selected_ind nor gene regex is not available")
      }

    combined_index = rep(TRUE,nrow(df_clonotype2))
    if(!is.null(chain1_v_gene)){
      grepl_chain1_v_gene = grepl(chain1_v_gene,df_clonotype2$chain1_v_gene)
      combined_index = combined_index & grepl_chain1_v_gene
    }
    if(!is.null(chain1_d_gene)){
      grepl_chain1_d_gene = grepl(chain1_d_gene,df_clonotype2$chain1_d_gene)
      combined_index = combined_index & grepl_chain1_d_gene
    }
    if(!is.null(chain1_j_gene)){
      grepl_chain1_j_gene = grepl(chain1_j_gene,df_clonotype2$chain1_j_gene)
      combined_index = combined_index & grepl_chain1_j_gene
    }
    if(!is.null(chain2_v_gene)){
      grepl_chain2_v_gene = grepl(chain2_v_gene,df_clonotype2$chain2_v_gene)
      combined_index = combined_index & grepl_chain2_v_gene
    }
    if(!is.null(chain2_d_gene)){
      grepl_chain2_d_gene = grepl(chain2_d_gene,df_clonotype2$chain2_d_gene)
      combined_index = combined_index & grepl_chain2_d_gene
    }
    if(!is.null(chain2_j_gene)){
      grepl_chain2_j_gene = grepl(chain2_j_gene,df_clonotype2$chain2_j_gene)
      combined_index = combined_index & grepl_chain2_j_gene
    }
    selected_ind = rownames(df_clonotype2)[combined_index]
  }

  df_clonotype2[selected_ind,]$selected = TRUE

  df_freq = ComputeFreq(mat_cdr3_pair,
                        df_clonotype2,
                        ident,
  chain1_v_gene = chain1_v_gene, chain1_d_gene = chain1_d_gene, chain1_j_gene = chain1_j_gene, 
  chain2_v_gene = chain2_v_gene, chain2_d_gene = chain2_d_gene, chain2_j_gene = chain2_j_gene)
  df_freq$cluster = factor(rownames(df_freq))
  if(!reverse){
    plt = ggplot(df_freq,aes(x=cluster,y=selected_pct*100)) +
      geom_point(size=3) +
      geom_line(group = 0) +
      xlab("") +
      ylab("% of non selected TCRs") +
      ylim((min(df_freq$selected_pct*100)-5),(max(df_freq$selected_pct*100)+5))
  } else {
    plt = ggplot(df_freq,aes(x=cluster,y=nonselected_pct*100)) +
      geom_point(size=3) +
      geom_line(group = 0) +
      xlab("") +
      ylab("% of selected TCRs") +
      ylim((min(df_freq$nonselected_pct*100)-5),(max(df_freq$nonselected_pct*100)+5))
  }
  print(plt)
}

#' Bar plot of frequencies by cluster
#'
#' Canonical/(Non)canoncial Frequency plot
#' 
#'
#' @param CDR3PairMat
#' @param ident It can be cluster or celltype
#' @param selected_ind a vector of canonical cell index
#' @param chain1_v_gene ex) c("TRAV11","TRAJ18"...) used to assign canonical and noncanonical
#' @param chain1_j_gene regex string for chain1_j_gene. default value is NULL. ex "TRAV11" 
#' @param chain1_d_gene regex string for chain1_d_gene. default value is NULL. ex "TRAJ18"
#' @param chain2_v_gene regex string for chain2_v_gene. default value is NULL. ex "(TRBV13|TRBV29|TRBV1\\*)"
#' @param chain2_j_gene regex string for chain2_j_gene. default value is NULL.
#' @param chain2_d_gene regex string for chain2_d_gene. default value is NULL.
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' sset_NKT = PlotTCR(sset, group.by = )
#' }
PlotProportionsByGene <- function(CDR3PairMat, ident = NULL, selected_ind = NULL,
  chain1_v_gene = NULL, chain1_d_gene = NULL, chain1_j_gene = NULL, 
  chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL){
  
  mat_cdr3_pair = CDR3PairMat[[1]]
  df_clonotype = CDR3PairMat[[2]]
  ident = as.factor(ident)
  df_mat_cdr3_pair = data.frame(mat_cdr3_pair)
  df_mat_cdr3_pair$ident = ""
  for(i in levels(ident)){
    cell_ind = names(ident[ident %in% c(i)])
    df_mat_cdr3_pair$ident[rownames(df_mat_cdr3_pair) %in% cell_ind] = i
  }

  for(i in 1:nrow(df_clonotype)){
    df_clonotype$ident[i] =  paste( unique(df_mat_cdr3_pair[mat_cdr3_pair[,i] > 0,'ident']), collapse = ",")
  }


  df_clonotype2 = subset(df_clonotype, n > 0)
  df_clonotype2$selected = FALSE
  if(is.null(selected_ind)){
    if(is.null(chain1_v_gene) & is.null(chain1_d_gene) & is.null(chain1_j_gene) & 
        is.null(chain2_v_gene) & is.null(chain2_d_gene) & is.null(chain2_j_gene)){
      stop("neither selected_ind nor gene regex is not available")
      }

    # selected
    combined_index = rep(TRUE,nrow(df_clonotype2))
    if(!is.null(chain1_v_gene)){
      grepl_chain1_v_gene = grepl(chain1_v_gene,df_clonotype2$chain1_v_gene)
      combined_index = combined_index & grepl_chain1_v_gene
    }
    if(!is.null(chain1_d_gene)){
      grepl_chain1_d_gene = grepl(chain1_d_gene,df_clonotype2$chain1_d_gene)
      combined_index = combined_index & grepl_chain1_d_gene
    }
    if(!is.null(chain1_j_gene)){
      grepl_chain1_j_gene = grepl(chain1_j_gene,df_clonotype2$chain1_j_gene)
      combined_index = combined_index & grepl_chain1_j_gene
    }
    if(!is.null(chain1_v_gene)){
      grepl_chain2_v_gene = grepl(chain2_v_gene,df_clonotype2$chain2_v_gene)
      combined_index = combined_index & grepl_chain2_v_gene
    }
    if(!is.null(chain2_d_gene)){
      grepl_chain2_d_gene = grepl(chain2_d_gene,df_clonotype2$chain2_d_gene)
      combined_index = combined_index & grepl_chain2_d_gene
    }
    if(!is.null(chain2_j_gene)){
      grepl_chain2_j_gene = grepl(chain2_j_gene,df_clonotype2$chain2_j_gene)
      combined_index = combined_index & grepl_chain2_j_gene
    }
    selected_ind = rownames(df_clonotype2)[combined_index]
  } 
  df_clonotype2[selected_ind,]$selected = TRUE

  df_freq = ComputeFreq(mat_cdr3_pair,
                        df_clonotype2,
                        ident,
  chain1_v_gene = chain1_v_gene, chain1_d_gene = chain1_d_gene, chain1_j_gene = chain1_j_gene, 
  chain2_v_gene = chain2_v_gene, chain2_d_gene = chain2_d_gene, chain2_j_gene = chain2_j_gene)
  df_freq$cluster = factor(rownames(df_freq))

  long_df_freq = data.frame(
    selected = melt(df_freq[,grepl("^selected_chain",colnames(df_freq))])$value,
    nonselected =  melt(df_freq[,grepl("^nonselected_chain",colnames(df_freq))])$value,
    cluster = rep(rownames(df_freq),ncol(df_freq[,grepl("^selected_chain",colnames(df_freq))])))
  
  long_long_df_freq = melt(long_df_freq, by = cluster)
  
  long_long_df_freq$TCR = rep(gsub("^selected_(chain.*)","\\1",colnames(df_freq)[grepl("^selected_chain",colnames(df_freq))]),2)
  names(long_long_df_freq)[names(long_long_df_freq) == "variable"] = "selected"
  long_long_long_df_freq = melt(long_long_df_freq, by=TCR)

  plt = ggplot(long_long_long_df_freq, aes(x=TCR, y=value, fill=selected)) +
    scale_fill_manual(values=c("#4d4d4d","#b2182b")) +
    geom_bar(stat="identity") +
    facet_wrap(~factor(cluster), nrow = 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(plt)
}

#' Plot of frequencies by cluster
#'
#' 
#'
#' @param CDR3PairMat
#' @param ident It can be either clusters or celltypes
#' @param cellLimit minimum paired cell number
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' sset_NKT = PlotTCR(sset, group.by = )
#' }
PlotClonotypes <- function(CDR3PairMat, ident, cellLimit = 1){
  mat_cdr3_pair = CDR3PairMat[[1]]
  df_clonotype = CDR3PairMat[[2]]
  ident = as.factor(ident)
  df_mat_cdr3_pair = data.frame(mat_cdr3_pair)

  cellcount_clonotypes = colSums(mat_cdr3_pair)
  df = data.frame(cellcount = cellcount_clonotypes, clonotypes = names(cellcount_clonotypes), stringsAsFactors = F)
  sub_df = subset(df, cellcount > cellLimit)
  sub_df = sub_df[rev(order(sub_df$cellcount)),]

    
  long_cellcount_clonotypes =  melt(mat_cdr3_pair[,rownames(sub_df)])
  colnames(long_cellcount_clonotypes) = c("cell","clonotype","cellcounts")
  long_cellcount_clonotypes = long_cellcount_clonotypes[long_cellcount_clonotypes$cellcounts > 0,]
  
  long_cellcount_clonotypes$cluster = ""
  for(i in 1:length(long_cellcount_clonotypes$cell)){
    long_cellcount_clonotypes$cluster[i] = ident[as.character(long_cellcount_clonotypes$cell[i])]
  }
  
  summary_df = long_cellcount_clonotypes %>% group_by(clonotype,cluster) %>% dplyr::count(cellcounts) %>% as.data.frame()


  plt = ggplot(summary_df, aes(x = factor(clonotype, levels=unique(summary_df$clonotype)), y =n, fill = factor(cluster))) +
  #scale_fill_manual(values=c("darkgreen","blue","darkred","black")) +
  geom_bar(stat="identity") + 
  ylab("cell counts") +
  xlab("clonotype") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

  print(plt)
}

#' Plot of frequencies by cluster
#'
#' 
#'
#' @param CDR3PairMat
#' @param ident It can be either clusters or celltypes
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' sset_NKT = PlotTCR(sset, group.by = )
#' }
PlotDiversity <- function(CDR3PairMat, ident = NULL){
  mat_cdr3_pair = CDR3PairMat[[1]]
  df_clonotype = CDR3PairMat[[2]]
  ident = as.factor(ident)
    
  clonotype_used = colnames(mat_cdr3_pair)
  mat_pseudotime_clonotype = matrix(data=0, ncol=nlevels(ident),
                                    nrow = length(clonotype_used))
  rownames(mat_pseudotime_clonotype) = clonotype_used
  colnames(mat_pseudotime_clonotype) = levels(ident)

  for(i in 1:nlevels(ident)){
    cluster_cell_ind = names(ident)[ident == levels(ident)[i]]
    
    c_mat_cdr3_pair = mat_cdr3_pair[rownames(mat_cdr3_pair) %in% cluster_cell_ind,]
    c_mat_cdr3_pair = c_mat_cdr3_pair[,clonotype_used]
    mat_pseudotime_clonotype[,levels(ident)[i]] = colSums(c_mat_cdr3_pair)
  }

  cluster_shannon_indices = apply(mat_pseudotime_clonotype, 2,ComputeShannonIndex)

  df_si = data.frame(cluster_shannon_indices)
  df_si$cluster = as.factor(rownames(df_si))
  plt = ggplot(df_si, aes(x=cluster,
                    y=cluster_shannon_indices,
                    group=1)) + geom_line() + geom_point(size=3)  +
    xlab("") +
    ylab("Normalized diversity Index") +
    ylim(low=0.4, high=1) +
    theme_bw() +
    theme(
      axis.text.x=element_text(size=20),
      axis.text.y=element_text(size=20),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=20),
      plot.background=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.x =element_blank()
    )
  print(plt)
}


#pair log2FC TCR 
#' Plot of PairUsage 
#'
#' 
#'
#' @param mat_1
#' @param mat_2
#' @param cell_ind
#' @param pseudocount
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' sset_NKT = PlotTCR(sset, group.by = )
#' }
PlotPairUsageHeatmap <- function(mat_1 = NULL, mat_2 = NULL, cell_ind = NULL, pseudocount = 0.1){
  if(is.null(cell_ind)){
    stop("cell_ind is null")
  }
  row_cells = rownames(mat_1)
  
  mat_12= ((t(mat_1) %*% mat_2))/length(unique(row_cells))

  cluster_cell_ind = cell_ind
  num_cluster_cell = length(cell_ind)
  c_mat_1 = mat_1[rownames(mat_1) %in% cluster_cell_ind,]
  c_mat_2 = mat_2[rownames(mat_2) %in% cluster_cell_ind,]
  c_mat_12 = ((t(c_mat_1) %*% c_mat_2))/num_cluster_cell

  log2FC_c_mat_12 = log2((c_mat_12+pseudocount)/(mat_12+pseudocount))
  paletteLength <- 100
  if(max(log2FC_c_mat_12) > 0 ){
    myBreaks <- c(seq(min(log2FC_c_mat_12), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(log2FC_c_mat_12)/paletteLength, max(log2FC_c_mat_12), length.out=floor(paletteLength/2)))
    ph =pheatmap(log2FC_c_mat_12,breaks = myBreaks,cluster_cols = F,cluster_rows = F, color =colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),)
  } else {
    myBreaks <- c(seq(min(log2FC_c_mat_12), 0, length.out=ceiling(paletteLength/2) + 1))
    ph = pheatmap(log2FC_c_mat_12,breaks = myBreaks,cluster_cols = F,cluster_rows = F, color =colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)[1:50],)
  }
  print(ph)
}

#geomean, logFC by cluster
#' Plot of frequencies by cluster
#'
#'
#'
#' @param CDR3PairMat
#' @param ident It can be cluster or celltype
#' @param compute geomean or logFC
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' sset_NKT = PlotTCR(sset, group.by = )
#' }
PlotOverlapHeatmap <- function(CDR3PairMat, ident, compute = "geomean"){

  if(compute == "geomean"){
    geo_mean = ComputeGeomean(CDR3PairMat, ident)
    mat_breaks <- seq(0, (max(geo_mean)+0.02), length.out = 50)
    plt = pheatmap(geo_mean, color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name =
        "RdBu")))(100)[51:100], breaks = mat_breaks,
        cluster_rows = F, cluster_cols = F, angle_col = c("0"))
  } else if(compute == "logFC"){
    mat_log2FC = ComputeLogFC(CDR3PairMat, ident)
    paletteLength <- 100
    myColor <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100)
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(mat_log2FC), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(mat_log2FC)/paletteLength, max(mat_log2FC), length.out=floor(paletteLength/2)))
    plt = pheatmap(mat_log2FC,  color = myColor, breaks = myBreaks,
        cluster_rows = F, cluster_cols = F, angle_col = c("0"))
  }
  print(plt)
}

#' Plot of frequencies by cluster
#'
#' plot was inspired by vdjtools(https://github.com/mikessh/vdjtools)
#'
#' @param CDR3PairMat
#' @param ident It can be cluster or celltype
#' @param cellLimit minimum paired cell number
#'
#' @export PlotFrequency
#'
#' @examples
#' \dontrun{
#' sset_NKT = PlotTCR(sset, group.by = )
#' }
PlotFancyChainUsage <- function(CDR3Mat, cellLimit = 3){
mat = CDR3Mat[[1]]
mat = mat[rowSums(mat) > cellLimit,colSums(mat) > cellLimit]


circos.par(gap.degree = c(rep(1, nrow(mat)-1), 10, rep(1, ncol(mat)-1), 15), start.degree = 5)

rcols <- rep(brewer.pal(12, "Paired"), nrow(mat)/12 + 1)[1:nrow(mat)]
ccols <- rep(brewer.pal(12, "Paired"), ncol(mat)/12 + 1)[1:ncol(mat)]

names(rcols) <- sort(rownames(mat))
names(ccols) <- sort(colnames(mat))

chordDiagram(mat, annotationTrack = "grid",
             grid.col = c(rcols, ccols),
             preAllocateTracks = list(track.height = 0.2), transparency = 0.5)

circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         sector.name = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                       }
)

}

