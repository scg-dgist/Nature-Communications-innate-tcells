
ComputeFreq <- function(mat_cdr3_pair,df_clonotype_labeled,ident,
  chain1_v_gene = NULL, chain1_d_gene = NULL, chain1_j_gene = NULL, 
  chain2_v_gene = NULL, chain2_d_gene = NULL, chain2_j_gene = NULL){

  nonselected_clonotypes = as.vector(df_clonotype_labeled$clonotype)[df_clonotype_labeled$selected == FALSE]
  selected_clonotypes = as.vector(df_clonotype_labeled$clonotype)[df_clonotype_labeled$selected == TRUE]
  
  nonselected_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_clonotypes]) > 0]
  selected_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_clonotypes]) > 0]
  

  df_freq = data.frame(row.names = levels(ident))

  
  
  if(!is.null(chain1_v_gene)){
    selected_chain1_v_clonotypes = as.vector(df_clonotype_labeled$clonotype)[grepl(chain1_v_gene,df_clonotype_labeled$chain1_v_gene)]
    selected_chain1_v_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_chain1_v_clonotypes]) > 0]
    selected_chain1_v_counts = rep(0,nlevels(ident))
    names(selected_chain1_v_counts) = levels(ident)
    for(i in rownames(df_freq)){
      selected_chain1_v_counts[i] = sum(selected_chain1_v_ind %in% names(ident)[ident == i])
    }
    nonselected_chain1_v_clonotypes = as.vector(df_clonotype_labeled$clonotype)[!grepl(chain1_v_gene,df_clonotype_labeled$chain1_v_gene)]
    nonselected_chain1_v_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_chain1_v_clonotypes]) > 0]
    nonselected_chain1_v_counts = rep(0,nlevels(ident))
    names(nonselected_chain1_v_counts) = levels(ident)
    for(i in rownames(df_freq)){
      nonselected_chain1_v_counts[i] = sum(nonselected_chain1_v_ind %in% names(ident)[ident == i])
    }
    chain1_v_counts = selected_chain1_v_counts + nonselected_chain1_v_counts
    df_freq$selected_chain1_v_pct = selected_chain1_v_counts/chain1_v_counts
    df_freq$nonselected_chain1_v_pct = 1 - df_freq$selected_chain1_v_pct
  }
  if(!is.null(chain1_d_gene)){
    selected_chain1_d_clonotypes = as.vector(df_clonotype_labeled$clonotype)[grepl(chain1_d_gene,df_clonotype_labeled$chain1_d_gene)]
    selected_chain1_d_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_chain1_d_clonotypes]) > 0]
    selected_chain1_d_counts = rep(0,nlevels(ident))
    names(selected_chain1_d_counts) = levels(ident)
    for(i in rownames(df_freq)){
      selected_chain1_d_counts[i] = sum(selected_chain1_d_ind %in% names(ident)[ident == i])
    }
    nonselected_chain1_d_clonotypes = as.vector(df_clonotype_labeled$clonotype)[!grepl(chain1_d_gene,df_clonotype_labeled$chain1_d_gene)]
    nonselected_chain1_d_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_chain1_d_clonotypes]) > 0]
    nonselected_chain1_d_counts = rep(0,nlevels(ident))
    names(nonselected_chain1_d_counts) = levels(ident)
    for(i in rownames(df_freq)){
      nonselected_chain1_d_counts[i] = sum(nonselected_chain1_d_ind %in% names(ident)[ident == i])
    }
    chain1_d_counts = selected_chain1_d_counts + nonselected_chain1_d_counts
    df_freq$selected_chain1_d_pct = selected_chain1_d_counts/chain1_d_counts
    df_freq$nonselected_chain1_d_pct = 1 - df_freq$selected_chain1_d_pct
  }
  if(!is.null(chain1_j_gene)){
    selected_chain1_j_clonotypes = as.vector(df_clonotype_labeled$clonotype)[grepl(chain1_j_gene,df_clonotype_labeled$chain1_j_gene)]
    selected_chain1_j_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_chain1_j_clonotypes]) > 0]
    selected_chain1_j_counts = rep(0,nlevels(ident))
    names(selected_chain1_j_counts) = levels(ident)
    for(i in rownames(df_freq)){
      selected_chain1_j_counts[i] = sum(selected_chain1_j_ind %in% names(ident)[ident == i])
    }
    nonselected_chain1_j_clonotypes = as.vector(df_clonotype_labeled$clonotype)[!grepl(chain1_j_gene,df_clonotype_labeled$chain1_j_gene)]
    nonselected_chain1_j_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_chain1_j_clonotypes]) > 0]
    nonselected_chain1_j_counts = rep(0,nlevels(ident))
    names(nonselected_chain1_j_counts) = levels(ident)
    for(i in rownames(df_freq)){
      nonselected_chain1_j_counts[i] = sum(nonselected_chain1_j_ind %in% names(ident)[ident == i])
    }
    chain1_j_counts = selected_chain1_j_counts + nonselected_chain1_j_counts
    df_freq$selected_chain1_j_pct = selected_chain1_j_counts/chain1_j_counts
    df_freq$nonselected_chain1_j_pct = 1 - df_freq$selected_chain1_j_pct
  }
  if(!is.null(chain2_v_gene)){
    selected_chain2_v_clonotypes = as.vector(df_clonotype_labeled$clonotype)[grepl(chain2_v_gene,df_clonotype_labeled$chain2_v_gene)]
    selected_chain2_v_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_chain2_v_clonotypes]) > 0]
    selected_chain2_v_counts = rep(0,nlevels(ident))
    names(selected_chain2_v_counts) = levels(ident)
    for(i in rownames(df_freq)){
      selected_chain2_v_counts[i] = sum(selected_chain2_v_ind %in% names(ident)[ident == i])
    }
    nonselected_chain2_v_clonotypes = as.vector(df_clonotype_labeled$clonotype)[!grepl(chain2_v_gene,df_clonotype_labeled$chain2_v_gene)]
    nonselected_chain2_v_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_chain2_v_clonotypes]) > 0]
    nonselected_chain2_v_counts = rep(0,nlevels(ident))
    names(nonselected_chain2_v_counts) = levels(ident)
    for(i in rownames(df_freq)){
      nonselected_chain2_v_counts[i] = sum(nonselected_chain2_v_ind %in% names(ident)[ident == i])
    }
    chain2_v_counts = selected_chain2_v_counts + nonselected_chain2_v_counts
    df_freq$selected_chain2_v_pct = selected_chain2_v_counts/chain2_v_counts
    df_freq$nonselected_chain2_v_pct = 1 - df_freq$selected_chain2_v_pct
  }
  if(!is.null(chain2_d_gene)){
    selected_chain2_d_clonotypes = as.vector(df_clonotype_labeled$clonotype)[grepl(chain2_d_gene,df_clonotype_labeled$chain2_d_gene)]
    selected_chain2_d_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_chain2_d_clonotypes]) > 0]
    selected_chain2_d_counts = rep(0,nlevels(ident))
    names(selected_chain2_d_counts) = levels(ident)
    for(i in rownames(df_freq)){
      selected_chain2_d_counts[i] = sum(selected_chain2_d_ind %in% names(ident)[ident == i])
    }
    nonselected_chain2_d_clonotypes = as.vector(df_clonotype_labeled$clonotype)[!grepl(chain2_d_gene,df_clonotype_labeled$chain2_d_gene)]
    nonselected_chain2_d_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_chain2_d_clonotypes]) > 0]
    nonselected_chain2_d_counts = rep(0,nlevels(ident))
    names(nonselected_chain2_d_counts) = levels(ident)
    for(i in rownames(df_freq)){
      nonselected_chain2_d_counts[i] = sum(nonselected_chain2_d_ind %in% names(ident)[ident == i])
    }
    chain2_d_counts = selected_chain2_d_counts + nonselected_chain2_d_counts
    df_freq$selected_chain2_d_pct = selected_chain2_d_counts/chain2_d_counts
    df_freq$nonselected_chain2_d_pct = 1 - df_freq$selected_chain2_d_pct
  }
  if(!is.null(chain2_j_gene)){
    selected_chain2_j_clonotypes = as.vector(df_clonotype_labeled$clonotype)[grepl(chain2_j_gene,df_clonotype_labeled$chain2_j_gene)]
    selected_chain2_j_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,selected_chain2_j_clonotypes]) > 0]
    selected_chain2_j_counts = rep(0,nlevels(ident))
    names(selected_chain2_j_counts) = levels(ident)
    for(i in rownames(df_freq)){
      selected_chain2_j_counts[i] = sum(selected_chain2_j_ind %in% names(ident)[ident == i])
    }
    nonselected_chain2_j_clonotypes = as.vector(df_clonotype_labeled$clonotype)[!grepl(chain2_j_gene,df_clonotype_labeled$chain2_j_gene)]
    nonselected_chain2_j_ind = rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair[,nonselected_chain2_j_clonotypes]) > 0]
    nonselected_chain2_j_counts = rep(0,nlevels(ident))
    names(nonselected_chain2_j_counts) = levels(ident)
    for(i in rownames(df_freq)){
      nonselected_chain2_j_counts[i] = sum(nonselected_chain2_j_ind %in% names(ident)[ident == i])
    }
    chain2_j_counts = selected_chain2_j_counts + nonselected_chain2_j_counts
    df_freq$selected_chain2_j_pct = selected_chain2_j_counts/chain2_j_counts
    df_freq$nonselected_chain2_j_pct = 1 - df_freq$selected_chain2_j_pct
  }

  selected_counts = rep(0,nlevels(ident))
  names(selected_counts) = levels(ident)
  for(i in rownames(df_freq)){
    selected_counts[i] = sum(selected_ind %in% names(ident)[ident == i])
  }
  nonselected_counts = rep(0,nlevels(ident))
  names(nonselected_counts) = levels(ident)
  for(i in rownames(df_freq)){
    nonselected_counts[i] = sum(nonselected_ind %in% names(ident)[ident == i])
  }

  clono_counts = selected_counts + nonselected_counts
  df_freq$selected_pct = selected_counts/clono_counts
  df_freq$nonselected_pct = 1 - df_freq$selected_pct
  
  return(df_freq)

}

ComputeShannonIndex <- function(vector){
  # k = length(vector)
  k = sum(vector > 0)
  return(-sum((vector/sum(vector))*log(vector/sum(vector)), na.rm=TRUE)/log(k))
}

ComputeGeomean <- function(CDR3PairMat, ident){
  mat_cdr3_pair = CDR3PairMat[[1]]
  df_clonotype = CDR3PairMat[[2]]
  mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(ident),
                                  nrow =ncol(mat_cdr3_pair))
  colnames(mat_cluster_clonotypes) = levels(ident)
  rownames(mat_cluster_clonotypes) = colnames(mat_cdr3_pair)
  for(i in 1:nlevels(ident)){
    cluster_cell_ind = names(ident)[ident == levels(ident)[i]]
    mat_cdr3_pair_cluster = mat_cdr3_pair[rownames(mat_cdr3_pair) %in% cluster_cell_ind,]
    mat_cluster_clonotypes[,levels(ident)[i]] = colSums(mat_cdr3_pair_cluster)
  }
  rownames(mat_cluster_clonotypes) = df_clonotype$clonotype

  colnames(mat_cluster_clonotypes) = levels(ident)
  geo_mean = matrix(data = 0, ncol = nlevels(ident),
                  nrow = nlevels(ident))
  rownames(geo_mean) = levels(ident)
  colnames(geo_mean) = levels(ident)
  for(i in 1:nlevels(ident)){
    for(j in 1:nlevels(ident)){
      if (i != j){
        sumA = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),i]/sum(mat_cluster_clonotypes[,i]))
        sumB = sum(mat_cluster_clonotypes[(mat_cluster_clonotypes[,i]*mat_cluster_clonotypes[,j] > 0),j]/sum(mat_cluster_clonotypes[,j]))
        geo_mean[i,j] = sqrt(sumA*sumB)
      }
    }
  }
  return(geo_mean)
}

ComputeLogFC <- function(CDR3PairMat, ident){
  mat_cdr3_pair = CDR3PairMat[[1]]
  df_clonotype = CDR3PairMat[[2]]
  geo_mean = ComputeGeomean(CDR3PairMat, ident)
  set.seed(42)
  mat_pvalue = matrix(data = 0, ncol = nlevels(ident),
                      nrow = nlevels(ident))
  rownames(mat_pvalue) = levels(ident)
  colnames(mat_pvalue) = levels(ident)
  mat_accum = matrix(data = 0, ncol = nlevels(ident),
                    nrow = nlevels(ident))
  rownames(mat_accum) = levels(ident)
  colnames(mat_accum) = levels(ident)

  totaliter =10000
  for(iter in 1:totaliter){
    sampled_cluster = names(ident)
    names(sampled_cluster) = names(ident)
    sampled_cluster =sampled_cluster[sample(names(ident))]
    
    cluster_counts = table(ident)
    sampled_cluster_value = c()
    for(i in 1:(nlevels(ident))){
      sampled_cluster_value = c(sampled_cluster_value, rep(names(cluster_counts)[i],cluster_counts[[i]]))
    }  
    names(sampled_cluster_value) = names(sampled_cluster)
    sampled_cluster = sampled_cluster_value
    
    sampled_mat_cluster_clonotypes = matrix(data=0, ncol=nlevels(ident),
                                            nrow =ncol(mat_cdr3_pair))
    colnames(sampled_mat_cluster_clonotypes) = levels(ident)
    rownames(sampled_mat_cluster_clonotypes) = colnames(mat_cdr3_pair)
    
    for(i in 1:nlevels(ident)){
      cluster_cell_ind = names(sampled_cluster)[sampled_cluster == levels(ident)[i]]
      mat_cdr3_pair_cluster = mat_cdr3_pair[rownames(mat_cdr3_pair) %in% cluster_cell_ind,]
      sampled_mat_cluster_clonotypes[,levels(ident)[i]] = colSums(mat_cdr3_pair_cluster)
    }
    
    rownames(sampled_mat_cluster_clonotypes) = df_clonotype$clonotype
    
    colnames(sampled_mat_cluster_clonotypes) = levels(ident)
    
    sampled_geo_mean = matrix(data = 0, ncol = nlevels(ident),
                              nrow = nlevels(ident))
    rownames(sampled_geo_mean) = levels(ident)
    colnames(sampled_geo_mean) = levels(ident)
    
    for(i in 1:nlevels(ident)){
      for(j in 1:nlevels(ident)){
        if (i != j){
          sumA = sum(sampled_mat_cluster_clonotypes[(sampled_mat_cluster_clonotypes[,i]*sampled_mat_cluster_clonotypes[,j] > 0),i]/sum(sampled_mat_cluster_clonotypes[,i]))
          sumB = sum(sampled_mat_cluster_clonotypes[(sampled_mat_cluster_clonotypes[,i]*sampled_mat_cluster_clonotypes[,j] > 0),j]/sum(sampled_mat_cluster_clonotypes[,j]))
          sampled_geo_mean[i,j] = sqrt(sumA*sumB)
        }
      }
    }
    for(i in 1:nlevels(ident)){
      for(j in 1:nlevels(ident)){
        if (i != j){
          mat_accum[i,j] = mat_accum[i,j] + sampled_geo_mean[i,j]
          if(sampled_geo_mean[i,j] >= geo_mean[i,j])
            mat_pvalue[i,j] = mat_pvalue[i,j] + 1
        }
      }
    }
  }
  P = (mat_pvalue +1) /(totaliter+1)
  mat_pvalue_twotail = matrix(data = 0, ncol = nlevels(ident),
                              nrow = nlevels(ident))
  rownames(mat_pvalue_twotail) = levels(ident)
  colnames(mat_pvalue_twotail) = levels(ident)
  for(i in 1:nlevels(ident)){
    for(j in 1:nlevels(ident)){
      mat_pvalue_twotail[i,j] = 2*min(P[i,j],(1-P[i,j]))
    }
  }
  print("Pvalue 2-tail ")
  print(P)
  mat_avg_geo_mean = mat_accum/10000

  mat_log2FC = log2((geo_mean+1)/(mat_avg_geo_mean+1))
  return(mat_log2FC)
}