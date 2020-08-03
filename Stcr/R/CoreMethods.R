#' Load contig annotation file from 10X(>2.2.0)
#'
#' Opens vdj result file provided by 10X genomics.
#' returns data.frame
#'
#' @param dataDir Directory containing contig_annotations.csv
#' @param filtered boolean option to choose filtered(TRUE)/all(FALSE)
#' @param highConf boolean option to set highConfidence
#' @param productive boolean option to set productive
#'
#' @export Load10xVDJ
#'
#' @examples
#' \dontrun{
#' data_dir = 'path/to/data/directory'
#' vdj_contig = StcrLoad(dataDir = data_dir, filtered = FALSE)
#' }
Load10xVDJ <- function(dataDir = NULL, filtered = TRUE, highConf = TRUE, productive = TRUE ){
  if(!dir.exists(paths = dataDir)){
    stop('Directory does not exist!')
  }
  if(filtered){
    contigLoc = file.path(dataDir, 'filtered_contig_annotations.csv')
  } else {
    contigLoc = file.path(dataDir, 'all_contig_annotations.csv')
  }
  data = read.csv(contigLoc, row.names = 'contig_id', stringsAsFactors = F)
  if(highConf){
    data = subset(data, high_confidence=='True')
  }
  if(productive){
    data = subset(data, productive=='True')
  }
  
  sset <- new( 'StcrSet',
               contig = data,
               cellName = vector(),
               reductions = list(),
               ident = factor(),
               vdjMat = list(),
               meta.data = data.frame())
  return(sset)
}

#' Assign sample ID to a sset
#'
#' @param sset StcrSet 
#' @param sampleID string 
#'
#' @export SetSample
#'
#' @examples
#' \dontrun{
#' ssetF3 = SetSample(ssetF3,"F3")
#' }
SetSample <- function(sset,sampleID = NULL, sep = "-"){
  if(is.null(sampleID)){
    stop('sampleID is NULL')
  }
  if(length(sampleID) == 1){
    sset@contig$cellName = ""
    print('Setting All data with given sample ID')
    sset@contig$cellName = paste0(sampleID,sep,gsub('([ATCG]*)-.*','\\1',rownames(sset@contig)))
  } else {
    if(length(sampleID) == nrow(sset@contig)){
      sset@contig$cellName = ""
      sset@contig$cellName = paste0(sampleID,sep,gsub('([ATCG]*)-.*','\\1',rownames(sset@contig)))
    } else {
      stop('in order to assign sampleID using vector input, # of contig has to match with the length of sampleID')
    }
  }
  return(sset)

}

#' merge two ssets
#'
#' @param sset1 First StcrSet which contains cellName
#' @param sset2 Second StcrSet which contains cellName 
#'
#' @export MergeSamples
#'
#' @examples
#' \dontrun{
#' sset = MergeSamples(ssetF3,ssetF5) 
#' }
MergeSamples <- function(sset1,sset2){
  # can't merge without assigne sample(HasCellName?)
  if(class(sset1)[1] != "StcrSet" | class(sset2)[1] != "StcrSet"){
    stop("your inputs are not StcrSet classes")
  }
  if(is.null(sset1@contig$cellName) |is.null(sset2@contig$cellName) ){
    message("we recommand to run SetSample first")
  }

  if(length(sset1@reductions) > 0 | length(sset1@ident) > 0  | length(sset1@vdjMat) > 0 |
     length(sset2@reductions) > 0 | length(sset2@ident) > 0  | length(sset2@vdjMat) > 0 ){
    message("resetting reductions, ident, vdjMat")
  }
  sset <- new( 'StcrSet',
               contig = rbind(sset1@contig, sset2@contig),
               cellName = vector(),
               reductions = list(),
               ident = factor(),
               vdjMat = list(),
               meta.data = data.frame())
  
  return(sset)
}

#' Import Seurat object to StcrSet
#'
#' @param sset StcrSet object which contains cellName
#' @param seuratObject Second StcrSet which contains cellName 
#' @param version Seurat object version (default v2). v3 is also available.
#'
#' @export ImportSeurat
#'
#' @examples
#' \dontrun{
#' sset_NKT = ImportSeurat(sset,seurat.NKT)
#' }
ImportSeurat <- function(sset,seuratObject, version = "v2"){

  if(version == "v3"){
    if(class(sset)[1] != "StcrSet" | class(seuratObject)[1] != "Seurat"){
      stop("your inputs are not correct classes")
    }
    if(!any(colnames(seuratObject) %in% sset@contig$cellName)){
      stop("sset and seurat cellnames do not match")
    }
  } else {
    if(class(sset)[1] != "StcrSet" | class(seuratObject)[1] != "seurat"){
      stop("your inputs are not correct classes")
    }
    if(!any(colnames(seuratObject@data) %in% sset@contig$cellName)){
      stop("sset and seurat cellnames do not match")
    }
  }
  if(is.null(sset@contig$cellName)){
    stop("cellName doesn't exist, use SetSample first")
  }

  if(version =="v3"){
    sset@cellName = colnames(seuratObject@assays$RNA@data)
    sset@ident = seuratObject@active.ident
    sset@reductions = seuratObject@reductions
    sset@contig = subset(sset@contig, cellName %in% sset@cellName)
    sset@meta.data = seuratObject@meta.data
  }
  else{
    sset@cellName = colnames(seuratObject@data)
    sset@ident = seuratObject@ident
    sset@reductions = seuratObject@dr
    sset@contig = subset(sset@contig, cellName %in% sset@cellName)
    sset@meta.data = seuratObject@meta.data
  }
    
  return(sset)
}

#' Create a matrix using StcrSet and one of vdj gene.
#' 
#' return a matrix with specific chain and gene
#'
#' @param sset StcrSet object which contains cellName
#' @param chain a set of string contains more than one chain ex) c("TRA") or "TRA"
#' @param ref default "IMGT"
#' @param v_gene only one gene can be selected among v_gene, d_gene, j_gene ex)"AV"
#' @param d_gene only one gene can be selected among v_gene, d_gene, j_gene ex)"AD"
#' @param j_gene only one gene can be selected among v_gene, d_gene, j_gene ex)"AJ"
#'
#' @export SetVDJcountMat
#'
#' @examples
#' \dontrun{
#' mat_AV = SetVDJcountMat(sset_NKT, chain = "TRA", v_gene = "AV")
#' }
SetVDJcountMat <- function(sset, chain = NULL, 
                           ref= "IMGT", v_gene = NULL, d_gene = NULL, j_gene = NULL){
  if(is.null(chain)){
    stop("need chain information")
  }
  subdf = subset(sset@contig, chain %in% chain)
  if(!is.null(v_gene)){
    v_gene_list = unique(as.vector(subdf$v_gene))
    v_gene_list = v_gene_list[grepl(v_gene, v_gene_list)]
    
    mat_TCR = matrix(data=0, ncol=length(v_gene_list), nrow=length(sset@cellName))
    rownames(mat_TCR) = sset@cellName
    colnames(mat_TCR) = v_gene_list

    subdf_v_gene = subset(subdf, v_gene %in% v_gene_list)
    ct1 = as.data.frame(plyr::count(subdf_v_gene, c("cellName", "v_gene")))
    ct2 = dcast(ct1, cellName ~ v_gene,value.var = "freq", fun.aggregate = length)
    rownames(ct2) = ct2[,1]
    ct3 = as.matrix(ct2[,2:ncol(ct2)])
    colnames(ct3) = colnames(ct2)[2:ncol(ct2)]
    rownames(ct3) = rownames(ct2)
    mat_TCR[rownames(ct3),colnames(ct3)] = ct3
    

    
  } else if(!is.null(d_gene)){
    d_gene_list = unique(as.vector(subdf$d_gene))
    d_gene_list = d_gene_list[grepl(d_gene, d_gene_list)]
    
    mat_TCR = matrix(data=0, ncol=length(d_gene_list), nrow=length(sset@cellName))
    rownames(mat_TCR) = sset@cellName
    colnames(mat_TCR) = d_gene_list

    subdf_d_gene = subset(subdf, d_gene %in% d_gene_list)
    ct1 = as.data.frame(plyr::count(subdf_d_gene, c("cellName", "d_gene")))
    ct2 = dcast(ct1, cellName ~ d_gene,value.var = "freq", fun.aggregate = length)
    rownames(ct2) = ct2[,1]
    ct3 = as.matrix(ct2[,2:ncol(ct2)])
    colnames(ct3) = colnames(ct2)[2:ncol(ct2)]
    rownames(ct3) = rownames(ct2)
    mat_TCR[rownames(ct3),colnames(ct3)] = ct3
    
  } else if(!is.null(j_gene)){
    j_gene_list = unique(as.vector(subdf$j_gene))
    j_gene_list = j_gene_list[grepl(j_gene, j_gene_list)]
    
    mat_TCR = matrix(data=0, ncol=length(j_gene_list), nrow=length(sset@cellName))
    rownames(mat_TCR) = sset@cellName
    colnames(mat_TCR) = j_gene_list

    subdf_j_gene = subset(subdf, j_gene %in% j_gene_list)
    ct1 = as.data.frame(plyr::count(subdf_j_gene, c("cellName", "j_gene")))
    ct2 = dcast(ct1, cellName ~ j_gene,value.var = "freq", fun.aggregate = length)
    rownames(ct2) = ct2[,1]
    ct3 = as.matrix(ct2[,2:ncol(ct2)])
    colnames(ct3) = colnames(ct2)[2:ncol(ct2)]
    rownames(ct3) = rownames(ct2)
    mat_TCR[rownames(ct3),colnames(ct3)] = ct3
  } else {
    stop("gene information required")
  }

  return(mat_TCR)

}

#' Create a matrix of CDR3 of each chains
#' 
#' chain1 -> row, chian2 -> column
#' VJ
#' return a matrix
#'
#' @param sset StcrSet object which contains cellName
#' @param cell_ind
#' @param chain1 chains(more than one can also be possible)
#' @param chain2 chains(more than one can also be possible)
#' @param v_gene1 additional gene info, if you use more than one chain
#' @param v_gene2 additional gene info, if you use more than one chain
#' @param j_gene1 additional gene info, if you use more than one chain
#' @param j_gene2 additional gene info, if you use more than one chain
#'
#' @export SetCDR3Mat
#'
#' @examples
#' \dontrun{
#' gdT_cdr3_mat = SetCDR3Mat(sset_gdT, cell_ind = colnames(seurat.gdT@scale.data), 
#' chain1 = c("TRG"), chain2 = c("TRD","Multi"), v_gene2 = "DV", j_gene2="DJ")
#' }
SetCDR3Mat <- function(sset, cell_ind = NULL, chain1 = c("TRA"), chain2 = c("TRB"), v_gene1 = NULL, v_gene2 = NULL
                       , j_gene1 = NULL, j_gene2 = NULL){
  
  if(is.null(v_gene1)){
    v_gene1 = paste0(substr(chain1,3,3),'V')
  }
  if(is.null(v_gene2)){
    v_gene2 = paste0(substr(chain2,3,3),'V')
  }
  
  if(is.null(j_gene1)){
    j_gene1 = paste0(substr(chain1,3,3),'J')
  }
  if(is.null(j_gene2)){
    j_gene2 = paste0(substr(chain2,3,3),'J')
  }
  
  df = sset@contig
  if(is.null(cell_ind)){
    cell_ind = sset@cellName
    sub_df = subset(df, cellName %in% cell_ind)
  } else {
    sub_df = subset(df, cellName %in% cell_ind)
  }
  
  sub_df_chain1 = subset(sub_df, chain %in% chain1)
  sub_df_chain1 = sub_df_chain1[grepl(paste0(".*(",v_gene1,"[0-9]+).*"),sub_df_chain1$v_gene) &
                                  grepl(paste0(".*(",j_gene1,"[0-9]+).*"),sub_df_chain1$j_gene) ,]
  sub_df_chain1$index = rownames(sub_df_chain1)
  sub_df_chain1$cdr3_gene = paste0(sub_df_chain1$cdr3,
                                   "-",
                                   gsub(paste0(".*(",v_gene1,"[0-9]+).*"),"\\1",sub_df_chain1$v_gene),
                                   "-",
                                   gsub(paste0(".*(",j_gene1,"[0-9]+).*"),"\\1",sub_df_chain1$j_gene))
  
  major_sub_df_chain1_cell_reads = aggregate(reads~cellName,sub_df_chain1,function(x) x[which.max(x)])
  major_sub_df_chain1_ind =c()
  for(i in 1:nrow(major_sub_df_chain1_cell_reads)){
    major_sub_df_chain1_ind = c(major_sub_df_chain1_ind,
                             subset(sub_df_chain1,cellName == major_sub_df_chain1_cell_reads$cellName[i] & reads == major_sub_df_chain1_cell_reads$reads[i])$index[1])
  }
  major_sub_df_chain1 = sub_df_chain1[major_sub_df_chain1_ind, ]
  
  sub_df_chain2 = subset(sub_df, chain %in% chain2)
  sub_df_chain2 = sub_df_chain2[grepl(paste0(".*(",v_gene2,"[0-9]+).*"),sub_df_chain2$v_gene) &
                                  grepl(paste0(".*(",j_gene2,"[0-9]+).*"),sub_df_chain2$j_gene) ,]
  sub_df_chain2$index = rownames(sub_df_chain2)
  sub_df_chain2$cdr3_gene = paste0(sub_df_chain2$cdr3,
                                   "-",
                                   gsub(paste0(".*(",v_gene2,"[0-9]+).*"),"\\1",sub_df_chain2$v_gene),
                                   "-",
                                   gsub(paste0(".*(",j_gene2,"[0-9]+).*"),"\\1",sub_df_chain2$j_gene))
  major_sub_df_chain2_cell_reads = aggregate(reads~cellName,sub_df_chain2,function(x) x[which.max(x)])
  major_sub_df_chain2_ind =c()
  for(i in 1:nrow(major_sub_df_chain2_cell_reads)){
    major_sub_df_chain2_ind = c(major_sub_df_chain2_ind,
                             subset(sub_df_chain2,cellName == major_sub_df_chain2_cell_reads$cellName[i] & reads == major_sub_df_chain2_cell_reads$reads[i])$index[1])
  }
  major_sub_df_chain2 = sub_df_chain2[major_sub_df_chain2_ind, ]
  major_df = rbind(major_sub_df_chain1,major_sub_df_chain2)


  
  
  v_gene_list = unique(major_df$v_gene)
  
  v_gene_list_chain1 = v_gene_list[grepl(v_gene1,v_gene_list)]
  chain1_cdr3_gene = unique(subset(major_df, v_gene %in% v_gene_list_chain1)$cdr3_gene)
  
  v_gene_list_chain2 = v_gene_list[grepl(v_gene2,v_gene_list)]
  chain2_cdr3_gene = unique(subset(major_df, v_gene %in% v_gene_list_chain2)$cdr3_gene)
  
  mat_chain1 = matrix(data = 0 , nrow = length(cell_ind), ncol = length(chain1_cdr3_gene))
  rownames(mat_chain1) = cell_ind
  colnames(mat_chain1) = chain1_cdr3_gene
  
  major_df_v_gene_chain1 = subset(major_df, v_gene %in% v_gene_list_chain1)
  ct1 = as.data.frame(plyr::count(major_df_v_gene_chain1, c("cellName", "cdr3_gene")))
  ct2 = dcast(ct1, cellName ~ cdr3_gene,value.var = "freq", fun.aggregate = length)
  rownames(ct2) = ct2[,1]
  ct3 = as.matrix(ct2[,2:ncol(ct2)])
  colnames(ct3) = colnames(ct2)[2:ncol(ct2)]
  rownames(ct3) = rownames(ct2)
  mat_chain1[rownames(ct3),colnames(ct3)] = ct3

  
  mat_chain2 = matrix(data = 0 , nrow = length(cell_ind), ncol = length(chain2_cdr3_gene))
  rownames(mat_chain2) = cell_ind
  colnames(mat_chain2) = chain2_cdr3_gene

  major_df_v_gene_chain2 = subset(major_df, v_gene %in% v_gene_list_chain2)
  ct1 = as.data.frame(plyr::count(major_df_v_gene_chain2, c("cellName", "cdr3_gene")))
  ct2 = dcast(ct1, cellName ~ cdr3_gene,value.var = "freq", fun.aggregate = length)
  rownames(ct2) = ct2[,1]
  ct3 = as.matrix(ct2[,2:ncol(ct2)])
  colnames(ct3) = colnames(ct2)[2:ncol(ct2)]
  rownames(ct3) = rownames(ct2)
  mat_chain2[rownames(ct3),colnames(ct3)] = ct3
  
  mat_cdr3= t(mat_chain1) %*% mat_chain2
  return(list(mat_cdr3,major_df))
}

#' Create a matrix of CDR3 of each chains
#' 
#' chain1 -> row, chian2 -> column
#' return a list with a matrix of two CDR3 chains and a matrix of its CDR3 pair and clonotypes
#'
#' @param sset StcrSet object which contains cellName
#' @param CDR3Mat a list output of SetCDR3Mat
#' @param chain1 default TRA
#' @param chain2 default TRB
#'
#' @export SetCDR3PairMat
#'
#' @examples
#' \dontrun{
#' gdT_cdr3_mat = SetCDR3Mat(sset_gdT, cell_ind = colnames(seurat.gdT@scale.data), 
#' chain1 = c("TRG"), chain2 = c("TRD","Multi"), v_gene2 = "DV", j_gene2="DJ")
#' gdT_cdr3_pair_mat = SetCDR3PairMat(sset_gdT, gdT_cdr3_mat, chain1 = c("TRG"), chain2 = c("TRD","Multi"))
#' }
SetCDR3PairMat <- function(sset, CDR3Mat, chain1 = c("TRA"), chain2 = c("TRB")){
  
  chain1_cdr3_ind = names(rowSums(CDR3Mat[[1]]))[rowSums(CDR3Mat[[1]]) > 0]
  chain2_cdr3_ind = names(colSums(CDR3Mat[[1]]))[colSums(CDR3Mat[[1]]) > 0]
  
  CDR3Mat_2 =CDR3Mat[[1]][chain1_cdr3_ind,chain2_cdr3_ind]
  major_df = CDR3Mat[[2]]
  
  cdr3_pairs = c()
  for(i in 1:nrow(CDR3Mat_2)){
    for(j in 1:ncol(CDR3Mat_2)){
      if(CDR3Mat_2[i,j] > 0){
        cdr3_pairs = c(cdr3_pairs, paste0(rownames(CDR3Mat_2)[i],"_",colnames(CDR3Mat_2)[j]))
      }
    }
  }
  
  mat_cdr3_pair =  matrix(data = 0 , nrow = length(sset@cellName), ncol = length(cdr3_pairs))
  rownames(mat_cdr3_pair) = sset@cellName
  colnames(mat_cdr3_pair) = cdr3_pairs
  for(j in rownames(mat_cdr3_pair)){
    cdr3s = subset(major_df, cellName %in% j)$cdr3_gene
    cdr3_chain1 = subset(major_df, cellName %in% j & chain %in% chain1)$cdr3_gene
    cdr3_chain2 = subset(major_df, cellName %in% j & chain %in% chain2)$cdr3_gene
    if(!is.null(cdr3_chain1) & !is.null(cdr3_chain2) & length(cdr3s) == 2){
      mat_cdr3_pair[j,paste0(cdr3_chain1,"_",cdr3_chain2)] = 1
    }
    
  }
  
  df_clonotype = data.frame(pair =colnames(mat_cdr3_pair), clonotype = paste0("clonotype",seq(1,ncol(mat_cdr3_pair))))
  colnames(mat_cdr3_pair) = df_clonotype$clonotype
  for(i in 1:ncol(mat_cdr3_pair)){
    df_clonotype$chain1_cdr3[i]= unlist(strsplit(as.vector(df_clonotype$pair[i]),"_"))[1]
    df_clonotype$chain2_cdr3[i]= unlist(strsplit(as.vector(df_clonotype$pair[i]),"_"))[2]
  }
  
  for(i in 1:ncol(mat_cdr3_pair)){
    # print(unique(as.vector(subset(major_df, chain %in% chain1 & cdr3_gene == df_clonotype$chain1_cdr3[i])$v_gene)))
    df_clonotype$chain1_v_gene[i]= paste(unique(as.vector(subset(major_df, chain %in% chain1 & cdr3_gene == df_clonotype$chain1_cdr3[i])$v_gene)), collapse=",")
    # print(unique(as.vector(subset(major_df, chain %in% chain1 & cdr3_gene == df_clonotype$chain1_cdr3[i])$j_gene)))
    df_clonotype$chain1_j_gene[i]= paste(unique(as.vector(subset(major_df, chain %in% chain1 & cdr3_gene == df_clonotype$chain1_cdr3[i])$j_gene)), collapse=",")
    # print(unique(as.vector(subset(major_df, chain %in% chain2 & cdr3_gene == df_clonotype$chain2_cdr3[i])$v_gene)))
    df_clonotype$chain2_v_gene[i]= paste(unique(as.vector(subset(major_df, chain %in% chain2 & cdr3_gene == df_clonotype$chain2_cdr3[i])$v_gene)), collapse=",")
    # print(unique(as.vector(subset(major_df, chain %in% chain2 & cdr3_gene == df_clonotype$chain2_cdr3[i])$j_gene)))
    df_clonotype$chain2_j_gene[i]= paste(unique(as.vector(subset(major_df, chain %in% chain2 & cdr3_gene == df_clonotype$chain2_cdr3[i])$j_gene)), collapse=",")
    
  }
  df_clonotype$n = colSums(mat_cdr3_pair)
  return(list(mat_cdr3_pair,df_clonotype))
}
