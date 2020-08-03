#' The stcrSet class
#' 
#' The stcrSet object is the main class object that contains vdj results of 10x vdj.
#' We recommands to use IMGT reference and Seurat package to analyze the single cell RNA sequence data.
#' 
#' 
#' @slot contig a data frame of 10x contig 
#' @slot cellName a vector of cell names
#' @slot ident a factor of cell types or clusters
#' @slot reductions a list of reductions, it has the same format as Seurat object 
#' @slot vdjMat a list of 
#' @slot meta.data a data frame of meta.data from Seurat object
#' 
#' @name StcrSet
#' @rdname StcrSet
#' @aliases StcrSet-class
#' @exportClass StcrSet
StcrSet <- setClass( 'StcrSet',
          slots = c(contig = 'data.frame',
                    cellName = 'vector',
                    ident = 'factor',
                    reductions = 'list',
                    vdjMat = 'list',
                    meta.data = 'data.frame'))
