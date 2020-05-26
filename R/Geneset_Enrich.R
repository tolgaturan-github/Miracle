#' Geneset Enrichment (GSE)
#' Calculates single sample GSE scores using "yaGST" method
#' @param matrix Numeric matrix of normalized gene expression values with unique column and rownames. 
#' 	This matrix holds the "Gene/Transcript Features" on the rows and the "Samples" on the columns. 
#'	Normalization is typically carried out by methods like "EDASeq" for RNAseq data and "RMA" or "Quantile" methods for Microarrays.
#' @param geneset_list a List of gene/Transcript features. 
#' 	Preferably all or at least 1 feature in each geneset must be present in the rownames of the expression matrix.
#' @param verbose Logical of length 1. User can choose to display messages. Defaults to TRUE. 
#'	If there are large numbers of genesets it can be set to FALSE
#' @param n_cores Integer of length 1. Specifies the number of parallel cores running the function.
#'	Defaults to 1
#' @param center Logical of length 1. User can choose whether to mean center the expression matrix. Defaults to TRUE.
#' @return data.frame holding geneset enrichment values for each geneset (rows) and for each sample (columns). 
#'	If the geneset is of length 1, then a named vector of enrichment values.
#'
#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{http://github.com/miccec/yaGST}
#' @examples 
#' normalized_enrichment_scores<-Geneset_Enrich(matrix1, geneset_list1)
#' @import yaGST
#' @export
#'

Geneset_Enrich<-function(matrix1, geneset_list, verbose=TRUE, center=TRUE, n_cores=1){
        for (i in 1:length(geneset_list)){
                        if (verbose) {if(any(!geneset_list[[i]] %in% rownames(matrix1))){
                                print(paste0("IDs", geneset_list[[i]][!geneset_list[[i]] %in% rownames(matrix1)], " not present in the expression matrix for the geneset ",names(geneset_list)[i]))}
                        else{print(paste0("All IDs are found in the matrix for the geneset: ", names(geneset_list)[i]))}}}
        library(yaGST)
        library(parallel)
	if(center){
		matrix1z<-matrix1
		for (i in 1:nrow(matrix1z)){
			matrix1z[i,] <-(matrix1[i,] - mean(matrix1[i,]))
			}
		matrix1<-matrix1z
		rm(matrix1z)
		}

        m1t<-t(matrix1)
        m1t_l1<-split(m1t, rownames(m1t))
        m1t_l2<-lapply(m1t_l1, function(x) setNames(x, rownames(matrix1)))
	if(.Platform$OS.type == "unix") {
        yaGST1_list1<-mclapply(m1t_l2, function(x) {rL<-sort(x, decreasing=TRUE)
                                                   sapply(geneset_list, function(y) mwwGST(rL, y, minLenGeneSet=1)[7:11])}, mc.cores=n_cores)
	}
	else{
	yaGST1_list1<-lapply(m1t_l2, function(x) {rL<-sort(x, decreasing=TRUE)
                                                   sapply(geneset_list, function(y) mwwGST(rL, y, minLenGeneSet=1)[7:11])})
	}
	
        yaGST1_list2<-lapply(yaGST1_list1, function(x) data.frame(x))
        yaGST1_matrix1<-sapply(yaGST1_list2, function(x) sapply(x, function(y) y[[2]]))
        if (!is.null(dim(yaGST1_matrix1))){
                yaGST1_matrix2<-yaGST1_matrix1[,order(order(colnames(matrix1)))]}
        else {yaGST1_matrix2<-yaGST1_matrix1[order(order(colnames(matrix1)))]
                 names(yaGST1_matrix2)<-gsub(paste0("\\.", names(geneset_list)), "", names(yaGST1_matrix2))
                 }
        yaGST1_matrix2
        }

