#' Mediators of Immune Response Against Cancer in soLid microEnvironments (MIRACLE)
#' 
#' Calculates Geneset Enrichment for ICR geneset and Miracle scores
#' @param matrix Numeric matrix of normalized gene expression values.
#'      This matrix holds the "Gene/Transcript Features" on the rows and the "Samples" on the columns.
#'      Normalization is typically carried out by methods like "EDASeq" for RNAseq data and "RMA" or "Quantile" methods for Microarrays.
#' @param platform Character string that shows the platform for the expression matrix. Can take values of c("ens", "u133p2", "entrez", "gene"). 
#'      The default is "ens" for Ensembl Gene Ids (ie. ENSG000..). 
#'      If the expression matrix is based on Affymetrix U133p2 or U133A platforms, please use "u133p2"
#'      If the expression matrix is based on Entrez Gene Ids or Gene Symbols, please use "entrez" or "gene", respectively.
#'      If the expression matrix is based on IlluminaHT12v3 or IlluminaHT12v4 platforms, please use "ilmn"
#'      When using platforms other than "ens", there will be some loss of information due to features not matching or multi-matching.
#' @param center Logical of length 1. User can choose whether to mean center the expression matrix. Defaults to FALSE.
#' @return data.frame of ICR (measure of immune infiltrate) and Miracle scores (Columns) for each sample (Rows).

#' @author Tolga Turan, \email{tolga.turan@abbvie.com}
#' @references \url{https://doi.org/10.1186/s40425-018-0355-5}
#' @examples
#' miracle_scores<-Calculate_Miracle(matrix)
#' @import yaGST
#' @export
#'

Calculate_Miracle<-function(matrix, platform="ens", center=FALSE){
        if (platform=="ens"){
                geneset_list1<-geneset_list}
        else if (platform=="u133p2"){
                geneset_list1<-geneset_u133p2}
        else if (platform=="entrez"){
                geneset_list1<-geneset_entrez}
        else if (platform=="gene"){
                geneset_list1<-geneset_gene}
        else if (platform=="ilmn"){
                geneset_list1<-geneset_ilmn}

      
	yaGST1_matrix2<-rbind(ICR=Geneset_Enrich(matrix, geneset_list1[1], center=TRUE),Geneset_Enrich(matrix, geneset_list1[2:3], center=center))
        yaGST1_df1<-data.frame(t(yaGST1_matrix2))
        yaGST1_df1$Miracle<-yaGST1_df1$Pos/yaGST1_df1$Neg
        #yaGST1_df1<-yaGST1_df1[,-c(2,3)]
	yaGST1_df1
        }
