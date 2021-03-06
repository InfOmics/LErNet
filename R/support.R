#' Creation of a coordinates data.frame from a GTF file
#'
#' Creates the coordinates data.frame by reading the data from a GTF file
#' having 9 columns (whihc is the typical format of GTF files from GENCODE).
#'
#' @param gtf_file the path tot he GTF file
#'
#' @return a data.frame with columns: \code{id type seqname start end}
#'
#' @examples
#' gtf_file <- system.file("extdata", "gencode.v34.annotation.gtf.gz", package = "LErNet")
#' complete_positions <- LErNet::load_gtf(gtf_file)
#'
#' @export
load_gtf <- function(
  gtf_file
)
{
  gtf<-read.table(gtf_file, header = FALSE, sep = "\t")
  colnames(gtf)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
  gtf[,c(1,2,3,6,7,8,9)]<-lapply(gtf[,c(1,2,3,6,7,8,9)], as.character)
  gtf <- gtf[gtf$feature == "gene",]

  ExtractAttributes<-function(attributes,c){
    res<-lapply(attributes,function(x){
      unlist(strsplit(x,";",fixed=TRUE))[c]
    })
    m<-matrix(unlist(res),ncol=length(c),byrow=TRUE)
    m<-apply(m,1,strsplit," ") #1 = rows
    m<-(unlist(m)[unlist(m)!=""])
    m<-matrix(m,ncol = length(c)*2,byrow=TRUE)
    m<-m[,(1:ncol(m))%%2==0]
  }

  #gtf$id <- ExtractAttributes(gtf$attribute,1)
  gtf$id <- sapply(strsplit( ExtractAttributes(gtf$attribute,1), "[.]"), `[`, 1)
  gtf$type <- ExtractAttributes(gtf$attribute,2)

  complete_positions <- gtf[, c("id","type","seqname","start","end")]

  return(complete_positions)
}

#' Retrieving of information from the STRING database
#'
#' Retrieves the PPI network form the STRING database via the STRINGdb.
#' STRINGdb often only associates a primary product to a gene,
#' thus other products are not reported.
#' The function also returns the proteins associated to each gene within the STRING database.
#'
#' @param stringdb_tax taxa of the species. Default human (9606)
#' @param stringdb_thr threshold to be applied to the score on the edges of the PPI. Default threshold value 900
#'
#' @return a list
#' \describe{
#'   \item{\code{ppi}}{a two columns data.frame representing the PPI network by listing its edges.}
#' }
#'
#' @examples
#' library(STRINGdb)
#' library(igraph)
#' stringdb_tax = 9606
#' stringdb_thr = 900
#' ppi_network <- LErNet::get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr)
#'
#' @export

get_stringdb <- function(
  stringdb_tax = 9606,
  stringdb_thr = 900
)
{
  ss<-STRINGdb$new( version="11", species=stringdb_tax, score_threshold=stringdb_thr)
  g<-ss$get_graph()
  ppi<-as.data.frame(get.edgelist(g))
  ppi[,c(1,2)]<-lapply(ppi[,c(1,2)], as.character)
  colnames(ppi)<-c("protein1", "protein2")
  #head(ppi)

  #rimuovo il prefisso (primi 6 caratteri)
  nchar(ppi$protein1[1]) #24 caratteri totali
  nchar(stringdb_tax) #5 (aggiungere +2 per rimuovere anche il punto)
  ppi$protein1<-substr(ppi$protein1,nchar(stringdb_tax)+2,nchar(ppi$protein1[1]))
  ppi$protein2<-substr(ppi$protein2,nchar(stringdb_tax)+2,nchar(ppi$protein2[1]))
  #head(ppi)

  ppi<-ppi[with(ppi, order(ppi$protein1, ppi$protein2)),]
  ppi<-ppi[!duplicated(ppi),]

  return( ppi )
}



#' Mapping from Ensembl to Entrez
#'
#' Maps a list of protein IDs in the Ensmbl format to the Entrez naming system
#'
#' @param ens_proteins list of Ensembl IDs
#' @paramt mart a biomaRt object for the given species
#'
#' @return a data.frame representing the mapping
#'
#'
#' @export
enps_to_entrez <-function(
  ens_proteins,
  mart
)
{
  #tryCatch(
  #  {mseeds <- getBM(attributes = c("ensembl_peptide_id","ensembl_gene_id","entrezgene"),
  #                filters = "ensembl_peptide_id", values = ens_proteins, mart = mart)},
  #  error= function(err){
  mseeds <- getBM(attributes = c("ensembl_peptide_id","ensembl_gene_id","entrezgene_id"),
                  filters = "ensembl_peptide_id", values = ens_proteins, mart = mart)
  #   }
  #)
  #mseeds <- getBM(attributes = c("ensembl_peptide_id","ensembl_gene_id","entrezgene"),
  #                filters = "ensembl_peptide_id", values = ens_proteins, mart = mart)
  return(mseeds)
}


