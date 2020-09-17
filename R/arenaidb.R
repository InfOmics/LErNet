
#' Retriving the aliases of lncrnas from Arena-idb
#' 
#' @export

get_arenaidb_ncrnas_aliases <- function(){
  ncrna_aliases_file <- system.file("extdata", "arenaidb", "alias.tsv.gz", package = "LErNet")
  file_i = readLines(ncrna_aliases_file)

  ll <- list()
  li <- 1
  for(row in file_i[2:length(file_i)]){
    cc <- unlist(strsplit(row, "\t"))
    nm <- c(cc[1],cc[2])
    if(length(cc) > 3){
      for(i in seq(4,length(cc)-1, by=2 )){
        ll[[li]] <- c(nm, cc[i], cc[i+1])
        li <- li + 1
      }
    }
    else if(length(cc) == 2){
      ll[[li]] <- c(nm, nm)
      li <- li + 1
    }
  }
  ncrna_aliases <- do.call(rbind, ll)
  colnames(ncrna_aliases) <- c('principal_name', 'principal_space', 'alias_name', 'alias_space')
  as.data.frame(ncrna_aliases)
}

#' Retriving the aliases of proteins from Arena-idb
#'
#' @export

get_arenaidb_protein_aliases <- function(){
  aliases_file <- system.file("extdata", "arenaidb", "proteins.aliases.tsv.gz", package = "LErNet")
  file_i = readLines(aliases_file)

  ll <- list()
  li <- 1
  for(row in file_i[2:length(file_i)]){
    cc <- unlist(strsplit(row, "\t"))
    nm <- c(cc[1],cc[2])
    if(length(cc) > 2){
      for(i in seq(3,length(cc)-1, by=2 )){
        ll[[li]] <- c(nm, cc[i], cc[i+1])
        li <- li + 1
      }
    }
    else if(length(cc) == 2){
      ll[[li]] <- c(nm, nm)
      li <- li + 1
    }
  }
  aliases <- do.call(rbind, ll)
  colnames(aliases) <- c('principal_name', 'principal_space', 'alias_name', 'alias_space')
  as.data.frame(aliases)
}

#' Retriving the alias of genes from Arena-idb
#'
#' @export

get_arenaidb_gene_aliases <- function(){
  aliases_file <- system.file("extdata", "arenaidb", "gene.aliases.tsv.gz", package = "LErNet")
  file_i = readLines(aliases_file)

  ll <- list()
  li <- 1
  for(row in file_i[2:length(file_i)]){
    cc <- unlist(strsplit(row, "\t"))
    nm <- c(cc[1],cc[2])
    if(length(cc) > 2){
      for(i in seq(3,length(cc)-1, by=2 )){
        ll[[li]] <- c(nm, cc[i], cc[i+1])
        li <- li + 1
      }
    }
    else if(length(cc) == 2){
      ll[[li]] <- c(nm, nm)
      li <- li + 1
    }
  }
  aliases <- do.call(rbind, ll)
  colnames(aliases) <- c('principal_name', 'principal_space', 'alias_name', 'alias_space')
  as.data.frame(aliases)
}

#' Retriving the interactions involving lncrnas
#' 
#' @return arenaidb a data.frame containing in the first two columns the informations about the lncrnas and in the following those about the interactions 
#'
#' @export

get_arenaidb_interactions <- function(){
  arenaidb_file <- system.file("extdata", "arenaidb", "interactions.tsv.gz", package = "LErNet")
  arenaidb <- read.csv(arenaidb_file , stringsAsFactors = FALSE, sep='\t')
  c <- colnames(arenaidb)
  c[1] <- 'principal_name'
  c[2] <- 'principal_space'
  colnames(arenaidb) <- c
  return(arenaidb)
}

#' Retriving the disease's section of Arena-idb
#'
#' @return arenaidb a data.frame where for each lncrna is specified the diseases in which it can be involved
#'
#' @export

get_arenaidb_diseases <- function(){
  arenaidb_file <- system.file("extdata", "arenaidb", "diseases.tsv.gz", package = "LErNet")
  arenaidb <- read.csv(arenaidb_file , stringsAsFactors = FALSE, sep='\t')
  return(arenaidb)
}

#' Retriving the lncrna's section of Arena-idb
#'
#' @return arenaidb a data.frame where for each ncrna are specified all the informations
#'
#' @export
get_arenaidb_ncrnas <- function(){
  arenaidb_file <- system.file("extdata", "arenaidb", "ncRNA.tsv.gz", package = "LErNet")
  arenaidb <- read.csv(arenaidb_file , stringsAsFactors = FALSE, sep='\t', rownames=NULL)
  return(arenaidb)
}



#' Retriving the informations about the aliases and the interactions of lncrna filtering by genes and proteins
#' 
#' @param lncrna_ensembl_ids a list of lncrnas in ensembl id format
#' @param pcgene_ensembl_ids a list of protein coding genes in ensembl id format
#' @param protein_ensembl_ids a list of protein in ensembl id format
#' @param protein_symbols a list of protein in symbol id format
#'
#'
#' @return return_value a list of three data.frame: alias_mapping, gene_icontext and protein_icontext                      
#' 
#' @export


search_in_arenaidb <- function(
  lncrna_ensembl_ids,
  pcgene_ensembl_ids = NULL,
  protein_ensembl_ids = NULL,
  protein_symbols = NULL
){
  # in Arena-idb each element has a set of aliases related to it and only principal name
  # thus, ensembl ids and others must be searched in the aliases and then the related principal name can be retrieved
  # once one has the principal name, it can be used to scan for interactions
  # interactions are reported as relationships between two principal names, that are the names of one ncrna and its mate

  return_value <- list()


  # we first search for the input ncrnas, for which ensembl ids are given
  arenaidb_ncrnas_aliases <- LErNet::get_arenaidb_ncrnas_aliases()

  de_lncrna_aliases <- data.frame(
    #'alias_name' = de_lncrnas$ensembl_transcript_id,
    'alias_name' = lncrna_ensembl_ids,
    'alias_space' = rep('ensembl', nrow(lncrna_ensembl_ids)) )


  de_lncrna_aliases <- merge(de_lncrna_aliases, arenaidb_ncrnas_aliases)

  message(c( 'number of mapped ncRNA aliases: ', length(unique(de_lncrna_aliases$alias_name)) ))
  message(c( 'linked to ', length(unique(de_lncrna_aliases$principal_name)), ' entities' ))


  return_value$alias_mapping <- de_lncrna_aliases
  return_value$alias_mapping$entity_type <- rep('ncrna', times=nrow(de_lncrna_aliases))



  # now we can retrieve the interactions from those ncrnas
  arenaidb <- LErNet::get_arenaidb_interactions()

  icontext <- merge(de_lncrna_aliases[,c('principal_name','principal_space')], arenaidb)

  message(c('retireved ',nrow(icontext),' interactions form them' ))

  message(c( nrow(icontext[icontext$mate.class == 'ncrna', ]), ' are ncRNA mates' ))
  message(c( nrow(icontext[icontext$mate.class == 'gene', ]), ' are gene mates' ))
  message(c( nrow(icontext[icontext$mate.class == 'protein', ]), ' are protein mates'))

  gene_icontext = icontext[icontext$mate.class == 'gene', ]
  protein_icontext = icontext[icontext$mate.class == 'protein', ]


  if(! is.null(pcgene_ensembl_ids)){
    # we filter for genes
    arenaidb_pcgene_aliases <- LErNet::get_arenaidb_gene_aliases()

    de_pcgene_aliases <- data.frame(
      'alias_name' = pcgene_ensembl_ids,
      'alias_space' = rep('ensembl', length(pcgene_ensembl_ids)) )

    de_pcgene_aliases <- merge(de_pcgene_aliases, arenaidb_pcgene_aliases)

    message(c('number of mapped gene aliases: ', length(unique(de_pcgene_aliases$alias_name)) ))
    message(c( 'linked to ', length(unique(de_pcgene_aliases$principal_name)), ' entities'))

    filtering_genes = de_pcgene_aliases[, c('principal_name','principal_space')]
    colnames(filtering_genes)  <- c('mate.name', 'mate.naming.resource')
    de_pcgene_icontext <- merge(filtering_genes, gene_icontext)

    message(c( nrow( unique(de_pcgene_icontext)), ' unique gene interactions after filtering' ))
    message(c( length( unique(de_pcgene_icontext$mate.name)), ' filtering genes are involved' ))

    x <- de_pcgene_aliases
    x$entity_type <- rep('gene', times=nrow(de_pcgene_aliases))
    return_value$alias_mapping <- rbind(return_value$alias_mapping, x)
    return_value$gene_icontext <- de_pcgene_icontext

  } else {
    return_value$gene_icontext <- gene_icontext
  }


  if( (!is.null(protein_ensembl_ids)) | (!is.null(protein_symbols)) ){
    # we filter for proteins
    arenaidb_protein_aliases <- LErNet::get_arenaidb_protein_aliases()


    de_protein_aliases <- data.frame(
      'alias_name' = protein_ensembl_ids,
      'alias_space' = rep('ensembl', length(protein_ensembl_ids)) )

    de_protein_aliases <- merge(de_protein_aliases, arenaidb_protein_aliases)

    message(c( length(unique(de_protein_aliases$alias_name)), ' mapped ensembl protein ids' ))


    de_protein_symbol_aliases <- data.frame(
      'alias_name' = protein_symbols,
      'alias_space' = rep('symbol', length(protein_symbols)) ) 
    de_protein_symbol_aliases <- merge(de_protein_symbol_aliases, arenaidb_protein_aliases)

    de_protein_aliases <- rbind(de_protein_aliases, de_protein_symbol_aliases )

    message(c( length(unique(de_protein_aliases$alias_name)), ' mapped protein symbols' ))

    message(c( length(unique(de_protein_aliases$principal_name)), ' resultant mapped protein entites'  ))


    filtering_proteins = de_protein_aliases[, c('principal_name','principal_space')]
    colnames(filtering_proteins)  <- c('mate.name', 'mate.naming.resource')
    de_protein_icontext <- merge(filtering_proteins, protein_icontext)

    message(c( nrow(de_protein_icontext), ' unique protein interactions after filtering' ))
    message(c( length( unique(de_protein_icontext$mate.name)), ' filtering proteins are involved' ))


    x <- de_protein_aliases
    x$entity_type <- rep('protein', times=nrow(de_protein_aliases))
    return_value$alias_mapping <- rbind(return_value$alias_mapping, x)
    return_value$protein_icontext <- de_protein_icontext
  } else {
    return_value$protein_icontext <- protein_icontext
  }

  return(return_value)

}
