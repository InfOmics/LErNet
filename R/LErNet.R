#if (!require(STRINGdb)) install.packages('STRINGdb')
#library(STRINGdb)
#if (!require(igraph)) install.packages('igraph')
#library(igraph)


#' Extraction of the genomic context
#'
#' Retrieves the genomic context of input lncRNAs.
#' The genomic context is defined as the set of protin coding genes
#' that resides within a given range.
#'
#' @param positions a data.frame reporting genomic positions. Columns are \code{id type seqname start end}. It may contain features into listed in \code{lncgenes} and \code{pcgenes}
#' @param lncgenes a list of lncRNA genes
#' @param pcgenes a list of protein-coding genes that are of interest for the study. If empty, the complete set of neighbors is extracted otherwise only pcgenes in this set are taken into account.
#' @param max_window the maximum size of the genomic range
#'
#' @return a two column data.frame reporting neighborhood information. The first column gives lncRNAs and the second column gives their associated neighbors.
#'
#' @examples
#' genomic_context <- LErNet::get_genomic_context(
#'positions = complete_positions,
#'lncgenes = de_lncrnas$ensembl_gene_id,
#'pcgenes = de_pcgenes,
#'max_window = 100000)
#' @export
get_genomic_context <- function(
  positions,
  lncgenes,
  pcgenes = NULL,
  max_window = 100000
)
{
  # return data.frame  lnc ensg -> neighbor ensg
  coding_gtf<-positions[positions$type == "protein_coding",]

  closest<-data.frame()
  for(k in 1:length(lncgenes)) {
    idx <- which(positions$id == lncgenes[k])
    chr = positions$seqname[idx]
    start = as.numeric(positions$start[idx])
    end = as.numeric(positions$end[idx])
    tmp<-coding_gtf[which(coding_gtf$seqname == chr),]
    if(!is.null(pcgenes)) {
      # Only DE-neighbors (protein coding)
      gid_tmp<-which(tmp$id %in% pcgenes)
      if(length(gid_tmp)<1) {
        next
      }
      tmp<-tmp[gid_tmp,]
    }
    tmp<-makeGRangesFromDataFrame(tmp,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field="seqname",
                                  start.field="start",
                                  end.field="end",
                                  starts.in.df.are.0based=FALSE)
    query<-GRanges(seqnames = chr,
                   ranges = IRanges(start = start-max_window, end = end+max_window))
    #, strand = s)

    neighbors<-nearest(query, tmp, select = "all", ignore.strand=TRUE)
    items<-subjectHits(neighbors)
    if(length(items)==0) {
      next
    }
    # values = to get the metadata (other functions: elementMetadata or mcols)
    res<-values(tmp[items])
    final<-data.frame(lnc_known = rep(lncgenes[k], times=length(res$id)),
                      partner_coding_gene = res$id)
    closest<-rbind(closest,final)
    remove(final)
    remove(tmp)
  }


  closest<-closest[!duplicated(closest),]
  closest[,c(1,2)]<-lapply(closest[,c(1,2)], as.character)
 # return data.frame  lnc ensg -> neighbor ensg

  return(closest)
}

#' Extraction of the connected components of the network
#'
#' @param net a data.frame reporting the edges of the network
#'
#' @return components a list of all connected proteins
#'
#'
#' @export
get_connected_components <- function(
  net
){
  components <- list()

  tovisit = unique(union(net$protein1, net$protein2))
  #print(tovisit)

  icomponent <- 1
  while( length(tovisit) > 0){
    components[[icomponent]] <- c(tovisit[1])
    tovisit <- tovisit[-1]

    iqueue <- 1
    while( iqueue <= length(components[[icomponent]]) ){
      node <- components[[icomponent]][iqueue]
      neighs <- intersect( unique(union(net[net$protein1 == node,]$protein2, net[net$protein2 == node,]$protein1)), tovisit)
      components[[icomponent]] <- c(components[[icomponent]], neighs)
      tovisit <- setdiff(tovisit, neighs)
      iqueue <- iqueue + 1
    }

    icomponent <- icomponent + 1
  }

  return(components)
}


#' Expansion of the seed network
#'
#' Expands with connectors the network formed by seed proteins,
#' that are the produtcs of the genes in the genomic context and in the interaction context,
#' by the expansion algorithm.
#' Connectors are neighbors of selected proteins in the input PPI network.
#'
#' @param seeds a list produced by \code{\link[LErNet]{get_genomic_context}} and \code{\link[LErNet]{search_in_arenaidb}}
#' @param ppi_network a two column data.frame representing PPI network edges (see also \code{\link[LErNet]{get_stringdb}} )
#' @param de_proteins a vector of stricted proteins. If empty, the complete set of neighbors is extracted otherwise only pcgenes in this set are taken into account.
#'
#' @return a list
#' \describe{
#'   \item{out_components}{a list of connected components of the resultant expanded network. Each compoent is a list of proteins.}
#' }
#'
#' @examples
#' components <- LErNet::expand_seeds(seeds,  ppi_network,  de_proteins)
#' @export

expand_seeds<- function(
  seeds,
  ppi_network,
  de_proteins = NULL
)
{
  #print(nrow(ppi))
  ppi <- ppi_network
  if(!is.null(de_proteins)){
    ppi <- ppi[ ppi$protein1 %in% de_proteins, ]
    ppi <- ppi[ ppi$protein2 %in% de_proteins, ]
  }
  #print(nrow(ppi))

  seed_ppi <- ppi[ ppi$protein1 %in% seeds, ]
  seed_ppi <- seed_ppi[ seed_ppi$protein2 %in% seeds, ]
  #print(c('@', nrow(seed_ppi)))

  init_components <- get_connected_components(seed_ppi)

  for(seed in setdiff( seeds, unlist(init_components) )){
    init_components <- append(init_components, list(seed) )
  }

  #selected <- unlist(init_components, recursive = FALSE)
  selected_n_co <- data.frame( matrix(ncol=2, nrow=length( unlist(init_components, recursive = FALSE) )) )
  colnames(selected_n_co) <- c('node', 'coco')
  selected_n_co$node = unlist(init_components, recursive = FALSE)
  for(i in 1:length(init_components)){
    selected_n_co$coco[ selected_n_co$node %in% init_components[[i]]] <- i
  }


  out_components <- list()


  while(nrow(selected_n_co) > 0){

    cocos <- unique(selected_n_co$coco)
    coco_sizes <- data.frame(matrix(ncol=2, nrow = length(cocos), dimnames=list(NULL, c('coco','size')) ))
    coco_sizes$coco <- cocos
    for(c in cocos){
      coco_sizes$size[coco_sizes$coco == c] <- nrow( selected_n_co[selected_n_co$coco == c, ] )
    }
    maxcoco <- coco_sizes[which.max(coco_sizes$size),'coco']


    coco_neighs <- ppi[ ppi$protein1 %in% selected_n_co$node[selected_n_co$coco == maxcoco] ,]
    coco_neighs <- ppi[ ppi$protein2 %in% selected_n_co$node[selected_n_co$coco == maxcoco],]
    coco_neighs <- unique( union( coco_neighs$protein1, coco_neighs$protein2 ) )
    coco_neighs <- setdiff(coco_neighs, selected_n_co$node)

    #print(coco_neighs)

    while(length(coco_neighs) > 0){


      scores <- data.frame(matrix(ncol=4, nrow=length(coco_neighs)))
      colnames(scores) <- c('node','lcc','triangles','seed_edges')
      scores$node <- coco_neighs

      for(irow in 1:nrow(scores)){
        node <- scores[irow,'node']
        node_neighs <-   intersect(selected_n_co$node,
                                   unique(union(ppi[ ppi$protein1 == node, ]$protein2,
                                                ppi[ ppi$protein2 == node, ]$protein1))
        )
        coco_neighs <- unique(selected_n_co[ selected_n_co$node %in% node_neighs, 'coco' ])

        scores[irow,'lcc'] <- nrow( selected_n_co[selected_n_co$coco %in% coco_neighs, ] )

        scores[irow, 'seed_edges'] <- length(coco_neighs)
      }
      scores <- scores[scores$lcc > 1,]


      if( max(scores$lcc) > max(coco_sizes$size) ){
        print(c('@',max(scores$lcc),max(coco_sizes$size)))


        maxneighs <- which( scores$lcc == max(scores$lcc) )
        if(length(maxneighs) == 1){
          imaxneigh <- maxneighs[1]
          new_connector  <- scores[imaxneigh,'node']

        } else {
          for(irow in maxneighs){
            node <- scores[irow,'node']
            node_neighs <- intersect(selected_n_co$node,
                                     unique(union(ppi[ ppi$protein1 == node, ]$protein2,
                                                  ppi[ ppi$protein2 == node, ]$protein1))
            )

            triangles <- seed_ppi[seed_ppi$protein1 %in% node_neighs, ]
            triangles <- triangles[triangles$protein2 %in% node_neighs, ]

            scores[irow,'triangles'] <- nrow(triangles)
          }

          maxneighs <- which( scores$triangles == max(scores$triangles) )
          if(length(maxneighs) == 1){
            imaxneigh <- maxneighs[1]
            new_connector  <- scores[imaxneigh,'node']

          } else {
            new_connector <- scores[which.max( scores$seed_edges ), 'node']
          }
        }


        # merge connected components by putting the same identifier for all of them
        node_neighs <- intersect(selected_n_co$node,
                                 unique(union(ppi[ ppi$protein1 == new_connector, ]$protein2,
                                              ppi[ ppi$protein2 == new_connector, ]$protein1))
        )
        neighs_coco <- unique(selected_n_co[ selected_n_co$node %in% node_neighs, 'coco' ])
        one_coco <- min(neighs_coco)
        selected_n_co$coco[selected_n_co$coco %in% neighs_coco] <- one_coco

        # add the connector to the merged connected components
        selected_n_co <- rbind(selected_n_co, c(new_connector, one_coco))

        # add all the edges between the connector and the seeds to the current seed network
        # it will be used for calculating triangles in future iterations
        seed_ppi <- rbind(seed_ppi, ppi[ ppi$protein1 == new_connector, ])
        seed_ppi <- rbind(seed_ppi, ppi[ ppi$protein2 == new_connector, ])
        ppi <- subset(ppi, protein1 != new_connector, protein2 != new_connector)


        # since the previous largest connected component is now merged with at least one other component
        # the list of nodes that are neighbors of the new largest component needs to be updated
        coco_neighs <- ppi[ ppi$protein1 %in% selected_n_co$node[selected_n_co$coco == one_coco] ,]
        coco_neighs <- ppi[ ppi$protein2 %in% selected_n_co$node[selected_n_co$coco == one_coco],]
        coco_neighs <- unique( union( coco_neighs$protein1, coco_neighs$protein2 ) )
        coco_neighs <- setdiff(coco_neighs, selected_n_co$node)



        cocos <- unique(selected_n_co$coco)
        coco_sizes <- data.frame(matrix(ncol=2, nrow = length(cocos), dimnames=list(NULL, c('coco','size')) ))
        coco_sizes$coco <- cocos
        for(c in cocos){
          coco_sizes$size[coco_sizes$coco == c] <- nrow( selected_n_co[selected_n_co$coco == c, ] )
        }
        maxcoco <- coco_sizes[which.max(coco_sizes$size),'coco']


      } else {
        coco_neighs <- list()
      }
    }

    cocos <- unique(selected_n_co$coco)
    coco_sizes <- data.frame(matrix(ncol=2, nrow = length(cocos), dimnames=list(NULL, c('coco','size')) ))
    coco_sizes$coco <- cocos
    for(c in cocos){
      coco_sizes$size[coco_sizes$coco == c] <- nrow( selected_n_co[selected_n_co$coco == c, ] )
    }
    maxcoco <- coco_sizes[which.max(coco_sizes$size),'coco']


    out_components <- append(out_components, list(  selected_n_co[ selected_n_co$coco == maxcoco, 'node']  ))
    selected_n_co <- subset(selected_n_co, coco != maxcoco)
  }


  for(seed in setdiff(seeds, unlist(out_components))){
    out_components <- append(out_components, list(seed))
  }
  return(out_components)
}


#' Visualization of the expanded network
#'
#' Visualizes the expanded network,
#' composed by seed proteins and connectors.
#' LncRNA are added together with extra edges in order to report the lncrna context.
#' LncRNAs are linked to the proteins that are products of their lncrna context.
#'
#'
#' @param lncrna_context the genomic and interaction context of the lncRNAs (see \code{\link[LErNet]{get_genomic_context}}) and \code{\link[LErNet]{search_in_arenaidb}}
#' @param de_proteins the complete vector of input proteins,that are of interest for the study
#' @param network_seeds the list of seed proteins
#' @param expanded_elements the list of proteins (see also \code{\link[LErNet]{expand_seeds}} ) that must be visualized
#' @param bg_ppi_network a two column data.frame representing PPI network edges (see also \code{\link[LErNet]{get_stringdb}} )
#'
#' @return visualizes the network in the Viewer window
#'
#' @examples
#' LErNet::visualize(
#'unique(lncrna_context),
#'unique(de_proteins),
#'unique(network_seeds),
#'unique(expanded_elements),
#'unique(bg_ppi_network),
#'labels)
#'
#' @export
visualize <- function(
  lncrna_context,
  de_proteins,
  network_seeds,
  expanded_elements,
  bg_ppi_network,
  labels
){
  #nodes <- data.frame( matrix(ncol=8, nrow=0, dimnames=list(NULL, c('id','name','group','flag','shape','values','label','font.size'))) )

  lncrnas = unique(lncrna_context$lncrna)

  x <- merge(lncrnas, labels, by.x=1, by.y='id')

  lnc_nodes <- data.frame(
    id = seq(from=1, to=nrow(x)),
    name = x$label,
    group = rep('lncRNA', times=nrow(x)),
    flag = rep('DE', times=nrow(x)),
    shape = rep('square', times=nrow(x)),
    values = rep(1, times=nrow(x)),
    label = x$label,
    font.size = as.integer(rep(20,times=nrow(x))),
    origin = x$x
  )


  y <- merge(  data.frame('protein'=unique( c(lncrna_context$mate, network_seeds, expanded_elements) )), labels, by.x='protein', by.y='id'  )

  protein_nodes <- data.frame(
    id = seq(from=nrow(x)+1, to=nrow(x)+nrow(y)),
    name = y$label,
    group = rep('Connector', times=nrow(y)),
    flag = rep(NA, times=nrow(y)),
    shape = rep('dot', times=nrow(y)),
    values = rep(50, times=nrow(y)),
    label = y$label,
    #label = y$protein,
    font.size = as.integer(rep(50,times=nrow(y))),
    origin = y$protein
  )

  protein_nodes$group[ protein_nodes$origin %in% network_seeds ] <- 'Seed'
  protein_nodes$group[ protein_nodes$origin %in% lncrna_context$mate ] <- 'Seed'
  protein_nodes$flag[ protein_nodes$origin %in% de_proteins ] <- 'DE'
  protein_nodes$shape[protein_nodes$flag != 'DE'] <- 'diamond'
  protein_nodes$shape[is.na(protein_nodes$flag)] <- 'diamond'


  context_edges <- (merge(lncrna_context, lnc_nodes, by.x='lncrna', by.y='origin'))[,c('id','mate')]
  colnames(context_edges) <- c('from', 'mate')
  context_edges <- merge(context_edges, protein_nodes, by.x='mate',by.y='origin')[,c('from','id')]
  colnames(context_edges) <- c('from', 'to')


  ppi_edges <- bg_ppi_network[bg_ppi_network$protein1 %in% protein_nodes$origin,]
  ppi_edges <- ppi_edges[ppi_edges$protein2 %in% protein_nodes$origin,]
  nrow(ppi_edges)
  ppi_edges <- (merge(ppi_edges, protein_nodes, by.x='protein1', by.y='origin'))[, c('id', 'protein2')]
  colnames(ppi_edges) <- c('from', 'protein2')
  nrow(ppi_edges)
  ppi_edges <- (merge(ppi_edges, protein_nodes, by.x='protein2', by.y='origin'))[, c('from', 'id')]
  colnames(ppi_edges) <- c('from', 'to')
  nrow(ppi_edges)




  nodes <- rbind(lnc_nodes, protein_nodes)
  rownames(nodes) <- seq(from=1, to=nrow(nodes))

  edges <- rbind(context_edges, ppi_edges)


  library(visNetwork)
  visNetwork(nodes, edges, width = "100%", height="1000px")%>%
    #visIgraphLayout(layout = "layout_with_fr") %>%
    #visIgraphLayout(layout = "layout_in_circle") %>%
    #visIgraphLayout(layout = "layout_with_sugiyama") %>%
    visIgraphLayout(layout = "layout_nicely" ,physics = TRUE, smooth = TRUE) %>%
    visEdges(color = "black") %>%
    visGroups(groupname = "lncRNA", color = "orange") %>%
    visGroups(groupname = "Seed", color = "purple") %>%
    visGroups(groupname = "Connector", color = "pink") %>%
    visOptions(highlightNearest = list(enabled = T)) %>%
    visOptions(selectedBy = "group")
}


#' Functional enrichment
#'
#' Computes functional enrichment of a given set of proteins via the ReactomePA package.
#'
#' @param entrez_genes list of proteins, in Entrez format, for which to compute the enrichment
#' @param oganism oganism name (see \code{\link[ReactomePA]{enrichPathway}})
#'
#' @return a ReactomePA result object.
#'
#' @examples
#'entrez_ids <- unique(annot[annot$ensembl_peptide %in% comp,]$entrezgene_id)
#'enrichment <- LErNet::enrich(entrez_ids, 'human')
#' @export

enrich <-function(
  entrez_genes,
  organism
)
{
  # library(ReactomePA)
  # #for human
  # BiocManager::install("org.Hs.eg.db")
  # library(org.Hs.eg.db)

  w<-enrichPathway(gene = unique(entrez_genes), organism = organism,
                   pvalueCutoff=0.05, pAdjustMethod = "BH",
                   qvalueCutoff = 0.2, readable=T,
                   minGSSize = 1, maxGSSize = 1000)

  enrichMappedSeeds<-as.data.frame(w)

  return(w)
}
