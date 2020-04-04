#if (!require(STRINGdb)) install.packages('STRINGdb')
#library(STRINGdb)
#if (!require(igraph)) install.packages('igraph')
#library(igraph)


#' Creation of a choordinates data.frame from a GTF file
#'
#' Creates the choordinates data.frame by reading the data from a GTF file
#' having 9 columns (whihc is the typical format of GTF files from GENCODE).
#'
#' @param gtf_file the path tot he GTF file
#'
#' @return a data.frame with columns: \code{id type seqname start end}
#'
#' @examples
#' library(LErNet)
#' gtf_file <- system.file("extdata", "gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz", package = "LErNet")
#' complete_positions <- LErNet.load_gtf(gtf_file)
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


#' Retrieving of infomraiton fron the STRING database
#'
#' Retrieves the PPI network formt he STRING database via the STRINGdb.
#' STRINGdb often onyl associates a primary product to a gene,
#' thus other products are not reported.
#' The function also returns the proteins associated to each gene within the STRING database.
#'
#' @param stringdb_tax taxa of the species
#' @param stringdb_thr threshold to be applied to the score on the edges of the PPI
#' @param mart a biomaRt object for mapping proteins to producer genes (Ensembl IDs).
#'
#' @return a list
#' \describe{
#'   \item{\code{ppi_network}}{a two columns data.frame representing the PPI network by listing its edges.}
#'   \item{\code{ensp_to_ensg}}{a two columns data.frame reporting for each protein the corresponding gene (Ensembl IDs)}
#' }
#'
#' @examples
#' library(biomaRt)
#' stringdb_tax = 9606
#' stringdb_thr = 900
#' mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#' ret <- LErNet.get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr, mart = mart)
#' ppi_network <- ret[["ppi_network"]]
#' ensp_to_ensg <- ret[["ensp_to_ensg"]]
#'
#' @export
get_stringdb <- function(
  stringdb_tax = 9606,
  stringdb_thr = 900,
  mart
)
{
  ss<-STRINGdb$new( version="10", species=stringdb_tax, score_threshold=stringdb_thr)
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

  #return complete_ensg_to_ensp, # dataframe  ENSP (stringdb) -> ENSG
  string_prot<-unique(union(ppi$protein1, ppi$protein2))
  string_genes<-getBM(attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
                      filters = "ensembl_peptide_id", values = string_prot, mart = mart)

  print(class(string_genes))
  return( list( "ppi_network" = ppi, "ensp_to_ensg" = string_genes) )
}


#' Extraction fo the genomic context
#'
#' Retrieves the genomic context of input lncRNAs.
#' The genomic context is defined as the set of protin coding genes
#' that resides within a given range.
#'
#' @param poitions a data.frame reporting gemoic positions. Columns are \code{id type seqname start end}. It may contain features nto lited in \code{lncgenes} and \code{pcgenes}
#' @param lncgenes a list of lncRNA genes
#' @param pngenes a list of protein-coding genes that are of interest for the study.
#' @param max_window the maximum size of the genomic range
#' @param stricg_genomics if \code{FALSE}, it allows the genomic context to be forme dby p.c. genes in the \code{pcgenes} list.
#'
#' @return a two column data.frame reporting neighborhood information. The first column gives lncRNAs and the second column gives their associated neighbors.
#'
#' @examples
#' library(R.utils)
#' library(xlsx)
#' lncrna_file <- system.file("extdata", "41598_2018_30359_MOESM2_ESM.xlsx", package = "LErNet")
#' pcrna_file <- system.file("extdata", "41598_2018_30359_MOESM3_ESM.xlsx", package = "LErNet")
#' gtf_file <- system.file("extdata", "gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz", package = "LErNet")
#' pcgenes<-read.xlsx(pcrna_file,sheetIndex = 1)
#' pcgenes<-as.character(pcgenes$gene_id)
#' lncrnaInfo<-read.xlsx(lncrna_file, sheetIndex = 1)
#' lncrnaInfo <- data.frame(lapply(lncrnaInfo, as.character), stringsAsFactors=FALSE)
#' last<-which(lncrnaInfo$significant == 'FALSE')[1]
#' lncrnaInfo<-lncrnaInfo[1:last-1,]
#' lncrnaAll<-as.character(lncrnaInfo$gene_id)
#' complete_positions <- LErNet.load_gtf(gtf_file)
#' novel<-lncrnaInfo
#' novel<-novel[novel$isoform_status == "lncRNA_Novel", ]
#' chrs <- paste0("chr",sapply(strsplit(sapply(strsplit( novel$locus, "-"), `[`, 1), ":"), `[`, 1))
#' starts <- sapply(strsplit(sapply(strsplit(novel$locus, "-"), `[`, 1), ":"), `[`, 2)
#' ends <- sapply(strsplit(novel$locus, "-"), `[`, 2)
#' novel_gtf <- data.frame( "id" = novel$gene_id, "type" = rep('novel lncRNA', times = nrow(novel)),
#' "seqname" = chrs, "start" = starts, "end" = ends )
#' complete_positions <- rbind(complete_positions, novel_gtf)
#' rownames(complete_positions) <- seq(1:nrow(complete_positions))
#'
#' genomic_context <- LErNet.get_genomic_context(positions = complete_positions, lncgenes = lncrnaAll, pcgenes = pcgenes, max_window = 100000, strict_genomics = TRUE)
#'
#' @export
get_genomic_context <- function(
  positions,
  lncgenes,
  pcgenes,
  max_window = 100000,
  strict_genomics = TRUE
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
    if(strict_genomics) {
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


#' Expansion of the seed network
#'
#' Expands with connectors the network formed by seed proteins,
#' that are the producs fo the genes int he genomic context,
#' by the expasion algorithm.
#' Connectors are neighbors of selected proteins in the input PPI network.
#'
#' @param genomic_context a two column data.fram produced by \code{\link[LErNet]{get_genomic_context}}
#' @param ppi_network a two column data.frame representing PPI network edges (see also \code{\link[LErNet]{get_stringdb}} )
#' @param ensp_to_ensg a two column data.frame for mapping proteins to their producer genes (see also \code{\link[LErNet]{get_stringdb}} )
#' @param strict_proteins a list of proteins
#' @param strict_connectors if \code{TRUE} connectors can onyl be choosen from the \code{strict_proteins} list
#'
#' @return a list
#' \describe{
#'   \item{network_components}{a list of connected components of the resultant expanded network. Each compoent is a list of proteins.}
#'   \item{network_seeds}{list of seed proteins that have succefully been mapped to the PPI network.}
#' }
#'
#' @examples
#'
#' @export
expand_seeds <- function(
  genomic_context,
  ppi_network,
  ensp_to_ensg,
  strict_proteins,
  strict_connectors = TRUE
)
{
  # map seeds to ppi nodes
  closestGene<-genomic_context$partner_coding_gene
  matching<-ensp_to_ensg[ensp_to_ensg$ensembl_gene_id %in% closestGene, ]
  res_prot<-matching$ensembl_peptide_id

  #pcgenes <- strinct_list
  #mrna_annot<-getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
  #                  filters = "ensembl_gene_id", values = unique(pcgenes), mart = mart)
  #strict_proteins<-mrna_annot$ensembl_peptide_id

  #empty<-which(strict_proteins == "")
  #strict_proteins<-strict_proteins[-empty]
  strict_proteins<-strict_proteins[ strict_proteins != ""  ]

  seedprot <- res_prot
  subprot <- seedprot
  ppi <- ppi_network

  G<-make_empty_graph(n = 0, directed = FALSE)
  G<-graph_from_data_frame(d = ppi, directed = FALSE, vertices = NULL)


  resConn<-list()
  resSub<-list()
  iter=1
  cond1=FALSE
  cond2=FALSE
  isSCA = FALSE

  while(TRUE) {
    connectors<-vector("character", length=0)
    #print(paste0("***ITERATION ", iter))
    c=1
    CandGene=vector(mode = "character", length = length(subprot))
    for(i in 1:length(subprot)) {
      if(subprot[i] %in% V(G)$name) {
        N = neighbors(graph = G, v = subprot[i])
        N = as.character(N$name)
        for(j in 1:length(N)) {
          if((!(N[j] %in% CandGene)) && (!(N[j] %in% subprot)) && !is.na(N[j])) {
            CandGene[c]<-N[j]
            c=c+1
          }
        }
      }
    }
    remove(N)
    remove(i)
    remove(j)

    length(CandGene)

    while(TRUE) {

      SubNet<-make_empty_graph(n = 0, directed = FALSE)
      keepV<-which(V(G)$name %in% subprot)
      SubNet<-induced_subgraph(graph = G, vids = keepV)

      subComp<-components(graph = SubNet)
      if(max(subComp$csize)<2) {
        cond1=TRUE
        #print("esci in subComp?")
        break
      }

      maxSubComp<-which(subComp$csize==max(subComp$csize))###############[1]
      #print(paste0("How many sub connected component? ", length(maxSubComp)))

      LCC<-vector(mode = "integer", length = length(maxSubComp))
      for(c in 1:length(maxSubComp)) {
        H<-names(subComp$membership[subComp$membership==maxSubComp[c] ])
        LCC[c]<-length(intersect(H, seedprot))
      }
      start_triangles<-length(triangles(SubNet))/3
      RankCand<-vector(mode="numeric", length=length(CandGene))
      noTriangles<-vector(mode="numeric", length=length(CandGene))

      for(k in 1:length(CandGene)) {
        if(strict_connectors && !(CandGene[k] %in% strict_proteins)) {
          next
        }
        TmpGene=subprot
        TmpGene[length(TmpGene)+1]<-CandGene[k]
        TmpNet<-induced_subgraph(graph = G, vids = which(V(G)$name %in% TmpGene))

        comp<-components(graph = TmpNet)
        if(max(comp$csize)<2) {
          cond2=TRUE
          #print("esci in comp?")
          break
        }
        maxComp<-which(comp$csize==max(comp$csize))############################[1]

        LCC1<-vector(mode = "integer", length = length(maxComp))
        for(c1 in 1:length(maxComp)) {
          H1<-names(comp$membership[comp$membership==maxComp[c1] ])
          LCC1[c1]<-length(intersect(H1, seedprot))
        }
        RankCand[k]<-max(LCC1)-max(LCC)

        triangles<-length(triangles(TmpNet))/3 #numero di triangoli totali nella TmpNet

        if(isSCA) {
          noTriangles[k]=0
        }
        else {
          noTriangles[k]<-triangles-start_triangles
        }

        remove(TmpGene)

      }

      # Get the IDs of CandGene that maximize the LCC value
      idx<-which(RankCand==max(RankCand))
      if(max(RankCand)==0) {
        #print("esci al rank?")
        break
      }

      # Massimizzo il numero di triangoli
      choice=-1 #quale idx scelgo
      maxTr=-1 #max numero di triangoli

      if(length(idx)>1) {
        #print("ho più massimi uguali")
      }

      idxSameTr<-c(rep(NA, length(idx)))

      for(y in 1:length(idx)) {
        if(maxTr < noTriangles[idx[y]]) {
          maxTr=noTriangles[idx[y]]
          choice=y
          idxSameTr[y]<-idx[y]
        }
        else if(maxTr == noTriangles[idx[y]]) {
          idxSameTr[y]<-idx[y]
        }
      }

      # Scelgo in base alla ratio maggiore, ovvero scelgo il nodo i cui vicini trovano
      # maggiore corrispondenza nella lista di subprot

      if(length(idxSameTr)>1) {

        #print("scelgo in base alla ratio ")
        remove(choice)

        choice=-1
        r=-1

        for(t in 1:length(idxSameTr)) {
          if(!is.na(idxSameTr[t])) {
            maxNeigh<-neighbors(graph = G, v = CandGene[idxSameTr[t]])
            maxNeigh<-as.character(maxNeigh$name)
            ratio<-length(intersect(subprot, unique(maxNeigh)))/length(unique(maxNeigh))

            if(r<ratio) {
              r=ratio
              choice=t
            }
          }
        }
      }

      #aggiorno le mie seed aggiungendo il nuovo gene
      subprot[length(subprot)+1]<-CandGene[idx[choice]]

      #aggiungo il gene alla lista dei connectors
      connectors[length(connectors)+1]<-CandGene[idx[choice]]


      #print(paste0("Candidate gene: ", CandGene[idx[choice]]))
      # print(paste0("Rank: ", RankCand[idx[choice]]))
      #print(paste0("LCC: ", LCC))

      #print('-----------------------------')

    }
    if(cond1 | cond2) break
    if(length(connectors)>1) {
      resConn[[iter]]<-connectors
    }
    else {
      break
    }

    # Update the starting graph to extract the second biggest connected component
    new_graph=induced_subgraph(graph = G, vids = which(V(G)$name %in% subprot))

    cc<-components(graph = new_graph)
    cc_idx<-which(cc$csize==max(cc$csize))
    cc_names<-names(cc$membership[cc$membership==cc_idx])

    if(all(subprot %in% cc_names)) {
      break
    }

    firstCC<-subprot[subprot %in% cc_names]
    resSub[[iter]]<-firstCC
    nodes<-V(G)$name[V(G)$name %in% firstCC]####
    ##  G <- G - nodes
    toRemove<-which(seedprot %in% firstCC)
    seedprot<-seedprot[-toRemove]
    subprot <- seedprot

    iter=iter+1
    remove(nodes)
    remove(toRemove)
  }

  return(  list(  "network_components" = resSub, "network_seeds" = res_prot) )

 # list of lists
}


#' Expansion of the seed network
#'
#' Expands with connectors the network formed by seed proteins,
#' that are the producs fo the genes int he genomic context,
#' by the expasion algorithm.
#' Connectors are neighbors of selected proteins in the input PPI network.
#'
#' @param genomic_context a two column data.fram produced by \code{\link[LErNet]{get_genomic_context}}
#' @param ppi_network a two column data.frame representing PPI network edges (see also \code{\link[LErNet]{get_stringdb}} )
#' @param ensp_to_ensg a two column data.frame for mapping proteins to their producer genes (see also \code{\link[LErNet]{get_stringdb}} )
#' @param strict_proteins a list of proteins
#' @param strict_connectors if \code{TRUE} connectors can onyl be choosen from the \code{strict_proteins} list
#'
#' @return a list
#' \describe{
#'   \item{network_components}{a list of connected components of the resultant expanded network. Each compoent is a list of proteins.}
#'   \item{network_seeds}{list of seed proteins that have succefully been mapped to the PPI network.}
#' }
#'
#' @examples
#'
#' @export
expand_seeds_refactor <- function(
  genomic_context,
  ppi_network,
  ensp_to_ensg,
  strict_proteins,
  strict_connectors=TRUE
)
{
  # map seeds to ppi nodes
  closest_gene <- genomic_context$partner_coding_gene
  matching <- ensp_to_ensg[ensp_to_ensg$ensembl_gene_id %in% closest_gene, ]
  res_prot <- matching$ensembl_peptide_id

  strict_proteins <- strict_proteins[ strict_proteins != ""  ]

  seedprot <- res_prot
  sub_prot <- seedprot
  ppi <- ppi_network

  graph_from_df <- graph_from_data_frame(d=ppi, directed=FALSE, vertices=NULL)

  res_conn <- list()
  res_sub <- list()
  iter <- 1
  sub_net_components_invalid_csize <- FALSE
  temp_components_invalid_csize <- FALSE
  max_component_csize_valid <- TRUE
  max_rank_cand_valid <- TRUE
  # is_sca <- FALSE

  while(TRUE) {
    connectors <- vector("character", length=0)
    c <- 1
    cand_gene <- vector(mode="character", length=length(sub_prot))
    for(i in 1:length(sub_prot)) {
      if(sub_prot[i] %in% V(graph_from_df)$name) {
        neighbors_graph <- neighbors(graph=graph_from_df, v=sub_prot[i])
        neighbors_graph <- as.character(neighbors_graph$name)
        for(j in 1:length(neighbors_graph)) {
          if((!(neighbors_graph[j] %in% cand_gene)) && (!(neighbors_graph[j] %in% sub_prot)) && !is.na(neighbors_graph[j])) {
            cand_gene[c] <- neighbors_graph[j]
            c <- c+1
          }
        }
      }
    }
    remove(neighbors_graph)
    remove(i)
    remove(j)

    while(max_component_csize_valid && max_rank_cand_valid) {

      keep_v <- which(V(graph_from_df)$name %in% sub_prot)
      sub_net <- induced_subgraph(graph=graph_from_df, vids=keep_v)

      sub_net_components <- components(graph=sub_net)
      max_component_cszie <- max(sub_net_components$csize)
      if( max_component_cszie < 2) {
        sub_net_components_invalid_csize <- TRUE
        max_component_csize_valid <- FALSE
      }
      else {

        max_sub_comp <- which(sub_net_components$csize== max_component_cszie)

        LCC <- vector(mode="integer", length=length(max_sub_comp))
        for(c in 1:length(max_sub_comp)) {
          H <- names(sub_net_components$membership[sub_net_components$membership==max_sub_comp[c] ])
          LCC[c] <- length(intersect(H, seedprot))
        }
        start_triangles <- length(triangles(sub_net))/3

        cand_gene_length <- length(cand_gene)
        rank_cand <- vector(mode="numeric", length=cand_gene_length)
        no_triangles <- vector(mode="numeric", length=cand_gene_length)

        k <- 1
        while (k <= cand_gene_length) {
          if(!strict_connectors || (cand_gene[k] %in% strict_proteins)) {

            tmp_gene <- sub_prot
            tmp_gene[length(tmp_gene)+1] <- cand_gene[k]

            vertices_in_temp_net <- which(V(graph_from_df)$name %in% tmp_gene)
            tmp_net <- induced_subgraph(graph=graph_from_df, vids= vertices_in_temp_net)
            temp_net_components <- components(graph=tmp_net)

            if(max(temp_net_components$csize) < 2) {
              temp_components_invalid_csize <- TRUE
              k <- cand_gene_length + 1
            }

            max_comp <- which(temp_net_components$csize==max(temp_net_components$csize))

            max_comp_length <- length(max_comp)
            LCC1 <- vector(mode="integer", length=max_comp_length)

            for(c1 in 1:max_comp_length) {
              H1 <- names(temp_net_components$membership[temp_net_components$membership==max_comp[c1] ])
              LCC1[c1] <- length(intersect(H1, seedprot))
            }

            rank_cand[k] <- max(LCC1)-max(LCC)

            tmp_net_tringles <- length(triangles(tmp_net))/3

            no_triangles[k] <- tmp_net_tringles - start_triangles

            remove(tmp_gene)
          }
          k <- k + 1
        }

        # Get the IDs of cand_gene that maximize the LCC value
        max_rank_cand <- max(rank_cand)
        idx <- which(rank_cand==max_rank_cand)
        if(max_rank_cand == 0) {
          max_rank_cand_valid <- FALSE
        }
        else {
          # Massimizzo il numero di triangoli
          chosen_idx <- -1
          max_triangles_number <- -1

          idx_length <- length(idx)
          idx_same_tr <- c(rep(NA, idx_length))

          for(y in 1:idx_length) {
            if(max_triangles_number < no_triangles[idx[y]]) {
              max_triangles_number <- no_triangles[idx[y]]
              chosen_idx <- y
              idx_same_tr[y] <- idx[y]
            }
            else if(max_triangles_number == no_triangles[idx[y]]) {
              idx_same_tr[y] <- idx[y]
            }
          }

          # Scelgo in base alla ratio maggiore, ovvero scelgo il nodo i cui vicini trovano
          # maggiore corrispondenza nella lista di sub_prot

          if(length(idx_same_tr)>1) {

            remove(chosen_idx)

            chosen_idx <- -1
            r <- -1

            for(t in 1:length(idx_same_tr)) {
              if(!is.na(idx_same_tr[t])) {
                max_neigh <- neighbors(graph=graph_from_df, v=cand_gene[idx_same_tr[t]])
                max_neigh <- as.character(max_neigh$name)
                ratio <- length(intersect(sub_prot, unique(max_neigh)))/length(unique(max_neigh))

                if(r<ratio) {
                  r <- ratio
                  chosen_idx <- t
                }
              }
            }
          }

          #aggiorno le mie seed aggiungendo il nuovo gene
          sub_prot[length(sub_prot)+1] <- cand_gene[idx[chosen_idx]]

          #aggiungo il gene alla lista dei connectors
          connectors[length(connectors)+1] <- cand_gene[idx[chosen_idx]]
        }
      }
    }
    if(sub_net_components_invalid_csize || temp_components_invalid_csize) break
    if(length(connectors)>1) {
      res_conn[[iter]] <- connectors
    }
    else {
      break
    }

    # Update the starting graph to extract the second biggest connected component
    new_graph <- induced_subgraph(graph=graph_from_df, vids=which(V(graph_from_df)$name %in% sub_prot))

    cc <- components(graph=new_graph)
    cc_idx <- which(cc$csize==max(cc$csize))
    cc_names <- names(cc$membership[cc$membership==cc_idx])

    if(all(sub_prot %in% cc_names)) {
      break
    }

    first_cc <- sub_prot[sub_prot %in% cc_names]
    res_sub[[iter]] <- first_cc
    nodes <- V(graph_from_df)$name[V(graph_from_df)$name %in% first_cc]####

    to_remove <- which(seedprot %in% first_cc)
    seedprot <- seedprot[-to_remove]
    sub_prot <- seedprot

    iter <- iter+1
    remove(nodes)
    remove(to_remove)
  }

  return(list(  "network_components"=res_sub, "network_seeds"=res_prot))

 # list of lists
}



#' Visualization of the expanded network
#'
#' Visualizes the expanded network,
#' composed by seed proteins and connectors.
#' LncRNA are added together with extra edges in order to report the genomic context.
#' LncRNAs are linked to the proteins that are products of their genomic context.
#'
#'
#' @param lncgenes a list of lncRNA genes
#' @param genomic_context the genomc context of the lncRNAs (see \code{\link[LErNet]{get_genomic_context}})
#' @param ppi_network a two column data.frame representing PPI network edges (see also \code{\link[LErNet]{get_stringdb}} )
#' @param ensp_to_ensg a two column data.frame for mapping proteins to their producer genes (see also \code{\link[LErNet]{get_stringdb}} )
#' @param input_proteins the complete list of input proteins,that are of interest for the study
#' @param network_seeds the lsit of seed protieins
#' @param expanded_elements the list of proteins that must be visualized
#' @param mart a biomaRt ojbect used to visualize symbols instead of Ensembl IDs
#' @param mart_symbol_column column of the biomaRt ojbect fomr wich symbols are retrieved
#'
#' @return visualizes the network in the Viewer window
#'
#' @examples
#'
#' @export
visualize <-function(
  lncgenes, #strict_lncgenes <- lncgenes
  genomic_context, #lnc_neighbors <- closest
  ensp_to_ensg, #complete_ensp_to_ensg <- string_genes
  input_proteins, #strict_proteins <- DEproteins
  network_seeds, #network_seeds <- res_prot
  ppi_network, #ppi <- ppi
  expanded_elements, # expanded_elements <- unlist(resSub)
  mart, #mart <- mart
  mart_symbol_column = NULL# mart_symbol_column <- "mgi_symbol"    #def <- NULL
)
{

  strict_lncgenes <- lncgenes
  closest <- genomic_context
  complete_ensp_to_ensg <- ensp_to_ensg
  strict_proteins <- input_proteins
  ppi <- ppi_network

  starting_seed_prot <- network_seeds
  resConn <- setdiff(expanded_elements,network_seeds )
  pred<-c(network_seeds, resConn)
  netFile<-subset(ppi, ppi$protein1 %in% pred & ppi$protein2 %in% pred)
  allNodes<-union(netFile$protein1, netFile$protein2)
  nodes<-data.frame(id = c(seq(1,length(allNodes))), name = c(allNodes))
  sourceNodes<-netFile$protein1
  IDsourceNodes<-vector(mode = "integer", length = length(sourceNodes))
  for(i in 1:length(sourceNodes)) {
    if(sourceNodes[i] %in% nodes$name) {
      IDsourceNodes[i]<-nodes$id[which(nodes$name == sourceNodes[i])]
    }
  }
  targetNodes<-netFile$protein2
  IDtargetNodes<-vector(mode = "integer", length = length(targetNodes))
  remove(i)
  for(i in 1:length(targetNodes)) {
    if(targetNodes[i] %in% nodes$name) {
      IDtargetNodes[i]<-nodes$id[which(nodes$name == targetNodes[i])]
    }
  }
  edges<-data.frame(from = IDsourceNodes, to = IDtargetNodes)
  seedProt<-starting_seed_prot
  seedConn<-unlist(resConn) # quindi ho 14 connectors
  seedPred<-vector(mode = "character", length = length(nodes$name))
  remove(i)
  for(i in 1:length(seedPred)) {
    if(nodes$name[i] %in% seedProt) {
      seedPred[i]<-"Seed Protein"
    }
    else if(nodes$name[i] %in% seedConn) {
      seedPred[i]<-"Seed Connector"
    }
  }
  nodes$group<-seedPred
  flag<-rep(NA, length = nrow(nodes))
  for(d in 1:nrow(nodes)) {
    if(nodes$name[d] %in% strict_proteins) {
      flag[d]<-"DE"
    }
  }
  nodes$flag<-flag
  shape<-vector(mode = "character", length = nrow(nodes))
  for(s in 1:length(shape)) {
    if(!is.na(flag[s])) {
      shape[s]<-"star"
    }
    else {
      shape[s]<-"dot"
    }
  }
  nodes$shape<-shape
  nodes$name<-as.character(nodes$name)
  starting_nodes<-nodes

  if(! is.null(mart_symbol_column)){
    label<-getBM(attributes = c("ensembl_peptide_id","ensembl_gene_id", mart_symbol_column),
                 filters = "ensembl_peptide_id", values = nodes$name, mart = mart)

    if(nrow(label) != nrow(nodes)) {
      for(p in 1:nrow(nodes)) {
        if((nodes$name[p] %in% label$ensembl_peptide_id)) {
        }else{
          addRow<-c(rep(nodes$name[p], times = 3))
          label<-rbind(label, addRow)
        }
      }
    }
    nodes$label<-label$mgi_symbol
  }
  font.size<-c(50)
  nodes$font.size<-font.size
  nodes$values <- rep(10, nrow(nodes) )

  visNetwork(nodes, edges, width = "100%", height="1000px") %>%
    visIgraphLayout(layout = "layout_with_fr") %>%
    visEdges(color = "black") %>%
    visGroups(groupname = "Seed Protein", color = "purple") %>%
    visGroups(groupname = "Seed Connector", color = "pink") %>%
    visOptions(highlightNearest = list(enabled = T))
  starting_nodes<-nodes
  starting_edges<-edges

  colnames(closest) <- c("lnc_known", "ensembl_gene_id")
  new_closest <- merge(closest, complete_ensp_to_ensg)
  new_closest<-new_closest[,c(2,1,3)]
  colnames(new_closest) <- c("lnc_known","ensembl_gene_id","partner_coding_peptide")
  lnc<-vector(mode = "character", length = 0)
  IDtarget_lnc<-vector(mode = "integer", length = 0)

  for(k in 1:nrow(new_closest)) {
    if(new_closest$partner_coding_peptide[k] %in% nodes$name) {
      IDtarget_lnc[length(IDtarget_lnc)+1]<-nodes$id[which(nodes$name %in% new_closest$partner_coding_peptide[k])]
      lnc[length(lnc)+1]<-new_closest$lnc[k]
    }
  }
  remove(k)
  tmp_lnc<-data.frame(id = seq(from = nrow(nodes)+1, to = nrow(nodes)+length(lnc)),
                      name = lnc,
                      group = rep("lncRNA", times = length(lnc)),
                      flag = rep("DE", length(lnc)),
                      shape = rep("square", times = length(lnc)),
                      values = rep(1, times = length(lnc)),
                      label = lnc,
                      font.size = as.integer(rep(10, times = length(lnc))))

  nodes<-rbind(nodes, tmp_lnc)
  # sappiamo già qual è il target quindi ricaviamo il sorgente
  IDsource_lnc<-nodes[nodes$name %in% tmp_lnc$name, 1]
  edges_lnc<-data.frame(from = IDsource_lnc, to = IDtarget_lnc)
  edges<-rbind(edges, edges_lnc)
  visNetwork(nodes, edges, width = "100%", height="1000px") %>%
    visIgraphLayout(layout = "layout_with_fr") %>%
    #visOptions(highlightNearest = list(enabled = T)) %>%
    visEdges(color = "black") %>%
    visGroups(groupname = "Seed Protein", color = "purple") %>%
    visGroups(groupname = "Seed Connector", color = "pink") %>%
    visGroups(groupname = "lncRNA", color = rgb(237,125,49,maxColorValue=255)) %>%
    visOptions(selectedBy = "group")


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


#' Functional enrichment
#'
#' Computes functional enrichment of a given set of proteins via the ReactomePA package.
#'
#' @param ens_proteins list of proteins, in Ensembl format, for which to compute the enrichment
#' @param oganism oganism name (see \code{\link[ReactomePA]{enrichPathway}})
#' @param mart a biomaRt object of the given species need for the conversion from  Ensembl to Entrez IDs
#'
#' @return a ReactomePA result object.
#'
#' @examples
#'
#' @export
enrich <-function(
  ens_proteins,
  organism,
  mart
)
{
  mseeds <- enps_to_entrez(ens_proteins, mart)

  w<-enrichPathway(gene = unique(mseeds$entrezgene), organism = organism,
                   pvalueCutoff=0.05, pAdjustMethod = "BH",
                   qvalueCutoff = 0.2, readable=T,
                   minGSSize = 1, maxGSSize = 1000)

  enrichMappedSeeds<-as.data.frame(w)

  return(w)
}


# LErNet.characterize <- function(
#   lncgenes,  # list of lncRNAs ENSG
#   pcgenes,	# list of pcnRNAs ENSG
#   complete_positions,	# dataframe  ENSG type[pc,lnc] chr strand start end # contains also non DE elements
#   max_window = 100000,
#   strict_genomics = TRUE,
#   complete_PPI, # dataframe of edges ENSP1 - ENSP2,
#   complete_ensg_to_ensp, # dataframe  ENSG -> ENSP (stringdb)
#   strict_connectors = TRUE
# )
# {
#   genomic_context <- LErNet.get_genomic_context(
#     positions = complete_positions,
#     lncgenes = lncgenes,
#     pcgenes = pcgenes,
#     max_window = max_window,
#     strict_genomics = strict_genomics)
#
#   annot<-getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
#                filters = "ensembl_gene_id", values = unique(pcgenes), mart = mart)
#   strict_proteins<-annot$ensembl_peptide_id
#   empty<-which(strict_proteins == "")
#   strict_proteins<-strict_proteins[-empty]
#
#   ret <- LErNet.expand(
#     genomic_context = genomic_context,
#     ppi_network = ppi_network,
#     ensp_to_ensg = ensp_to_ensg,
#     strict_proteins = strict_proteins,
#     strict_connectors = TRUE)
# }


