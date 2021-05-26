library(LErNet)
library(biomaRt)
library(GenomicRanges)
library(STRINGdb)
library(igraph)
library(Organism.dplyr)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationHub)
library(enrichplot)
library(ReactomePA)

#load("MART.RData")
table4 <- read.csv(
  system.file("extdata","558557_Yang", "Table4.csv", package = "LErNet"),
  sep='\t'
)

filtered_table4 <- subset(table4,  (pvalue<0.05) & (padj<0.05) & (abs(log2FoldChange)>2))

#annot <- read.csv("mart_export.csv", sep=',', stringsAsFactors = FALSE)
annot <- read.csv(system.file("extdata", "mart_export.csv", package = "LErNet"), 
                  sep=',', stringsAsFactors = FALSE)
mart_snapshot <- annot
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

de_lncrnas <- subset(annot,
                     ( transcript_biotype=='lncRNA' | transcript_biotype=='processed_transcript' ) &
                       (transcript_length>=200) &
                       (ensembl_gene_id_version %in% filtered_table4$Gid)
)
de_pc_elements <- subset(annot, (transcript_biotype == 'protein_coding') & 
                           (ensembl_gene_id_version %in% filtered_table4$Gid)  )
de_pcgenes <- unique(de_pc_elements$ensembl_gene_id)
de_proteins <- unique(de_pc_elements$ensembl_peptide_id)
de_proteins_symbols <- unique(de_pc_elements$hgnc_symbol)

#gtf_file <- "gencode.v34.annotation.gtf.gz"
gtf_file <- system.file("extdata", "gencode.v34.annotation.gtf.gz", package = "LErNet")
complete_positions <- LErNet::load_gtf(gtf_file)

genomic_context <- LErNet::get_genomic_context(
  positions = complete_positions,
  lncgenes = unique(de_lncrnas$ensembl_gene_id),
  pcgenes = de_pcgenes,
  max_window = 100000)

interaction_context <- LErNet::search_in_arenaidb(
  unique(de_lncrnas$ensembl_transcript_id),
  unique(de_pc_elements$ensembl_gene_id),
  unique(de_pc_elements$ensembl_peptide_id),
  unique(de_pc_elements$hgnc_symbol)
)
ncrnas_info <- get_arenaidb_ncrnas()[, c('name','naming.resource','biotype')]
colnames(ncrnas_info) <- c('principal_name', 'principal_space', 'biotype')
x <- merge(interaction_context$alias_mapping, ncrnas_info)
table(x$biotype)
y1 <- merge(interaction_context$gene_icontext, x)
table(y1$biotype)
y2 <- merge(interaction_context$protein_icontext, x)
table(y2$biotype)
interaction_context$gene_icontext <- y1[y1$biotype != "mirna_primary_transcript",]
interaction_context$protein_icontext <- y2[y2$biotype != "mirna_primary_transcript",]

all <- read.csv("converted_ppi.txt")
ppi_disease <- list()
disease <- c(2,4,166,167,169,170,171,174,175,176,177,181,184,186,188,189,190,192,193,194,195,196,202,203,204,205,206,207,208,210,211,
             214,221,222,223,224,234,235,236,237,239,240,241,242,243,244)
subsetdisease <- all[,disease]
k <- 1
for(l in 3:ncol(subsetdisease)){
  ppi_disease[[k]] <- subsetdisease[subsetdisease[,l]=='1',c(1,2)]
  k <- k+1
}
labels <- colnames(subsetdisease)[3:ncol(subsetdisease)]
names(ppi_disease) <- labels

results <- as.data.frame(matrix(,ncol = 2))
colnames(results) <- c("pathway", "rateo")
results <- na.omit(results)

de_pcgene_icontext <- interaction_context$gene_icontext
de_pcgene_aliases <- LErNet::get_arenaidb_gene_aliases()
#de_pcgene_aliases <- subset(interaction_context$alias_mapping[interaction_context$alias_mapping$entity_type == 'gene' ,], select=-c(entity_type))
nrow(de_pcgene_icontext)


de_protein_icontext <- interaction_context$protein_icontext##
de_protein_aliases <- LErNet::get_arenaidb_protein_aliases()
#de_protein_aliases <- subset(interaction_context$alias_mapping[interaction_context$alias_mapping$entity_type == 'protein' ,], select=-c(entity_type))
nrow(de_protein_icontext)






for(i in 1:length(ppi_disease)){
  ppi_network <- ppi_disease[[i]]
  ###########################
  print(names(ppi_disease[i]))
  ###########################
  colnames(ppi_network) <- c("protein1", "protein2")
  ensp_to_ensg <- subset(annot, ensembl_peptide_id %in% unique(union(ppi_network$protein1, ppi_network$protein2)) )[, c('ensembl_peptide_id','ensembl_gene_id')]
  
  seeds <- list() #ENSP format
  # genomic seeds
  # interaction seeds - genes
  # interaction seeds - proteins
  
  # genomic seeds
  genomic_seeds <- unique((merge(genomic_context, ensp_to_ensg, by.x='partner_coding_gene', by.y='ensembl_gene_id'))$ensembl_peptide_id)
  print("genomic seeds")
  print(length(genomic_seeds))
  
  # interaction seeds - genes
  ipcgene_seeds <- unique(de_pcgene_icontext[, c('mate.name','mate.naming.resource')])
  nrow(ipcgene_seeds)
  colnames(ipcgene_seeds) <- c('principal_name', 'principal_space')
  ipcgene_seeds <- (merge(ipcgene_seeds, de_pcgene_aliases))$alias_name
  length(ipcgene_seeds)
  ipcgene_seeds <- ensp_to_ensg[ ensp_to_ensg$ensembl_gene_id %in% ipcgene_seeds, ]$ensembl_peptide_id
  print("ipcgene_seeds")
  print(length(ipcgene_seeds))
  
  # interaction seeds - proteins
  iprotein_seeds <- unique(de_protein_icontext[, c('mate.name','mate.naming.resource')])##33
  nrow(iprotein_seeds)
  colnames(iprotein_seeds) <- c('principal_name', 'principal_space')
  all_iprotein_seeds <- (merge(iprotein_seeds, de_protein_aliases))
  nrow(all_iprotein_seeds)##33
  iprotein_seeds <- ensp_to_ensg[ ensp_to_ensg$ensembl_gene_id %in% all_iprotein_seeds$alias_name, ]$ensembl_peptide_id
  nrow(iprotein_seeds)##NULL -> BECAUSE THEY AREN'T IN ENSEMBL_GENE BUT IN REFSEQ,UNIGENE,UNIPROT,GENERIC AND SYMBOL.
  
  #SYMBOL
  s <- unique(subset(all_iprotein_seeds,alias_space=="symbol")$alias_name)
  x <-  unique(mart_snapshot[ mart_snapshot$hgnc_symbol %in% s, ]$ensembl_peptide_id)
  x <- x[x != ""]
  length(x)##225
  
  #UNIGENE
  conv_protein <- matrix(,nrow = 0, ncol = 2)
  
  n <- unique(subset(all_iprotein_seeds,alias_space=="unigene")$alias_name)
  src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
  #mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  #' select
  for(alias_name in n){
    b <- select(src, columns = c("ensembl"),
                keytype = "unigene", keys = alias_name)
    
    if(nrow(b)!= 0){
      conv_protein <- rbind(conv_protein,b)
      break
    }
  }
  colnames(conv_protein)<-c("alias_name","ensembl_gene_id")
  
  #REFSEQ-UNIPROT-GENERIC
  fi <- c("external_gene_name","external_synonym", "refseq_mrna","uniprot_gn_id")
  names(fi) <- fi
  m <- unique(subset(all_iprotein_seeds,alias_space!="unigene" & alias_space!="symbol")$alias_name)
  
  for(alias_name in m){
    print(alias_name)
    for(i in fi){
      b <- getBM(attributes = c("ensembl_gene_id"),
                 filters = i, values = alias_name, mart = mart)
      
      if(nrow(b)!= 0){
        b <- cbind(alias_name,b)
        conv_protein <- rbind(conv_protein,b)
        break
      }
    }
  }
  
  conv_protein <- as.data.frame(unique(conv_protein))
  ens_prot <- unique(mart_snapshot[mart_snapshot$ensembl_gene_id %in% conv_protein$ensembl_gene_id,]$ensembl_peptide_id)
  
  length(ens_prot)##201
  
  
  iprotein_seeds <- unique(c(ens_prot, x))##4251conSOLOTRANSCRIPTS
  print("iprotein_seeds")
  print(length(iprotein_seeds))
  
  seeds <- unique(c(genomic_seeds, ipcgene_seeds, iprotein_seeds))
  print("unique seeds")
  print(length(seeds))
  
  ## lncrna
  lncrna_context <- unique((merge(genomic_context, mart_snapshot, by.x='partner_coding_gene', by.y='ensembl_gene_id'))[,c('lnc_known', 'ensembl_peptide_id')])
  colnames(lncrna_context) <- c('lncrna','mate')
  length(unique(lncrna_context$lncrna))##
  nrow(lncrna_context)##
  
  #de_pcgene_icontext
  ipcgene_seeds <- unique(de_pcgene_icontext[, c('mate.name','mate.naming.resource')])
  colnames(ipcgene_seeds) <- c('principal_name', 'principal_space')
  ipcgene_seeds <- (merge(ipcgene_seeds, de_pcgene_aliases))
  ipcgene_seeds <- merge(ipcgene_seeds, ensp_to_ensg, by.x='alias_name', by.y='ensembl_gene_id')
  
  de_pcgene_context <- ipcgene_seeds[,c('principal_name', 'ensembl_peptide_id')]
  nrow(de_pcgene_context)##26
  colnames(de_pcgene_context) <- c('lncrna','mate')
  length(unique(de_pcgene_context$lncrna))##25
  lncrna_context <- rbind(lncrna_context, de_pcgene_context)##9806
  nrow(lncrna_context)##
  
  # de protein icontext
  iprotein_seeds <- unique(de_protein_icontext[, c('mate.name','mate.naming.resource')])
  colnames(iprotein_seeds) <- c('principal_name', 'principal_space')
  all_iprotein_seeds <- (merge(iprotein_seeds, de_protein_aliases))
  y <- unique((merge( all_iprotein_seeds, mart_snapshot, by.x='alias_name', by.y='hgnc_symbol' ))[, c('alias_name', 'ensembl_peptide_id')])
  y <- y[y$ensembl_peptide_id != "",]##225
  colnames(y) <- c('lncrna','mate')
  length(unique(y$lncrna))##33
  lncrna_context <- unique(rbind(lncrna_context, y))##10024
  
  
  y2 <- unique((merge(conv_protein,mart_snapshot, by.x='ensembl_gene_id', by.y='ensembl_gene_id'))[,c('alias_name','ensembl_peptide_id')])
  y2 <- y2[y2$ensembl_peptide_id != "",]##200
  colnames(y2) <- c('lncrna','mate')
  length(unique(y2$lncrna))##29
  lncrna_context <- unique(rbind(lncrna_context, y2))##10224
  
  
  lncrna_context <- unique(lncrna_context[ lncrna_context$mate %in% seeds,])
  lncrna.context <- lncrna_context[lncrna_context$mate != "",]
  de.proteins <- de_proteins
  network_seeds <- seeds
  bg_ppi_network <- ppi_network
  
  
  labels <- data.frame( matrix( ncol=2, nrow=0, dimnames=list(NULL, c('id','label')) ) )
  x <- data.frame(
    'id' = unique(lncrna.context$lncrna),
    'label' = unique(lncrna.context$lncrna)
  )
  labels <- rbind(labels, x)
  
  # ## WITHOUT EXPANSION
  # ppi <- ppi_network##324152
  # if(!is.null(de_proteins)){
  #   ppi <- ppi[ ppi$protein1 %in% de_proteins, ]
  #   ppi <- ppi[ ppi$protein2 %in% de_proteins, ]
  # }
  # nrow(ppi)##12981
  # 
  # seed_ppi <- ppi[ ppi$protein1 %in% seeds, ]
  # seed_ppi <- seed_ppi[ seed_ppi$protein2 %in% seeds, ]
  # nrow(seed_ppi)##1094
  # 
  # init_components <- get_connected_components(seed_ppi)##24
  # 
  # for(seed in setdiff( seeds, unlist(init_components) )){
  #   init_components <- append(init_components, list(seed) )
  # }##455
  # 
  # print("init_components")
  # print(length(init_components))
  # com_distribution <- as.data.frame(summary(init_components))
  # com_distribution <- as.numeric(com_distribution[1:length(init_components),3])
  # 
  # coDistribution <- as.data.frame(table(com_distribution))
  # g <- as.numeric(coDistribution$Freq)
  # 
  # print(barplot(g,data=coDistribution, names.arg = coDistribution$com_distribution,
  #               main="Histogram of the distribution of the components sizes", xlab = "Components sizes"))
  # 
  # print(coDistribution)
  # 
  # expanded_elements <- unlist(init_components)
  
  # x <- unique(mart_snapshot[ mart_snapshot$ensembl_peptide_id %in% unlist(init_components), c('ensembl_peptide_id','hgnc_symbol')])
  # colnames(x) <- c('id','label')
  # x <- x[x$id != "",]
  # labels <- rbind(labels, x)
  # 
  # Vis <- LErNet::visualize(
  #   unique(lncrna.context),
  #   unique(de.proteins),
  #   unique(network_seeds),
  #   unique(expanded_elements),
  #   unique(bg_ppi_network),
  #   labels)
  # 
  # edges <- (Vis[["x"]][["edges"]])
  # print("edges")
  # print(nrow(edges))
  # node <- (Vis[["x"]][["nodes"]])
  # node1 <- (node[node$group=="Connector","name"])
  # node2 <- (node[node$group=="Seed","name"])
  # nodes <- unique(c(node1,node2))
  # print("nodes")
  # print(length(nodes))
  # node3 <- (node[node$group=="lncRNA","name"])
  # nodes <- unique(c(nodes,node3))
  # print("nodes with lncRNA")
  # print(length(nodes))
  
  components <- LErNet::expand_seeds(seeds,  ppi_network,  de_proteins)
  # print("components")
  # print(length(components))
  # connected component distribution
  # expanded_distribution <- as.data.frame(summary(components))
  # expanded_distribution <- expanded_distribution[1:length(components),3]
  # exDistribution<- as.data.frame(table(expanded_distribution))
  # l <- as.numeric(exDistribution$Freq)
  # 
  # print(barplot(l,data=exDistribution,names.arg = exDistribution$expanded_distribution, 
  #               main="Expanded network histogram of the distribution of the components sizes",xlab = "Components sizes"))
  # print(exDistribution)
  # 
  labels <- data.frame( matrix( ncol=2, nrow=0, dimnames=list(NULL, c('id','label')) ) )
  x <- data.frame(
    'id' = unique(lncrna.context$lncrna),
    'label' = unique(lncrna.context$lncrna)
  )
  labels <- rbind(labels, x)
  x <- unique(mart_snapshot[ mart_snapshot$ensembl_peptide_id %in% unlist(components), c('ensembl_peptide_id','hgnc_symbol')])
  colnames(x) <- c('id','label')
  x <- x[x$id != "",]
  labels <- rbind(labels, x)
  
  expanded_elements <- unlist(components)
  
  Vis <- LErNet::visualize(
    unique(lncrna.context),
    unique(de.proteins),
    unique(network_seeds),
    unique(expanded_elements),
    unique(bg_ppi_network),
    labels)
  
  # edges <- (Vis[["x"]][["edges"]])
  # print("edges")
  # print(nrow(edges))
  # node <- (Vis[["x"]][["nodes"]])
  # node1 <- (node[node$group=="Connector","name"])
  # node2 <- (node[node$group=="Seed","name"])
  # nodes <- unique(c(node1,node2))
  # print("nodes")
  # print(length(nodes))
  # node3 <- (node[node$group=="lncRNA","name"])
  # nodes <- unique(c(nodes,node3))
  # print("nodes with lncRNA")
  # print(length(nodes))
  
  
  
  node <- (Vis[["x"]][["nodes"]])
  node1 <- (node[node$group=="Connector","name"])
  node2 <- (node[node$group=="Seed","name"])
  nodes <- length(unique(c(node1,node2)))
  
  for(comp in components){
    entrez_ids <- unique(annot[annot$ensembl_peptide %in% comp,]$entrezgene_id)
    enrichment <- LErNet::enrich(entrez_ids, 'human')
    if(!(typeof(enrichment) == "NULL")){
      rateo <- length(comp)/nodes
      path <- (enrichment@result[1,2])
      print(c(path, rateo))
      if(path %in% results[,1]){
        results[results$pathway==path,2] <- results[results$pathway==path,2] + rateo
      }else{
        add <- as.data.frame(matrix(nrow = 1,ncol=2))
        add[1,1] <- path
        add[1,2] <- rateo
        colnames(add) <- colnames(results)
        results <- rbind(results, add)
      }
      
    }
  }
  rm(list = c("ppi_network","ensp_to_ensg", "seeds", "genomic_seeds","ipcgene_seeds", "iprotein_seeds", "all_iprotein_seeds", "x", "init_components","components", "comp", "entrez_ids", "enrichment", "lncrna_context", "de_pcgene_context", "y1", "y2", "y", "bg_ppi_network", "expanded_elements", "labels", "Vis", "network_seeds")
  )
}

write.csv(results, "results_diseasewithnotEX.csv", row.names = F)











