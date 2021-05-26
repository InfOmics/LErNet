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

gtf_file <- system.file("extdata","gencode.v34.annotation.gtf.gz", package = "LErNet")

#from a total of 30863 samples, we filtered by Genetype ("lincRNA", "processed_transcript","lncRNA") and length >= 200, obtaining
#a total of 15541. Moreover we filtered even by coloumns of the counts thus resulting 5994 samples. After performing DEseq2
#we subsetted the "results" by abs(log2FC) > 1.5 and pvalue and padj <= 0.05 getting 679 samples
all.lncrna <- read.csv("all_lncrna_filterd.csv", stringsAsFactors = F)

#from a total of 45006 samples, we filtered by Genetype ("protein_coding") obtaining 44979. Moreover we filtered even by coloumns of
#the counts thus resulting 38386. After performing DEseq2, we subsetted the "results" by abs(log2FC) > 1.5 and pvalue and padj <= 0.05
#getting 6070 samples
all.pcrna <- read.csv("all_pcrna_filterd.csv", stringsAsFactors = F)

#from a total of 7418 samples, we filtered by length >= 200, obtaining a total of 7404. Moreover we filtered even by coloumns of the
#counts thus resulting 4921 samples. After performing DEseq2, we subsetted the "results" by abs(log2FC) > 1.5 and pvalue and padj <= 0.05
#getting 934 samples
all.lncnovel <- read.csv("all_lncnovel_filterd.csv", stringsAsFactors = F)

all.lncrna$transcript_id <- substring(all.lncrna$X[], 1, nchar(all.lncrna$X[])-2)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#mart = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
#mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest")
mart_snapshot <- read.csv(system.file("extdata", "mart_export.csv", package = "LErNet"), sep=',', stringsAsFactors = FALSE)

annot <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "refseq_ncrna"), 
               filters = "ensembl_transcript_id", 
               mart, 
               values = all.lncrna$transcript_id)
annot2 <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "refseq_ncrna"), 
                filters = "refseq_ncrna", 
                mart, 
                values = all.lncrna$transcript_id)
annot <- unique(rbind(annot,annot2)) 

##all information about novel_lncrnas
novel <- all.lncnovel
chrs <- sapply(strsplit(sapply(strsplit( as.character(novel$Pos), "-"), `[`, 1), ":"), `[`, 1)
starts <- sapply(strsplit(sapply(strsplit(as.character(novel$Pos), "-"), `[`, 1), ":"), `[`, 2)
ends <- sapply(strsplit(as.character(novel$Pos), "-"), `[`, 2)
novel_gtf <- data.frame("id" = novel$X, "type" = rep('novel lncRNA', times = nrow(novel)),
                        "seqname" = chrs, "start" = starts, "end" = ends)


nrow(novel_gtf)

##all information about lncrnas
de_lncrnas <- unique(rbind(merge(annot, all.lncrna, by.x = 'ensembl_transcript_id', by.y = 'transcript_id'),
                    merge(annot, all.lncrna, by.x = 'refseq_ncrna', by.y = 'transcript_id')))

#mart_snapshot <- read.csv(system.file("extdata", "mart_export.csv", package = "LErNet"), sep=',', stringsAsFactors = FALSE)

all.pcrna$transcript_id <- substring(all.pcrna$X[], 1, nchar(all.pcrna$X[])-2)

transcript.pcgenes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna"), 
                           filters = "refseq_mrna", 
                           mart, 
                           values = all.pcrna$transcript_id)
##all information about pcrnas
pcrnaInfo <- unique(merge(transcript.pcgenes, all.pcrna, by.x = 'refseq_mrna', by.y = 'transcript_id'))

de_pc_elements <- subset(mart_snapshot, ensembl_gene_id %in% pcrnaInfo$ensembl_gene_id )
de_pcgenes <- unique(de_pc_elements$ensembl_gene_id)
de_proteins <- unique(de_pc_elements$ensembl_peptide_id)
de_proteins_symbols <- unique(de_pc_elements$hgnc_symbol)

complete_positions <- LErNet::load_gtf(gtf_file)
complete_positions <- rbind(complete_positions, novel_gtf)
rownames(complete_positions) <- seq(1:nrow(complete_positions))

de_lncrnas_total <- unique(c(de_lncrnas$ensembl_gene_id, novel$X))
##get_genomic_context
genomic_context <- LErNet::get_genomic_context(
  positions = complete_positions,
  lncgenes = de_lncrnas_total,
  pcgenes = de_pcgenes,
  max_window = 100000)##1492 objs

# Number of genomic seeds
length(unique(genomic_context$partner_coding_gene))
# Number of lncRNAs
length(unique(genomic_context$lnc_known))
# Mean number of genomic seeds for each lncRNA
mean(table(genomic_context$lnc_known))
# Distribution
tab_distribution <- as.data.frame(table(genomic_context$lnc_known))
head(tab_distribution)
t <- (tab_distribution[,2])
t <- as.data.frame(table(t))
d <- as.numeric(t$Freq)

barplot(d,data=t,names.arg = t$t, 
        main=" Histogram of the distribution of the seeds for each lncRNAs",xlab = "Quantity of the genomic seeds for lncRNAs") 

#search_in_arenaidb
interaction_context <- LErNet::search_in_arenaidb(
  de_lncrnas$ensembl_transcript_id,
  unique(de_pc_elements$ensembl_gene_id),
  unique(de_pc_elements$ensembl_peptide_id),
  unique(de_pc_elements$hgnc_symbol)
)

# interaction_context <- LErNet::search_in_arenaidb(
#   de_lncrnas$ensembl_transcript_id
# 
# )


stringdb_tax = 9606
stringdb_thr = 900

ppi_network <- LErNet::get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr)


ensp_to_ensg <- subset(mart_snapshot, ensembl_peptide_id %in% unique(union(ppi_network$protein1, ppi_network$protein2)) )[, c('ensembl_peptide_id','ensembl_gene_id')]


nrow(ppi_network)
nrow(ensp_to_ensg)


seeds <- c() #ENSP format
# genomic seeds
# interaction seeds - genes
# interaction seeds - proteins

# genomic seeds
genomic_seeds <- unique((merge(genomic_context, ensp_to_ensg, by.x='partner_coding_gene', by.y='ensembl_gene_id'))$ensembl_peptide_id)
length(genomic_seeds)


de_pcgene_icontext <- interaction_context$gene_icontext
de_pcgene_aliases <- LErNet::get_arenaidb_gene_aliases()
#de_pcgene_aliases <- subset(interaction_context$alias_mapping[interaction_context$alias_mapping$entity_type == 'gene' ,], select=-c(entity_type))
nrow(de_pcgene_icontext)


de_protein_icontext <- interaction_context$protein_icontext##
de_protein_aliases <- LErNet::get_arenaidb_protein_aliases()
#de_protein_aliases <- subset(interaction_context$alias_mapping[interaction_context$alias_mapping$entity_type == 'protein' ,], select=-c(entity_type))
nrow(de_protein_icontext)

# interaction seeds - genes
ipcgene_seeds <- unique(de_pcgene_icontext[, c('mate.name','mate.naming.resource')])
nrow(ipcgene_seeds)
colnames(ipcgene_seeds) <- c('principal_name', 'principal_space')
ipcgene_seeds <- (merge(ipcgene_seeds, de_pcgene_aliases))$alias_name
length(ipcgene_seeds)
ipcgene_seeds <- ensp_to_ensg[ ensp_to_ensg$ensembl_gene_id %in% ipcgene_seeds, ]$ensembl_peptide_id
length(ipcgene_seeds)


# interaction seeds - proteins
iprotein_seeds <- unique(de_protein_icontext[, c('mate.name','mate.naming.resource')])##33
nrow(iprotein_seeds)
colnames(iprotein_seeds) <- c('principal_name', 'principal_space')
all_iprotein_seeds <- (merge(iprotein_seeds, de_protein_aliases))
nrow(all_iprotein_seeds)

iprotein_seeds <- ensp_to_ensg[ ensp_to_ensg$ensembl_gene_id %in% all_iprotein_seeds$alias_name, ]$ensembl_peptide_id
nrow(iprotein_seeds)##NULL -> BECAUSE THEY AREN'T IN ENSEMBL_GENE BUT IN REFSEQ,UNIGENE,UNIPROT,GENERIC AND SYMBOL.

#SYMBOL
s <- unique(subset(all_iprotein_seeds,alias_space=="symbol")$alias_name)
x <-  unique(mart_snapshot[ mart_snapshot$hgnc_symbol %in% s, ]$ensembl_peptide_id)
x <- x[x != ""]
length(x)

#UNIGENE
conv_protein <- matrix(,nrow = 0, ncol = 2)

n <- unique(subset(all_iprotein_seeds,alias_space=="unigene")$alias_name)
src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

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


iprotein_seeds <- unique(c(ens_prot, x))
length(iprotein_seeds)


seeds <- unique(c(genomic_seeds, ipcgene_seeds, iprotein_seeds))

length(seeds)

length(de_proteins)

# genomic context
lncrna_context <- unique((merge(genomic_context, mart_snapshot, by.x='partner_coding_gene', by.y='ensembl_gene_id'))[,c('lnc_known', 'ensembl_peptide_id')])
colnames(lncrna_context) <- c('lncrna','mate')
length(unique(lncrna_context$lncrna))
nrow(lncrna_context)

#de_pcgene_icontext
ipcgene_seeds <- unique(de_pcgene_icontext[, c('mate.name','mate.naming.resource')])
colnames(ipcgene_seeds) <- c('principal_name', 'principal_space')
ipcgene_seeds <- (merge(ipcgene_seeds, de_pcgene_aliases))
ipcgene_seeds <- merge(ipcgene_seeds, ensp_to_ensg, by.x='alias_name', by.y='ensembl_gene_id')

de_pcgene_context <- ipcgene_seeds[,c('principal_name', 'ensembl_peptide_id')]
nrow(de_pcgene_context)
colnames(de_pcgene_context) <- c('lncrna','mate')
length(unique(de_pcgene_context$lncrna))
lncrna_context <- rbind(lncrna_context, de_pcgene_context)
nrow(lncrna_context)

# de protein icontext
iprotein_seeds <- unique(de_protein_icontext[, c('mate.name','mate.naming.resource')])
colnames(iprotein_seeds) <- c('principal_name', 'principal_space')
all_iprotein_seeds <- (merge(iprotein_seeds, de_protein_aliases))
y <- unique((merge( all_iprotein_seeds, mart_snapshot, by.x='alias_name', by.y='hgnc_symbol' ))[, c('alias_name', 'ensembl_peptide_id')])
y <- y[y$ensembl_peptide_id != "",]
colnames(y) <- c('lncrna','mate')
length(unique(y$lncrna))
lncrna_context <- unique(rbind(lncrna_context, y))


y2 <- unique((merge(conv_protein,mart_snapshot, by.x='ensembl_gene_id', by.y='ensembl_gene_id'))[,c('alias_name','ensembl_peptide_id')])
y2 <- y2[y2$ensembl_peptide_id != "",]
colnames(y2) <- c('lncrna','mate')
length(unique(y2$lncrna))
lncrna_context <- unique(rbind(lncrna_context, y2))


sumlncRNA <- length(unique(y$lncrna)) + length(unique(y2$lncrna))


## WITHOUT EXPANSION
ppi <- ppi_network##324152
if(!is.null(de_proteins)){
  ppi <- ppi[ ppi$protein1 %in% de_proteins, ]
  ppi <- ppi[ ppi$protein2 %in% de_proteins, ]
}
nrow(ppi)##12981

seed_ppi <- ppi[ ppi$protein1 %in% seeds, ]
seed_ppi <- seed_ppi[ seed_ppi$protein2 %in% seeds, ]
nrow(seed_ppi)##1094

init_components <- get_connected_components(seed_ppi)##24

for(seed in setdiff( seeds, unlist(init_components) )){
  init_components <- append(init_components, list(seed) )
}


com_distribution <- as.data.frame(summary(init_components))
com_distribution <- as.numeric(com_distribution[1:3837,3])

coDistribution <- as.data.frame(table(com_distribution))
g <- as.numeric(coDistribution$Freq)

barplot(g,data=coDistribution, names.arg = coDistribution$com_distribution,
        main="Histogram of the distribution of the components sizes", xlab = "Components sizes")

##EXPANSION
components <- LErNet::expand_seeds(seeds,  ppi_network,  de_proteins)
length(unlist(components))
length(components)
# connected component distribution
expanded_distribution <- as.data.frame(summary(components))
expanded_distribution <- expanded_distribution[1:3769,3]
exDistribution<- as.data.frame(table(expanded_distribution))
l <- as.numeric(exDistribution$Freq)

barplot(l,data=exDistribution,names.arg = exDistribution$expanded_distribution, 
                     main="Expanded network histogram of the distribution of the components sizes",xlab = "Components sizes") 

###only DE
lncrna_context <- lncrna_context[ lncrna_context$mate %in% seeds,]

lncrna.context <- unique(lncrna_context)
lncrna.context <- lncrna.context[lncrna.context$mate != "",]
de.proteins <- de_proteins
network_seeds <- seeds
bg_ppi_network <- ppi_network
expanded_elements <- unlist(components)


labels <- data.frame( matrix( ncol=2, nrow=0, dimnames=list(NULL, c('id','label')) ) )
x <- data.frame(
  'id' = unique(lncrna.context$lncrna),
  'label' = unique(lncrna.context$lncrna)
)
labels <- rbind(labels, x)

x <- unique(mart_snapshot[ mart_snapshot$ensembl_peptide_id %in% unlist(components), c("ensembl_peptide_id","hgnc_symbol")])
colnames(x) <- c('id','label')
x <- x[x$id != "",]
labels <- unique(rbind(labels, x))


Vis <- LErNet::visualize(
  lncrna.context,
  de.proteins,
  network_seeds,
  expanded_elements,
  unique(bg_ppi_network),
  labels)

edges <- (Vis[["x"]][["edges"]])
node <- (Vis[["x"]][["nodes"]])
node1 <- (node[node$group=="Connector","name"])
node2 <- (node[node$group=="Seed","name"])
nodes <- unique(c(node1,node2))
node3 <- (node[node$group=="lncRNA","name"])
nodes <- unique(c(nodes,node3))

library(ReactomePA)
#for human
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(enrichplot)

## rateo
node <- (Vis[["x"]][["nodes"]])
node1 <- (node[node$group=="Connector","name"])
node2 <- (node[node$group=="Seed","name"])
nodes <- length(unique(c(node1,node2)))

results <- as.data.frame(matrix(,ncol = 2))
colnames(results) <- c("pathway", "rateo")
results <- na.omit(results)

#h <- 1
for(comp in components){
  #print(h)
  entrez_ids <- unique((merge( data.frame('ensembl_peptide_id' = comp), mart_snapshot ))$entrezgene_id)
  enrichment <- LErNet::enrich(entrez_ids, 'human')
  if(!(typeof(enrichment) == "NULL")){
    rateo <- length(comp)/nodes
    path <- (enrichment@result[1,2])
    #print(path)
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
  #h <- h + 1
}
write.csv(results, "resultsAllcompnsDE_hantaan.csv", row.names = F)

## enrichments
#com <- list()
#x <- 1
#for(comp in components){
#  df <- data.frame()
#  entrez_ids <- unique(annot[annot$ensembl_peptide %in% comp,]$entrezgene_id)
#  enrichment <- LErNet::enrich(entrez_ids, 'human')
  #View(enrichment@result[1:20,])
#  if(!is.null(enrichment)){
#    #print(enrichment)
#    com[[x]]<-barplot(enrichment, title = x, font.size = 8)  
#    plot(com[[x]])
#    print(paste(c(x,":"), collapse = " "))
#    print(enrichment@result[1:20,c(2,6)])
#  }else{
#   print("_______________________________________________________")
#    print(paste(c(x,": NULL"), collapse = " "))
#    print("_______________________________________________________")
#  }
#  x <- x+1
#}

## first component's enrichment
comp <- components[[1]]
entrez_ids <- unique((merge( data.frame('ensembl_peptide_id' = comp), mart_snapshot ))$entrezgene_id)
enrichment <- LErNet::enrich(entrez_ids, 'human')
print(enrichment)
barplot(enrichment)

## DE protein's enrichment
entrez_ids <- unique((merge( data.frame('ensembl_peptide_id' = de_proteins), mart_snapshot ))$entrezgene_id)
enrichment <- LErNet::enrich(entrez_ids, 'human')
print(enrichment)
barplot(enrichment)

