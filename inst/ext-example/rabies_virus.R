
lncrna_file <- system.file("extdata", "mouse_rabies_virus", "41598_2018_30359_MOESM2_ESM.csv", package = "LErNet")
pcrna_file <- system.file("extdata", "mouse_rabies_virus","41598_2018_30359_MOESM3_ESM.csv", package = "LErNet")
gtf_file <- system.file("extdata", "mouse_rabies_virus","gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz", package = "LErNet")

#lncrna_file <- "extdata/mouse_rabies_virus/41598_2018_30359_MOESM2_ESM.csv"
##pcrna_file <- "extdata/mouse_rabies_virus/41598_2018_30359_MOESM3_ESM.csv"
#gtf_file <- "extdata/mouse_rabies_virus/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz"


pcgenes<-(read.csv(pcrna_file, stringsAsFactors=FALSE))$gene_id

length(pcgenes)



lncrnaInfo<-read.csv(lncrna_file, stringsAsFactors=FALSE)
lncrnaInfo<-lncrnaInfo[lncrnaInfo$significant != 'FALSE',  ]
lncrnaAll<-as.character(lncrnaInfo$gene_id)

nrow(lncrnaInfo)
length(lncrnaAll)


complete_positions <- LErNet::load_gtf(gtf_file)
nrow(complete_positions)


# Extract the novel lncRNAs from the dataframe "lncrnaInfo"
novel<-lncrnaInfo[lncrnaInfo$isoform_status == "lncRNA_Novel", ]
nrow(novel)

# Some elaboration to extract the necessary information about lncRNAs
chrs <- paste0("chr",sapply(strsplit(sapply(strsplit( novel$locus, "-"), `[`, 1), ":"), `[`, 1))
starts <- sapply(strsplit(sapply(strsplit(novel$locus, "-"), `[`, 1), ":"), `[`, 2)
ends <- sapply(strsplit(novel$locus, "-"), `[`, 2)
novel_gtf <- data.frame( "id" = novel$gene_id, "type" = rep('novel lncRNA', times = nrow(novel)),
                         "seqname" = chrs, "start" = starts, "end" = ends )


nrow(novel_gtf)

# Add information of novel lncRNAs into the dataframe containing information about known genes/lncRNAs
complete_positions <- rbind(complete_positions, novel_gtf)
rownames(complete_positions) <- seq(1:nrow(complete_positions))


nrow(complete_positions)


library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
annot<-getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_peptide_id",
                            "hgnc_symbol", 'transcript_biotype', 'transcript_length', 'entrezgene_id', "mgi_symbol"),  mart = mart)
#write.csv(annot, system.file("extdata", "mouse_mart_export.csv", package = "LErNet"), row.names = FALSE)
#annot <- read.csv(system.file("extdata", "mouse_mart_export.csv", package = "LErNet"), sep=',', stringsAsFactors = FALSE)



library(GenomicRanges)
genomic_context <- LErNet::get_genomic_context(
  positions = complete_positions,
  lncgenes = lncrnaAll,
  pcgenes = pcgenes,
  max_window = 100000)


print(nrow(genomic_context))

# Number of genomic seeds
length(unique(genomic_context$partner_coding_gene))

# Mean number of genomic seeds for each lncRNA
mean(table(genomic_context$lnc_known))



stringdb_tax = 10090
stringdb_thr = 900

library(STRINGdb)
library(igraph)


ppi_network <- LErNet::get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr)
ensp_to_ensg <- subset(annot, ensembl_peptide_id %in% unique(union(ppi_network$protein1, ppi_network$protein2)) )[, c('ensembl_peptide_id','ensembl_gene_id')]

print(nrow(ppi_network))
print(nrow(ensp_to_ensg))


# string10file <- system.file("extdata", "mouse_rabies_virus","10090.protein.links.v10.txt", package = "LErNet")
# string10 <- read.csv(string10file, stringsAsFactors=FALSE, sep=' ')
# string10$combined_score <- as.numeric(string10$combined_score)
# ppi10 <- string10[string10$combined_score >= 900, c('protein1','protein2')]
# ppi10$protein1<-substr(ppi10$protein1,nchar(stringdb_tax)+2,nchar(ppi10$protein1[1]))
# ppi10$protein2<-substr(ppi10$protein2,nchar(stringdb_tax)+2,nchar(ppi10$protein2[1]))
# ppi_network <- ppi10
# ensp_to_ensg <- subset(annot, ensembl_peptide_id %in% unique(union(ppi_network$protein1, ppi_network$protein2)) )[, c('ensembl_peptide_id','ensembl_gene_id')]


length(unique(genomic_context$partner_coding_gene))

genomic_seeds <- unique((merge(genomic_context, ensp_to_ensg, by.x='partner_coding_gene', by.y='ensembl_gene_id'))$ensembl_peptide_id)
length(genomic_seeds)

length(pcgenes)
de_proteins <- merge(annot[annot$ensembl_peptide_id != '', ], data.frame('ensembl_gene_id' = unique(pcgenes)))$ensembl_peptide_id
length(de_proteins)

components <- LErNet::expand_seeds(genomic_seeds,  ppi_network,  de_proteins)

length(unlist(components))
length(components[[1]])
components[[1]]

lncrna_context <- unique( (merge(genomic_context, annot[annot$ensembl_peptide_id != '',], by.x='partner_coding_gene', by.y='ensembl_gene_id'))[ , c('lnc_known','ensembl_peptide_id')])
colnames(lncrna_context) <- c('lncrna','mate')
lncrna_context <- unique(lncrna_context[ lncrna_context$mate %in% genomic_seeds,])

labels <- data.frame(
  'id' = unique(lncrna_context$lncrna),
  'label' = unique(lncrna_context$lncrna)
)

x <- (merge(data.frame('ensembl_peptide_id' = as.character(unique(unlist(components)))), annot))[, c('ensembl_peptide_id','mgi_symbol')]
colnames(x) <- c('id','label')

labels <- rbind(labels, x)

setdiff(unlist(components), de_proteins)


LErNet::visualize(
  lncrna_context,
  de_proteins,
  genomic_seeds,
  unlist(components),
  ppi_network,
  labels)

library(ReactomePA)

#BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

# for human
# BiocManager::install("org.Hs.eg.db")
# library(org.Hs.eg.db)

comp <- components[[1]]
entrez_ids <- unique((merge( data.frame('ensembl_peptide_id' = comp), annot ))$entrezgene_id)
enrichment <- LErNet::enrich(entrez_ids, 'mouse')
#print(enrichment)

# for human: organism = "human"
barplot(enrichment)
