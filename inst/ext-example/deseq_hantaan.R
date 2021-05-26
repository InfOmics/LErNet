library(LErNet)
library(GenomicRanges)
library(DESeq2)

########### LNCRNA_KNOWN
lncrna_known_file <- system.file("extdata","homo_hantaan", "original", "gse133751_lncrna-exp.txt", package = "LErNet")

lncrna_known <- read.csv(lncrna_known_file, sep = '\t', stringsAsFactors = FALSE)

unique(lncrna_known$RNA_type) 

lncrna_known_subset <- subset(lncrna_known, (Length >= 200) &
                                (Gene_type=='lincRNA' |
                                Gene_type=='lncRNA' |
                                Gene_type=='processed_transcript'))

unique(lncrna_known_subset$Gene_type) 


dati_lncrna <- lncrna_known_subset[,6:11]
rownames(dati_lncrna)<- lncrna_known_subset[,1]
exp <- matrix(,ncol = 6)
#' we take only the samples with at least two columns higher than 0 for untreated and treated
for (i in 1:nrow(dati_lncrna)) {
  d <- 0
  if(dati_lncrna[i,1]>= 1){
    d <- d+1
  }
  if(dati_lncrna[i,2]>= 1){
    d <- d+1
  }
  if(dati_lncrna[i,3]>= 1){
    d <- d+1
  }

  d2 <- 0
  if(dati_lncrna[i,4]>= 1){
    d2 <- d2+1
  }
  if(dati_lncrna[i,5]>= 1){
    d2 <- d2+1
  }
  if(dati_lncrna[i,6]>= 1){
    d2 <- d2+1
  }

  if(d >= 2 ){
    if(d2 >= 2){
      exp <- rbind(exp, as.matrix(dati_lncrna[i,]))
    }
  }
}

exp_lncrna <- exp[2:nrow(exp),]
df <- DataFrame(rep(c("exp","con"),c(3,3)))
colnames(df) <- c("status")
rownames(df) <- colnames(lncrna_known_subset)[6:11]
df
ddsMat <- DESeqDataSetFromMatrix(countData = exp_lncrna,
                                 colData = df,
                                 design = ~ status)

dds <- DESeq(ddsMat)

res <- results(dds)
res_lncrna_filterd <- subset(
  res,
  (abs(log2FoldChange) > 1.5 |
  pvalue <= 0.05 |
  padj <= 0.05)  
)

########## DE_PCGENES/MIRNA

pcrna_file <- system.file("extdata","homo_hantaan", "original", "gse133751_mrna-exp.txt", package = "LErNet")

pcrna <- read.csv(pcrna_file, sep = '\t', stringsAsFactors = FALSE)

unique(pcrna$RNA_type)
unique(pcrna$Gene_type)

pcrna_subset <- subset(pcrna, Gene_type == 'protein_coding')

dati_pcgenes <- pcrna_subset[,6:11]
rownames(dati_pcgenes)<- pcrna_subset[,1]
exp <- matrix(,ncol = 6)
#' we take only the samples with at least two columns higher than 0 for untreated and treated
for (i in 1:nrow(dati_pcgenes)) {
  d <- 0
  if(dati_pcgenes[i,1]>= 1){
    d <- d+1
  }
  if(dati_pcgenes[i,2]>= 1){
    d <- d+1
  }
  if(dati_pcgenes[i,3]>= 1){
    d <- d+1
  }

  d2 <- 0
  if(dati_pcgenes[i,4]>= 1){
    d2 <- d2+1
  }
  if(dati_pcgenes[i,5]>= 1){
    d2 <- d2+1
  }
  if(dati_pcgenes[i,6]>= 1){
    d2 <- d2+1
  }

  if(d >= 2 ){
    if(d2 >= 2){
      exp <- rbind(exp, as.matrix(dati_pcgenes[i,]))
    }
  }
}

exp_pcrna <- exp[2:nrow(exp),]
df <- DataFrame(rep(c("exp","con"),c(3,3)))
colnames(df) <- c("status")
rownames(df) <- colnames(pcrna_subset)[6:11]
df
ddsMat <- DESeqDataSetFromMatrix(countData = exp_pcrna,
                                 colData = df,
                                 design = ~ status)

dds <- DESeq(ddsMat)

res <- results(dds)
res_pcrna_filterd <- subset(
  res,
  (abs(log2FoldChange) > 1.5 |
  pvalue <= 0.05 |
  padj <= 0.05)  
)

################# NEWLNCRNA

lncrna_novel_file <- system.file("extdata","homo_hantaan", "original", "gse133751_new_lncrna-exp.txt", package = "LErNet")

lncrna_novel <- read.csv(lncrna_novel_file, sep = '\t', stringsAsFactors = FALSE)

lncrna_novel_subset <- subset(lncrna_novel, Length >= 200)


dati_lncnovel <- lncrna_novel_subset[,5:10]
rownames(dati_lncnovel)<- lncrna_novel_subset[,1]
exp <- matrix(,ncol = 6)
#' we take only the samples with at least two columns higher than 0 for untreated and treated
for (i in 1:nrow(dati_lncnovel)) {
  d <- 0
  if(dati_lncnovel[i,1]>= 1){
    d <- d+1
  }
  if(dati_lncnovel[i,2]>= 1){
    d <- d+1
  }
  if(dati_lncnovel[i,3]>= 1){
    d <- d+1
  }

  d2 <- 0
  if(dati_lncnovel[i,4]>= 1){
    d2 <- d2+1
  }
  if(dati_lncnovel[i,5]>= 1){
    d2 <- d2+1
  }
  if(dati_lncnovel[i,6]>= 1){
    d2 <- d2+1
  }

  if(d >= 2 ){
    if(d2 >= 2){
      exp <- rbind(exp, as.matrix(dati_lncnovel[i,]))
    }
  }
}

exp_lncnovel <- exp[2:nrow(exp),]
df <- DataFrame(rep(c("exp","con"),c(3,3)))
colnames(df) <- c("status")
rownames(df) <- colnames(lncrna_novel_subset)[5:10]
df
ddsMat <- DESeqDataSetFromMatrix(countData = exp_lncnovel,
                                 colData = df,
                                 design = ~ status)

dds <- DESeq(ddsMat)

res <- results(dds)
res_lncnovel_filterd <- subset(
  res,
  (abs(log2FoldChange) > 1.5 |
  pvalue <= 0.05 |
  padj <= 0.05)  
)

#################
write.csv(res_lncrna_filterd, "lncrna_filtered.csv", row.names = T)
write.csv(res_pcrna_filterd, "pcrna_filtered.csv", row.names = T)
write.csv(res_lncnovel_filterd, "lncnovel_filtered.csv", row.names = T)
#################

lncnovel_filterd <- read.csv("lncnovel_filtered.csv")
All_lncnovel_filterd <- merge(lncnovel_filterd, lncrna_novel_subset, by.x = 'X', by.y = 'tracking_id')
write.csv(All_lncnovel_filterd, "all_lncnovel_filterd.csv", row.names = F)

################

pcrna_filterd <- read.csv("pcrna_filtered.csv")
All_pcrna_filterd <- merge(pcrna_filterd, pcrna_subset, by.x = 'X', by.y = 'tracking_id')
write.csv(All_pcrna_filterd, "all_pcrna_filterd.csv", row.names = F)

###############

lncrna_filterd <- read.csv("lncrna_filtered.csv")
All_lncrna_filterd <- merge(lncrna_filterd, lncrna_known_subset, by.x = 'X', by.y = 'tracking_id')
write.csv(All_lncrna_filterd, "all_lncrna_filterd.csv", row.names = F)




