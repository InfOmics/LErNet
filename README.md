### LErNet
LErNet: characterization of lncRNAs via context-aware network expansion and enrichment analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [](#lang-en)

<hr />


## License
LErNet is distributed under the MIT license. This means that it is free for both academic and commercial use. Note however that some third party components in LErNet require that you reference certain works in scientific publications.
You are free to link or use LErNet inside source code of your own program. If do so, please reference (cite) LErNet and this website. We appreciate bug fixes and would be happy to collaborate for improvements. 
[License](https://raw.githubusercontent.com/InfOmics/LErNet/master/LICENSE)

<hr />

## Citation
If you have used the package LErNet in your project, please cite the following paper:

     Bonnici V., Caligola S., Fiorini G., Giudice L., Giugno R.
     LErNet: characterization of lncRNAs via context-aware network expansion and enrichment analysis.
     
<hr />

## Installation

```
install.packages("devtools") # if you have not installed "devtools" package
library(devtools)
install_github("InfOmics/LErNet")
```

<hr />

## Running example

We report a step-by-step example to execute LErNet on the data provided by Zhao et al. The dataset is composed by a list of differentially expressed genes and long non-coding RNA (lncRNA). Original excel files are provided with the LErNet package in order to correctly execute the analysis. Further, a GTF file is provided to retrieve genomic context of genes and lncRNAs.

To run the example..

```
library(R.utils)
library(xlsx)
library(biomaRt)

```


```
lncrna_file <- system.file("extdata", "41598_2018_30359_MOESM2_ESM.xlsx", package = "LErNet")
pcrna_file <- system.file("extdata", "41598_2018_30359_MOESM3_ESM.xlsx", package = "LErNet")
gtf_file <- system.file("extdata", "gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz", package = "LErNet")
```



```
pcgenes<-read.xlsx(pcrna_file,sheetIndex = 1)
pcgenes<-as.character(pcgenes$gene_id)
# <- list of strict pcgenes

lncrnaInfo<-read.xlsx(lncrna_file, sheetIndex = 1)
lncrnaInfo <- data.frame(lapply(lncrnaInfo, as.character), stringsAsFactors=FALSE)
last<-which(lncrnaInfo$significant == 'FALSE')[1]
lncrnaInfo<-lncrnaInfo[1:last-1,]
lncrnaAll<-as.character(lncrnaInfo$gene_id)
# <- list of strict lncgenes

complete_positions <- LErNet.load_gtf(gtf_file)

# mi ricavo i novel per ottenere le loro coordinate sul gtf
novel<-lncrnaInfo
#novel<-novel[grep(pattern = "XLOC", x = novel$gene_id), ]
novel<-novel[novel$isoform_status == "lncRNA_Novel", ]

chrs <- paste0("chr",sapply(strsplit(sapply(strsplit( novel$locus, "-"), `[`, 1), ":"), `[`, 1))
starts <- sapply(strsplit(sapply(strsplit(novel$locus, "-"), `[`, 1), ":"), `[`, 2)
ends <- sapply(strsplit(novel$locus, "-"), `[`, 2)
novel_gtf <- data.frame( "id" = novel$gene_id, "type" = rep('novel lncRNA', times = nrow(novel)),
                         "seqname" = chrs, "start" = starts, "end" = ends )

complete_positions <- rbind(complete_positions, novel_gtf)
rownames(complete_positions) <- seq(1:nrow(complete_positions))

# complete_positions




mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
stringdb_tax = 10090
stringdb_thr = 900
ret <- LErNet.get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr, mart = mart)

ppi_network <- ret[["ppi_network"]]
ensp_to_ensg <- ret[["ensp_to_ensg"]]


genomic_context <- LErNet.get_genomic_context(positions = complete_positions, lncgenes = lncrnaAll, pcgenes = pcgenes, max_window = 100000, strict_genomics = TRUE)

# Number of genomic seeds
length(unique(genomic_context$partner_coding_gene))
# Mean
mean(table(genomic_context$partner_coding_gene))


annot<-getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
             filters = "ensembl_gene_id", values = unique(pcgenes), mart = mart)
strict_proteins<-annot$ensembl_peptide_id
empty<-which(strict_proteins == "")
strict_proteins<-strict_proteins[-empty]

ret <- LErNet.expand(
                genomic_context = genomic_context,
                ppi_network = ppi_network,
                ensp_to_ensg = ensp_to_ensg,
                strict_proteins = strict_proteins,
                strict_connectors = TRUE)

network_components <- ret[["network_components"]]
input_proteins <- ret[["input_proteins"]]
network_seeds <- ret[["network_seeds"]]

LErNet.visualize(
  lncgenes = lncrnaAll,
  genomic_context = genomic_context,
  ensp_to_ensg = ensp_to_ensg,
  input_proteins = input_proteins,
  network_seeds = network_seeds,
  ppi_network = ppi_network,
  expanded_elements = unlist(network_components) ,
  mart = mart,
  mart_symbol_column = "mgi_symbol"  # "hgnc_symbol" for human
)

enrichment <- LErNet.enrich(  ens_proteins = unlist(network_components),  organism = "mouse",  mart = mart)
#LErNet.enrich(  ens_proteins = unlist(network_components),  organism = "mouse",  mart = mart, max_to_show =2)
```
