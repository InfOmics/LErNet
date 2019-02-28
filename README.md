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

We report a step-by-step example to execute LErNet on the data provided by Zhao et al. The dataset is composed by a list of differentially expressed genes and long non-coding RNA (lncRNA). Original excel files are provided with the LErNet package in order to correctly execute the analysis. Further, a GTF file (from ENCODE database) is provided to retrieve genomic context of genes and lncRNAs.

To run the example it's necessary to install and load the following libraries:

```
library(R.utils)
library(xlsx)
library(biomaRt)

```
The first step of the analysis is to retrieve a set of genes and lncRNAs of interest and the information of the genomic context. In the folowing lines of code DE genes and lncRNAs are obtained directly from the excel files provided by LErNet and loaded after several preprocess operations. 
 

```
lncrna_file <- system.file("extdata", "41598_2018_30359_MOESM2_ESM.xlsx", package = "LErNet")
pcrna_file <- system.file("extdata", "41598_2018_30359_MOESM3_ESM.xlsx", package = "LErNet")
gtf_file <- system.file("extdata", "gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz", package = "LErNet")

pcgenes<-read.xlsx(pcrna_file,sheetIndex = 1)
pcgenes<-as.character(pcgenes$gene_id)

lncrnaInfo<-read.xlsx(lncrna_file, sheetIndex = 1)
lncrnaInfo <- data.frame(lapply(lncrnaInfo, as.character), stringsAsFactors=FALSE)
last<-which(lncrnaInfo$significant == 'FALSE')[1]
lncrnaInfo<-lncrnaInfo[1:last-1,]
lncrnaAll<-as.character(lncrnaInfo$gene_id)

```
LErNet provides the function "load_gtf" to load into a dataframe the necessary information from a GTF file.

```
complete_positions <- LErNet::load_gtf(gtf_file)
```

It is necessary that the dataframe must contain the information for all genes and lncRNAs in input. In this example data come with information about novel lncRNAs, these information must be added to the dataframe "complete_positions":


```
# Sxtract the novel lncRNAs from the dataframe "lncrnaInfo"
novel<-lncrnaInfo
novel<-novel[novel$isoform_status == "lncRNA_Novel", ]

# Some elaboration to extract the necessary information about lncRNAs 
chrs <- paste0("chr",sapply(strsplit(sapply(strsplit( novel$locus, "-"), `[`, 1), ":"), `[`, 1))
starts <- sapply(strsplit(sapply(strsplit(novel$locus, "-"), `[`, 1), ":"), `[`, 2)
ends <- sapply(strsplit(novel$locus, "-"), `[`, 2)
novel_gtf <- data.frame( "id" = novel$gene_id, "type" = rep('novel lncRNA', times = nrow(novel)),
                         "seqname" = chrs, "start" = starts, "end" = ends )

# Add information of novel lncRNAs into the dataframe containing information about known genes/lncRNAs
complete_positions <- rbind(complete_positions, novel_gtf)
rownames(complete_positions) <- seq(1:nrow(complete_positions))
```

LErNet exploits PPI network to expand a set of protein coding genes associated with lncRNAs. In this example the database STRING is exploited to build the PPI network, however LErNet can take as input a dataframe with 2 columns containing the edges of the network. Each element of the netwoek must be identified with its ENSEMBL id. To build the network with STRING is necessary to specfy a threshold
of significance for protein interactions, the taxonomy id of the organism of interest:


```
mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
stringdb_tax = 10090
stringdb_thr = 900
```

After the function "get_stringdb" must be executed to map ENSEMBL protein ids into ENSEMBL gene ids:

```
ret <- LErNet::get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr, mart = mart)
```


```
ppi_network <- ret[["ppi_network"]]
ensp_to_ensg <- ret[["ensp_to_ensg"]]


genomic_context <- LErNet::get_genomic_context(positions = complete_positions, lncgenes = lncrnaAll, pcgenes = pcgenes, max_window = 100000, strict_genomics = TRUE)

# Number of genomic seeds
length(unique(genomic_context$partner_coding_gene))
# Mean
mean(table(genomic_context$partner_coding_gene))


annot<-getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
             filters = "ensembl_gene_id", values = unique(pcgenes), mart = mart)
strict_proteins<-annot$ensembl_peptide_id
empty<-which(strict_proteins == "")
strict_proteins<-strict_proteins[-empty]

ret <- LErNet::expand_seeds(
                genomic_context = genomic_context,
                ppi_network = ppi_network,
                ensp_to_ensg = ensp_to_ensg,
                strict_proteins = strict_proteins,
                strict_connectors = TRUE)

network_components <- ret[["network_components"]]
input_proteins <- ret[["input_proteins"]]
network_seeds <- ret[["network_seeds"]]

LErNet::visualize(
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

enrichment <- LErNet::enrich(  ens_proteins = unlist(network_components),  organism = "mouse",  mart = mart)
#LErNet::enrich(  ens_proteins = unlist(network_components),  organism = "mouse",  mart = mart, max_to_show =2)
```
