# LErNet 1.0
*LErNet*: characterization of lncRNAs via context-aware network expansion and enrichment analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [](#lang-en)

<hr />

# WARNING!!!
# THIS IS THE NEW VERSION RELEASED ON 28th OCTOBER 2020 
[the old verion can be found here](https://github.com/InfOmics/LErNet/tree/lernet.0.1)

<hr />


## License
*LErNet* is distributed under the MIT license. This means that it is free for both academic and commercial use. Note however that some third party components in *LErNet* require that you reference certain works in scientific publications.
You are free to link or use *LErNet* inside source code of your own program. If do so, please reference (cite) *LErNet* and this website. We appreciate bug fixes and would be happy to collaborate for improvements. 
[License](https://raw.githubusercontent.com/InfOmics/LErNet/master/LICENSE)

<hr />

## Citation
If you have used the package LErNet in your project, please cite the following paper:

     Bonnici V., Caligola S., Fiorini G., Giudice L., Giugno R.
     LErNet: characterization of lncRNAs via context-aware network expansion and enrichment analysis.
     In 2019 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB) (pp. 1-8). IEEE.
<hr />

## Installation
Before to install the *LErNet* package make sure that all the required packages are already installed on your computer.

```R
# if you have not installed "devtools" package
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
BiocManager::install("STRINGdb")
BiocManager::install("biomaRt")
BiocManager::install("ReactomePA")
install.packages("visNetwork")
install.packages("igraph")
install.packages("R.utils")
install.packages("rmarkdown")
```

Then, you can install the *LErNet* package.

```R
install_github("InfOmics/LErNet")
```
or 

```R
install_github("InfOmics/LErNet", ref="lernet.1.0")
```
<hr />

## Running example

We report a step-by-step example to execute LErNet on published data (Zhao et al., *Scientific reports*, 2018). The dataset is composed by a list of differentially expressed genes and long non-coding RNA (lncRNA). Original excel files are provided with the *LErNet* package in order to correctly execute the analysis. Further, a GTF file (from GENCODE database) is provided to retrieve genomic context of genes and lncRNAs.

It's necessary to install and load the following libraries to run the example:

```R
library(R.utils)
library(biomaRt)

```
The first step of the analysis is to retrieve a set of genes and lncRNAs of interest and the information of the genomic context. In the following lines of code DE genes and DE lncRNAs are obtained directly from the excel files provided within the *LErNet* package and loaded after several preprocess operations. 
 
```R
lncrna_file <- system.file("extdata", "41598_2018_30359_MOESM2_ESM.csv", package = "LErNet")
pcrna_file <- system.file("extdata", "41598_2018_30359_MOESM3_ESM.csv", package = "LErNet")
gtf_file <- system.file("extdata", "gencode.vM20.chr_patch_hapl_scaff.annotation.gtf.gz", package = "LErNet")

pcgenes<-(read.csv(pcrna_file, stringsAsFactors=FALSE))$gene_id

length(pcgenes)

lncrnaInfo<-read.csv(lncrna_file, stringsAsFactors=FALSE)
lncrnaInfo<-lncrnaInfo[lncrnaInfo$significant != 'FALSE',  ]
lncrnaAll<-as.character(lncrnaInfo$gene_id)

nrow(lncrnaInfo)
length(lncrnaAll)

```
*LErNet* provides the function `load_gtf` to load into a dataframe the necessary information from a GTF file.

```R
complete_positions <- LErNet::load_gtf(gtf_file)
nrow(complete_positions)
```

It is necessary that the dataframe contains the information for all genes and lncRNAs in input. In this example, data comes with information about novel lncRNAs. These information must be added to the dataframe `complete_positions`.


```R
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
```

*LErNet* exploits PPI network to expand a set of protein coding genes associated with lncRNAs. In this example the database STRING is exploited to build the PPI network, however *LErNet* can take as input a dataframe with 2 columns containing the edges of the network. Each element of the network must be identified with its ENSEMBL id. 
To build the network with STRING is necessary to specify a threshold of significance for protein interactions and the taxonomy id of the organism of interest.

```R
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#WARNING: on R <= 3.4 this may cause mutiple errors. Please, run it until no errors are arised.
# or, artenatively, use
#mart <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
stringdb_tax = 10090
stringdb_thr = 900

# alteratively, for human
# mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# #WARNING: on R <= 3.4 this may cause mutiple errors. Please, run it until no errors are arised.
# stringdb_tax = 9606

```

The function `get_stringdb` returns a list with a dataframe named `ppi` containing the PPI network. This dataframe can be provided to *LErNEt* without the use of STRING.

```R
library(STRINGdb)
library(igraph)

ppi_network <- LErNet::get_stringdb( stringdb_tax = stringdb_tax, stringdb_thr = stringdb_thr)

annot<-getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_peptide_id",
                            "hgnc_symbol", 'transcript_biotype', 'transcript_length', 'entrezgene_id', "mgi_symbol"),  mart = mart)

# write.csv(annot, system.file("extdata", "mouse_mart_export.csv", package = "LErNet"), row.names = FALSE)
# annot <- read.csv(system.file("extdata", "mouse_mart_export.csv", package = "LErNet"), sep=',', stringsAsFactors = FALSE)

ensp_to_ensg <- subset(annot, ensembl_peptide_id %in% unique(union(ppi_network$protein1, ppi_network$protein2)) )[, c('ensembl_peptide_id','ensembl_gene_id')]

print(nrow(ppi_network))
print(nrow(ensp_to_ensg))
```

This step is used to generate the genomic context, i.e. to find the genomic seeds necessary to run the expansion phase through th PPI network. The function to perform accomplish this task is `get_genomic_context`. The function takes in input the information retrieved from the GTF file (`complete_positions`), the list of protein coding genes and lncRNAs and a window in which to search for genomic neighbors. The function returns a dataframe containing for each lncRNA one or more partner coding genes. 

```R
library(GenomicRanges)
genomic_context <- LErNet::get_genomic_context(
  positions = complete_positions,
  lncgenes = lncrnaAll,
  pcgenes = pcgenes,
  max_window = 100000)
  
nrow(genomic_context)
```
It can be useful to show some basic statistics on the generated seeds:

```R
# Number of genomic seeds
length(unique(genomic_context$partner_coding_gene))

# Mean number of genomic seeds for each lncRNA
mean(table(genomic_context$partner_coding_gene))

```
The following lines of codes are necessary to match the seeds

```R
length(unique(genomic_context$partner_coding_gene))

genomic_seeds <- unique((merge(genomic_context, ensp_to_ensg, by.x='partner_coding_gene', by.y='ensembl_gene_id'))$ensembl_peptide_id)
length(genomic_seeds)

length(pcgenes)
de_proteins <- merge(annot[annot$ensembl_peptide_id != '', ], data.frame('ensembl_gene_id' = unique(pcgenes)))$ensembl_peptide_id
length(de_proteins)
```

The next step is the core phase of *LErNEet*, i.e. the expansion phase with the function `expand_seeds`. Expansion takes as input the genomic context, the PPI network and the list of starting proteins.
The function `expand_seeds` returns a list containing a dataframe with the network components.

```R
components <- LErNet::expand_seeds(genomic_seeds,  ppi_network,  de_proteins=NULL)

length(unlist(components))
length(components[[1]])
components[[1]]
```

This step is necessary to retrive the list of the protein mates which interact with lncRNAs of interest.

```R
lncrna_context <- unique( (merge(genomic_context, annot[annot$ensembl_peptide_id != '',], by.x='partner_coding_gene', by.y='ensembl_gene_id'))[ , c('lnc_known','ensembl_peptide_id')])

colnames(lncrna_context) <- c('lncrna','mate')

```
To go through with the next step and then to facilitate the reading of the network, we need to compute the labels.

```R
labels <- data.frame(
  'id' = unique(genomic_context$lnc_known),
  'label' = unique(genomic_context$lnc_known)
)

x <- (merge(data.frame('ensembl_peptide_id' = as.character(unique(unlist(components)))), annot))[, c('ensembl_peptide_id','mgi_symbol')]
colnames(x) <- c('id','label')

labels <- rbind(labels, x)
```

*LErNet* allows to visualize the results of the analysis through the use of the package `visNetwork`. The function `visualize` takes in input the lncrna context, the list of strict starting proteins, the network seeds, the PPI network, one or more network components extracted by *LErNet* and the labels.

```R
LErNet::visualize(
  lncrna_context,
  de_proteins,
  genomic_seeds,
  unlist(components),
  ppi_network,
  labels)
```

![This is the image returned by the function. In the left upper box it is possible to select only a group to be viewed in the plot (lncRNA, Seed Connector and Seed Protein).](https://i.imgur.com/lS6Gzxo.png)

The last step is the functional enrichment of the results. Basically *LErNet* exploits the package ReactomePA to retrieve significant pathways through the function `enrich`: 

```R
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
```

![These are the most significant pathways retrieved by LErNet for the example with mouse genes.](https://i.imgur.com/ivBV1S1.png)

However, the user can use the preferred tool to make functional enrichment.

