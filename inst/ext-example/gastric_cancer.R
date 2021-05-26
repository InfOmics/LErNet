library(LErNet)


arenaidb <- LErNet::get_arenaidb_diseases()
View(arenaidb)
gc_interactions <- arenaidb[arenaidb$phenotype.name == 'GASTRIC CANCER', ]
View(gc_interactions)

unique(gc_interactions$ncrna.naming.resource)


ncrnas <- LErNet::get_arenaidb_ncrnas()
unique(ncrnas$biotype)
View(ncrnas)

lncrnas <- ncrnas[(nchar(ncrnas$sequence) >= 200) | (ncrnas$biotype == 'lncrna') , ]
View(lncrnas)
View( lncrnas[lncrnas$sequence == '',])


x <- lncrnas[!is.na(lncrnas$location),]
View(x)

gc_genes <- read.csv(system.file("extdata", "gastric_cancer", "gene_result.txt", package = "LErNet"), sep='\t', stringsAsFactors = FALSE)
View(gc_genes)
