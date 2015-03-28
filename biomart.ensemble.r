# http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf

library(biomart)

f = '/Users/benjaminkellman/Documents/Grad/Labs/Couchesne/ASDanyl/Diagnosis_Biomarker/Classifier/testset.125samples/MFseperated/lasso.gaussian.male1.coef.s80.csv'

#ensembl=useMart("ensembl")
#listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

#attributes = listAttributes(ensembl)

results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                 filters = "ensembl_transcript_id", values = "ENST00000296026",
                 mart = ensembl)


##

#kegg pathway -> genes <- transcripts
library(KEGGREST)
# get all genes in each pathway
p = keggLink("pathway","hsa")
map = split(names(p), unname(p))

hsa2ens = list()
ensg2enst = list()
#iterate over pathways
for(m in 1:length(map)){
  #iterate over genes in pathways
  for(h in 1:length(map[[m]])){
	gene.hsa = map[[m]][h]
	if(!(gene.hsa %in% names(hsa2ens))){
	  query <- keggGet(gene.hsa)
	  hsa2ens[gene.hsa] <- strsplit(query[[1]]$DBLINKS[6],' ')[[1]][2]
	}
	ensg <- hsa2ens[gene.hsa]
	if(!(ensg %in% names(ensg2enst))){
	  
	}
  }
}
