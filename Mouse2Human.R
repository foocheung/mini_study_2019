## check the genes of Human and Mouse Genes are just capitals for the pathways of interest.
## https://www.biostars.org/p/147484/
## host="https://dec2021.archive.ensembl.org/

library(biomaRt)

Mouse2Human <- function(MouseGenes){

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  genesMousetoHuman = getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), 
                             filters = "mgi_symbol", 
                             values = MouseGenes , 
                             mart = mouse, 
                             attributesL = c("ensembl_gene_id", "hgnc_symbol"), 
                             martL = human, 
                             uniqueRows = TRUE)

  colnames(genesMousetoHuman) <- c("Mouse.Gene_ID", "MGI", "Human.Gene_ID", "HGNC")

  return(genesMousetoHuman) 

}

## Get mouse genes
mmusculus_genes <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),  
                         mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                         useCache = FALSE)

## create the conversion table
Mouse2HumanTable <- Mouse2Human(MouseGenes = mmusculus_genes$mgi_symbol)
