#########################################################################################
## The R script to retrieve SNP rs IDs for genes exported from GeneCards
## Get the genes list first: http://www.genecards.org/Search/Keyword?queryString=PUFA
#########################################################################################

# install.packages("XML")
# install.packages("scholar", dependencies=T)
# install.packages("stringi", repos="http://cran.rstudio.com/", dependencies=TRUE)
library(stringi)
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
library(biomaRt)
# listMarts()
library(httr)

rm(list = ls())
version <- "grch37."
species <- "homo_sapiens"

# load the human variation data
variation = useEnsembl(biomart="snp", dataset="hsapiens_snp")

wd <- "D:/Dropbox/GeneCard/"
setwd(normalizePath(wd))
d <- read.csv("GeneCards-SearchResults.csv", sep=",", header = T, fill = T)
number_of_records <- nrow(d)

for (gene_to_read in 1:number_of_records) {
  gene <- d$Gene.Symbol[gene_to_read]
  gene <- as.character(gene)
  # show some progress
  cat("Doing:", gene, "\tTotal progress:", 100*gene_to_read%/%number_of_records, "%\n")
  # gene symbol to Ensembl ID
  response <- GET(sprintf("http://%srest.ensembl.org/xrefs/symbol/%s/%s?content-type=application/json",
                          version, species, gene))
  parsed <- content(response, "parsed", encoding = "UTF8")
  # ignore missing genes
  EnsemblID <- tryCatch({parsed[[1]]$id}, error = function(cond)"")
  # look up a single gene and get SNP data
  if (EnsemblID != "") {
    a <- getBM(attributes = c(
      "ensembl_gene_stable_id",
      "refsnp_id",
      "chr_name",
      "chrom_start",
      "chrom_end"),
      filters = "ensembl_gene",
      values = EnsemblID,
      mart = variation)
    # ignore missing SNP rs IDs
    if (nrow(a) != 0) {
      # replace Ensembl ID by gene name
      a$ensembl_gene_stable_id <- gene
      # fix the column name
      colnames(a)[colnames(a)=="ensembl_gene_stable_id"] <- "gene_name"
      # write the data
      write.table(a, file = "all_SNPs_from_GeneCards.csv", sep = ",", append = T, row.names = FALSE)
    }
  }
}
