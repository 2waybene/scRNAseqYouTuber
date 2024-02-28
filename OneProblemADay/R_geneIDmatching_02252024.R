##  need to give him the credit: https://www.gungorbudak.com/blog/2018/08/07/convert-gene-symbols-to-entrez-ids-in-r/
library('org.Hs.eg.db')
columns(org.Hs.eg.db)
symbols <- c('AHNAK', 'BOD1L1', 'HSPB1', 'SMARCA4', 'TRIM28')

# use mapIds method to obtain Entrez IDs
mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')

availableCols <- columns(org.Hs.eg.db)
for (i in 1:length(availableCols)){
  results <-mapIds(org.Hs.eg.db, symbols, availableCols[i], 'SYMBOL')
  print (results)
}


i = 10
availableCols[i]
results <-mapIds(org.Hs.eg.db, symbols, availableCols[i], 'SYMBOL')
print (results)

# input list of Ensembl ID's
ensembl.ids <- mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL')
ensembl.ids <- read.delim('gene_ids.txt', header = F)

##==========================================================
##  credit from Ms. Patel
##  https://www.youtube.com/watch?v=cWe359VnfaY&t=501s
# script to convert ensembl gene Id to gene symbols
# setwd("~/Desktop/demo/convert_geneID_to_gene_symbols")

library(biomaRt)
library(annotables)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

# input list of Ensembl ID's
ensembl.ids <- read.delim('gene_ids.txt', header = F)


# method 1: biomaRt
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)

ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

getBM(attributes = c('ensembl_gene_id','external_gene_name'),
      filters = "ensembl_gene_id",
      values = ensembl.ids,
      mart = ensembl.con)

# method 2: annotables
grch38 %>%
  dplyr::filter(ensgene %in% ensembl.ids)


# method 3: annotation DBs

keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

mapIds(org.Hs.eg.db,
       keys = ensembl.ids,
       keytype = 'ENSEMBL',
       column = 'SYMBOL')

keytypes(EnsDb.Hsapiens.v86)
columns(EnsDb.Hsapiens.v86)

mapIds(EnsDb.Hsapiens.v86,
       keys = ensembl.ids,
       keytype = 'GENEID',
       column = 'SYMBOL')






