
##=====================================================
##  Now, let's work on Mouse related ID questions
##=====================================================


##======================================================
##  Goal 1:  ...convert gene symbol to EntrezID -----------
##======================================================

## An extremely user's friendly package
## https://biit.cs.ut.ee/gprofiler/page/r
library(gprofiler2)
setwd("x:/myGit/scRNAseqYouTuber/myLearningDir/geneCentricScripts/")
RNAseq_datamatrix <- read.csv("sample_mouse_rnaseq_datamatrix.csv", header = TRUE)

head(RNAseq_datamatrix)

mouseEntrezSymbol <- RNAseq_datamatrix$X
head(mouseEntrezSymbol)

##==============================
##  use "gconvert" function
##==============================
geneDesc <- gconvert(query = normdata_comparison_1$X, organism = "mmusculus",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)

head(geneDesc)

##  retain two columns: entrezSymbol and Description 
annot <- as.data.frame(cbind(geneDesc$input,geneDesc$description))
dim(annot)
head(annot)

colnames(annot) <- c("EntrezID", "Description")
head(annot)

##  clean up the gene descriptions
geneDescription = gsub(" \\[Source:MGI Symbol;Acc:MGI:\\d+\\]", "", geneDesc$description)
length(geneDescription)
dim(annot)
annot$Description <- geneDescription
head(annot)

##   Finally update the original file with gene decription
df1 <- merge(x = as.data.frame(mouseEntrezSymbol), 
                        y = geneDesc, 
                        by.x="mouseEntrezSymbol",
                        by.y="input")
head(df1)

df <- merge(x = RNAseq_datamatrix, 
            y = geneDesc, 
            by.x="X",
            by.y="input")
head(df)
dim(df)

##====================================
##  use org.Mm.eg.db database query
##====================================

library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggridges)



mm <- org.Mm.eg.db
my.symbols <- c(as.character(RNAseq_datamatrix$X))
df <- AnnotationDbi::select(mm, 
                            keys = my.symbols,
                            columns = c("ENTREZID", "GENENAME"),
                            keytype = "SYMBOL")

columns(mm)
head(df)
dim(df)
length(which(!is.na(df$ENTREZID)))
df <- df [which(!is.na(df$ENTREZID)),]

df1 <- merge(x = RNAseq_datamatrix, 
             y = df, 
             by.x="X",
             by.y="SYMBOL")
head(df1)
dim(df1)

##======================================================
##  Goal 2:  ...convert ensemble ID to gene symbol
##======================================================

library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(biomaRt)


dm <- read.csv ("mm39STAR_mouse_gene_count.csv" ,header = TRUE)

head(dm)

##========================================
##  remove ensemble id verion: .number
##==========================================
ensembl.ids <- dm[,1]
ensembl.ids <- gsub('\\..+$', '', ensembl.ids)
dm.mod <- cbind (ensembl.ids, dm)
head(dm.mod)

##  biomaRt function
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
m_ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)


symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","mgi_symbol"),
                values = ensembl.ids, 
                mart = m_ensembl)
dim(symbol)
# [1] 3000    2
head(symbol)

df <- merge(x = symbol, 
            y = dm.mod, 
            by.x="ensembl_gene_id",
            by.y="ensembl.ids")
head(df)





