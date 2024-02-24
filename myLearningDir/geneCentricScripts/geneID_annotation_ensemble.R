normalizedReadCountsAll <- read.table(
               "x:/project2023/CidlowskiRNAseq/DavidRNAseq/RNAseqAnalysis/analysisDir/results/DESeq2/mouse_STAR_normed.txt",
               header = TRUE, sep = "\t")

dim(normalizedReadCountsAll)
head(normalizedReadCountsAll)

#Comparison 1: dFlox M0 Vehicle vs dFlox M0 Dex

## wt(\d)VehM0 wt(\d)DexM0
colnames(normalizedReadCountsAll) [c( which(grepl("wt.VehM0", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM0", colnames(normalizedReadCountsAll))))]
# [1] "wt2VehM0" "wt3VehM0" "wt4VehM0" "wt2DexM0" "wt3DexM0" "wt4DexM0"   

normdata_comparison_1 <- normalizedReadCountsAll [, c(1,  which(grepl("wt.VehM0", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM0", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_1)


# write.table (normdata_comparison_1, file = "x:/project2024/DavidRNAseq/GSEAfiles/comp_1_data.txt", sep="\t", row.names = TRUE, col.names = NA)


# ...convert gene symbol to EntrezID -----------

library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggridges)

## An extremely user's friendly package
## https://biit.cs.ut.ee/gprofiler/page/r

library(gprofiler2)

mouseEntrezSymbol <- normdata_comparison_1$X

head(mouseEntrezSymbol)
geneDesc <- gconvert(query = normdata_comparison_1$X, organism = "mmusculus",
         target="ENSG", mthreshold = Inf, filter_na = TRUE)

head(geneDesc)

annot <- as.data.frame(cbind(geneDesc$input,geneDesc$description))

dim(annot)

colnames(annot) <- c("EntrezID", "Description")
head(annot)

gsub("[*", "", annot$Description[1:3])


gsub(" \\[Source", "", head(geneDesc$description))
gsub(" \\[Source", "", head(geneDesc$description))
gsub(" \\[Source:MGI Symbol;Acc:MGI:", "", head(geneDesc$description))
gsub(" \\[Source:MGI Symbol;Acc:MGI:\\d+\\]", "", head(geneDesc$description))

geneDescription = gsub(" \\[Source:MGI Symbol;Acc:MGI:\\d+\\]", "", geneDesc$description)


length(normdata_comparison_1$X)


merge(head(mouseEntrezSymbol), geneDesc, by='input' )

mm <- org.Mm.eg.db
my.symbols <- c(as.character(normdata_comparison_1$X))
df <- AnnotationDbi::select(mm, 
                            keys = my.symbols,
                            columns = c("ENTREZID", "GENENAME"),
                            keytype = "SYMBOL")

columns(mm)
head(df)
is.na(df$ENTREZID)
