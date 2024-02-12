## DESeq2 workflow for bulk RNA-seq data analysis

library(DESeq2)
library(tidyverse)

### Example: Load count matrix input from bulk RNA-seq data (*.csv)

## 1. mouse sciatic nerve injury bulk RNA-seq data
counts <- read.delim("../DESeq2/Sciatic_nerve.csv", header = T, sep = ",")
counts <- column_to_rownames(counts, var = "ENSEMBL_ID")

## 2. Load coldata for bulk RNA-seq
coldata <- read.delim("../DESeq2/Coldata_nerve.txt", header = T, sep = "\t")
coldata <- column_to_rownames(coldata, var = "Samples")

all(rownames(coldata) == colnames (counts)) # same name and same order

## 3. create a DESeq Data Set (dds)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Condition) # design formula

dds

## 4. Pre-filtering (optional)
# specify the smallest group size, the minimal number of samples
smallestGroupSize <- 4 # 4 Intact, 4 Injury

# keep only rows that have a count of at least 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep, ]

## 5. Set reference level for differential analysis
# factor
dds$Condition <- factor(dds$Condition, levels = c("Intact", "Injury"))

# or relevel
dds$Condition <- relevel(dds$Condition, ref = "Intact")

## 6. DESeq2 differential expression analysis
dds <- DESeq(dds) # three functions wrapped into a single function

## 7. Generate the results data frame and explore the results
resultsNames (dds) # lists the coefficients
res <- results(dds) # generate results table using the results function 
res

summary(res)
sum(res$padj < 0.1, na.rm = TRUE) # How many adjusted p-values were < 0.1?

res05 <- results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm = TRUE)

# plot the results
plotMA(res, ylim = c(-3, 3)) # blue padj <0.1
plotCounts(dds, gene = "ENSMUSG00000020218", intgroup = "Condition") # Wif1

library(org.Mm.eg.db) # Genome wide annotation database for Mouse
Symbol <- mapIds(org.Mm.eg.db, keys = rownames(res), keytype= "ENSEMBL",
                 column = "SYMBOL")
res$Symbol <- Symbol
res

## 8. Export and save the results
write.csv(as.data.frame(res), "../DESeq2/Sciatic_nerve_res.csv")

# subset by padj
resSig <- subset(res, padj < 0.05)
write.csv(as.data.frame(resSig), "../DESeq2/Sciatic_nerve__resSig.csv")