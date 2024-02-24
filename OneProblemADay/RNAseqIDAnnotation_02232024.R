##  preparing GSEA files for the analysis
##   with gprofiler2 package

prepareGSEAdata <- function (dmIN, dmOUT)
{
  ## dmIN$X has the entrez gene ID
  ## dmOUT is the file to write to
  ## Use gprofiler2 to get gene annotation
  library(gprofiler2)

  mouseEntrezSymbol <- dmIN$X
  colnames(dmIN)[1] <- "EntrezID"  
  geneDesc <- gconvert(query = normdata_comparison_1$X, organism = "mmusculus",
                     target="ENSG", mthreshold = Inf, filter_na = TRUE)


  annot <- as.data.frame(cbind(geneDesc$input,geneDesc$description))
  geneDescription = gsub(" \\[Source:MGI Symbol;Acc:MGI:\\d+\\]", "", geneDesc$description)
  annot$V2 =  geneDescription
  colnames(annot) <- c("EntrezID", "Description")

  dataWithAnnot <- merge(annot, dmIN, by = c("EntrezID"))

  dataOut <- dataWithAnnot [-which(rowSums (dataWithAnnot[,-c(1,2)]) ==0),] 
  
  write.table (dataOut, dmOUT, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
##==============================================================================


normalizedReadCountsAll <- read.table(
  "x:/project2023/CidlowskiRNAseq/DavidRNAseq/RNAseqAnalysis/analysisDir/results/DESeq2/mouse_STAR_normed.txt",
  header = TRUE, sep = "\t")

dim(normalizedReadCountsAll)
colnames(normalizedReadCountsAll)

##=====================================================
##  Comparison 1: dFlox M0 Vehicle vs dFlox M0 Dex
##
## wt(\d)VehM0 wt(\d)DexM0
##=====================================================
colnames(normalizedReadCountsAll) [c( which(grepl("wt.VehM0", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM0", colnames(normalizedReadCountsAll))))]
# [1] "wt2VehM0" "wt3VehM0" "wt4VehM0" "wt2DexM0" "wt3DexM0" "wt4DexM0"   

normdata_comparison_1 <- normalizedReadCountsAll [, c(1,  which(grepl("wt.VehM0", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM0", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_1)


prepareGSEAdata (normdata_comparison_1, "x:/project2024/DavidRNAseq/GSEAfiles/comp_1_data.gct")


##===================================================
#Comparison 10: dFlox M0 Dex vs dKO M0 Dex
##  ## wt(\d)DexM0 ko(\d)DexM0
##====================================================
colnames(normalizedReadCountsAll)

colnames(normalizedReadCountsAll) [c( which(grepl("ko.DexM0", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM0", colnames(normalizedReadCountsAll))))]
# [1] "ko1DexM0" "ko2DexM0" "ko3DexM0" "ko4DexM0" "ko5DexM0" "wt2DexM0" "wt3DexM0" "wt4DexM0" 

normdata_comparison_10 <- normalizedReadCountsAll [, c(1, which(grepl("ko.DexM0", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM0", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_10)

prepareGSEAdata (normdata_comparison_10, "x:/project2024/DavidRNAseq/GSEAfiles/comp_10_data.gct")




##===================================================
#Comparison 2: dFlox M1 Vehicle vs dFlox M1 Dex 
##
## wt(\d)VehM1 wt(\d)DexM1
##====================================================
colnames(normalizedReadCountsAll) [c( which(grepl("wt.VehM1", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM1", colnames(normalizedReadCountsAll))))]
# [1] "wt2VehM1" "wt3VehM1" "wt4VehM1" "wt2DexM1" "wt3DexM1" "wt4DexM1"   

normdata_comparison_2 <- normalizedReadCountsAll [, c(1,  which(grepl("wt.VehM1", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM1", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_2)
dim(normdata_comparison_2)

prepareGSEAdata (normdata_comparison_2, "x:/project2024/DavidRNAseq/GSEAfiles/comp_2_data.gct")


##=================================

#Comparison 11: dFlox M1 Dex vs dKO M1 Dex 
##  ## wt(\d)DexM1 ko(\d)DexM1
##====================================================
colnames(normalizedReadCountsAll)

colnames(normalizedReadCountsAll) [c( which(grepl("ko.DexM1", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM1", colnames(normalizedReadCountsAll))))]
# [1] "ko1DexM1" "ko2DexM1" "ko3DexM1" "ko4DexM1" "ko5DexM1" "wt2DexM1" "wt3DexM1" "wt4DexM1"

normdata_comparison_11 <- normalizedReadCountsAll [, c(1, which(grepl("ko.DexM1", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM1", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_11)

prepareGSEAdata (normdata_comparison_11, "x:/project2024/DavidRNAseq/GSEAfiles/comp_11_data.gct")



##==================================================
#Comparison 3: dFlox M2 Vehicle vs dFlox M2 Dex 
##
## wt(\d)VehM2 wt(\d)DexM2
##====================================================
colnames(normalizedReadCountsAll) [c( which(grepl("wt.VehM2", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM2", colnames(normalizedReadCountsAll))))]
# [1] "wt2VehM2" "wt3VehM2" "wt4VehM2" "wt2DexM2" "wt3DexM2" "wt4DexM2" 

normdata_comparison_3 <- normalizedReadCountsAll [, c(1,  which(grepl("wt.VehM2", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM2", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_3)
dim(normdata_comparison_3)

prepareGSEAdata (normdata_comparison_3, "x:/project2024/DavidRNAseq/GSEAfiles/comp_3_data.gct")


##===================================================
#Comparison 12: dFlox M2 Dex vs dKO M2 Dex
##  ## wt(\d)DexM2 ko(\d)DexM2
##====================================================
colnames(normalizedReadCountsAll)

colnames(normalizedReadCountsAll) [c( which(grepl("ko.DexM2", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM2", colnames(normalizedReadCountsAll))))]
# [1] "ko1DexM2" "ko2DexM2" "ko3DexM2" "ko4DexM2" "ko5DexM2" "wt2DexM2" "wt3DexM2" "wt4DexM2"

normdata_comparison_12 <- normalizedReadCountsAll [, c(1, which(grepl("ko.DexM2", colnames(normalizedReadCountsAll))), which(grepl("wt.DexM2", colnames(normalizedReadCountsAll))))]
head(normdata_comparison_12)

prepareGSEAdata (normdata_comparison_12, "x:/project2024/DavidRNAseq/GSEAfiles/comp_12_data.gct")



