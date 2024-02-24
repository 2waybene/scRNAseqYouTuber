##=================================================================
##  Problem: I have an RNAseq normalized datamatrix
##  the first column is entrez gene ID, and I want to 
##  get the "gene description" to further prepere an expression 
##  data set for GSEA analysis
##================================================================
##  Answer: to get this done, I need 
##        a database that I can map the entrez ID with description
##        only retain those genes with description
##  
##        gprofiler2 package provide a function "gconvert" that 
##        will do the trick
##=================================================================
prepareGSEAdata <- function (dmIN)
{
  ## dmIN$X has the entrez gene ID
  ## dmOUT is the file to write to
  ## Use gprofiler2 to get gene annotation
  library(gprofiler2)

  ## get entrez symbol, in my case the first column
  mouseEntrezSymbol <- dmIN$X
  colnames(dmIN)[1] <- "EntrezSymbol"  
  
  ##  Main function "gconvert"
  geneDesc <- gconvert(query =  mouseEntrezSymbol, organism = "mmusculus",
                     target="ENSG", mthreshold = Inf, filter_na = TRUE)

  # I ONLY need to columns
  annot <- as.data.frame(cbind(geneDesc$input,geneDesc$description))
  
  # clean up the description
  geneDescription = gsub(" \\[Source:MGI Symbol;Acc:MGI:\\d+\\]", "", geneDesc$description)
  annot$V2 =  geneDescription
  
  colnames(annot) <- c("EntrezSymbol", "Description")

  dataWithAnnot <- merge(annot, dmIN, by = c("EntrezSymbol"))

  ## remove all zero rows
  dataOut <- dataWithAnnot [-which(rowSums (dataWithAnnot[,-c(1,2)]) ==0),] 
  
  return(dataOut)
}


##==============================================================================




##=====================================================
##  Comparison 1: dFlox M0 Vehicle vs dFlox M0 Dex
##
## wt(\d)VehM0 wt(\d)DexM0
##=====================================================

dt <- readRDS("/Users/jyli/myGit/scRNAseqYT/OneProblemADay/data.rds")
dim(dt)
dt.w.annot <- prepareGSEAdata (dt)
dim(dt.w.annot)

