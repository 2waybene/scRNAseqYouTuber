##  credit: https://jmw86069.github.io/splicejam/articles/create-a-sashimi-plot.html#overview-1

options("warn"=-1);
suppressPackageStartupMessages(library(splicejam));
suppressPackageStartupMessages(library(jamba));
suppressPackageStartupMessages(library(kableExtra));


if (suppressPackageStartupMessages(require(TxDb.Mmusculus.UCSC.mm10.knownGene))) {
  
  # First obtain exons by transcript
  exonsByTxMm10 <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,
                           by="tx",
                           use.names=TRUE);
  values(exonsByTxMm10@unlistData)$feature_type <- "exon";
  values(exonsByTxMm10@unlistData)$subclass <- "exon";
  
  # For added insight, obtain CDS exons by transcript (optional)
  cdsByTxMm10 <- cdsBy(TxDb.Mmusculus.UCSC.mm10.knownGene,
                       by="tx",
                       use.names=TRUE);
  values(cdsByTxMm10@unlistData)$feature_type <- "cds";
  values(cdsByTxMm10@unlistData)$subclass <- "cds";
}



if (require(TxDb.Mmusculus.UCSC.mm10.knownGene)) {
  
  # Now prepare tx_name, gene_id, gene_name data.frame,
  # surprisingly difficult
  tx2geneMm10 <- suppressMessages(
    AnnotationDbi::select(TxDb.Mmusculus.UCSC.mm10.knownGene,
                          keys(TxDb.Mmusculus.UCSC.mm10.knownGene, "GENEID"),
                          columns=c("GENEID","TXNAME"),
                          keytype="GENEID"
    )
  );
  tx2geneMm10 <- renameColumn(tx2geneMm10,
                              from=c("GENEID", "TXNAME", "TXTYPE"),
                              to=c("gene_id", "transcript_id", "transcript_type"));
  
  # add gene_name using org.Mm.eg.db
  if (suppressPackageStartupMessages(require(org.Mm.eg.db))) {
    gene_ids <- values(genes(TxDb.Mmusculus.UCSC.mm10.knownGene))$gene_id;
    gene_namesL <- mget(gene_ids,
                        org.Mm.egSYMBOL,
                        ifnotfound=NA);
    ## Convert list to vector taking the first gene_name each
    ## (All genes should only have one SYMBOL but there is no
    ## hard constraint so we should make absolutely sure to
    ## use only one value per gene.)
    gene_names <- unlist(heads(S4Vectors::List(gene_namesL), 1));
    ## Replace NA with LOC# format
    ## Note we use gsub() to ensure the data fits the expected format
    if (any(is.na(gene_names))) {
      gene_na <- which(is.na(gene_names));
      gene_names[gene_na] <- gsub("^([0-9]+)$", "LOC\\1",
                                  names(gene_names[gene_na]));
    }
    tx2geneMm10$gene_name <- gene_names[as.character(tx2geneMm10$gene_id)];
  } else {
    ## If we have no gene annotations, use the gene_id values
    tx2geneMm10$gene_name <- as.character(tx2geneMm10$gene_id);
  }
  # print the first 20 rows to show the content
  print(head(tx2geneMm10, 20));
}
#> Loading required package: TxDb.Mmusculus.UCSC.mm10.knownGene
#> 


if (require(TxDb.Mmusculus.UCSC.mm10.knownGene)) {
  # flatten exons to the gene level
  # for speed, we will only process "Gria1", and "Ntrk2"
  flatExonsByGeneMm10 <- flattenExonsBy(exonsByTx=exonsByTxMm10,
                                        cdsByTx=cdsByTxMm10,
                                        by="gene",
                                        genes=c("Gria1", "Ntrk2"),
                                        tx2geneDF=tx2geneMm10,
                                        verbose=FALSE);
  
  # to be fancy, also flatten transcripts, to include CDS ranges   
  flatExonsByTxMm10 <- flattenExonsBy(exonsByTx=exonsByTxMm10,
                                      cdsByTx=cdsByTxMm10,
                                      tx2geneDF=tx2geneMm10,
                                      by="tx",
                                      genes=c("Gria1", "Ntrk2"));
}
#> Loading required package: TxDb.Mmusculus.UCSC.mm10.knownGene


if (require(TxDb.Mmusculus.UCSC.mm10.knownGene)) {
  # Pull out Gria1
  grlGria1 <- flatExonsByGeneMm10[["Gria1"]];
  
 # grlGria1 <- flatExonsByGeneMm10[["Atg9b"]];
  
  
  ## Plot a basic gene-exon structure
  ggGria1exons <- gene2gg(gene="Gria1",
                          flatExonsByGene=flatExonsByGeneMm10,
                          exonLabelSize=6);
  print(ggGria1exons + ggtitle("Gria1 exons"));
  
  ## Compare to the gene structure without compressing introns
  gg1full <- gene2gg(gene="Gria1",
                     flatExonsByGene=flatExonsByGeneMm10,
                     compressGaps=FALSE)
  print(gg1full);
  
  ## Plot a slightly more detailed gene-transcript-exon structure
  ggGria1exonsTx <- gene2gg(gene="Gria1",
                            flatExonsByGene=flatExonsByGeneMm10,
                            flatExonsByTx=flatExonsByTxMm10,
                            tx2geneDF=tx2geneMm10);
  print(ggGria1exonsTx + ggtitle("Gria1 (compressed introns)"));
  
  ## Notice how difficult it is to see exon15 and exon16 are
  ## mutually exclusive exons
  gg2full <- gene2gg(gene="Gria1",
                     flatExonsByGene=flatExonsByGeneMm10,
                     flatExonsByTx=flatExonsByTxMm10,
                     tx2geneDF=tx2geneMm10,
                     compressGaps=FALSE);
  print(gg2full + ggtitle("Gria1 (uncompressed introns)"));
  
}
#> Loading required package: TxDb.Mmusculus.UCSC.mm10.knownGene

## assemble a data.frame
baseurl <- "https://orio.niehs.nih.gov/ucscview/farrisHub/mm10/";

# BED files with junction reads
bedext <- ".STAR_mm10.combinedJunctions.bed";
bwext <- c("492_1.sickle.merged.cutadapt.STAR_mm10.pos.bw",
           "492_1.sickle.merged.cutadapt.STAR_mm10.neg.bw");
c1 <- c("CA1", "CA2");
r1 <- c("CB", "DE");
bedsamples <- paste0(rep(c1, each=2), "_", r1);
bedurls <- paste0(baseurl,
                  bedsamples,
                  bedext);

# bigWig files with strand-specific read coverage
bwsamples <- paste0(rep(c1, each=4),
                    rep(r1, each=2));
bwsamples1 <- paste0(rep(c1, each=4),
                     "_",
                     rep(r1, each=2));
bwurls <- paste0(baseurl, "NS50211",
                 bwsamples,
                 bwext);

# Assemble into a data.frame
filesDF <- data.frame(stringsAsFactors=FALSE,
                      check.names=FALSE,
                      url=c(bedurls, bwurls),
                      type=rep(c("junction", "bw"), c(4,8)),
                      sample_id=c(bedsamples, bwsamples1),
                      scale_factor=rep(c(1,3), c(length(bedurls), length(bwurls))));
filesDF;

data(test_exon_wide_gr);
test_exon_wide_gr;
#> GRanges object with 4 ranges and 1 metadata column:
#>         seqnames      ranges strand |   gene_name
#>            <Rle>   <IRanges>  <Rle> | <character>
#>   wide1     chr1     100-200      + |   TestGene1
#>   wide2     chr1 10300-10400      + |   TestGene1
#>   wide3     chr1 20500-20750      + |   TestGene1
#>   wide4     chr1 39900-40000      + |   TestGene1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
data(test_cov_wide_gr);
test_cov_wide_gr;
#> GRanges object with 4 ranges and 1 metadata column:
#>         seqnames      ranges strand |        sample_A
#>            <Rle>   <IRanges>  <Rle> |   <NumericList>
#>   wide1     chr1     100-200      + | 246,248,261,...
#>   wide2     chr1 10300-10400      + | 195,202,198,...
#>   wide3     chr1 20500-20750      + | 195,189,178,...
#>   wide4     chr1 39900-40000      + | 260,257,253,...
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#>   

# To plot a simple GRanges object
widecovdf <- exoncov2polygon(test_cov_wide_gr, covNames="sample_A");
suppressPackageStartupMessages(require(ggplot2));
ggWide3 <- ggplot(widecovdf,
                  aes(x=x, y=y, group=gr, fill=gr, color=gr)) +
  ggforce::geom_shape(alpha=0.7) +
  colorjam::theme_jam() +
  colorjam::scale_fill_jam() +
  colorjam::scale_color_jam();
print(ggWide3);

# Now compress the introns keeping axis labels
ref2c <- make_ref2compressed(test_cov_wide_gr,
                             nBreaks=10);

# Repeat ggplot using ref2c$trans_grc
ggWide3c <- ggWide3 +
  scale_x_continuous(trans=ref2c$trans_grc) +
  xlab("chr1 (compressed introns)") +
  ggtitle("exons (compressed introns)");
print(ggWide3c);


data(test_junc_wide_gr);
test_junc_wide_gr;

# To plot junctions, use grl2df(..., shape="junction")
junc_wide_df <-grl2df(test_junc_wide_gr,
                      shape="junction");

ggWide1 <- ggplot(junc_wide_df,
                  aes(x=x, y=y, group=gr_name, fill=gr_name, color=gr_name)) +
  splicejam::geom_diagonal_wide_arc() +
  colorjam::theme_jam() +
  colorjam::scale_fill_jam(alpha=0.7) +
  colorjam::scale_color_jam() +
  xlab("chr1") +
  ggtitle("junctions (full intron width)")
print(ggWide1);


## problem in preparing sashimi file for shGria1

if (exists("flatExonsByGeneMm10")) {
  shGria1 <- prepareSashimi(gene="Gria1",
                            flatExonsByGene=flatExonsByGeneMm10,
                            minJunctionScore=100,
                            sample_id=c("CA1_CB", "CA2_CB"),
                            filesDF=filesDF);
}

plotSashimi(shGrial)


data(test_exon_wide_gr);
data(test_junc_wide_gr);
data(test_cov_wide_gr);
testfilesDF <- data.frame(url="sample_A",
                          type="coverage_gr",
                          sample_id="sample_A",
                          scale_factor=1);
testfilesDF;
#>        url        type sample_id scale_factor
#> 1 sample_A coverage_gr  sample_A            1

# Now prepare sashimi data
sh1 <- prepareSashimi(GRangesList(TestGene1=test_exon_gr),
                      filesDF=testfilesDF,
                      gene="TestGene1",
                      sample_id="sample_A",
                      covGR=test_cov_gr,
                      juncGR=test_junc_gr);


plotSashimi(sh1)

plotSashimi(sh1) + 
  ggforce::facet_zoom(xlim=sh1$ref2c$transform(c(370, 410)));

