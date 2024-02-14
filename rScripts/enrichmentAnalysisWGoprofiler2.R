##  https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html#enrichment-analysis
##  this is the online document from CRAN

library(gprofiler2)


gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)


names(gostres)
head(gostres$result, 3)


gostres2 <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL, highlight = TRUE)

head(gostres2$result, 3)

gostres_link <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                     as_short_link = TRUE)

gostplot(gostres, capped = TRUE, interactive = TRUE)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p

pp <- publish_gostplot(p, highlight_terms = c("GO:0048013", "REAC:R-HSA-3928663"), 
                       width = NA, height = NA, filename = NULL )

pp


publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,120),],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)


gorth(query = c("Klf4", "Sox2", "71950"), source_organism = "mmusculus", 
      target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE,
      numeric_ns = "ENTREZGENE_ACC")


