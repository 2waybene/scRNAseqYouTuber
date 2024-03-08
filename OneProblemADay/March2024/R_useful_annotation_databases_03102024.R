library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
m_ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)
h_ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

today()
#[1] "2024-03-08"
length(columns(m_ensembl))
# [1] 2988
length(columns(h_ensembl))
# [1] 3172