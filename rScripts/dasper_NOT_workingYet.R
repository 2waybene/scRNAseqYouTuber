library("dasper")

ref <- GenomicState::GenomicStateHub(
  version = "31",
  genome = "hg38",
  filetype = "TxDb"
)[[1]]


junctions_processed <- junction_process(
  junctions_example,
  ref,
  types = c("ambig_gene", "unannotated")
)

sashimi_plot <- plot_sashimi(
  junctions = junction_filter(junctions_processed),
  ref = ref,
  gene_tx_id = "ENSG00000142156.14",
  gene_tx_col = "gene_id",
  sum_func = NULL
)