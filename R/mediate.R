# Mediation

pos_Mbp <- 83.81371

## Mediate `r physioshort` across all genes.

gene_info <- readRDS("data/geno_info.rds") |>
  dplyr::filter(chromosome_name == chrname) |>
  dplyr::rename(chr = "chromosome_name",
                pos = "start_position") |>
  dplyr::mutate(
    ensembl_gene_id =  stringr::str_remove(ensembl_gene_id, "ENSMUSG0*")) |>
  tidyr::unite(id, external_gene_name, ensembl_gene_id)

# Limit to genes on chr

gene <- gene[, colnames(gene) %in% gene_info$id]
gene <- apply(gene, 2, qtl2mediate::nqrank)

# Sort `gene` to match Physio

m <- match(rownames(Physio), rownames(gene))
gene <- gene[m,, drop = FALSE]

peak_mar <- qtl2::find_marker(pmap, chrname, pos_Mbp)
geno_max <- subset(aprobs_chr, chr = chrname, mar = peak_mar)[[1]][,,1]

if(exists("medscan") && medscan) {
  med_scan <- intermediate::mediation_scan(
    target = physio_trait,
    mediator = gene,
    driver = geno_max,
    annotation = gene_info,
    covar = addcovarSD,
    kinship = kinship_chr,
    fitFunction =  qtl2mediate::fitQtl2)
  
  summary(med_scan, minimal = TRUE)
}

med_test <- intermediate::mediation_test(
  target = physio_trait,
  mediator = gene_trait,
  driver = geno_max,
  annotation = gene_info,
  covar_tar = addcovarSD,
  covar_med = addcovar,
  kinship = kinship_chr,
  fitFunction =  qtl2mediate::fitQtl2)
summary(med_test)
