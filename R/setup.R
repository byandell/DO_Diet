# Dataset Setup

datadir <- "data"

## 8 allele probabilities

fst_dir <- file.path(datadir, "fst_genoprob")
aprobs_chr <- readRDS(file.path(fst_dir, "aprobs_fstindex.rds"))
aprobs_chr <- qtl2fst::subset_fst_genoprob(aprobs_chr, chr = chrname)

## Physical Map and Kinship

# data/gigamuga_map_test_v2.RDS
pmap <- readRDS(file.path(datadir, "pmap.rds"))

# data/K_1176_DietDO.rds
kinship_chr <- readRDS(file.path(datadir, "kinship.rds"))[chrname]

## Covariates

# data/Phenotype files/clinical data sets/Final in vivo DO data for analysis
# (raw and computed).xlsx
Physio <- readxl::read_excel(file.path(datadir, "Physio.xlsx"),
                             sheet = "Data") |>
  tibble::column_to_rownames("mouse")
addcovarSD <- stats::model.matrix(
  ~ Sex * Diet.name,
  Physio |> dplyr::select(Sex, Diet.name))[,-1]
addcovar <- stats::model.matrix(
  ~ Gen_Lit + Sex * Diet.name,
  Physio |> dplyr::select(Gen_Lit, Sex, Diet.name))[,-1]
intcovar <- addcovar[,"Diet.nameHF", drop = FALSE]

physio_trait <- Physio[,
  grep(paste(physioname, collapse = "|"), colnames(Physio)),
  drop = FALSE]

## RNA Data

# Read gene data and combine gene name and ENSEMBL ID.
gene <- readRDS("data/gene.rds")

gene_trait <- gene[,
  grep(paste(paste0(genename, "_"), collapse = "|"), colnames(gene)),
  drop = FALSE]

## Isoform data

isoform_trait <- readRDS("data/isoform.rds")
isoformname <- grep(paste(paste0(genename, "_"), collapse = "|"),
                    colnames(isoform_trait))
if(isoformid != "")
  isoformname <- isoformname[
    grep(isoformid, colnames(isoform_trait)[isoformname])]

isoform_trait <- isoform_trait[, isoformname, drop = FALSE]
isoformid <- stringr::str_remove(colnames(isoform_trait), "^.*_")

## Normal Score

if(normscore) {
  gene_trait <- qtl2mediate::nqrank(gene_trait)
  if(isoformid != "")
    isoform_trait <- qtl2mediate::nqrank(isoform_trait)
}

## Sort `gene_trait` to match Physio

m <- match(rownames(Physio), rownames(gene_trait))
gene_trait <- gene_trait[m,, drop = FALSE]
m <- match(rownames(Physio), rownames(isoform_trait))
isoform_trait <- isoform_trait[m,, drop = FALSE]
