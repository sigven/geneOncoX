
source("data-raw/utils_driver_annotation.R")
source("data-raw/utils_predisposition_annotation.R")
source("data-raw/utils_gencode_annotation.R")
source("data-raw/utils_other.R")

## get metadata from metadata.xlsx
metadata <- list()
for (elem in c("basic", "predisposition", "panels", "alias", "gencode")) {
  metadata[[elem]] <- as.data.frame(openxlsx::read.xlsx(
    "data-raw/metadata.xlsx",
    sheet = elem, colNames = TRUE
  ) |>
    dplyr::mutate(version = dplyr::if_else(
      is.na(version) &
        (abbreviation == "ncbi" |
          abbreviation == "other" |
          abbreviation == "appris"),
      as.character(Sys.Date()),
      as.character(version)
    )))
}

## set logging layout
lgr::lgr$appenders$console$set_layout(
  lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T")
)

gene_info <- get_gene_info_ncbi() |>
  dplyr::select(
    entrezgene, symbol, gene_biotype,
    synonyms, name, other_genename_designations,
    ensembl_gene_id, hgnc_id
  ) |>
  dplyr::distinct()

gene_panels <- list()
gene_panels$metadata <- metadata[["panels"]]
gene_panels$records <- dplyr::bind_rows(
  get_panel_app_genes(gene_info = gene_info, build = "grch37"),
  get_panel_app_genes(gene_info = gene_info, build = "grch38")
)

gene_alias <- list()
gene_alias$metadata <- metadata[["alias"]]
gene_alias$records <- get_gene_aliases_ncbi(
  gene_info = gene_info,
  path_data_raw = file.path(here::here(), "data-raw")
)

gene_gencode <- list()
gene_gencode$records <- list()
gene_gencode$metadata <- metadata[["gencode"]]

gene_gencode$records[["grch38"]] <- gencode_get_transcripts(
  build = "grch38",
  gene_info = gene_info,
  gencode_version = as.integer(metadata$gencode[1, ]$version),
  ensembl_version = as.integer(metadata$gencode[2, ]$version),
  uniprot_version = "2022_04",
  gene_alias = gene_alias
)
gene_gencode$records[["grch37"]] <- gencode_get_transcripts(
  build = "grch37",
  gene_info = gene_info,
  gencode_version = as.integer(19),
  ensembl_version = as.integer(metadata$gencode[2, ]$version),
  uniprot_version = "2022_04",
  gene_alias = gene_alias
)

## "Rescue" some UniProt identifiers from
## grch38 - missing/not found for grch37

up_xref_grch38 <- gene_gencode[["records"]][["grch38"]] |>
  dplyr::select(
    entrezgene, uniprot_acc,
    uniprot_id
  ) |>
  dplyr::filter(!is.na(entrezgene) &
    !is.na(uniprot_acc) &
    !is.na(uniprot_id)) |>
  # dplyr::select(-uniprot_acc) |>
  dplyr::distinct()

up_xref <- list()
up_xref[["found"]] <-
  gene_gencode[["records"]][["grch37"]] |>
  dplyr::select(
    entrezgene, uniprot_acc,
    ensembl_transcript_id, uniprot_id
  ) |>
  dplyr::filter(!is.na(entrezgene) &
    !is.na(uniprot_acc) &
    !is.na(uniprot_id)) |>
  dplyr::distinct()

up_xref[["rescued_from_grch38"]] <-
  gene_gencode[["records"]][["grch37"]] |>
  dplyr::select(
    entrezgene, uniprot_acc,
    ensembl_transcript_id, uniprot_id
  ) |>
  dplyr::filter(!is.na(entrezgene) &
    !is.na(uniprot_acc) &
    is.na(uniprot_id)) |>
  dplyr::select(-c(uniprot_id)) |>
  dplyr::left_join(
    up_xref_grch38,
    by = c("entrezgene", "uniprot_acc")
  ) |>
  dplyr::filter(!is.na(uniprot_id)) |>
  dplyr::distinct()

up_xref[["found"]] <- as.data.frame(
  up_xref[["found"]] |>
    dplyr::bind_rows(up_xref[["rescued_from_grch38"]]) |>
    dplyr::distinct() |>
    dplyr::group_by(
      ensembl_transcript_id,
      uniprot_acc, uniprot_id
    ) |>
    dplyr::summarise(
      entrezgene = paste(entrezgene, collapse = "&"),
      .groups = "drop"
    ) |>
    dplyr::filter(!stringr::str_detect(
      entrezgene, "&"
    )) |>
    dplyr::group_by(
      entrezgene
    ) |>
    dplyr::mutate(
      entrezgene = as.integer(entrezgene)
    )
)


gene_gencode$records[["grch37"]] <-
  gene_gencode$records[["grch37"]] |>
  dplyr::select(-c(uniprot_id, uniprot_acc)) |>
  dplyr::left_join(
    up_xref$found,
    by = c(
      "entrezgene",
      "ensembl_transcript_id"
    )
  ) |>
  dplyr::distinct()

gene_summary <- get_function_summary_ncbi(gene_df = gene_info)

cgc_all <- get_cancer_gene_census(origin = "all")
cgc_som <- get_cancer_gene_census(origin = "somatic")
cgc_gl <- get_cancer_gene_census(origin = "germline")

cgc_som_gl <- cgc_som |>
  dplyr::full_join(
    cgc_gl,
    by = c(
      "cgc_moi",
      "cgc_tsg",
      "cgc_hallmark",
      "cgc_tier",
      "entrezgene",
      "cgc_oncogene"
    )
  ) |>
  dplyr::mutate(cgc_somatic = dplyr::if_else(
    is.na(cgc_somatic),
    as.logical(FALSE),
    as.logical(cgc_somatic)
  )) |>
  dplyr::mutate(cgc_germline = dplyr::if_else(
    is.na(cgc_germline),
    as.logical(FALSE),
    as.logical(cgc_germline)
  )) |>
  dplyr::select(-cgc_moi)

cgc <- cgc_all |>
  dplyr::left_join(cgc_som_gl)


intogen_drivers <- get_intogen_driver_genes(gene_info = gene_info)
fp_drivers <- get_curated_fp_cancer_genes(gene_info = gene_info)
ncg <- get_network_of_cancer_genes()
tcga_drivers <- get_tcga_driver_genes()
f1cdx <- get_f1cdx(gene_info = gene_info)
tso500 <- get_tso500(gene_info = gene_info, gene_alias = gene_alias)
dna_repair <- get_dna_repair_genes(gene_info = gene_info)
cancermine_genes <- get_cancermine_genes(
  cancermine_version = "48"
)
signaling_genes <- get_signaling_pathway_genes(gene_info = gene_info)
dbnsfp_annotations <- get_dbnsfp_gene_annotations()

gene_basic <- list()
gene_basic$metadata <- metadata[["basic"]]

gene_basic$records <- gene_info |>
  dplyr::select(-ensembl_gene_id) |>
  dplyr::left_join(gene_summary) |>
  dplyr::left_join(cgc) |>
  dplyr::left_join(ncg) |>
  dplyr::left_join(intogen_drivers) |>
  dplyr::left_join(fp_drivers) |>
  dplyr::left_join(dna_repair) |>
  dplyr::left_join(tso500) |>
  dplyr::left_join(f1cdx) |>
  dplyr::left_join(signaling_genes) |>
  dplyr::left_join(cancermine_genes) |>
  dplyr::left_join(tcga_drivers) |>
  dplyr::left_join(dbnsfp_annotations) |>
  dplyr::mutate(cancermine_n_cit_oncogene = dplyr::if_else(
    is.na(cancermine_n_cit_oncogene),
    as.integer(0),
    as.integer(cancermine_n_cit_oncogene)
  )) |>
  dplyr::mutate(cancermine_n_cit_tsg = dplyr::if_else(
    is.na(cancermine_n_cit_tsg),
    as.integer(0),
    as.integer(cancermine_n_cit_tsg)
  )) |>
  dplyr::mutate(cancermine_n_cit_driver = dplyr::if_else(
    is.na(cancermine_n_cit_driver),
    as.integer(0),
    as.integer(cancermine_n_cit_driver)
  )) |>
  dplyr::mutate(cgc_somatic = dplyr::if_else(
    is.na(cgc_somatic),
    as.logical(FALSE),
    as.logical(cgc_somatic)
  )) |>
  dplyr::mutate(cgc_driver_tier1 = dplyr::if_else(
    is.na(cgc_driver_tier1),
    as.logical(FALSE),
    as.logical(cgc_driver_tier1)
  )) |>
  dplyr::mutate(cgc_driver_tier2 = dplyr::if_else(
    is.na(cgc_driver_tier2),
    as.logical(FALSE),
    as.logical(cgc_driver_tier2)
  )) |>
  dplyr::mutate(cgc_germline = dplyr::if_else(
    is.na(cgc_germline),
    as.logical(FALSE),
    as.logical(cgc_germline)
  )) |>
  dplyr::mutate(cgc_tsg = dplyr::if_else(
    is.na(cgc_tsg),
    as.logical(FALSE),
    as.logical(cgc_tsg)
  )) |>
  dplyr::mutate(cgc_oncogene = dplyr::if_else(
    is.na(cgc_oncogene),
    as.logical(FALSE),
    as.logical(cgc_oncogene)
  )) |>
  dplyr::mutate(ncg_tsg = dplyr::if_else(
    is.na(ncg_tsg),
    as.logical(FALSE),
    as.logical(ncg_tsg)
  )) |>
  dplyr::mutate(ncg_driver = dplyr::if_else(
    is.na(ncg_driver),
    as.logical(FALSE),
    as.logical(ncg_driver)
  )) |>
  dplyr::mutate(ncg_oncogene = dplyr::if_else(
    is.na(ncg_oncogene),
    as.logical(FALSE),
    as.logical(ncg_oncogene)
  )) |>
  dplyr::mutate(bailey2018_fp_driver = dplyr::if_else(
    is.na(bailey2018_fp_driver),
    as.logical(FALSE),
    as.logical(bailey2018_fp_driver)
  )) |>
  dplyr::mutate(intogen_driver = dplyr::if_else(
    is.na(intogen_driver),
    as.logical(FALSE),
    as.logical(intogen_driver)
  )) |>
  dplyr::mutate(tcga_driver = dplyr::if_else(
    is.na(tcga_driver),
    as.logical(FALSE),
    as.logical(tcga_driver)
  )) |>
  dplyr::select(-synonyms)

gene_predisposition <- list()
gene_predisposition[["metadata"]] <- metadata$predisposition
gene_predisposition[["records"]] <-
  get_predisposition_genes(
    gene_info = gene_info,
    gene_panels = gene_panels,
    cache_dir = file.path(
      here::here(), "data-raw"
    )
  )

## clean up
rm(cgc_gl)
rm(cgc)
rm(cgc_som)
rm(intogen_drivers)
rm(fp_drivers)
rm(dna_repair)
rm(signaling_genes)
rm(cancermine_genes)
rm(tso500)
rm(gene_summary)
rm(gene_info)
rm(dbnsfp_annotations)
rm(metadata)
rm(ncg)
rm(cgc_all)
rm(tcga_drivers)
rm(cgc_som_gl)
rm(f1cdx)
rm(up_xref)
rm(up_xref_grch38)


## upload to Google Drive
version_minor_bumped <- paste0(
  "0.",
  as.character(as.integer(substr(
    as.character(
      packageVersion("geneOncoX")
    ), 3, 3
  )) + 1),
  ".0"
)

gd_records <- list()
db_id_ref <- data.frame()

db <- list()
db[["basic"]] <- gene_basic
db[["gencode"]] <- gene_gencode
db[["alias"]] <- gene_alias
db[["predisposition"]] <- gene_predisposition
db[["panels"]] <- gene_panels

rm(gene_alias)
rm(gene_basic)
rm(gene_panels)
rm(gene_predisposition)

# googledrive::drive_auth_configure(api_key = Sys.getenv("GD_KEY"))

for (elem in c(
  "basic", "predisposition",
  "panels", "alias", "gencode"
)) {
  saveRDS(db[[elem]],
    file = paste0(
      "data-raw/gd_local/gene_", elem, "_v",
      version_minor_bumped, ".rds"
    )
  )

  (gd_records[[elem]] <- googledrive::drive_upload(
    paste0(
      "data-raw/gd_local/gene_",
      elem, "_v", version_minor_bumped, ".rds"
    ),
    paste0("geneOncoX/gene_", elem, "_v", version_minor_bumped, ".rds")
  ))

  google_rec_df <-
    dplyr::select(
      as.data.frame(gd_records[[elem]]), name, id
    ) |>
    dplyr::rename(
      gid = id,
      filename = name
    ) |>
    dplyr::mutate(
      name = stringr::str_replace(
        stringr::str_replace(filename, "_v\\S+$", ""),
        "gene_", ""
      ),
      date = as.character(Sys.Date()),
      pVersion = version_minor_bumped
    ) |>
    dplyr::mutate(
      md5Checksum =
        gd_records[[elem]]$drive_resource[[1]]$md5Checksum
    )

  db_id_ref <- db_id_ref |>
    dplyr::bind_rows(google_rec_df)
}

usethis::use_data(db_id_ref, internal = TRUE, overwrite = TRUE)
