#' Function to write virtual panel BED files (CPSR)
#'
#' @param gwas_bed_fpath path to assembly-specific GWAS BED track
#' @param gox_gencode object returned by `geneOncoX::get_gencode()`
#' @param gox_panels object returned by `geneOncoX::get_panels()`
#' @param gox_predisposition object returned by
#'   `geneOncoX::get_predisposition()`
#' @param super_panel logical indicating if a BED files for a CPSR super panel
#'   are to be generated
#' @param build genome assembly build (grch37/grch38)
#' @param dest_dir destination directory for BED file
#'
#' @return integer indicating success (0) or failure (-1)
#'
#' @keywords internal
#'
create_virtual_panel_bed <- function(gwas_bed_fpath = NA,
                                     gox_gencode = NULL,
                                     gox_panels = NULL,
                                     gox_predisposition = NULL,
                                     super_panel = FALSE,
                                     build = "grch37",
                                     dest_dir = NA) {
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T")
  )

  if (is.na(dest_dir)) {
    lgr::lgr$error(paste0(
      "Argument dest_dir = '",
      dest_dir, "' is not defined"
    ))
  }

  if (!dir.exists(dest_dir)) {
    lgr::lgr$error(paste0(
      "Argument dest_dir = '",
      dest_dir, "' does not exist"
    ))
  }

  if (is.na(gwas_bed_fpath)) {
    lgr::lgr$error(paste0(
      "Argument gwas_bed_fpath = '",
      gwas_bed_fpath, "' is undefined/NA"
    ))
  }

  if (nchar(gwas_bed_fpath) > 0) {
    if (!(file.exists(gwas_bed_fpath) &&
      endsWith(gwas_bed_fpath, ".bed.gz"))) {
      lgr::lgr$error(paste0(
        "File gwas_bed_file = '",
        gwas_bed_fpath, "' does not exist"
      ))
    }
  }

  ## read build-specific GWAS BED file (argument)
  gwas_bed_fpath_uncomp <- file.path(
    dest_dir, stringr::str_replace(basename(gwas_bed_fpath), ".gz", "")
  )
  system(paste0("bgzip -dc ", gwas_bed_fpath, " > ", gwas_bed_fpath_uncomp))
  gwas_bed_data <- readr::read_tsv(
    gwas_bed_fpath_uncomp,
    col_names = FALSE,
    show_col_types = FALSE
  )
  colnames(gwas_bed_data) <- c("chrom", "start", "end", "name")
  system(paste0("rm -f ", gwas_bed_fpath_uncomp))

  ## Get all build-specific Genomics England PanelApp records
  all_panel_records <- gox_panels$records |>
    dplyr::filter(.data$genome_build == build)

  ## Get build-specific GENCODE transcripts
  gencode_transcripts <- gox_gencode$records[[build]] |>
    dplyr::mutate(chrom = as.character(
      stringr::str_replace(.data$chrom, "chr", "")
    ))

  ## Get ACMG secondary-findings genes and intersect
  ## them with GENCODE tracks -> BED data (chrom, start, end, name)
  acmg_sf_bed_data <- gox_predisposition$records |>
    dplyr::filter(stringr::str_detect(.data$predisp_source, "ACMG_SF")) |>
    dplyr::select(.data$entrezgene, .data$moi) |>
    dplyr::inner_join(
      dplyr::select(
        gencode_transcripts,
        .data$chrom,
        .data$start,
        .data$end,
        .data$ensembl_gene_id,
        .data$ensembl_transcript_id,
        .data$entrezgene,
        .data$symbol
      ),
      by = c("entrezgene")
    ) |>
    dplyr::mutate(name = paste(
      .data$symbol,
      .data$ensembl_gene_id,
      .data$entrezgene,
      .data$ensembl_transcript_id,
      "ACMG_SF", .data$moi,
      sep = "|"
    )) |>
    dplyr::select(
      .data$chrom, .data$start,
      .data$end, .data$name
    ) |>
    dplyr::mutate(chrom = as.character(
      stringr::str_replace(.data$chrom, "chr", "")
    ))


  unique_panel_ids <- unique(all_panel_records$id)
  i <- 1
  while (i <= length(unique_panel_ids)) {
    panel_id <- unique_panel_ids[i]
    panel_name <- unique(
      all_panel_records[all_panel_records$id == i, "gepa_panel_name"]
    )

    ## Get genes for the specific panel, and intersect them
    ## with GENCODE transcripts  --> BED format (chrom, start, end, name)
    panel_bed_data <- dplyr::filter(
      all_panel_records, .data$id == i
    ) |>
      dplyr::select(
        c("id", 
        "entrezgene",
        "ensembl_gene_id",
        "gepa_panel_name",
        "gepa_panel_id",
        "gepa_confidence_level",
        "gepa_moi",
        "gepa_panel_version")
      ) |>
      dplyr::inner_join(
        dplyr::select(
          gencode_transcripts,
          c("chrom",
          "start",
          "end",
          "ensembl_gene_id",
          "ensembl_transcript_id",
          "entrezgene", 
          "symbol")
        ),
        by = c("ensembl_gene_id", "entrezgene")
      ) |>
      dplyr::mutate(name = paste(
        .data$symbol,
        .data$ensembl_gene_id,
        .data$entrezgene,
        .data$ensembl_transcript_id,
        stringr::str_replace_all(
          .data$gepa_panel_name, " ", "_"
        ),
        paste0("GEPA_PANEL_", .data$id),
        .data$gepa_confidence_level,
        .data$gepa_panel_id,
        .data$gepa_panel_version,
        .data$gepa_moi,
        sep = "|"
      )) |>
      dplyr::select(
        c("chrom",
        "start",
        "end",
        "name",
        "gepa_confidence_level")
      ) |>
      dplyr::mutate(chrom = as.character(
        stringr::str_replace(.data$chrom, "chr", "")
      ))


    ## Add BED records with ACMG genes and GWAS variants
    ## - For all and high-confidence genes only (green)

    virtual_panels <- list()
    virtual_panels[["green"]] <- panel_bed_data |>
      dplyr::filter(.data$gepa_confidence_level == 3) |>
      dplyr::select(-.data$gepa_confidence_level)
    virtual_panels[["all"]] <- panel_bed_data |>
      dplyr::select(-.data$gepa_confidence_level)

    if (nrow(virtual_panels[["green"]]) == 0) {
      lgr::lgr$warn(paste0(
        "Panel - ",
        panel_name, " has zero GREEN entries - using all"
      ))
      virtual_panels[["green"]] <- virtual_panels[["all"]]
    }

    ## Green
    acmg_sf_bed_data_missing <- acmg_sf_bed_data |>
      dplyr::anti_join(virtual_panels[["green"]],
        by = c("chrom", "start", "end")
      )

    virtual_panels[["green"]] <- virtual_panels[["green"]] |>
      dplyr::bind_rows(acmg_sf_bed_data_missing) |>
      dplyr::bind_rows(gwas_bed_data)

    f <- write_bed_file(
      bed_data = virtual_panels[["green"]],
      bed_fname = paste0(panel_id, ".", build, ".GREEN.bed"),
      dest_dir = dest_dir
    )

    lgr::lgr$info(
      paste0("Wrote BED track: ", f)
    )

    ## All
    acmg_sf_bed_data_missing <- acmg_sf_bed_data |>
      dplyr::anti_join(virtual_panels[["all"]],
        by = c("chrom", "start", "end")
      )

    virtual_panels[["all"]] <- virtual_panels[["all"]] |>
      dplyr::bind_rows(acmg_sf_bed_data_missing) |>
      dplyr::bind_rows(gwas_bed_data)

    f <- write_bed_file(
      bed_data = virtual_panels[["all"]],
      bed_fname = paste0(
        panel_id, ".", build, ".bed"
      ),
      dest_dir = dest_dir
    )

    lgr::lgr$info(
      paste0("Wrote BED track: ", f)
    )

    i <- i + 1
  }
}

#' Function to sort data frame with genomic regions (simple BED format)
#'
#' @param unsorted_regions data frame with unsorted regions (BED)
#' @param chr_names logical indicating if chromosome names contain 'chr' or not
#'
#' @return data frame with sorted BED regions
#'
#' @keywords internal
#'
sort_bed_regions <- function(unsorted_regions, chr_names = FALSE) {
  sorted_regions <- NULL
  assertable::assert_colnames(
    unsorted_regions, c("start", "end"),
    only_colnames = FALSE, quiet = TRUE
  )
  if ("chrom" %in% colnames(unsorted_regions) &&
    "start" %in% colnames(unsorted_regions) &&
    "end" %in% colnames(unsorted_regions)) {
    chrOrder <- c(
      as.character(seq_len(22)),
      "X", "Y", "M"
    )
    if (chr_names == TRUE) {
      chrOrder <- paste0("chr", c(
        as.character(seq_len(22)),
        "X", "Y", "M"
      ))
    }
    unsorted_regions$chrom <-
      factor(unsorted_regions$chrom, levels = chrOrder)
    unsorted_regions <- unsorted_regions[order(unsorted_regions$chrom), ]

    sorted_regions <- data.frame()
    for (chrom in chrOrder) {
      if (nrow(unsorted_regions[unsorted_regions$chrom == chrom, ]) > 0) {
        chrom_regions <-
          unsorted_regions[unsorted_regions$chrom == chrom, ]
        chrom_regions_sorted <-
          chrom_regions[with(chrom_regions, order(start, end)), ]
        sorted_regions <-
          dplyr::bind_rows(sorted_regions, chrom_regions_sorted)
      }
    }

    sorted_regions$start <- as.integer(sorted_regions$start)
    sorted_regions$end <- as.integer(sorted_regions$end)
  }
  return(sorted_regions)
}

#' Function that writes BED data to file
#'
#' @param bed_data data frame with BED data
#' @param bed_fname file name for BED file
#' @param dest_dir directory to write BED file
#'
#' @return file name of bgzipped BED file with BED track data
#'
#' @keywords internal
#'
write_bed_file <- function(bed_data,
                           bed_fname,
                           dest_dir = NA) {
  if (!dir.exists(dest_dir)) {
    lgr::lgr$error(paste0(
      "Argument dest_dir = '",
      dest_dir, "' does not exist"
    ))
  }

  bed_fname_full <- file.path(dest_dir, bed_fname)
  bed_data_sorted <- sort_bed_regions(bed_data)

  utils::write.table(bed_data_sorted,
    file = bed_fname_full, sep = "\t",
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  system(paste0("bgzip -f ", bed_fname_full))
  system(paste0("tabix -p bed ", bed_fname_full, ".gz"))

  bed_fname_bgzipped <-
    paste0(bed_fname_full, ".gz")

  return(bed_fname_full)
}

#' Function to classify genes as tumor suppressors/proto-oncogenes,
#' driver etc. based on multiple lines of evidence.
#'
#' @param gox_basic output from geneOncoX::get_basic()
#' @param min_citations_tsg minimum citations (CancerMine) for tsg support
#' @param min_citations_oncogene minimum citations
#'   (CancerMine) for oncogene support
#'
#' @return data frame with cancer gene annotations
#' @keywords internal
#'
assign_cancer_gene_evidence <- function(gox_basic = NULL,
                                        min_citations_tsg = 15,
                                        min_citations_oncogene = 15) {
  tsg_oncogene_driver <- list()

  #### CancerMine: Tumor suppressor genes and oncogene annotation
  ## Considering
  ## A) minimum number of supporting citations, and
  ## B) ratio of tsg/oncogene citations (ignore entries in which the support
  ##    is three times as many for the opposite role)

  tsg_oncogene_driver[["CancerMine"]] <- gox_basic$records |>
    dplyr::mutate(cm_oncogene_tsg_citratio = dplyr::if_else(
      !is.na(.data$cancermine_n_cit_tsg),
      as.numeric(.data$cancermine_n_cit_oncogene) /
        .data$cancermine_n_cit_tsg,
      as.numeric(NA)
    )) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      !is.na(.data$cancermine_n_cit_oncogene) &
        .data$cancermine_n_cit_oncogene >= min_citations_oncogene,
      TRUE, FALSE
    )) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      !is.na(.data$cm_oncogene_tsg_citratio) &
        .data$cm_oncogene_tsg_citratio <= 0.33,
      FALSE,
      as.logical(.data$oncogene)
    )) |>
    dplyr::mutate(tsg = dplyr::if_else(
      !is.na(.data$cancermine_n_cit_tsg) &
        .data$cancermine_n_cit_tsg >= min_citations_tsg,
      TRUE, FALSE
    )) |>
    dplyr::mutate(tsg = dplyr::if_else(
      !is.na(.data$cm_oncogene_tsg_citratio) &
        .data$cm_oncogene_tsg_citratio > 3,
      FALSE,
      as.logical(.data$tsg)
    )) |>
    dplyr::mutate(
      links_oncogene = dplyr::if_else(
        !is.na(.data$cancermine_cit_links_oncogene),
        paste0(
          "Oncogenic role (CancerMine): ",
          .data$cancermine_cit_links_oncogene
        ),
        ""
      ),
      links_tsg = dplyr::if_else(
        !is.na(.data$cancermine_cit_links_tsg),
        paste0(
          "Tumor suppressive role (CancerMine): ",
          .data$cancermine_cit_links_tsg
        ),
        ""
      ),
      links_driver = dplyr::if_else(
        !is.na(.data$cancermine_cit_links_driver),
        paste0(
          "Cancer driver role (CancerMine): ",
          .data$cancermine_cit_links_driver
        ),
        ""
      )
    ) |>
    dplyr::mutate(source = "CancerMine") |>
    dplyr::select(
      c("entrezgene",
        "oncogene",
        "tsg",
        "links_driver",
        "links_oncogene",
        "links_tsg",
        "source")
    ) |>
    dplyr::filter(.data$oncogene == TRUE |
      .data$tsg == TRUE |
      stringr::str_detect(
        .data$links_driver, ": <"
      )) |>
    dplyr::distinct()

  ### Network of cancer genes (NCG): Proto-oncogenes,
  ### tumor suppressor genes and cancer drivers

  tsg_oncogene_driver[["NCG"]] <- gox_basic$records |>
    dplyr::mutate(oncogene = dplyr::if_else(
      .data$ncg_oncogene == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(tsg = dplyr::if_else(
      .data$ncg_tsg == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(driver = dplyr::if_else(
      .data$ncg_driver == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(
      links_oncogene = dplyr::if_else(
        .data$oncogene == TRUE,
        paste0(
          "Oncogenic role (NCG): ",
          "<a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "\" target=\"_blank\">",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      ),
      links_tsg = dplyr::if_else(
        .data$tsg == TRUE,
        paste0(
          "Tumor suppressive role (NCG): ",
          "<a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "\" target=\"_blank\">",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      ),
      links_driver = dplyr::if_else(
        .data$driver == TRUE,
        paste0(
          "Cancer driver role (NCG): ",
          "<a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "\" target=\"_blank\">",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      )
    ) |>
    dplyr::mutate(source = "NCG") |>
    dplyr::select(
      c("entrezgene",
        "oncogene",
        "tsg",
        "driver",
        "links_driver",
        "links_oncogene",
        "links_tsg",
        "source")
    ) |>
    dplyr::filter(
      .data$oncogene == TRUE |
        .data$tsg == TRUE |
        stringr::str_detect(
          .data$links_driver, ": <"
        )
    ) |>
    dplyr::distinct()


  ### Cancer Gene Census: Proto-oncogenes,
  ### tumor suppressor genes and cancer drivers (all somatic + germline)

  tsg_oncogene_driver[["CGC"]] <- gox_basic$records |>
    dplyr::mutate(oncogene = dplyr::if_else(
      .data$cgc_oncogene == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(tsg = dplyr::if_else(
      .data$cgc_tsg == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(driver = dplyr::if_else(
      .data$cgc_driver_tier1 == TRUE |
        .data$cgc_driver_tier2 == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(hallmark = dplyr::if_else(
      .data$cgc_hallmark == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(
      links_oncogene = dplyr::if_else(
        .data$oncogene == TRUE,
        paste0(
          "Oncogenic role (CGC): ",
          "<a href=\"https://cancer.sanger.ac.uk/census",
          "\" target=\"_blank\">",
          "YES (TIER ",
          .data$cgc_tier, ")</a>"
        ),
        ""
      ),
      links_tsg = dplyr::if_else(
        .data$tsg == TRUE,
        paste0(
          "Tumor suppressive role (CGC): ",
          "<a href=\"https://cancer.sanger.ac.uk/census",
          "\" target=\"_blank\">",
          "YES (TIER ",
          .data$cgc_tier, ")</a>"
        ),
        ""
      ),
      links_driver = dplyr::if_else(
        .data$driver == TRUE,
        paste0(
          "Cancer driver role (CGC): ",
          "<a href=\"https://cancer.sanger.ac.uk/census",
          "\" target=\"_blank\">",
          "YES (TIER ",
          .data$cgc_tier, ")</a>"
        ),
        ""
      )
    ) |>
    dplyr::mutate(source = dplyr::if_else(
      .data$cgc_tier == 1,
      "CGC_TIER1",
      "CGC_TIER2"
    )) |>
    dplyr::select(
      c("entrezgene",
        "oncogene",
        "tsg",
        "driver",
        "hallmark",
        "links_driver",
        "links_oncogene",
        "links_tsg",
        "source")
    ) |>
    dplyr::filter(.data$oncogene == TRUE |
      .data$tsg == TRUE |
      .data$driver == TRUE) |>
    dplyr::distinct()

  ### IntOGen: Predicted cancer driver genes

  tsg_oncogene_driver[["IntOGen"]] <- gox_basic$records |>
    dplyr::mutate(driver = dplyr::if_else(
      .data$intogen_driver == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(
      links_driver = dplyr::if_else(
        .data$intogen_driver == TRUE,
        paste0(
          "Cancer driver role (IntOGen): ",
          "<a href=\"https://www.intogen.org/search?gene=",
          .data$symbol,
          "\" target=\"_blank\">",
          "YES</a>"
        ),
        ""
      )
    ) |>
    dplyr::mutate(source = "IntOGen") |>
    dplyr::filter(.data$intogen_driver == TRUE) |>
    dplyr::select(
      c("entrezgene",
        "driver",
        "links_driver",
        "source")
    ) |>
    dplyr::distinct()


  tsg_oncogene_driver[["TCGA"]] <- gox_basic$records |>
    dplyr::mutate(driver = dplyr::if_else(
      .data$tcga_driver == TRUE,
      TRUE, FALSE
    )) |>
    dplyr::mutate(
      links_driver = dplyr::if_else(
        .data$tcga_driver == TRUE,
        paste0(
          "Cancer driver role (TCGA PanCancer): ",
          paste0(
            "<a href=\"https://gdc.cancer.gov/about-data/",
            "publications/pancan-driver"
          ),
          "\" target=\"_blank\">YES</a>"
        ),
        ""
      )
    ) |>
    dplyr::mutate(source = "TCGA") |>
    dplyr::filter(.data$tcga_driver == TRUE) |>
    dplyr::select(
      c("entrezgene",
        "driver",
        "links_driver",
        "source")
    ) |>
    dplyr::distinct()


  driver_genes <- as.data.frame(dplyr::bind_rows(
    tsg_oncogene_driver$CGC,
    tsg_oncogene_driver$IntOGen
  ) |>
    dplyr::bind_rows(tsg_oncogene_driver$NCG) |>
    dplyr::bind_rows(
      dplyr::filter(
        tsg_oncogene_driver$CancerMine,
        nchar(.data$links_driver) > 0
      )
    ) |>
    dplyr::bind_rows(tsg_oncogene_driver$TCGA) |>
    dplyr::mutate(
      driver = dplyr::if_else(
        is.na(.data$driver),
        FALSE,
        as.logical(.data$driver)
      )
    ) |>
    dplyr::mutate(
      links_driver = dplyr::if_else(
        is.na(.data$links_driver),
        "",
        as.character(.data$links_driver)
      )
    ) |>
    dplyr::group_by(.data$entrezgene) |>
    dplyr::summarise(
      driver = paste(unique(.data$driver), collapse = "&"),
      driver_links = paste(unique(.data$links_driver),
        collapse = ", "
      ),
      driver_support = paste(unique(sort(.data$source)),
        collapse = "&"
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$driver != FALSE) |>
    dplyr::mutate(driver = dplyr::if_else(
      stringr::str_detect(.data$driver, "&"),
      TRUE,
      as.logical(.data$driver)
    )))


  oncogenes <- as.data.frame(dplyr::bind_rows(
    tsg_oncogene_driver$CGC,
    tsg_oncogene_driver$NCG
  ) |>
    dplyr::bind_rows(tsg_oncogene_driver$CancerMine) |>
    dplyr::mutate(
      oncogene = dplyr::if_else(
        is.na(.data$oncogene),
        FALSE,
        as.logical(.data$oncogene)
      )
    ) |>
    dplyr::mutate(
      links_oncogene = dplyr::if_else(
        is.na(.data$links_oncogene),
        "",
        as.character(.data$links_oncogene)
      )
    ) |>
    dplyr::group_by(.data$entrezgene) |>
    dplyr::summarise(
      oncogene = paste(unique(.data$oncogene), collapse = "&"),
      oncogene_links = paste(unique(.data$links_oncogene),
        collapse = ", "
      ),
      oncogene_support = paste(unique(sort(.data$source)),
        collapse = "&"
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$oncogene != FALSE) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      stringr::str_detect(.data$oncogene, "&"),
      TRUE,
      as.logical(.data$oncogene)
    )))

  tsgs <- as.data.frame(dplyr::bind_rows(
    tsg_oncogene_driver$CGC,
    tsg_oncogene_driver$NCG
  ) |>
    dplyr::bind_rows(tsg_oncogene_driver$CancerMine) |>
    dplyr::mutate(
      tsg = dplyr::if_else(
        is.na(.data$tsg),
        FALSE,
        as.logical(.data$tsg)
      )
    ) |>
    dplyr::mutate(
      links_tsg = dplyr::if_else(
        is.na(.data$links_tsg),
        "",
        as.character(.data$links_tsg)
      )
    ) |>
    dplyr::group_by(.data$entrezgene) |>
    dplyr::summarise(
      tsg = paste(unique(.data$tsg), collapse = "&"),
      tsg_links = paste(unique(.data$links_tsg),
        collapse = ", "
      ),
      tsg_support = paste(unique(sort(.data$source)),
        collapse = "&"
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$tsg != FALSE) |>
    dplyr::mutate(tsg = dplyr::if_else(
      stringr::str_detect(.data$tsg, "&"),
      TRUE,
      as.logical(.data$tsg)
    )))

  tsg_oncogene_driver_evidence <-
    dplyr::full_join(driver_genes, tsgs, by = "entrezgene") |>
    dplyr::full_join(oncogenes, by = "entrezgene") |>
    dplyr::mutate(oncogene = dplyr::if_else(
      is.na(.data$oncogene), FALSE, as.logical(.data$oncogene)
    )) |>
    dplyr::mutate(tsg = dplyr::if_else(
      is.na(.data$tsg), FALSE, as.logical(.data$tsg)
    )) |>
    dplyr::mutate(driver = dplyr::if_else(
      is.na(.data$driver), FALSE, as.logical(.data$driver)
    )) |>
    dplyr::mutate(oncogene_links = stringr::str_replace_all(
      stringr::str_replace(
        .data$oncogene_links, "^, ", ""
      ),
      ", , ", ", "
    )) |>
    dplyr::mutate(tsg_links = stringr::str_replace_all(
      stringr::str_replace_all(
        .data$tsg_links, "^, ", ""
      ),
      ", , ", ", "
    )) |>
    dplyr::mutate(driver_links = stringr::str_replace_all(
      stringr::str_replace_all(
        .data$driver_links, "^, ", ""
      ),
      ", , ", ", "
    ))

  num_oncogene_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$oncogene == TRUE
    )
  num_oncogene_multisupport <- num_oncogene_all |>
    dplyr::filter(
      stringr::str_detect(
        .data$oncogene_support, "&"
      )
    )
  num_tsg_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$tsg == TRUE
    )
  num_tsg_multisupport <- num_tsg_all |>
    dplyr::filter(
      stringr::str_detect(
        .data$tsg_support, "&"
      )
    )
  num_driver_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$driver == TRUE
    )
  num_driver_multisupport <- num_driver_all |>
    dplyr::filter(
      stringr::str_detect(
        .data$driver_support, "&"
      )
    )

  lgr::lgr$info(
    paste0(
      "Sources: Network of Cancer Genes (NCG), Cancer Gene ",
      "Census (CGC_TIER1, CGC_TIER2)"
    )
  )
  lgr::lgr$info(
    paste0(
      "Sources: CancerMine, IntOGen, TCGA (PanCancer driver ",
      "gene classification)"
    )
  )
  lgr::lgr$info(
    paste0(
      "Driver/proto-oncogene/tumor suppressor gene stats - ",
      "support from at least one source:"
    )
  )
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(num_oncogene_all),
      " classified proto-oncogenes"
    )
  )
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(num_tsg_all),
      " classified tumor suppressor genes"
    )
  )
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(num_driver_all),
      " genes were annotated as potential cancer drivers"
    )
  )
  lgr::lgr$info(
    "-----------------------------------------------------------------"
  )
  lgr::lgr$info(
    paste0(
      "Driver/proto-oncogene/tumor suppressor gene stats ",
      "- support from multiple sources:"
    )
  )
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(num_oncogene_multisupport),
      " classified proto-oncogenes"
    )
  )
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(num_tsg_multisupport),
      " classified tumor suppressor genes"
    )
  )
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(num_driver_multisupport),
      " genes were annotated as potential drivers"
    )
  )


  return(tsg_oncogene_driver_evidence)
}

#' Tidy eval helpers
#'
#' <https://cran.r-project.org/web/packages/dplyr/vignettes/programming.html>
#'
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang .data :=
NULL
