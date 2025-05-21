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
#' @param min_citation_support minimum citation count for support from CancerMine 
#' @param min_citation_support_cm_only minimum citation count for genes annotated 
#' with cancer gene roles from CancerMine only 
#' @param min_sources_driver minimum number of sources (NCG, CGC, CancerMine,  
#' TCGA, IntoGen) that must contribute to cancer driver gene status
#'
#' @return data frame with cancer gene annotation status
#' @keywords internal
#'
assign_cancer_gene_roles_old <- function(gox_basic = NULL,
                                     min_citation_support = 5,
                                     min_citation_support_cm_only = 20,
                                     min_sources_driver = 2) {
  tsg_oncogene_driver <- list()
  
  #### CancerMine: Tumor suppressor genes and oncogene annotation
  ## Considering
  ## A) minimum number of supporting citations, and
  ## B) ratio of tsg/oncogene citations (ignore entries in which the support
  ##    is three times as many for the opposite role)
  
  tsg_oncogene_driver[["CancerMine"]] <- gox_basic$records |>
    dplyr::mutate(cancermine_oncogene_tsg_citratio = dplyr::if_else(
      !is.na(.data$cancermine_n_cit_tsg) & .data$cancermine_n_cit_tsg > 0,
      round(as.numeric(.data$cancermine_n_cit_oncogene) /
        .data$cancermine_n_cit_tsg, digits = 1),
      as.numeric(NA)
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
        "cancermine_n_cit_oncogene",
        "cancermine_n_cit_tsg",
        "cancermine_n_cit_driver",
        "cancermine_oncogene_tsg_citratio",
        "links_driver",
        "links_oncogene",
        "links_tsg",
        "source")
    ) |>
    dplyr::filter(
      stringr::str_detect(
        .data$links_oncogene, ": <"
      ) |
        stringr::str_detect(
          .data$links_driver, ": <"
        ) |
        stringr::str_detect(
          .data$links_oncogene, ": <"
        )
      
      ) |>
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
          "<a href='http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "' target='_blank'>",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      ),
      links_tsg = dplyr::if_else(
        .data$tsg == TRUE,
        paste0(
          "Tumor suppressive role (NCG): ",
          "<a href='http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "' target='_blank'>",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      ),
      links_driver = dplyr::if_else(
        .data$driver == TRUE,
        paste0(
          "Cancer driver role (NCG): ",
          "<a href='http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "' target='_blank'>",
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
          "<a href='https://cancer.sanger.ac.uk/census",
          "' target='_blank'>",
          "YES (TIER ",
          .data$cgc_tier, ")</a>"
        ),
        ""
      ),
      links_tsg = dplyr::if_else(
        .data$tsg == TRUE,
        paste0(
          "Tumor suppressive role (CGC): ",
          "<a href='https://cancer.sanger.ac.uk/census",
          "' target='_blank'>",
          "YES (TIER ",
          .data$cgc_tier, ")</a>"
        ),
        ""
      ),
      links_driver = dplyr::if_else(
        .data$driver == TRUE,
        paste0(
          "Cancer driver role (CGC): ",
          "<a href='https://cancer.sanger.ac.uk/census",
          "' target='_blank'>",
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
    dplyr::filter(
      .data$oncogene == TRUE |
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
          "<a href='https://www.intogen.org/search?gene=",
          .data$symbol,
          "' target='_blank'>",
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
            "<a href='https://gdc.cancer.gov/about-data/",
            "publications/pancan-driver"
          ),
          "' target='_blank'>YES</a>"
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
    dplyr::filter(tsg_oncogene_driver$CGC, .data$driver == T),
    tsg_oncogene_driver$IntOGen) |>
      dplyr::select(-c("tsg","oncogene","links_oncogene",
                       "links_tsg","hallmark")) |>
      dplyr::bind_rows(
        dplyr::filter(tsg_oncogene_driver$NCG, .data$driver == T)) |>
      dplyr::select(-c("tsg","oncogene","links_oncogene",
                       "links_tsg")) |>
      dplyr::bind_rows(
        dplyr::filter(
          tsg_oncogene_driver$CancerMine,
          nchar(.data$links_driver) > 0 &
            .data$cancermine_n_cit_driver >= min_citation_support
        )
      ) |>
      dplyr::select(-c("cancermine_n_cit_oncogene",
                       "cancermine_n_cit_tsg",
                       "cancermine_n_cit_driver",
                       "cancermine_oncogene_tsg_citratio",
                       "links_oncogene",
                       "links_tsg")) |>
      dplyr::bind_rows(tsg_oncogene_driver$TCGA) |>
      dplyr::mutate(
        driver = dplyr::if_else(
          is.na(.data$driver),
          TRUE,
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
      dplyr::filter(
        stringr::str_count(.data$driver_support, "&") + 1 >= 
          min_sources_driver
      )
  )
  
  
  oncogenes <- as.data.frame(dplyr::bind_rows(
    dplyr::filter(tsg_oncogene_driver$CGC, .data$oncogene == T),
    dplyr::filter(tsg_oncogene_driver$NCG, .data$oncogene == T)) |>
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
        full_oncogene_links = paste(unique(.data$links_oncogene),
                               collapse = ", "
        ),
        full_oncogene_support = paste(unique(sort(.data$source)),
                                 collapse = "&"
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(oncogene = dplyr::if_else(
        stringr::str_detect(.data$oncogene, "&"),
        TRUE,
        as.logical(.data$oncogene)
      ))) |>
    dplyr::full_join(
      dplyr::select(
        tsg_oncogene_driver$CancerMine,
        "entrezgene",
        "cancermine_n_cit_oncogene",
        "cancermine_n_cit_tsg",
        "cancermine_oncogene_tsg_citratio",
        "links_oncogene"
      ), by = "entrezgene"
    ) |>
    dplyr::mutate(
      cancermine_n_cit_oncogene = dplyr::if_else(
        is.na(.data$cancermine_n_cit_oncogene),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_oncogene))) |> 
    dplyr::mutate(
      cancermine_n_cit_tsg = dplyr::if_else(
        is.na(.data$cancermine_n_cit_tsg),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_tsg))) |>
    dplyr::filter(
      !(is.na(.data$oncogene) & 
          .data$cancermine_n_cit_oncogene < min_citation_support)) |>
    dplyr::mutate(
      links_oncogene = dplyr::if_else(
        !is.na(.data$links_oncogene) &
          .data$cancermine_n_cit_oncogene < min_citation_support,
        "",
        as.character(.data$links_oncogene)
      )
    ) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      
      ## Text-mined oncogene support only
      ((is.na(oncogene) & 
        cancermine_oncogene_tsg_citratio >= 4 &
        cancermine_n_cit_oncogene >= min_citation_support_cm_only) |
         
         ## Text-mined dual role TSG/oncogene 
        (is.na(oncogene) & 
           cancermine_oncogene_tsg_citratio <= 2 &
           cancermine_oncogene_tsg_citratio >= 0.5 &
           cancermine_n_cit_oncogene >= min_citation_support_cm_only &
           cancermine_n_cit_tsg >= min_citation_support_cm_only)),
      TRUE,
      as.logical(oncogene)
    )) |>
    dplyr::mutate(
      oncogene_links = dplyr::case_when(
        is.na(full_oncogene_links) &
          oncogene == T &
          (!is.na(.data$links_oncogene) |
             nchar(.data$links_oncogene) > 0) ~ links_oncogene,
        !is.na(full_oncogene_links) &
          oncogene == T &
          !is.na(.data$links_oncogene) &
          nchar(.data$links_oncogene) > 0 ~
          paste(
            full_oncogene_links, links_oncogene, sep=", "
          ),
        TRUE ~ as.character(.data$full_oncogene_links)
      ) 
    ) |>
    dplyr::select(-c("links_oncogene","full_oncogene_links")) |>
    dplyr::mutate(oncogene_support = dplyr::case_when(
      oncogene == T & 
      !is.na(full_oncogene_support) &
        cancermine_n_cit_oncogene >= min_citation_support ~
        paste(full_oncogene_support, 
              paste0("CancerMine:",cancermine_n_cit_oncogene), 
              sep = "&"),
      oncogene == T &
      is.na(full_oncogene_support) ~
        paste0("CancerMine:",cancermine_n_cit_oncogene),
      TRUE ~ as.character(full_oncogene_support))) |>
    dplyr::mutate(oncogene_support = stringr::str_replace(
      .data$oncogene_support, "^NA&","")
    ) |>
    dplyr::filter(oncogene == TRUE)
  
  ranked_oncogenes <- list()
  ranked_oncogenes[['part1']] <- oncogenes |> 
    dplyr::filter(
      stringr::str_detect(
        oncogene_support, "&")) |> 
    dplyr::arrange(
      desc(stringr::str_count(oncogene_support,"&")),
      desc(nchar(oncogene_support)),
      desc(cancermine_n_cit_oncogene))
  
  ranked_oncogenes[['part2']] <- oncogenes |>
    dplyr::filter(
      !stringr::str_detect(
        oncogene_support, "&") &
        !stringr::str_detect(oncogene_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_oncogene))
  
  ranked_oncogenes[['part3']] <- oncogenes |>
    dplyr::filter(
     !stringr::str_detect(
        oncogene_support, "&") &
        stringr::str_detect(oncogene_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_oncogene))
  
  oncogenes <- 
    dplyr::bind_rows(
      ranked_oncogenes[['part1']], 
      ranked_oncogenes[['part2']],
      ranked_oncogenes[['part3']]) |>
    dplyr::mutate(ind = dplyr::row_number()) |>
    dplyr::mutate(
      oncogene_rank = round(
        100 - (100 * dplyr::percent_rank(ind)), digits = 0)) |>
    dplyr::select(
      -c("full_oncogene_support",
         "cancermine_n_cit_tsg",
         "ind",
         "cancermine_n_cit_oncogene",
         "cancermine_oncogene_tsg_citratio")) |>
    dplyr::mutate(
      oncogene_confidence_level = dplyr::case_when(
        stringr::str_count(oncogene_support,"&") == 2 ~ "Very strong",
        stringr::str_count(oncogene_support,"&") == 1 ~ "Strong",
        stringr::str_count(oncogene_support,"&") == 0 ~ "Moderate",
        TRUE ~ as.character(NA)
      )
    )
  
    
  tsgs <- as.data.frame(dplyr::bind_rows(
    dplyr::filter(tsg_oncogene_driver$CGC, .data$tsg == T),
    dplyr::filter(tsg_oncogene_driver$NCG, .data$tsg == T)) |>
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
        full_tsg_links = paste(unique(.data$links_tsg),
                                    collapse = ", "),
        full_tsg_support = paste(unique(sort(.data$source)),
                                      collapse = "&"),
        .groups = "drop"
      ) |>
      dplyr::mutate(tsg = dplyr::if_else(
        stringr::str_detect(.data$tsg, "&"),
        TRUE,
        as.logical(.data$tsg)
      ))) |>
    dplyr::full_join(
      dplyr::select(
        tsg_oncogene_driver$CancerMine,
        "entrezgene",
        "cancermine_n_cit_oncogene",
        "cancermine_n_cit_tsg",
        "cancermine_oncogene_tsg_citratio",
        "links_tsg"
      ), by = "entrezgene"
    ) |>
    dplyr::mutate(
      cancermine_n_cit_oncogene = dplyr::if_else(
        is.na(.data$cancermine_n_cit_oncogene),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_oncogene))) |> 
    dplyr::mutate(
      cancermine_n_cit_tsg = dplyr::if_else(
        is.na(.data$cancermine_n_cit_tsg),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_tsg))) |>
    
    dplyr::filter(
      !(is.na(.data$tsg) & 
          .data$cancermine_n_cit_tsg < min_citation_support)) |>
    dplyr::mutate(
      links_tsg = dplyr::if_else(
        !is.na(.data$links_tsg) &
          .data$cancermine_n_cit_tsg < min_citation_support,
        "",
        as.character(.data$links_tsg)
      )
    ) |>
    dplyr::mutate(tsg = dplyr::if_else(
      
      ## text-mined tumor-suppressive role only
      ((is.na(tsg) & 
        cancermine_oncogene_tsg_citratio <= 0.25 &
        cancermine_n_cit_tsg >= min_citation_support_cm_only) |
         
         ## text-mined dual role TSGs/Proto-oncogenes
        (is.na(tsg) & 
           cancermine_oncogene_tsg_citratio <= 2 &
           cancermine_oncogene_tsg_citratio >= 0.5 &
           cancermine_n_cit_oncogene >= min_citation_support_cm_only &
           cancermine_n_cit_tsg >= min_citation_support_cm_only)),
      TRUE,
      as.logical(tsg)
    )) |>
    dplyr::mutate(
      tsg_links = dplyr::case_when(
        is.na(.data$full_tsg_links) &
          tsg == T &
          (!is.na(.data$links_tsg) |
             nchar(.data$links_tsg) > 0) ~ links_tsg,
        !is.na(.data$full_tsg_links) &
          tsg == T &
          nchar(.data$links_tsg) > 0 ~ paste(
            .data$full_tsg_links, .data$links_tsg, sep=", "
          ),
        TRUE ~ as.character(
          .data$full_tsg_links)
      ) 
    ) |>
    dplyr::select(-c("links_tsg","full_tsg_links")) |>
    dplyr::mutate(tsg_support = dplyr::case_when(
      tsg == T & 
        !is.na(full_tsg_support) &
        cancermine_n_cit_tsg >= min_citation_support ~
        paste(full_tsg_support, 
              paste0("CancerMine:",cancermine_n_cit_tsg), 
              sep = "&"),
      tsg == T &
        is.na(full_tsg_support) ~
        paste0("CancerMine:",cancermine_n_cit_tsg),
      TRUE ~ as.character(full_tsg_support))) |>
    dplyr::mutate(tsg_support = stringr::str_replace(
      .data$tsg_support, "^NA&","")
    ) |>
    dplyr::filter(tsg == TRUE)
  
  
  ranked_tsgs <- list()
  ranked_tsgs[['part1']] <- tsgs |> 
    dplyr::filter(
      stringr::str_detect(
        tsg_support, "&")) |> 
    dplyr::arrange(
      desc(stringr::str_count(tsg_support,"&")),
      desc(nchar(tsg_support)),
      desc(cancermine_n_cit_tsg))
  
  ranked_tsgs[['part2']] <- tsgs |>
    dplyr::filter(
      !stringr::str_detect(
        tsg_support, "&") &
        !stringr::str_detect(tsg_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_tsg))
  
  ranked_tsgs[['part3']] <- tsgs |>
    dplyr::filter(
      !stringr::str_detect(
        tsg_support, "&") &
        stringr::str_detect(tsg_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_tsg))
  
  tsgs <- 
    dplyr::bind_rows(
      ranked_tsgs[['part1']], 
      ranked_tsgs[['part2']],
      ranked_tsgs[['part3']]) |>
    dplyr::mutate(ind = dplyr::row_number()) |>
    dplyr::mutate(
      tsg_rank = round(
        100 - (100 * dplyr::percent_rank(ind)), digits = 0)) |>
    dplyr::select(
      -c("full_tsg_support",
         "cancermine_n_cit_tsg",
         "ind",
         "cancermine_n_cit_oncogene",
         "cancermine_oncogene_tsg_citratio")) |>
    dplyr::mutate(
      tsg_confidence_level = dplyr::case_when(
        stringr::str_count(tsg_support,"&") == 2 ~ "Very strong",
        stringr::str_count(tsg_support,"&") == 1 ~ "Strong",
        stringr::str_count(tsg_support,"&") == 0 ~ "Moderate",
        TRUE ~ as.character(NA)
      )
    )
  
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
    )) |>
    dplyr::mutate(oncogene_links = dplyr::if_else(
      is.na(.data$oncogene_links),'', 
      as.character(.data$oncogene_links)
    )) |>
    dplyr::mutate(tsg_links = dplyr::if_else(
      is.na(.data$tsg_links),'', 
      as.character(.data$tsg_links)
    )) |>
    dplyr::mutate(driver_links = dplyr::if_else(
      is.na(.data$driver_links),'', 
      as.character(.data$driver_links)
    )) |>
    dplyr::mutate(cancergene_evidence = dplyr::if_else(
      !is.na(.data$tsg_links) |
        !is.na(.data$oncogene_links) |
        !is.na(.data$driver_links),
      paste0(  
        .data$oncogene_links, ", ",
        .data$tsg_links, ", ",
        .data$driver_links),
      as.character(NA)
    )
    ) |>
    dplyr::mutate(cancergene_evidence = stringr::str_replace_all(
      .data$cancergene_evidence, "^(, ){1,}", ""
    )) |>
    dplyr::mutate(cancergene_evidence = stringr::str_replace_all(
      .data$cancergene_evidence, "(, ){2,}", ", "
    )) |>
    dplyr::mutate(driver_links = dplyr::if_else(
      nchar(.data$driver_links) == 0,
      as.character(NA), 
      as.character(.data$driver_links)
    )) |>
    dplyr::mutate(tsg_links = dplyr::if_else(
      nchar(.data$tsg_links) == 0,
      as.character(NA), 
      as.character(.data$tsg_links)
    )) |>
    dplyr::mutate(oncogene_links = dplyr::if_else(
      nchar(.data$oncogene_links) == 0,
      as.character(NA), 
      as.character(.data$oncogene_links)
    ))
  
  oncogene_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$oncogene == TRUE
    )
  
  tsg_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$tsg == TRUE
    )
  
  driver_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$driver == TRUE
    )
  lgr::lgr$info(
    "-----------------------------------------------------------------"
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
      "A total of n = ", nrow(oncogene_all),
      " classified proto-oncogenes"
    )
  )
  
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(tsg_all),
      " classified tumor suppressor genes"
    )
  )
  
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(driver_all),
      " genes were annotated as potential cancer drivers"
    )
  )
  lgr::lgr$info(
    "-----------------------------------------------------------------"
  )

  
  
  return(tsg_oncogene_driver_evidence)
}

#' Function to classify genes as tumor suppressors/proto-oncogenes,
#' driver etc. based on multiple lines of evidence.
#'
#' @param gox_basic output from geneOncoX::get_basic()
#' @param min_citation_support minimum citation count for support from CancerMine 
#' @param min_citation_support_cm_only minimum citation count for genes annotated 
#' with cancer gene roles from CancerMine only 
#' @param min_sources_driver minimum number of sources (NCG, CanVar-UK, CancerMine,  
#' TCGA, IntoGen) that must contribute to cancer driver gene status
#'
#' @return data frame with cancer gene annotation status
#' @keywords internal
#'
assign_cancer_gene_roles <- function(gox_basic = NULL,
                                     min_citation_support = 5,
                                     min_citation_support_cm_only = 20,
                                     min_sources_driver = 2) {
  tsg_oncogene_driver <- list()
  
  #### CancerMine: Tumor suppressor genes and oncogene annotation
  ## Considering
  ## A) minimum number of supporting citations, and
  ## B) ratio of tsg/oncogene citations (ignore entries in which the support
  ##    is three times as many for the opposite role)
  
  tsg_oncogene_driver[["CancerMine"]] <- gox_basic$records |>
    dplyr::mutate(cancermine_oncogene_tsg_citratio = dplyr::if_else(
      !is.na(.data$cancermine_n_cit_tsg) & .data$cancermine_n_cit_tsg > 0,
      round(as.numeric(.data$cancermine_n_cit_oncogene) /
              .data$cancermine_n_cit_tsg, digits = 1),
      as.numeric(NA)
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
        "cancermine_n_cit_oncogene",
        "cancermine_n_cit_tsg",
        "cancermine_n_cit_driver",
        "cancermine_oncogene_tsg_citratio",
        "links_driver",
        "links_oncogene",
        "links_tsg",
        "source")
    ) |>
    dplyr::filter(
      stringr::str_detect(
        .data$links_oncogene, ": <"
      ) |
        stringr::str_detect(
          .data$links_driver, ": <"
        ) |
        stringr::str_detect(
          .data$links_oncogene, ": <"
        )
      
    ) |>
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
          "<a href='http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "' target='_blank'>",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      ),
      links_tsg = dplyr::if_else(
        .data$tsg == TRUE,
        paste0(
          "Tumor suppressive role (NCG): ",
          "<a href='http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "' target='_blank'>",
          .data$ncg_phenotype, "</a>"
        ),
        ""
      ),
      links_driver = dplyr::if_else(
        .data$driver == TRUE,
        paste0(
          "Cancer driver role (NCG): ",
          "<a href='http://ncg.kcl.ac.uk/query.php?gene_name=",
          .data$entrezgene, "' target='_blank'>",
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
  
  # tsg_oncogene_driver[["CGC"]] <- gox_basic$records |>
  #   dplyr::mutate(oncogene = dplyr::if_else(
  #     .data$cgc_oncogene == TRUE,
  #     TRUE, FALSE
  #   )) |>
  #   dplyr::mutate(tsg = dplyr::if_else(
  #     .data$cgc_tsg == TRUE,
  #     TRUE, FALSE
  #   )) |>
  #   dplyr::mutate(driver = dplyr::if_else(
  #     .data$cgc_driver_tier1 == TRUE |
  #       .data$cgc_driver_tier2 == TRUE,
  #     TRUE, FALSE
  #   )) |>
  #   dplyr::mutate(hallmark = dplyr::if_else(
  #     .data$cgc_hallmark == TRUE,
  #     TRUE, FALSE
  #   )) |>
  #   dplyr::mutate(
  #     links_oncogene = dplyr::if_else(
  #       .data$oncogene == TRUE,
  #       paste0(
  #         "Oncogenic role (CGC): ",
  #         "<a href='https://cancer.sanger.ac.uk/census",
  #         "' target='_blank'>",
  #         "YES (TIER ",
  #         .data$cgc_tier, ")</a>"
  #       ),
  #       ""
  #     ),
  #     links_tsg = dplyr::if_else(
  #       .data$tsg == TRUE,
  #       paste0(
  #         "Tumor suppressive role (CGC): ",
  #         "<a href='https://cancer.sanger.ac.uk/census",
  #         "' target='_blank'>",
  #         "YES (TIER ",
  #         .data$cgc_tier, ")</a>"
  #       ),
  #       ""
  #     ),
  #     links_driver = dplyr::if_else(
  #       .data$driver == TRUE,
  #       paste0(
  #         "Cancer driver role (CGC): ",
  #         "<a href='https://cancer.sanger.ac.uk/census",
  #         "' target='_blank'>",
  #         "YES (TIER ",
  #         .data$cgc_tier, ")</a>"
  #       ),
  #       ""
  #     )
  #   ) |>
  #   dplyr::mutate(source = dplyr::if_else(
  #     .data$cgc_tier == 1,
  #     "CGC_TIER1",
  #     "CGC_TIER2"
  #   )) |>
  #   dplyr::select(
  #     c("entrezgene",
  #       "oncogene",
  #       "tsg",
  #       "driver",
  #       "hallmark",
  #       "links_driver",
  #       "links_oncogene",
  #       "links_tsg",
  #       "source")
  #   ) |>
  #   dplyr::filter(
  #     .data$oncogene == TRUE |
  #       .data$tsg == TRUE |
  #       .data$driver == TRUE) |>
  #   dplyr::distinct()
  
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
          "<a href='https://www.intogen.org/search?gene=",
          .data$symbol,
          "' target='_blank'>",
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
            "<a href='https://gdc.cancer.gov/about-data/",
            "publications/pancan-driver"
          ),
          "' target='_blank'>YES</a>"
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
  
  
  driver_genes <- as.data.frame(
    tsg_oncogene_driver$IntOGen |>
    #driver_genes <- as.data.frame(dplyr::bind_rows(
    #dplyr::filter(tsg_oncogene_driver$CGC, .data$driver == T),
    #tsg_oncogene_driver$IntOGen) |>
    #dplyr::select(-c("tsg","oncogene","links_oncogene",
    #                 "links_tsg","hallmark")) |>
    dplyr::bind_rows(
      dplyr::filter(tsg_oncogene_driver$NCG, .data$driver == T)) |>
    dplyr::select(-c("tsg","oncogene","links_oncogene",
                     "links_tsg")) |>
      dplyr::bind_rows(
        dplyr::filter(
          tsg_oncogene_driver$CancerMine,
          nchar(.data$links_driver) > 0 &
            .data$cancermine_n_cit_driver >= min_citation_support
        )
      ) |>
      dplyr::select(-c("cancermine_n_cit_oncogene",
                       "cancermine_n_cit_tsg",
                       "cancermine_n_cit_driver",
                       "cancermine_oncogene_tsg_citratio",
                       "links_oncogene",
                       "links_tsg")) |>
      dplyr::bind_rows(tsg_oncogene_driver$TCGA) |>
      dplyr::mutate(
        driver = dplyr::if_else(
          is.na(.data$driver),
          TRUE,
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
      dplyr::filter(
        stringr::str_count(.data$driver_support, "&") + 1 >= 
          min_sources_driver
      )
  )
  
  
  oncogenes <- as.data.frame(
    dplyr::filter(tsg_oncogene_driver$NCG, .data$oncogene == T) |>
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
        full_oncogene_links = paste(unique(.data$links_oncogene),
                                    collapse = ", "
        ),
        full_oncogene_support = paste(unique(sort(.data$source)),
                                      collapse = "&"
        ),
        .groups = "drop"
      ) |>
      dplyr::mutate(oncogene = dplyr::if_else(
        stringr::str_detect(.data$oncogene, "&"),
        TRUE,
        as.logical(.data$oncogene)
      ))) |>
    dplyr::full_join(
      dplyr::select(
        tsg_oncogene_driver$CancerMine,
        "entrezgene",
        "cancermine_n_cit_oncogene",
        "cancermine_n_cit_tsg",
        "cancermine_oncogene_tsg_citratio",
        "links_oncogene"
      ), by = "entrezgene"
    ) |>
    dplyr::mutate(
      cancermine_n_cit_oncogene = dplyr::if_else(
        is.na(.data$cancermine_n_cit_oncogene),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_oncogene))) |> 
    dplyr::mutate(
      cancermine_n_cit_tsg = dplyr::if_else(
        is.na(.data$cancermine_n_cit_tsg),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_tsg))) |>
    dplyr::filter(
      !(is.na(.data$oncogene) & 
          .data$cancermine_n_cit_oncogene < min_citation_support)) |>
    dplyr::mutate(
      links_oncogene = dplyr::if_else(
        !is.na(.data$links_oncogene) &
          .data$cancermine_n_cit_oncogene < min_citation_support,
        "",
        as.character(.data$links_oncogene)
      )
    ) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      
      ## Text-mined oncogene support only
      ((is.na(oncogene) & 
          cancermine_oncogene_tsg_citratio >= 4 &
          cancermine_n_cit_oncogene >= min_citation_support_cm_only) |
         
         ## Text-mined dual role TSG/oncogene 
         (is.na(oncogene) & 
            cancermine_oncogene_tsg_citratio <= 2 &
            cancermine_oncogene_tsg_citratio >= 0.5 &
            cancermine_n_cit_oncogene >= min_citation_support_cm_only &
            cancermine_n_cit_tsg >= min_citation_support_cm_only)),
      TRUE,
      as.logical(oncogene)
    )) |>
    dplyr::mutate(
      oncogene_links = dplyr::case_when(
        is.na(full_oncogene_links) &
          oncogene == T &
          (!is.na(.data$links_oncogene) |
             nchar(.data$links_oncogene) > 0) ~ links_oncogene,
        !is.na(full_oncogene_links) &
          oncogene == T &
          !is.na(.data$links_oncogene) &
          nchar(.data$links_oncogene) > 0 ~
          paste(
            full_oncogene_links, links_oncogene, sep=", "
          ),
        TRUE ~ as.character(.data$full_oncogene_links)
      ) 
    ) |>
    dplyr::select(-c("links_oncogene","full_oncogene_links")) |>
    dplyr::mutate(oncogene_support = dplyr::case_when(
      oncogene == T & 
        !is.na(full_oncogene_support) &
        cancermine_n_cit_oncogene >= min_citation_support ~
        paste(full_oncogene_support, 
              paste0("CancerMine:",cancermine_n_cit_oncogene), 
              sep = "&"),
      oncogene == T &
        is.na(full_oncogene_support) ~
        paste0("CancerMine:",cancermine_n_cit_oncogene),
      TRUE ~ as.character(full_oncogene_support))) |>
    dplyr::mutate(oncogene_support = stringr::str_replace(
      .data$oncogene_support, "^NA&","")
    ) |>
    dplyr::filter(oncogene == TRUE)
  
  ranked_oncogenes <- list()
  ranked_oncogenes[['part1']] <- oncogenes |> 
    dplyr::filter(
      stringr::str_detect(
        oncogene_support, "&")) |> 
    dplyr::arrange(
      desc(stringr::str_count(oncogene_support,"&")),
      desc(nchar(oncogene_support)),
      desc(cancermine_n_cit_oncogene))
  
  ranked_oncogenes[['part2']] <- oncogenes |>
    dplyr::filter(
      !stringr::str_detect(
        oncogene_support, "&") &
        !stringr::str_detect(oncogene_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_oncogene))
  
  ranked_oncogenes[['part3']] <- oncogenes |>
    dplyr::filter(
      !stringr::str_detect(
        oncogene_support, "&") &
        stringr::str_detect(oncogene_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_oncogene))
  
  oncogenes <- 
    dplyr::bind_rows(
      ranked_oncogenes[['part1']], 
      ranked_oncogenes[['part2']],
      ranked_oncogenes[['part3']]) |>
    dplyr::mutate(ind = dplyr::row_number()) |>
    dplyr::mutate(
      oncogene_rank = round(
        100 - (100 * dplyr::percent_rank(ind)), digits = 0)) |>
    dplyr::select(
      -c("full_oncogene_support",
         "cancermine_n_cit_tsg",
         "ind",
         "cancermine_n_cit_oncogene",
         "cancermine_oncogene_tsg_citratio")) |>
    dplyr::mutate(
      oncogene_confidence_level = dplyr::case_when(
        stringr::str_count(oncogene_support,"&") == 2 ~ "Very strong",
        stringr::str_count(oncogene_support,"&") == 1 ~ "Strong",
        stringr::str_count(oncogene_support,"&") == 0 ~ "Moderate",
        TRUE ~ as.character(NA)
      )
    )
  
  
  tsgs <- as.data.frame(
    dplyr::filter(tsg_oncogene_driver$NCG, .data$tsg == T) |>
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
        full_tsg_links = paste(unique(.data$links_tsg),
                               collapse = ", "),
        full_tsg_support = paste(unique(sort(.data$source)),
                                 collapse = "&"),
        .groups = "drop"
      ) |>
      dplyr::mutate(tsg = dplyr::if_else(
        stringr::str_detect(.data$tsg, "&"),
        TRUE,
        as.logical(.data$tsg)
      ))) |>
    dplyr::full_join(
      dplyr::select(
        tsg_oncogene_driver$CancerMine,
        "entrezgene",
        "cancermine_n_cit_oncogene",
        "cancermine_n_cit_tsg",
        "cancermine_oncogene_tsg_citratio",
        "links_tsg"
      ), by = "entrezgene"
    ) |>
    dplyr::mutate(
      cancermine_n_cit_oncogene = dplyr::if_else(
        is.na(.data$cancermine_n_cit_oncogene),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_oncogene))) |> 
    dplyr::mutate(
      cancermine_n_cit_tsg = dplyr::if_else(
        is.na(.data$cancermine_n_cit_tsg),
        as.numeric(0),
        as.numeric(.data$cancermine_n_cit_tsg))) |>
    
    dplyr::filter(
      !(is.na(.data$tsg) & 
          .data$cancermine_n_cit_tsg < min_citation_support)) |>
    dplyr::mutate(
      links_tsg = dplyr::if_else(
        !is.na(.data$links_tsg) &
          .data$cancermine_n_cit_tsg < min_citation_support,
        "",
        as.character(.data$links_tsg)
      )
    ) |>
    dplyr::mutate(tsg = dplyr::if_else(
      
      ## text-mined tumor-suppressive role only
      ((is.na(tsg) & 
          cancermine_oncogene_tsg_citratio <= 0.25 &
          cancermine_n_cit_tsg >= min_citation_support_cm_only) |
         
         ## text-mined dual role TSGs/Proto-oncogenes
         (is.na(tsg) & 
            cancermine_oncogene_tsg_citratio <= 2 &
            cancermine_oncogene_tsg_citratio >= 0.5 &
            cancermine_n_cit_oncogene >= min_citation_support_cm_only &
            cancermine_n_cit_tsg >= min_citation_support_cm_only)),
      TRUE,
      as.logical(tsg)
    )) |>
    dplyr::mutate(
      tsg_links = dplyr::case_when(
        is.na(.data$full_tsg_links) &
          tsg == T &
          (!is.na(.data$links_tsg) |
             nchar(.data$links_tsg) > 0) ~ links_tsg,
        !is.na(.data$full_tsg_links) &
          tsg == T &
          nchar(.data$links_tsg) > 0 ~ paste(
            .data$full_tsg_links, .data$links_tsg, sep=", "
          ),
        TRUE ~ as.character(
          .data$full_tsg_links)
      ) 
    ) |>
    dplyr::select(-c("links_tsg","full_tsg_links")) |>
    dplyr::mutate(tsg_support = dplyr::case_when(
      tsg == T & 
        !is.na(full_tsg_support) &
        cancermine_n_cit_tsg >= min_citation_support ~
        paste(full_tsg_support, 
              paste0("CancerMine:",cancermine_n_cit_tsg), 
              sep = "&"),
      tsg == T &
        is.na(full_tsg_support) ~
        paste0("CancerMine:",cancermine_n_cit_tsg),
      TRUE ~ as.character(full_tsg_support))) |>
    dplyr::mutate(tsg_support = stringr::str_replace(
      .data$tsg_support, "^NA&","")
    ) |>
    dplyr::filter(tsg == TRUE)
  
  
  ranked_tsgs <- list()
  ranked_tsgs[['part1']] <- tsgs |> 
    dplyr::filter(
      stringr::str_detect(
        tsg_support, "&")) |> 
    dplyr::arrange(
      desc(stringr::str_count(tsg_support,"&")),
      desc(nchar(tsg_support)),
      desc(cancermine_n_cit_tsg))
  
  ranked_tsgs[['part2']] <- tsgs |>
    dplyr::filter(
      !stringr::str_detect(
        tsg_support, "&") &
        !stringr::str_detect(tsg_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_tsg))
  
  ranked_tsgs[['part3']] <- tsgs |>
    dplyr::filter(
      !stringr::str_detect(
        tsg_support, "&") &
        stringr::str_detect(tsg_support, "CancerMine")) |>
    dplyr::arrange(
      desc(cancermine_n_cit_tsg))
  
  tsgs <- 
    dplyr::bind_rows(
      ranked_tsgs[['part1']], 
      ranked_tsgs[['part2']],
      ranked_tsgs[['part3']]) |>
    dplyr::mutate(ind = dplyr::row_number()) |>
    dplyr::mutate(
      tsg_rank = round(
        100 - (100 * dplyr::percent_rank(ind)), digits = 0)) |>
    dplyr::select(
      -c("full_tsg_support",
         "cancermine_n_cit_tsg",
         "ind",
         "cancermine_n_cit_oncogene",
         "cancermine_oncogene_tsg_citratio")) |>
    dplyr::mutate(
      tsg_confidence_level = dplyr::case_when(
        stringr::str_count(tsg_support,"&") == 2 ~ "Very strong",
        stringr::str_count(tsg_support,"&") == 1 ~ "Strong",
        stringr::str_count(tsg_support,"&") == 0 ~ "Moderate",
        TRUE ~ as.character(NA)
      )
    )
  
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
    )) |>
    dplyr::mutate(oncogene_links = dplyr::if_else(
      is.na(.data$oncogene_links),'', 
      as.character(.data$oncogene_links)
    )) |>
    dplyr::mutate(tsg_links = dplyr::if_else(
      is.na(.data$tsg_links),'', 
      as.character(.data$tsg_links)
    )) |>
    dplyr::mutate(driver_links = dplyr::if_else(
      is.na(.data$driver_links),'', 
      as.character(.data$driver_links)
    )) |>
    dplyr::mutate(cancergene_evidence = dplyr::if_else(
      !is.na(.data$tsg_links) |
        !is.na(.data$oncogene_links) |
        !is.na(.data$driver_links),
      paste0(  
        .data$oncogene_links, ", ",
        .data$tsg_links, ", ",
        .data$driver_links),
      as.character(NA)
    )
    ) |>
    dplyr::mutate(cancergene_evidence = stringr::str_replace_all(
      .data$cancergene_evidence, "^(, ){1,}", ""
    )) |>
    dplyr::mutate(cancergene_evidence = stringr::str_replace_all(
      .data$cancergene_evidence, "(, ){2,}", ", "
    )) |>
    dplyr::mutate(driver_links = dplyr::if_else(
      nchar(.data$driver_links) == 0,
      as.character(NA), 
      as.character(.data$driver_links)
    )) |>
    dplyr::mutate(tsg_links = dplyr::if_else(
      nchar(.data$tsg_links) == 0,
      as.character(NA), 
      as.character(.data$tsg_links)
    )) |>
    dplyr::mutate(oncogene_links = dplyr::if_else(
      nchar(.data$oncogene_links) == 0,
      as.character(NA), 
      as.character(.data$oncogene_links)
    ))
  
  oncogene_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$oncogene == TRUE
    )
  
  tsg_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$tsg == TRUE
    )
  
  driver_all <-
    dplyr::filter(
      tsg_oncogene_driver_evidence,
      .data$driver == TRUE
    )
  lgr::lgr$info(
    "-----------------------------------------------------------------"
  )
  lgr::lgr$info(
    paste0(
      "Sources: Network of Cancer Genes (NCG)"
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
      "A total of n = ", nrow(oncogene_all),
      " classified proto-oncogenes"
    )
  )
  
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(tsg_all),
      " classified tumor suppressor genes"
    )
  )
  
  lgr::lgr$info(
    paste0(
      "A total of n = ", nrow(driver_all),
      " genes were annotated as potential cancer drivers"
    )
  )
  lgr::lgr$info(
    "-----------------------------------------------------------------"
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
