get_intogen_driver_genes <- function(gene_info = NULL) {
  intogen_drivers <- as.data.frame(
    read.table(
      file = file.path(
        "data-raw", "intogen",
        "Compendium_Cancer_Genes.tsv"
      ), sep = "\t",
      header = TRUE, quote = NULL, stringsAsFactors = FALSE
    ) |>
      janitor::clean_names() |>
      dplyr::group_by(symbol) |>
      dplyr::mutate(role = dplyr::if_else(
        role == "", as.character(NA), as.character(role)
      )) |>
      dplyr::mutate(role = dplyr::if_else(
        role == "Act", "Activating", as.character(role)
      )) |>
      dplyr::mutate(role = dplyr::if_else(
        role == "LoF", "Loss_of_Function", as.character(role)
      )) |>
      dplyr::summarise(
        intogen_role = paste(
          unique(Hmisc::capitalize(role)),
          collapse = "&"
        ),
        intogen_phenotype = paste(
          unique(cancer_type), collapse = "&"),
        .groups = "drop"
      ) |>
      dplyr::mutate(intogen_role = dplyr::if_else(
        is.na(intogen_role) | intogen_role == "NA",
        "Unknown",
        as.character(intogen_role)
      )) |>
      dplyr::left_join(
        dplyr::select(gene_info, symbol, entrezgene),
        by = c("symbol"), multiple = "all"
      ) |>
      dplyr::filter(!is.na(entrezgene)) |>
      dplyr::select(-symbol) |>
      dplyr::mutate(intogen_driver = T) |>
      dplyr::distinct()
  )
  lgr::lgr$info(paste0(
    "Parsed n = ", nrow(intogen_drivers),
    " predicted driver genes (Intogen)"
  ))

  return(intogen_drivers)
}


get_curated_fp_cancer_genes <- function(gene_info = NULL) {
  invisible(assertthat::assert_that(
    !is.null(gene_info),
    msg = "'gene_info' object is NULL"
  ))
  assertable::assert_colnames(gene_info, c("symbol", "entrezgene"),
    only_colnames = FALSE, quiet = TRUE
  )
  assertable::assert_coltypes(
    gene_info, list(
      symbol = character(),
      entrezgene = integer()
    ),
    quiet = TRUE
  )
  xlsx_fname <- file.path(
    "data-raw", "tcga_driver_genes", "bailey_2018_cell.xlsx"
  )
  invisible(assertthat::assert_that(
    file.exists(xlsx_fname),
    msg = paste0(
      "File ",
      xlsx_fname, " does not exist"
    )
  ))
  tmp <- openxlsx::read.xlsx(xlsx_fname, sheet = 8, startRow = 4)
  fp_cancer_drivers <- data.frame(
    "symbol" = tmp[, 2],
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(bailey2018_fp_driver = TRUE) |>
    dplyr::filter(!is.na(symbol)) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, entrezgene),
      by = c("symbol"), multiple = "all"
    ) |>
    dplyr::filter(!is.na(entrezgene)) |>
    dplyr::select(-symbol)

  lgr::lgr$info(
    paste0(
      "Parsed n = ", nrow(fp_cancer_drivers),
      " false positive driver genes (TCGA, Bailey et al., Cell, 2018"
    )
  )

  return(fp_cancer_drivers)
}

get_tcga_driver_genes <- function() {
  tcga_projects <- readRDS(
    paste0("data-raw/tcga_projects.rds")
  ) |>
    dplyr::bind_rows(
      data.frame(
        "tumor" = "PANCAN",
        "name" = "Pancancer",
        stringsAsFactors = FALSE
      )
    )

  tcga_drivers <- as.data.frame(
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "tcga_driver_genes",
        "bailey_2018_cell.xlsx"
      ),
      sheet = 2, startRow = 4, cols = c(1, 2, 4, 6, 7)
    ) |>
      janitor::clean_names() |>
      dplyr::mutate(
        tissue_frequency =
          dplyr::if_else(
            is.na(tissue_frequency) |
              tissue_frequency == "NA",
            as.character(pancan_frequency),
            as.character(tissue_frequency)
          )
      ) |>
      dplyr::rename(
        symbol = gene,
        tumor = cancer,
        frequency = tissue_frequency
      ) |>
      dplyr::mutate(frequency = as.numeric(frequency)) |>
      dplyr::left_join(tcga_projects, by = c("tumor"), multiple = "all") |>
      dplyr::select(symbol, tumor, name, frequency) |>
      dplyr::mutate(frequency = as.numeric(frequency)) |>
      dplyr::mutate(name_freq = paste0(
        tumor, ":", round(frequency, 3)
      )) |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(
        tcga_driver_support = paste(
          name_freq,
          collapse = "|"
        ),
        .groups = "drop"
      ) |>
      dplyr::filter(tcga_driver_support != "PANCAN:0")
  ) |>
    dplyr::mutate(tcga_driver = TRUE) |>
    dplyr::mutate(symbol = dplyr::case_when(
      symbol == "CBWD3" ~ "ZNG1C",
      symbol == "FAM46D" ~ "TENT5D",
      symbol == "H3F3A" ~ "H3-3A",
      symbol == "H3F3C" ~ "H3-5",
      symbol == "HIST1H1C" ~ "H1-2",
      symbol == "HIST1H1E" ~ "H1-4",
      symbol == "RQCD1" ~ "CNOT9",
      symbol == "TCEB1" ~ "ELOC",
      symbol == "WHSC1" ~ "NSD2",
      TRUE ~ as.character(symbol)
    ))

  lgr::lgr$info(paste0(
    "Parsed n = ", nrow(tcga_drivers),
    " predicted driver genes (TCGA, Bailey et al., Cell, 2018)"
  ))

  return(tcga_drivers)
}

get_signaling_pathway_genes <- function(gene_info) {
  signaling_genes_tcga_cell <-
    file.path(
      "data-raw",
      "oncogenic_signaling",
      "tcga_pancan_2018_sanchez_vega.xlsx"
    )

  cell_cycle_genes <-
    openxlsx::read.xlsx(signaling_genes_tcga_cell,
      sheet = 1, startRow = 3
    ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "Cell Cycle", signaling_pathway = "CC")

  hippo_genes <- openxlsx::read.xlsx(
    signaling_genes_tcga_cell,
    sheet = 2, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(
      pathway = "Hippo Signaling Cascade",
      signaling_pathway = "HIPPO"
    )

  myc_genes <- openxlsx::read.xlsx(
    signaling_genes_tcga_cell,
    sheet = 3, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "MYC", signaling_pathway = "MYC")

  notch_genes <- openxlsx::read.xlsx(
    signaling_genes_tcga_cell,
    sheet = 4, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "NOV", "CCN3", as.character(symbol)
    )) |>
    dplyr::mutate(
      pathway = "Notch Signaling Pathway",
      signaling_pathway = "NOTCH"
    )

  nrf2_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
    sheet = 5, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(
      pathway = "Oxidative Stress Response",
      signaling_pathway = "NRF2"
    )

  pi3k_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
    sheet = 6, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(
      pathway = "PI-3-Kinase signaling",
      signaling_pathway = "PI3K"
    )

  rtkras_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
    sheet = 8, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(
      pathway =
        "Receptor-Tyrosine Kinase (RTK)/RAS/MAP-Kinase Signaling",
      signaling_pathway = "RTKRAS"
    )

  tgfbeta_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
    sheet = 7, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(
      pathway = "TGFβ Signaling",
      signaling_pathway = "TFGBETA"
    )

  tp53_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
    sheet = 9, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "TP53 Signaling", signaling_pathway = "TP53")

  wnt_genes <- openxlsx::read.xlsx(
    signaling_genes_tcga_cell,
    sheet = 10, startRow = 1
  ) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(
      pathway = "β-catenin/Wnt Signaling", signaling_pathway = "WNT"
    )

  signaling_genes <- as.data.frame(
    dplyr::bind_rows(
      cell_cycle_genes, hippo_genes,
      nrf2_genes, myc_genes,
      wnt_genes, tp53_genes,
      rtkras_genes, tgfbeta_genes,
      pi3k_genes, notch_genes
    ) |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(
        sanchezvega2018_signaling_pathway_short =
          paste(unique(signaling_pathway), collapse = "&"),
        sanchezvega2018_signaling_pathway =
          stringi::stri_enc_toascii(
            paste(unique(pathway),
              collapse = "&"
            )
          ), .groups = "drop"
      )
  ) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, entrezgene),
      by = "symbol", multiple = "all"
    ) |>
    dplyr::select(-c(symbol, sanchezvega2018_signaling_pathway_short))

  lgr::lgr$info(
    paste0(
      "Parsed n = ", nrow(signaling_genes),
      " genes with pathway signalling annotation ",
      "(TCGA, Sanchez-Vega et al., Cell, 2018"
    )
  )

  return(signaling_genes)
}

get_cancer_gene_census <- function(
    origin = "somatic",
    opentargets_version = "2024.09",
    cgc_version = "v101") {
  
  cosmic_genes <- readr::read_tsv(
    file = file.path(
      "data-raw",
      "cancer_gene_census",
      paste0("Cosmic_Genes_",
             cgc_version,"_GRCh38.tsv.gz")
    ),
    show_col_types = FALSE) |>
    janitor::clean_names() |>
    dplyr::select(
      cosmic_gene_id,
      entrez_id
    ) |>
    dplyr::rename(
      entrez_gene_id = entrez_id
    ) |>
    dplyr::distinct()
  
  cosmic_cgc <- readr::read_tsv(
    file = file.path(
      "data-raw", 
      "cancer_gene_census",
      paste0("Cosmic_CancerGeneCensus_",
             cgc_version,"_GRCh38.tsv.gz")
    ),
    show_col_types = FALSE) |>
    janitor::clean_names() |>
    dplyr::left_join(
      cosmic_genes,
      by = "cosmic_gene_id") |>
    dplyr::select(
      -c("chromosome","genome_start","genome_stop",
         "chr_band","synonyms")
    )

  cgc_hallmark_genes <- readRDS(
    file = file.path(
      "data-raw", "opentargets",
      paste0(
        "opentargets_target_",
        opentargets_version,
        ".rds"
      )
    )
  ) |>
    dplyr::select(
      target_symbol,
      cancer_hallmark
    ) |>
    dplyr::filter(
      !is.na(cancer_hallmark)
    ) |>
    dplyr::rename(
      symbol = target_symbol
    ) |>
    dplyr::select(symbol) |>
    dplyr::distinct() |>
    dplyr::mutate(cgc_hallmark = TRUE)

  if (origin == "somatic") {
    cosmic_cgc <- cosmic_cgc |>
      dplyr::filter(somatic == "y")
  } else {
    if (origin == "germline") {
      cosmic_cgc <- cosmic_cgc |>
        dplyr::filter(germline == "y")
    }
  }

  cosmic_cgc <- cosmic_cgc |>
    dplyr::select(
      gene_symbol,
      entrez_gene_id,
      germline,
      somatic,
      tier,
      role_in_cancer,
      molecular_genetics,
      cancer_syndrome,
      tumour_types_somatic,
      tumour_types_germline
    ) |>
    dplyr::rename(
      symbol = gene_symbol,
      entrezgene = entrez_gene_id,
      cgc_tier = tier,
      moi = molecular_genetics
    ) |> 
    dplyr::left_join(
      cgc_hallmark_genes, by = "symbol", 
      relationship = "many-to-many") |>
    dplyr::mutate(cgc_hallmark = dplyr::if_else(
      is.na(cgc_hallmark),
      FALSE,
      as.logical(cgc_hallmark)
    )) |>
    dplyr::mutate(cgc_oncogene = dplyr::if_else(
      !is.na(role_in_cancer) &
        stringr::str_detect(role_in_cancer, "oncogene"),
      TRUE, FALSE
    )) |>
    dplyr::mutate(cgc_tsg = dplyr::if_else(
      !is.na(role_in_cancer) &
        stringr::str_detect(role_in_cancer, "TSG"),
      TRUE, FALSE
    )) |>
    dplyr::mutate(moi = stringr::str_trim(moi)) |>
    dplyr::mutate(
      moi = dplyr::if_else(
        stringr::str_detect(
          moi, "^Rec$"
        ), "AR", as.character(moi)
      )
    ) |>
    dplyr::mutate(moi = dplyr::if_else(
      stringr::str_detect(
        moi, "^Dom$"
      ), "AD", as.character(moi)
    )) |>
    dplyr::mutate(moi = dplyr::if_else(
      stringr::str_detect(
        moi, "^Dom/Rec$"
      ), "AD/AR", as.character(moi)
    )) |>
    dplyr::mutate(cgc_moi = dplyr::if_else(
      stringr::str_detect(
        moi, "Rec/X"
      ), as.character(NA), as.character(moi)
    )) |>
    dplyr::select(
      symbol, cgc_moi,
      entrezgene,
      cgc_tsg,
      cgc_oncogene,
      cgc_hallmark,
      cgc_tier,
      tumour_types_somatic,
      cancer_syndrome,
      tumour_types_germline
    ) |>
    
    ## bug in cancer gene census for ERCC5
    dplyr::mutate(entrezgene = dplyr::case_when(
      symbol == "MDS2" ~ as.integer(259283),
      symbol == "MALAT1" ~ as.integer(378938),
      symbol == "HMGN2P46" ~ as.integer(283651),
      TRUE ~ as.integer(entrezgene)
    )) |>
    dplyr::select(-symbol)

  if (origin == "somatic") {
    cosmic_cgc <- cosmic_cgc |>
      dplyr::mutate(
        cgc_somatic = TRUE
      ) |>
      dplyr::mutate(
        cgc_phenotype_somatic = tumour_types_somatic
      ) |>
      dplyr::select(
        -c(
          cancer_syndrome,
          tumour_types_germline,
          tumour_types_somatic
        )
      )
  }
  if (origin == "germline") {
    cosmic_cgc <- cosmic_cgc |>
      dplyr::mutate(cgc_germline = TRUE) |>
      dplyr::mutate(cgc_phenotype_germline = paste(
        tumour_types_germline,
        cancer_syndrome,
        sep = "|"
      )) |>
      dplyr::select(-c(
        cancer_syndrome,
        tumour_types_germline,
        tumour_types_somatic
      ))
  }
  if (origin == "all") {
    cosmic_cgc <- cosmic_cgc |>
      dplyr::select(entrezgene, cgc_hallmark, cgc_tier) |>
      dplyr::distinct() |>
      dplyr::mutate(cgc_driver_tier1 = dplyr::if_else(
        cgc_tier == 1,
        TRUE, FALSE
      )) |>
      dplyr::mutate(cgc_driver_tier2 = dplyr::if_else(
        cgc_tier == 2,
        TRUE, FALSE
      ))
  }

  lgr::lgr$info(paste0(
    "Parsed n = ", nrow(cosmic_cgc),
    " genes in COSMIC's Cancer Gene Census (version ",
    cgc_version, ") - (germline/somatic: ", origin, ")"
  ))

  return(cosmic_cgc)
}

get_network_of_cancer_genes <- function(ncg_version = "7.2") {
  ncg <- read.table(
    file = file.path(
      "data-raw", "ncg", "ncg.tsv"
    ), header = TRUE,
    stringsAsFactors = FALSE, sep = "\t",
    quote = "", comment.char = "") |>
    janitor::clean_names() |> 
    
    # dplyr::select(
    #   ncg_tsg, entrez, vogelstein_annotation,
    #   cgc_annotation, saito_annotation, cancer_type
    # ) |>
    dplyr::distinct() |>
    dplyr::mutate(ncg_driver = dplyr::if_else(
      !is.na(type) & type == "\"Canonical Cancer Driver\"",
      as.logical(TRUE), as.logical(FALSE)
    )) |>
    dplyr::mutate(cancer_type = stringr::str_replace_all(
      cancer_type, "\\\"", ""
    )) |>
    dplyr::mutate(
      ncg_tsg = dplyr::if_else(
        stringr::str_detect(cgc_annotation, "TSG") |
          ncg_tsg == 1,
        as.logical(TRUE), as.logical(FALSE)
      )
    ) |>
    dplyr::mutate(
      ncg_oncogene = dplyr::if_else(
        stringr::str_detect(cgc_annotation, "oncogene") |
          ncg_oncogene == 1,
        as.logical(TRUE), as.logical(FALSE)
      )
    ) |>
    dplyr::rename(entrezgene = entrez) |>
    dplyr::filter(!(ncg_tsg == FALSE &
      ncg_oncogene == FALSE &
      ncg_driver == FALSE)) |>
    dplyr::distinct() |>
    dplyr::group_by(entrezgene) |>
    dplyr::summarise(
      ncg_phenotype = paste(
        sort(unique(cancer_type)),
        collapse = ","
      ),
      ncg_pmid = paste(
        sort(unique(pubmed_id)),
        collapse = ";"
      ),
      ncg_tsg = paste(
        unique(ncg_tsg),
        collapse = ","
      ),
      ncg_oncogene = paste(
        unique(ncg_oncogene),
        collapse = ","
      ),
      ncg_driver = paste(
        unique(ncg_driver),
        collapse = ","
      )
    ) |>
    dplyr::mutate(ncg_driver = dplyr::if_else(
      stringr::str_detect(ncg_driver, ","),
      as.logical(TRUE),
      as.logical(ncg_driver)
    )) |>
    dplyr::mutate(
      ncg_tsg = as.logical(ncg_tsg),
      ncg_oncogene = as.logical(ncg_oncogene)
    ) |>
    dplyr::select(
      entrezgene, ncg_driver, ncg_tsg,
      ncg_oncogene, ncg_pmid, ncg_phenotype
    ) |>
    dplyr::mutate(ncg_phenotype = stringr::str_replace_all(
      ncg_phenotype, ",,", ","
    )) |>
    dplyr::mutate(ncg_phenotype = stringr::str_replace_all(
      ncg_phenotype, "^,", ""
    )) |>
    dplyr::distinct()


  ncg_phenotype_cleaned <- as.data.frame(
    ncg |>
      dplyr::select(entrezgene, ncg_phenotype) |>
      tidyr::separate_rows(ncg_phenotype, sep = ",") |>
      dplyr::mutate(ncg_phenotype = stringr::str_replace_all(
        stringr::str_to_title(ncg_phenotype), "_", " "
      )) |>
      dplyr::mutate(ncg_phenotype = stringr::str_trim(
        stringr::str_replace(
          ncg_phenotype, "b-Cell", "B-Cell"
        )
      )) |>
      dplyr::mutate(ncg_phenotype = stringr::str_replace(
        ncg_phenotype, "t-Cell|t cell", "T-Cell"
      )) |>
      dplyr::mutate(ncg_phenotype = stringr::str_replace(
        ncg_phenotype, "dlbcl|Dlblc", "DLBCL"
      )) |>
      dplyr::mutate(ncg_phenotype = stringr::str_replace(
        ncg_phenotype, "gianT-Cell", "giant cell"
      )) |>
      dplyr::mutate(ncg_phenotype = stringr::str_replace(
        ncg_phenotype, "high-Grade", "high grade"
      )) |>
      dplyr::mutate(ncg_phenotype = stringr::str_replace(
        ncg_phenotype, "Pan-Cancer", "Pan-cancer"
      )) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(ncg_phenotype = paste(
        unique(ncg_phenotype),
        collapse = ", "
      )) |>
      dplyr::mutate(ncg_phenotype = dplyr::if_else(
        is.na(ncg_phenotype) |
          nchar(ncg_phenotype) == 0,
        "Undefined tumor type(s)",
        as.character(ncg_phenotype)
      ))
  )

  ncg$ncg_phenotype <- NULL
  ncg <- ncg |>
    dplyr::left_join(
      ncg_phenotype_cleaned,
      by = "entrezgene", multiple = "all"
    )

  lgr::lgr$info(
    paste0(
      "Found n = ", nrow(ncg),
      " cancer-relvant genes (drivers, proto-oncogenes, tumor ",
      "suppressors) in Network of Cancer Genes (version ",
      ncg_version, ")"
    )
  )
  return(ncg)
}

get_cancermine_genes <- function(cancermine_version = "51",
                                 num_links_per_role = 15) {
  
  sjr_impact_fname <- paste0(
    "data-raw/cancermine/scimagojr_2024.csv"
  )
  
  cancermine_sentences_fname <-
    paste0(
      "data-raw/cancermine/cancermine_sentences.v",
      cancermine_version, ".tsv.gz"
    )
  cancermine_citations_fname <-
    paste0(
      "data-raw/cancermine/cancermine_citations.v",
      cancermine_version, ".tsv.gz"
    )
  cancermine_collated_fname <-
    paste0(
      "data-raw/cancermine/cancermine_collated.v",
      cancermine_version, ".tsv.gz"
    )

  journal_impact <- 
    suppressMessages(readr::read_delim(
      file=sjr_impact_fname, 
      show_col_types = F, 
      guess_max = 100000, delim = ";")) |> 
    dplyr::select(Title, SJR) |> 
    dplyr::mutate(SJR = stringr::str_replace(SJR,",",".")) |> 
    dplyr::rename(journal = "Title") |> 
    dplyr::mutate(journal_impact = as.numeric(SJR)) |>
    dplyr::mutate(journal_lc = tolower(journal)) |>
    dplyr::select(journal_impact, journal_lc) |>
    dplyr::filter(!is.na(journal_impact)) |>
    dplyr::group_by(journal_lc) |>
    dplyr::filter(dplyr::n() == 1) |>
    dplyr::ungroup()
  
  pmids <- as.data.frame(
    read.table(
      file = gzfile(cancermine_sentences_fname),
      header = TRUE, comment.char = "", quote = "",
      sep = "\t", stringsAsFactors = FALSE) |>
      #dplyr::filter(predictprob >= 0.8) |>
      dplyr::filter(
        subsection %in% 
          c("background", "conclusion", "introduction",
            "discussion", "results","conclusions",
            "None","summary")
      ) |>
      ## some entries wrongly captured, are in
      ## fact mentions of anti-sense non-coding genes
      dplyr::filter(
        !stringr::str_detect(
          formatted_sentence, "-<b>AS1|b>-AS1"
        )
      ) |>
      dplyr::mutate(journal = stringr::str_replace(
        .data$journal, "&", "and"
      )) |>
      dplyr::mutate(
        journal2 = 
          stringr::str_replace(
            .data$journal, " \\(.+\\)","")
      ) |> 
      dplyr::mutate(journal2 = stringr::str_replace(
        .data$journal2, "( )?: .+$","")) |>
      dplyr::mutate(journal2 = stringr::str_replace(
        .data$journal2, "\\. "," "
      )) |>
      dplyr::mutate(journal2 = stringr::str_replace(
        .data$journal2, "^The ",""
      )) |>
      dplyr::mutate(journal2 = stringr::str_replace(
        .data$journal2, "Genes, chromosomes",
        "Genes chromosomes"
      )) |>
      dplyr::mutate(journal_lc = tolower(journal2)) |>
      dplyr::select(-c("journal2")) |>
      dplyr::left_join(
        journal_impact, by = "journal_lc", 
        relationship = "many-to-one"
      ) |>
      dplyr::mutate(pmid = as.character(pmid)) |>
      dplyr::group_by(role, gene_entrez_id, pmid, 
                      journal, journal_impact, year) |>
      dplyr::summarise(doid = paste(
        unique(cancer_id),
        collapse = ","
      ), .groups = "drop") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::distinct()
  )

  all_citations <- data.frame()
  if (file.exists(cancermine_citations_fname)) {
    all_citations <- read.table(
      file = gzfile(cancermine_citations_fname),
      sep = "\t", header = FALSE, quote = "",
      comment.char = "", stringsAsFactors = FALSE) |>
      magrittr::set_colnames(c("pmid", "citation", "citation_link")) |>
      dplyr::mutate(pmid = as.character(pmid)) |>
      # dplyr::mutate(
      #   citation = stringi::stri_enc_toascii(citation)) |>
      # dplyr::mutate(
      #   citation_link = stringi::stri_enc_toascii(citation_link)) |>
      dplyr::distinct() |>
      dplyr::arrange(desc(pmid))
  }

  pmids <- pmids |>
    dplyr::inner_join(
      all_citations, by = c("pmid"), 
      relationship = "many-to-many") |>
    dplyr::mutate(
      pmid = as.integer(pmid)
    )

  
  citation_links_oncogene <- 
    as.data.frame(
      pmids |>
        dplyr::filter(role == "Oncogene") |>
        dplyr::arrange(entrezgene, 
                       dplyr::desc(journal_impact)) |>
        dplyr::group_by(entrezgene) |>
        dplyr::slice_head(n = num_links_per_role) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          entrezgene, 
          dplyr::desc(pmid)
        ) |>
        dplyr::group_by(entrezgene) |>
        dplyr::summarise(
          citations_oncogene = paste(
            citation,
            collapse = "; "
          ),
          citation_links_oncogene = paste(
            unique(citation_link),
            collapse = ", "
          ), .groups = "drop"
        )
    )

  pmids_oncogene <- as.data.frame(
    pmids |>
      dplyr::filter(role == "Oncogene") |>
      dplyr::arrange(entrezgene, 
                     dplyr::desc(pmid)) |> 
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        pmids_oncogene = paste(pmid, collapse = ";"),
        #citations_oncogene = paste(
        #  citation,
        #  collapse = "; "
        #),
        .groups = "drop"
      ) |>
      dplyr::mutate(n_citations_oncogene = as.integer(
        stringr::str_count(pmids_oncogene, ";")
      ) + 1) |>
      dplyr::filter(n_citations_oncogene > 1) |>
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(
        n_citations_oncogene >= 6,
        "MC",
        as.character("LC")
      )) |>
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(
        n_citations_oncogene >= 15,
        "HC",
        as.character(oncogene_cancermine)
      )) |>
      dplyr::left_join(
        citation_links_oncogene, by = "entrezgene", 
        relationship = "one-to-one"
      )
  )

  citation_links_tsgene <- 
    as.data.frame(
      pmids |>
        dplyr::filter(role == "Tumor_Suppressor") |>
        dplyr::arrange(entrezgene, 
                       dplyr::desc(journal_impact)) |>
        dplyr::group_by(entrezgene) |>
        dplyr::slice_head(n = num_links_per_role) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          entrezgene, 
          dplyr::desc(pmid)
        ) |>
        dplyr::group_by(entrezgene) |>
        dplyr::summarise(
          citations_tsgene = paste(
            citation,
            collapse = "; "
          ),
          citation_links_tsgene = paste(
            unique(citation_link),
            collapse = ", "
          ), .groups = "drop"
        )
    )
  
  pmids_tsgene <- as.data.frame(
    pmids |>
      dplyr::filter(role == "Tumor_Suppressor") |>
      dplyr::arrange(entrezgene, desc(pmid)) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        pmids_tsgene = paste(unique(pmid), collapse = ";"),
        #citations_tsgene = paste(
        #  citation,
        #  collapse = "; "
        #),
        .groups = "drop"
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(n_citations_tsgene = as.integer(
        stringr::str_count(pmids_tsgene, ";")
      ) + 1) |>
      dplyr::filter(n_citations_tsgene > 1) |>
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(
        n_citations_tsgene >= 6,
        "MC",
        as.character("LC")
      )) |>
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(
        n_citations_tsgene >= 15,
        "HC",
        as.character(tumor_suppressor_cancermine)
      )) |>
      dplyr::left_join(
        citation_links_tsgene, by = "entrezgene",
        relationship = "many-to-one"
      )
  )
  
  
  citation_links_cdriver <- 
    as.data.frame(
      pmids |>
        dplyr::filter(role == "Driver") |>
        dplyr::arrange(entrezgene, 
                       dplyr::desc(journal_impact)) |>
        dplyr::group_by(entrezgene) |>
        dplyr::slice_head(n = num_links_per_role) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          entrezgene, 
          dplyr::desc(pmid)
        ) |>
        dplyr::group_by(entrezgene) |>
        dplyr::summarise(
          citations_cdriver = paste(
            citation, collapse = "; "),
          citation_links_cdriver = paste(
            unique(citation_link),
            collapse = ", "
          ), .groups = "drop"
        )
    )

  pmids_cdriver <- as.data.frame(
    pmids |>
      dplyr::filter(role == "Driver") |>
      dplyr::arrange(entrezgene, desc(pmid)) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        pmids_cdriver = paste(pmid, collapse = ";"),
        #citations_cdriver = paste(citation, collapse = "; "),
        .groups = "drop"
      ) |>
      dplyr::mutate(n_citations_cdriver = as.integer(
        stringr::str_count(pmids_cdriver, ";")
      ) + 1) |>
      dplyr::filter(n_citations_cdriver > 1) |>
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(
        n_citations_cdriver >= 6,
        "MC",
        as.character("LC")
      )) |>
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(
        n_citations_cdriver >= 10,
        "HC",
        as.character(cancer_driver_cancermine)
      )) |>
      dplyr::left_join(
        citation_links_cdriver, by = "entrezgene",
        relationship = "many-to-one"
      )
  )


  lgr::lgr$info(
    paste0(
      "Retrieving known proto-oncogenes/tumor suppressor",
      "genes from CancerMine"
    )
  )
  oncogene <- as.data.frame(
    readr::read_tsv(cancermine_collated_fname,
      col_names = TRUE, na = "-",
      comment = "#", quote = "",
      show_col_types = FALSE
    ) |>
      dplyr::filter(role == "Oncogene") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        doid_oncogene = paste(
          unique(cancer_id), collapse = ","),
        .groups = "drop"
      ) |>
      dplyr::inner_join(
        pmids_oncogene, 
        by = "entrezgene", 
        multiple = "all") |>
      dplyr::distinct()
  )

  n_hc_oncogene <- oncogene |>
    dplyr::filter(oncogene_cancermine == "HC") |>
    nrow()


  tsgene <- as.data.frame(
    readr::read_tsv(cancermine_collated_fname,
      col_names = TRUE, na = "-",
      comment = "#", quote = "",
      show_col_types = FALSE
    ) |>
      dplyr::filter(role == "Tumor_Suppressor") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(doid_tsgene = paste(
        unique(cancer_id),
        collapse = ","
      ), .groups = "drop") |>
      dplyr::inner_join(
        pmids_tsgene, 
        by = "entrezgene", 
        multiple = "all") |>
      dplyr::distinct()
  )
  n_hc_tsgene <- tsgene |>
    dplyr::filter(tumor_suppressor_cancermine == "HC") |>
    nrow()


  cdriver <- as.data.frame(
    readr::read_tsv(cancermine_collated_fname,
      col_names = TRUE, na = "-",
      comment = "#", quote = "",
      show_col_types = FALSE
    ) |>
      dplyr::filter(role == "Driver") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(doid_cdriver = paste(
        unique(cancer_id),
        collapse = ","
      ), .groups = "drop") |>
      dplyr::inner_join(
        pmids_cdriver, 
        by = "entrezgene", 
        multiple = "all") |>
      dplyr::distinct()
  )
  n_hc_cdriver <- cdriver |>
    dplyr::filter(cancer_driver_cancermine == "HC") |>
    nrow()

  cancermine_full <- cdriver |>
    dplyr::full_join(tsgene, by = "entrezgene", multiple = "all") |>
    dplyr::full_join(oncogene, by = "entrezgene", multiple = "all") |>
    dplyr::rename(
      cancermine_pmid_driver = pmids_cdriver,
      cancermine_pmid_tsg = pmids_tsgene,
      cancermine_pmid_oncogene = pmids_oncogene,
      cancermine_doid_driver = doid_cdriver,
      cancermine_doid_tsg = doid_tsgene,
      cancermine_doid_oncogene = doid_oncogene,
      cancermine_n_cit_driver = n_citations_cdriver,
      cancermine_n_cit_tsg = n_citations_tsgene,
      cancermine_n_cit_oncogene = n_citations_oncogene,
      cancermine_cit_tsg = citations_tsgene,
      cancermine_cit_oncogene = citations_oncogene,
      cancermine_cit_driver = citations_cdriver,
      cancermine_cit_links_tsg = citation_links_tsgene,
      cancermine_cit_links_oncogene = citation_links_oncogene,
      cancermine_cit_links_driver = citation_links_cdriver
    ) |>
    dplyr::select(
      entrezgene,
      cancermine_pmid_driver,
      cancermine_pmid_oncogene,
      cancermine_pmid_tsg,
      cancermine_doid_driver,
      cancermine_doid_oncogene,
      cancermine_doid_tsg,
      cancermine_n_cit_driver,
      cancermine_n_cit_oncogene,
      cancermine_n_cit_tsg,
      cancermine_cit_tsg,
      cancermine_cit_oncogene,
      cancermine_cit_driver,
      cancermine_cit_links_driver,
      cancermine_cit_links_oncogene,
      cancermine_cit_links_tsg
    ) |>
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
    ))


  return(cancermine_full)
}
