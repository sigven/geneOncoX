#' A function that splits an array into chunks of equal size
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

#' A function that returns a citation with first author, journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @param raw_db_dir base directory
#' @param chunk_size Size of PMID chunks
#'
#' @export
#'
#' @return citation PubMed citation, with first author, journal and year
#'
get_citations_pubmed <- function(
    pmid,
    raw_db_dir = NULL,
    chunk_size = 100){
  
  ## set logging layout
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  
  literature_vault_fname <- file.path(
    raw_db_dir,
    "literature_vault",
    "literature_vault.rds")
  
  pmid_df <- data.frame('pmid' = pmid)
  all_citations <- data.frame()
  literature_vault <- data.frame()
  
  if(file.exists(literature_vault_fname)){
    lgr::lgr$info("Checking cached literature vault")
    literature_vault <- readRDS(literature_vault_fname) |>
      dplyr::filter(!is.na(.data$pmid)) |>
      dplyr::distinct()
    
    missing_pmids <- pmid_df |>
      dplyr::anti_join(literature_vault, by = "pmid")
    
    all_citations <- literature_vault |>
      dplyr::semi_join(pmid_df, by = "pmid") |>
      dplyr::distinct()
    
    lgr::lgr$info(
      paste0('Found ', nrow(all_citations),
             " PMIDs in cached literature vault - missing ", nrow(missing_pmids)))
    
    pmid <- missing_pmids$pmid
  }
  
  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(
    pmid, ceiling(length(pmid)/chunk_size))
  
  
  if(length(pmid_chunks) == 0){
    lgr::lgr$info("No PMIDs to process")
    return(all_citations)
  }
  
  j <- 0
  
  lgr::lgr$info(
    paste0('Retrieving PubMed citations for PMID list, total length: ',
           length(pmid)))
  while (j < length(pmid_chunks)) {
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    lgr::lgr$info(
      paste0('Processing chunk ',j,' with ',length(pmid_chunk),' PMIDS'))
    pmid_string <- paste(pmid_chunk,collapse = " ")
    res <- RISmed::EUtilsGet(
      RISmed::EUtilsSummary(
        pmid_string, type = "esearch", db = "pubmed", retmax = 5000)
    )
    year <- RISmed::YearPubmed(res)
    authorlist <- RISmed::Author(res)
    pmid_list <- RISmed::PMID(res)
    i <- 1
    first_author <- c()
    while (i <= length(authorlist)) {
      if (length(authorlist[[i]]) == 5) {
        first_author <- c(
          first_author,
          paste(authorlist[[i]][1,]$LastName," et al.",sep = ""))
      } else{
        first_author <- c(
          first_author, as.character("Unknown et al.")
        )
      }
      i <- i + 1
    }
    journal <- RISmed::ISOAbbreviation(res)
    citations <- data.frame(
      'pmid' = as.integer(pmid_list),
      'citation' = paste(
        first_author, year, journal, sep = ", "),
      stringsAsFactors = F)
    citations$citation_link <- paste0(
      '<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/',
      citations$pmid,'\' target=\'_blank\'>',
      citations$citation,'</a>')
    all_citations <- dplyr::bind_rows(
      all_citations, citations) |>
      dplyr::filter(!is.na(pmid)) |>
      dplyr::distinct()
    
    literature_vault <- dplyr::bind_rows(
      literature_vault, citations
    )
    saveRDS(literature_vault, file = literature_vault_fname)
    j <- j + 1
  }
  
  return(all_citations)
  
}


get_cpic_genes <- function(update = T,
                           gene_info = NULL){
  datestamp <- Sys.Date()
  remote_url <-
    paste0(
      "https://files.cpicpgx.org/data/report/current/",
      "pair/cpic_gene-drug_pairs.xlsx"
    )
  cpic_genes_oncology <- data.frame()
  if (RCurl::url.exists(remote_url)) {
    cpic_genes <- openxlsx2::read_xlsx(file=remote_url, col_names = T)
    colnames(cpic_genes) <- 
      c('symbol','drug_name','drug_rxnorm_id',
        'atc_id','guideline','cpic_level',
        'cpic_level_status',
        'pharmgkb_evidence_level', 
        'pgx_on_fda_label','pmids')
    cpic_genes_oncology <- cpic_genes |>
      dplyr::filter(
        stringr::str_detect(atc_id,"^L|, L")
      ) |>
      dplyr::mutate(cpic_entry = paste(
        cpic_level,
        stringr::str_replace_all(atc_id, ", ", "&"),
        sep=":")) |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(
        cpic_pgx_oncology = paste(sort(cpic_entry), collapse=";"), 
        .groups="drop")
  }
  cpic_genes_oncology <- cpic_genes_oncology |>
    dplyr::inner_join(
      dplyr::select(gene_info, symbol, entrezgene), 
      multiple = "all", by = "symbol"
    ) |>
    dplyr::select(-symbol) |>
    dplyr::distinct()
  
  return(cpic_genes_oncology)
  
}

get_gene_info_ncbi <- function(update = T) {
  datestamp <- Sys.Date()
  remote_url <-
    paste0(
      "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/",
      "Mammalia/Homo_sapiens.gene_info.gz"
    )
  if (RCurl::url.exists(remote_url)) {
    gene_info <- suppressWarnings(readr::read_tsv(
      remote_url,
      show_col_types = FALSE, skip = 1,
      comment = "#", quote = "",
      col_select = c(1, 2, 3, 5, 6, 9, 10, 14),
      progress = FALSE, col_names = FALSE
    ))
  } else {
    lgr::lgr$error(
      paste0(
        "Cannot read gene_info - download not available: ",
        remote_url
      )
    )
  }

  gene_info <- gene_info |>
    dplyr::filter(X1 == 9606) |>
    dplyr::rename(
      entrezgene = X2,
      synonyms = X5,
      symbol = X3,
      name = X9,
      gene_biotype = X10,
      other_genename_designations = X14
    ) |>
    dplyr::mutate(
      ensembl_gene_id = stringr::str_replace(
        stringr::str_match(X6, "Ensembl:ENSG[0-9]{1,}"), "Ensembl:", ""
      )
    ) |>
    dplyr::mutate(hgnc_id = stringr::str_replace(
      stringr::str_match(X6, "HGNC:HGNC:[0-9]{1,}"), "HGNC:HGNC:", ""
    )) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::mutate(gene_biotype = dplyr::if_else(
      gene_biotype == "protin-coding",
      "protein_coding", as.character(gene_biotype)
    )) |>
    dplyr::mutate(
      gene_biotype = stringr::str_replace(
        gene_biotype, "-","_"
      )) |>
    dplyr::mutate(gene_biotype = dplyr::case_when(
      gene_biotype == "ncRNA" ~ "lncRNA",
      gene_biotype == "pseudo" ~ "pseudogene",
      TRUE ~ as.character(gene_biotype)
    )) |>
    dplyr::filter(
      gene_biotype != "biological_region"
    ) |>
    dplyr::select(-c(X1, X6))

  ### for genes annotated with the same ensembl gene ids, ignore this

  ensgene_id_count <- as.data.frame(
    dplyr::filter(gene_info, !is.na(ensembl_gene_id)) |>
      dplyr::group_by(ensembl_gene_id) |>
      dplyr::summarise(
        n = dplyr::n(),
        .groups = "drop"
      )
  )

  gene_info <- gene_info |>
    dplyr::left_join(
      ensgene_id_count,
      by = c("ensembl_gene_id"), 
      multiple = "all"
    ) |>
    dplyr::mutate(ensembl_gene_id = dplyr::if_else(
      !is.na(n) & n > 1,
      as.character(NA),
      as.character(ensembl_gene_id)
    )) |>
    dplyr::select(-n)

  ## custom fixes
  gene_info[gene_info$symbol == "PTCSC3", ]$ensembl_gene_id <-
    "ENSG00000259104"
  gene_info[gene_info$symbol == "RUNX1", ]$ensembl_gene_id <-
    "ENSG00000159216"
  gene_info[gene_info$symbol == "TERF2IP", ]$ensembl_gene_id <-
    "ENSG00000166848"

  gene_info <- gene_info |>
    ## TEC
    dplyr::filter(entrezgene != 100124696) |>
    ## HBD
    dplyr::filter(entrezgene != 100187828) |>
    ## MMD2
    dplyr::filter(entrezgene != 100505381) |>
    ## MEMO1
    dplyr::filter(entrezgene != 7795) |>
    dplyr::mutate(other_genename_designations = dplyr::if_else(
      !is.na(other_genename_designations) &
        other_genename_designations == "-",
      as.character(NA),
      as.character(other_genename_designations)
    ))

  ## exclude entries with ambiguous gene symbols (multiple records)
  ambig_symbols <- gene_info |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::filter(n > 1)

  gene_info <- as.data.frame(gene_info |>
    dplyr::anti_join(ambig_symbols, by = "symbol"))

  attributes(gene_info)$spec <- NULL

  lgr::lgr$info(paste0(
    "Parsed n = ", nrow(gene_info),
    " genes from NCBI's Gene resource (datestamp = ",
    datestamp, ")"
  ))

  return(gene_info)
}

get_gene_aliases_ncbi <- function(gene_info,
                                  path_data_raw = NULL) {
  primary_to_primary_all <-
    gene_info |>
    janitor::clean_names() |>
    dplyr::mutate(alias = symbol) |>
    dplyr::select(entrezgene,symbol, alias) |>
    dplyr::distinct()

  gene_synonyms <-
    gene_info |>
    dplyr::select(symbol, synonyms, entrezgene) |>
    dplyr::rename(alias = synonyms) |>
    tidyr::separate_rows(alias, sep = "\\|") |>
    dplyr::filter(!(symbol == "H3P10" & alias == "p16")) |>
    dplyr::filter(!(symbol == "NRAS" & alias == "KRAS")) |>
    dplyr::filter(nchar(alias) > 2) |>
    dplyr::filter(!stringr::str_detect(alias, "^(-|[0-9])")) |>
    dplyr::anti_join(primary_to_primary_all, by = "alias")

  unique_aliases <- as.data.frame(
    gene_synonyms |>
      dplyr::group_by(alias) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n == 1) |>
      dplyr::select(-n)
  )

  ambiguous_aliases <- as.data.frame(
    gene_synonyms |>
      dplyr::group_by(alias) |>
      dplyr::summarise(n_primary_map = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n_primary_map > 1) |>
      dplyr::inner_join(gene_synonyms, by = "alias", multiple = "all") |>
      dplyr::select(
        alias, symbol, entrezgene,
        n_primary_map
      ) |>
      dplyr::mutate(
        ambiguous = TRUE,
        source = "NCBI"
      )
  )

  unambiguous_aliases <- gene_synonyms |>
    dplyr::inner_join(unique_aliases, by = "alias", multiple = "all") |>
    dplyr::bind_rows(primary_to_primary_all) |>
    dplyr::filter(
      !(alias == "PD-1" & symbol != "PDCD1")
    ) |>
    dplyr::filter(
      !(symbol == "HRAS" &
        (alias == "KRAS" |
          alias == "c-K-ras" |
          alias == "c-Ki-ras"))
    ) |>
    dplyr::arrange(symbol, alias) |>
    dplyr::mutate(symbol = dplyr::if_else(
      is.na(symbol) & !is.na(alias),
      as.character(alias),
      as.character(symbol)
    )) |>
    dplyr::mutate(
      ambiguous = FALSE,
      source = "NCBI",
      n_primary_map = 1
    ) |>
    dplyr::select(
      alias, symbol, entrezgene,
      n_primary_map, ambiguous, source
    )

  alias_custom <- as.data.frame(readr::read_tsv(
    file = file.path(path_data_raw, "custom_gene_aliases.tsv"),
    show_col_types = FALSE
  )) |>
    dplyr::mutate(
      ambiguous = FALSE,
      source = "custom",
      n_primary_map = 1
    )

  gene_alias_info <- dplyr::bind_rows(
    ambiguous_aliases,
    alias_custom,
    unambiguous_aliases
  ) |>
    dplyr::arrange(alias, symbol) |>
    dplyr::mutate(is_primary_symbol = dplyr::if_else(
      alias == symbol,
      TRUE,
      FALSE
    ))
  
  ## Find more aliases by looking at substrings
  ## (among other gene name designations) that are very
  ## similar to existing aliases - e.g. Trop-2
  
  valid_aliases <-
    gene_alias_info |> 
    dplyr::filter(ambiguous == FALSE) |> 
    dplyr::arrange(entrezgene)
  
  ambiguous_aliases <- 
    gene_alias_info |> 
    dplyr::filter(ambiguous == TRUE) |> 
    dplyr::arrange(entrezgene)
  
  other_designations <- 
    gene_info |> 
    dplyr::select(entrezgene, other_genename_designations) |> 
    tidyr::separate_rows(other_genename_designations, sep="\\|") |> 
    dplyr::mutate(other_index = dplyr::row_number()) |> 
    dplyr::rename(designation_substring = other_genename_designations) |> 
    tidyr::separate_rows(designation_substring, sep=" ") |> 
    dplyr::filter(
      stringr::str_detect(
        designation_substring, "[0-9]") & 
        nchar(designation_substring) > 4) |> 
    dplyr::mutate(
      designation_substring = stringr::str_replace_all(
        designation_substring,"^\\(|\\)$|,$","")) |>
    dplyr::mutate(
      designation_substring = stringr::str_replace_all(
        designation_substring,"\\)$","")) |>
    dplyr::filter(
        nchar(designation_substring) > 4)
 
 
  more_alias_candidates <- valid_aliases |> 
    dplyr::left_join(
      other_designations,
      by = "entrezgene",
      relationship = "many-to-many") |> 
    dplyr::mutate(
      substring_distance = stringdist::stringdist(
        tolower(alias),tolower(designation_substring))) |> 
    dplyr::filter(substring_distance == 1) |>
    dplyr::select(
      entrezgene, other_index, designation_substring) |>
    dplyr::rename(alias = designation_substring) |>
    dplyr::distinct() |>
    dplyr::group_by(alias) |>
    dplyr::summarise(
      other_index = paste(
        unique(other_index), collapse=", "),
      entrezgene = paste(
        unique(entrezgene), collapse=", "
      )) |>
    dplyr::filter(!stringr::str_detect(
      entrezgene,","
    )) |>
    dplyr::mutate(
      source = "NCBI_designation_substring",
      entrezgene = as.integer(entrezgene)) |>
    dplyr::anti_join(
      valid_aliases, by = c("alias", "entrezgene")) |>
    dplyr::anti_join(
      ambiguous_aliases, by = c("alias", "entrezgene")) |>
    dplyr::filter(
      stringr::str_detect(
        alias, "[A-Za-z]{3,}"
      )
    ) |>
    dplyr::mutate(
      ambiguous = FALSE,
      n_primary_map = 1,
      is_primary_symbol = FALSE
    ) |>
    dplyr::left_join(
      dplyr::select(
        gene_info, entrezgene, symbol
      ), by = "entrezgene"
    ) |>
    dplyr::distinct()
  
  gene_alias_info <- gene_alias_info |>
    dplyr::bind_rows(more_alias_candidates)
  
  return(gene_alias_info)
}



get_function_summary_ncbi <- function(gene_df = NULL,
                                      update = FALSE) {
  assertable::assert_colnames(
    gene_df, c("entrezgene", "gene_biotype"),
    only_colnames = FALSE,
    quiet = TRUE
  )

  pcg <- gene_df |>
    dplyr::filter(gene_biotype == "protein_coding" |
                    gene_biotype == "lincRNA") |>
    dplyr::select(entrezgene) |>
    dplyr::distinct()

  i <- 1
  ncbi_gene_summary <- data.frame()
  while (i < NROW(pcg) - 500) {
    queryset_stop <- i + 499
    queryset <- pcg[i:queryset_stop, "entrezgene"]

    summary_results <- suppressMessages(as.data.frame(
      mygene::queryMany(queryset,
        scopes = "entrezgene",
        fields = "summary",
        species = "human"
      )
    ))


    if ("summary" %in% colnames(summary_results)) {
      res <- summary_results |>
        dplyr::select(query, summary) |>
        dplyr::rename(
          entrezgene = query,
          ncbi_function_summary = summary
        ) |>
        dplyr::mutate(entrezgene = as.integer(entrezgene))

      ncbi_gene_summary <-
        dplyr::bind_rows(
          ncbi_gene_summary,
          res
        )
    }
    i <- i + 500
  }

  queryset <- pcg[i:nrow(pcg), "entrezgene"]

  summary_results <- suppressMessages(
    as.data.frame(
      mygene::queryMany(queryset,
        scopes = "entrezgene",
        fields = "summary",
        species = "human"
      )
    )
  )

  if ("summary" %in% colnames(summary_results)) {
    res <- summary_results |>
      dplyr::select(query, summary) |>
      dplyr::rename(
        entrezgene = query,
        ncbi_function_summary = summary
      ) |>
      dplyr::mutate(entrezgene = as.integer(entrezgene))

    ncbi_gene_summary <-
      dplyr::bind_rows(
        ncbi_gene_summary,
        res
      )
  }

  ncbi_gene_summary <- ncbi_gene_summary |>
    dplyr::mutate(
      ncbi_function_summary = stringr::str_replace(
        ncbi_function_summary,
        " \\[provided by RefSeq, [A-Za-z]{3} [0-9]{4}\\]\\.",
        ""
      )
    ) |>
    dplyr::mutate(
      ncbi_function_summary = stringr::str_replace(
        ncbi_function_summary,
        "\\[supplied by OMIM, [A-Za-z]{3} [0-9]{4}\\]\\.",
        ""
      )
    ) |>
    dplyr::mutate(
      ncbi_function_summary = stringr::str_replace(
        ncbi_function_summary,
        " \\[PubMed [0-9]{1,}\\]",
        ""
      )
    )
  # dplyr::mutate(
  #   gene_summary_ncbi =
  #     dplyr::if_else(
  #       !is.na(gene_summary_ncbi),
  #       paste0("<b>NCBI/RefSeq/OMIM:</b> ",
  #              gene_summary_ncbi),
  #       as.character(gene_summary_ncbi)
  #     )
  # )

  return(ncbi_gene_summary)
}

get_f1cdx <- function(gene_info = NULL) {
  lgr::lgr$info("Retrieving genes covered by Foundation One's F1CDx panel")
  f1cdx_raw <- readxl::read_xlsx(
    "data-raw/f1cdx/f1cdx.xlsx",
    col_names = TRUE
  )


  rearrangement_targets <-
    c(
      "ALK", "BCL2", "BCR", "BRAF", "BRCA1",
      "BRCA2", "CD74", "EGFR", "ETV4",
      "ETV5", "ETV6", "EWSR1", "EZR", "FGFR1",
      "FGFR2", "FGFR3", "KIT", "KMT2A",
      "MSH2", "MYB", "MYC", "NOTCH2", "NTRK1",
      "NTRK2", "NUTM1", "PDGFRA", "RAF1",
      "RARA", "RET", "ROS1", "RSPO2",
      "SDC4", "SLC34A2", "TERC",
      "TMPRSS2"
    )
  f1cdx <- dplyr::bind_rows(
    data.frame("symbol" = f1cdx_raw$sym1, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym2, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym3, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym4, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym5, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym6, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym7, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym8, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym9, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym10, stringsAsFactors = FALSE),
    data.frame("symbol" = f1cdx_raw$sym11, stringsAsFactors = FALSE)
  ) |>
    dplyr::filter(!is.na(symbol)) |>
    dplyr::mutate(symbol = dplyr::case_when(
      symbol == "C11orf30" ~ "EMSY",
      symbol == "FAM46C" ~ "TENT5C",
      symbol == "H3F3A" ~ "H3-3A",
      symbol == "MRE11A" ~ "MRE11",
      symbol == "PARK2" ~ "PRKN",
      symbol == "WHSC1" ~ "NSD2",
      symbol == "WHSC1L1" ~ "NSD3",
      TRUE ~ as.character(symbol)
    )) |>
    dplyr::mutate(F1CDx = "CNA,SNV_INDEL") |>
    dplyr::mutate(F1CDx = dplyr::if_else(
      symbol %in% rearrangement_targets,
      "CNA,SNV_INDEL,FUSION",
      as.character(F1CDx)
    )) |>
    dplyr::bind_rows(
      data.frame(
        symbol = "TERT",
        F1CDx = "PROMOTER",
        stringsAsFactors = FALSE
      )
    ) |>
    dplyr::arrange(symbol) |>
    dplyr::rename(foundation_one_f1cdx = F1CDx) |>
    dplyr::inner_join(
      dplyr::select(gene_info, symbol, entrezgene), multiple = "all"
    ) |>
    dplyr::select(-symbol)

  return(f1cdx)
}

get_tso500 <- function(gene_info = NULL,
                       gene_alias = NULL) {
  lgr::lgr$info("Retrieving genes covered by Illumina's TSO500 panel")
  tso500_all <- openxlsx::read.xlsx(
    "data-raw/tso500/journal.pone.0260089.s001.xlsx",
    sheet = "ST3"
  )

  tso500_snv_indel <-
    data.frame(
      "symbol" = stringr::str_split(tso500_all[1, ], pattern = ", "),
      stringsAsFactors = FALSE
    )
  colnames(tso500_snv_indel) <- c("symbol")
  tso500_snv_indel <- tso500_snv_indel |>
    dplyr::mutate(symbol = stringr::str_trim(symbol)) |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol),
      by = "symbol", multiple = "all"
    )
  tso500_snv_indel_complete <- tso500_snv_indel |>
    dplyr::filter(is.na(entrezgene)) |>
    dplyr::select(-c(entrezgene)) |>
    dplyr::rename(alias = symbol) |>
    dplyr::left_join(
      dplyr::filter(gene_alias$records, ambiguous == FALSE),
      by = "alias", multiple = "all"
    ) |>
    dplyr::select(-c(alias, n_primary_map, ambiguous)) |>
    dplyr::bind_rows(
      dplyr::filter(tso500_snv_indel, !is.na(entrezgene))
    ) |>
    dplyr::mutate(
      TSO500 = "SNV_INDEL"
    )

  tso500_cna <-
    data.frame(
      "symbol" = stringr::str_split(tso500_all[3, ], pattern = ", "),
      stringsAsFactors = FALSE
    )
  colnames(tso500_cna) <- c("symbol")
  tso500_cna <- tso500_cna |>
    dplyr::mutate(symbol = stringr::str_trim(symbol)) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "MYCL1", "MYCL", as.character(symbol)
    )) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        entrezgene,
        symbol
      ),
      by = "symbol", multiple = "all"
    ) |>
    dplyr::mutate(
      TSO500 = "CNA_GAIN"
    ) |>
    dplyr::mutate(TSO500 = dplyr::if_else(
      symbol == "BRCA1" |
        symbol == "BRCA2" |
        symbol == "PTEN" |
        symbol == "ATM",
      "CNA_LOSS",
      as.character(TSO500)
    ))


  tso500_snv_indel_complete <- tso500_snv_indel |>
    dplyr::filter(is.na(entrezgene)) |>
    dplyr::select(-c(entrezgene)) |>
    dplyr::rename(alias = symbol) |>
    dplyr::left_join(
      dplyr::filter(gene_alias$records, ambiguous == FALSE),
      by = "alias", multiple = "all"
    ) |>
    dplyr::select(-c(alias, n_primary_map, ambiguous)) |>
    dplyr::bind_rows(
      dplyr::filter(tso500_snv_indel, !is.na(entrezgene))
    ) |>
    dplyr::mutate(
      TSO500 = "SNV_INDEL"
    )

  tso500_rna <-
    data.frame(
      "symbol" = stringr::str_split(tso500_all[5, ], pattern = ", "),
      stringsAsFactors = FALSE
    )
  colnames(tso500_rna) <- c("symbol")
  tso500_rna <- tso500_rna |>
    dplyr::mutate(symbol = stringr::str_trim(symbol)) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "KMT2A(MLL)", "KMT2A", as.character(symbol)
    )) |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol),
      by = "symbol", multiple = "all"
    ) |>
    dplyr::mutate(
      TSO500 = "RNA_FUSION"
    )

  tso500_df <- tso500_snv_indel_complete |>
    dplyr::bind_rows(tso500_rna) |>
    dplyr::bind_rows(tso500_cna) |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(
      illumina_tso500 = paste(sort(unique(TSO500)), collapse = ","),
      .groups = "drop"
    )


  return(tso500_df)
}

get_dna_repair_genes <- function(gene_info = NULL) {
  lgr::lgr$info("Retrieving genes in DNA repair genes database")
  all_genes <- readr::read_tsv(
    file.path(
      "data-raw",
      "predisposition",
      "ge_panelapp",
      "dna_repair.tsv"
    ),
    show_col_types = FALSE
  ) |>
    janitor::clean_names() |>
    dplyr::rename(symbol = gene_symbol) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "H2AFX", "H2AX", as.character(symbol)
    )) |>
    dplyr::select(symbol, phenotypes) |>
    tidyr::separate_rows(phenotypes, sep = ";") |>
    dplyr::filter(stringr::str_detect(phenotypes, "^(Class|Activity)"))


  dnarepair_class <- all_genes |>
    dplyr::filter(stringr::str_detect(phenotypes, "^(Class)")) |>
    dplyr::rename(woods_dnarepair_class = phenotypes) |>
    dplyr::mutate(woods_dnarepair_class = stringr::str_replace(
      woods_dnarepair_class, "Class: ", ""
    )) |>
    dplyr::bind_rows(
      data.frame(symbol = "MNAT1", woods_dnarepair_class = "Unknown")
    )

  dnarepair_activity <- all_genes |>
    dplyr::filter(stringr::str_detect(phenotypes, "^(Activity)")) |>
    dplyr::rename(woods_dnarepair_activity = phenotypes) |>
    dplyr::mutate(woods_dnarepair_activity = stringr::str_replace(
      woods_dnarepair_activity, "Activity: ", ""
    ))

  dna_repair_all <-
    dnarepair_class |>
    dplyr::left_join(dnarepair_activity, 
                     by = "symbol", multiple = "all") |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol),
      by = "symbol", multiple = "all"
    ) |>
    dplyr::select(-symbol)

  return(dna_repair_all)
}

get_dbnsfp_gene_annotations <- function() {
  lgr::lgr$info(
    "Retrieving gene damage scores/OMIM annotation from dbNSFP_gene"
  )
  dbnsfp_gene <- read.table(
    file = gzfile(file.path("data-raw", "dbnsfp", "dbNSFP5.2_gene.gz")),
    sep = "\t",
    header = TRUE, stringsAsFactors = FALSE,
    na.strings = c(".", ""), comment.char = "",
    quote = NULL
  ) |>
    janitor::clean_names() |>
    dplyr::rename(
      prob_gnomad_lof_intolerant = gnom_ad_p_li,
      prob_gnomad_lof_intolerant_hom = gnom_ad_p_rec,
      prob_gnomad_lof_tolerant_null = gnom_ad_p_null,
      prob_exac_lof_intolerant = ex_ac_p_li,
      prob_haploinsuffiency = p_hi,
      prob_exac_lof_intolerant_hom = ex_ac_p_rec,
      prob_exac_lof_tolerant_null = ex_ac_p_null,
      prob_exac_nontcga_lof_intolerant = ex_ac_non_tcga_p_li,
      prob_exac_nontcga_lof_intolerant_hom = ex_ac_non_tcga_p_rec,
      prob_exac_nontcga_lof_tolerant_null =
        ex_ac_non_tcga_p_null
    ) |>
    dplyr::select(
      entrez_gene_id,
      mim_id, mim_phenotype_id,
      function_description,
      prob_haploinsuffiency,
      gene_indispensability_score,
      gene_indispensability_pred,
      essential_gene_crispr,
      essential_gene_crispr2,
      prob_gnomad_lof_intolerant,
      prob_gnomad_lof_intolerant_hom,
      prob_gnomad_lof_tolerant_null,
      prob_exac_lof_intolerant,
      prob_exac_lof_intolerant_hom,
      prob_exac_lof_tolerant_null,
      prob_exac_nontcga_lof_intolerant,
      prob_exac_nontcga_lof_intolerant_hom,
      prob_exac_nontcga_lof_tolerant_null
    )

  for (e in c(
    "gene_indispensability_score",
    "prob_gnomad_lof_intolerant",
    "prob_gnomad_lof_intolerant_hom",
    "prob_gnomad_lof_tolerant_null",
    "prob_exac_lof_intolerant",
    "prob_exac_lof_intolerant_hom",
    "prob_exac_lof_tolerant_null",
    "prob_exac_nontcga_lof_intolerant",
    "prob_exac_nontcga_lof_intolerant_hom",
    "prob_exac_nontcga_lof_tolerant_null"
  )) {
    dbnsfp_gene[, e] <- suppressWarnings(
      signif(as.numeric(dbnsfp_gene[, e]), digits = 4)
    )
  }

  dbnsfp_gene <- dbnsfp_gene |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description, "HAMAP- Rule", "HAMAP-Rule"
    )) |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description,
      paste0(
        "FUNCTION: |\\{ECO:[0-9]{1,}(\\|(PubMed|HAMAP-Rule|",
        "UniProtKB):\\S+){0,1}(, ECO:[0-9]{1,}(\\|(PubMed|UniProtKB",
        "|HAMAP-Rule):\\S+){0,}){0,}\\}\\."
      ), ""
    )) |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description, "( )?; $", ""
    )) |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description,
      "( )?\\(PubMed:[0-9]{1,}(, PubMed:[0-9]{1,}){0,}\\)", ""
    )) |>
    dplyr::mutate(dbnsfp_function_description = stringr::str_replace_all(
      function_description,
      "( )?\\(PubMed:[0-9]{1,}(, PubMed:[0-9]{1,}){0,}\\)", ""
    )) |>
    dplyr::rename(entrezgene = entrez_gene_id) |>
    dplyr::select(-function_description) |>
    # dplyr::bind_cols(pmid_all) |>
    dplyr::filter(
      !is.na(entrezgene) & stringr::str_detect(entrezgene, "^[0-9]{1,}$")
    ) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::mutate(mim_phenotype_id = stringr::str_replace_all(
      stringr::str_replace(mim_phenotype_id, ";$", ""), ";", "&"
    )) |>
    dplyr::mutate(mim_id = stringr::str_replace_all(
      stringr::str_replace(mim_id, ";$", ""), ";", "&"
    )) |>
    dplyr::mutate(prob_exac_nontcga_lof_intolerant = dplyr::if_else(
      !is.na(prob_exac_nontcga_lof_intolerant),
      format(prob_exac_nontcga_lof_intolerant,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_exac_nontcga_lof_intolerant)
    )) |>
    dplyr::mutate(prob_exac_nontcga_lof_intolerant_hom = dplyr::if_else(
      !is.na(prob_exac_nontcga_lof_intolerant_hom),
      format(prob_exac_nontcga_lof_intolerant_hom,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_exac_nontcga_lof_intolerant_hom)
    )) |>
    dplyr::mutate(prob_exac_lof_tolerant_null = dplyr::if_else(
      !is.na(prob_exac_lof_tolerant_null),
      format(prob_exac_lof_tolerant_null,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_exac_lof_tolerant_null)
    )) |>
    dplyr::mutate(prob_exac_lof_intolerant = dplyr::if_else(
      !is.na(prob_exac_lof_intolerant),
      format(prob_exac_lof_intolerant,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_exac_lof_intolerant)
    )) |>
    dplyr::mutate(prob_exac_lof_intolerant_hom = dplyr::if_else(
      !is.na(prob_exac_lof_intolerant_hom),
      format(prob_exac_lof_intolerant_hom,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_exac_lof_intolerant_hom)
    )) |>
    dplyr::mutate(prob_exac_nontcga_lof_tolerant_null = dplyr::if_else(
      !is.na(prob_exac_nontcga_lof_tolerant_null),
      format(prob_exac_nontcga_lof_tolerant_null,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_exac_nontcga_lof_tolerant_null)
    )) |>
    dplyr::mutate(prob_gnomad_lof_intolerant = dplyr::if_else(
      !is.na(prob_gnomad_lof_intolerant),
      format(prob_gnomad_lof_intolerant,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_gnomad_lof_intolerant)
    )) |>
    dplyr::mutate(prob_gnomad_lof_intolerant_hom = dplyr::if_else(
      !is.na(prob_gnomad_lof_intolerant_hom),
      format(prob_gnomad_lof_intolerant_hom,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_gnomad_lof_intolerant_hom)
    )) |>
    dplyr::mutate(prob_gnomad_lof_tolerant_null = dplyr::if_else(
      !is.na(prob_gnomad_lof_tolerant_null),
      format(prob_gnomad_lof_tolerant_null,
        scientific = TRUE, trim = TRUE, digits = 4
      ),
      as.character(prob_gnomad_lof_tolerant_null)
    ))


  return(dbnsfp_gene)
}

get_gene_signatures <- function(
    raw_db_dir = NULL,
    db_version = 'v2023.1.Hs (March 2023)'){
  
  ## get full dataset: Broad Institute's Molecular Signatures Database
  invisible(assertthat::assert_that(
    dir.exists(raw_db_dir),
    msg = paste0("Directory '",
                 raw_db_dir,"' does not exist")))
  msigdb_xml_fname <- file.path(
    raw_db_dir, "msigdb", "msigdb.xml")
  invisible(assertthat::assert_that(
    file.exists(msigdb_xml_fname),
    msg = paste0("File '",
                 msigdb_xml_fname,
                 "' does not exist")))
  
  msig_data_xml <- xml2::read_xml(msigdb_xml_fname)
  
  ## make data frame with signatures, one record pr. gene-signature association
  all_genesets <- msig_data_xml |> xml2::xml_find_all("//GENESET")
  category_code <- all_genesets |> xml2::xml_attr("CATEGORY_CODE")
  all_msigdb <- data.frame('category_code' = category_code, stringsAsFactors = F)
  all_msigdb$description <- all_genesets |> xml2::xml_attr("DESCRIPTION_BRIEF")
  all_msigdb$standard_name <- all_genesets |> xml2::xml_attr("STANDARD_NAME")
  all_msigdb$organism <- all_genesets |> xml2::xml_attr("ORGANISM")
  all_msigdb$pmid <- all_genesets |> xml2::xml_attr("PMID")
  all_msigdb$systematic_name <- all_genesets |> xml2::xml_attr("SYSTEMATIC_NAME")
  all_msigdb$subcategory_code <- all_genesets |> xml2::xml_attr("SUB_CATEGORY_CODE")
  all_msigdb$entrezgene <- all_genesets |> xml2::xml_attr("MEMBERS_EZID")
  all_msigdb$contributor <- all_genesets |> xml2::xml_attr("CONTRIBUTOR")
  all_msigdb$exact_source <- all_genesets |> xml2::xml_attr("EXACT_SOURCE")
  all_msigdb$external_url <- all_genesets |> xml2::xml_attr("EXTERNAL_DETAILS_URL")
  all_msigdb <- all_msigdb |>
    tidyr::separate_rows(entrezgene,sep=",") |>
    dplyr::filter(organism == "Homo sapiens") |>
    #dplyr::filter(subcategory_code != 'CP:KEGG') |>
    dplyr::arrange(category_code) |>
    dplyr::mutate(pmid = dplyr::if_else(
      nchar(pmid) == 0,
      as.character(NA),
      as.character(pmid))) |>
    dplyr::mutate(subcategory_code = dplyr::if_else(
      nchar(subcategory_code) == 0,
      as.character("ALL"),
      as.character(subcategory_code))) |>
    dplyr::filter(
      category_code != "ARCHIVED" &
        category_code != "C1") |>
    dplyr::mutate(external_url = dplyr::if_else(
      nchar(external_url) == 0,
      paste0("http://software.broadinstitute.org/gsea/msigdb/cards/",standard_name),
      as.character(external_url))) |>
    dplyr::mutate(description = stringr::str_replace_all(
      description,
      "( \\[ICI [0-9]{1,}(;[0-9]{1,})*\\]( )?)|( \\[GeneID=[0-9]{1,}(;[0-9]{1,})*\\]( )?)|( \\[PubChem=[0-9]{1,}(;[0-9]{1,})*\\]( )?)",""))
  
  msigdb_category_description <- read.table(
    file = file.path(
      raw_db_dir, "msigdb",
      "msigdb_collection_description.tsv"),
    sep = "\t", header = T, stringsAsFactors = F)
  
  msigdb_complete <- as.data.frame(
    all_msigdb |>
      dplyr::left_join(msigdb_category_description,
                       by = c("category_code", "subcategory_code"),
                       multiple = "all") |>
      dplyr::mutate(db = "MSigDB", db_version = db_version) |>
      dplyr::select(
        db, db_version, category_code, category_description,
        subcategory_code, subcategory_description,
        standard_name, description, organism,
        entrezgene) |>
      dplyr::rename(signature_description = description)) |>
    dplyr::filter(subcategory_code != "MIR:MIR_Legacy") |>
    dplyr::filter(subcategory_code != "TFT:TFT_Legacy") |>
    dplyr::filter(subcategory_code != "VAX") |>
    dplyr::distinct() |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^GO(BP|CC|MF)_"),
      subcategory_code,as.character(NA))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^HP_"),
      subcategory_code,as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^REACTOME_"),
      "REACTOME",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^BIOCARTA_"),
      "BIOCARTA",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^WP_"),
      "WIKIPATHWAYS",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^KEGG_"),
      "KEGG",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      stringr::str_detect(standard_name,"^PID_"),
      "PATHWAY_INTERACTION_DB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "H",
      "HALLMARK",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C2" & subcategory_code == "CGP",
      "CHEM_GEN_PERTURB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C2" & subcategory_code == "CP",
      "CANONICAL_PATHWAY",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C3" & subcategory_code == "MIR:MIRDB",
      "MICRORNA_TARGET_MIRDB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C3" & subcategory_code == "TFT:GTRD",
      "TF_TARGET_GTRD",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C4" & subcategory_code == "CGN",
      "CANCER_NEIGHBOURHOOD",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C4" & subcategory_code == "CM",
      "CANCER_MODULE",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C6","ONCOGENIC",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C7","IMMUNESIGDB",as.character(db))) |>
    dplyr::mutate(db = dplyr::if_else(
      category_code == "C8","CELLTYPE_SIGNATURES",as.character(db))) |>
    dplyr::rename(signature_name = standard_name)
  
  return(msigdb_complete)
}
  
  