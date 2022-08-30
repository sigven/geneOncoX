
get_gene_info_ncbi <- function(
    update = T){

  datestamp <- Sys.Date()
  remote_url <- "ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
  if(RCurl::url.exists(remote_url)){
    gene_info <- suppressWarnings(readr::read_tsv(
      remote_url, show_col_types = F, skip = 1,
      comment = "#", quote = "",
      col_select = c(1,2,3,5,6,9,10,14),
      progress = F, col_names = F))
  }else{
    lgr::lgr$error(
      paste0("Cannot read gene_info - download not available: ",
             remote_url))
  }

  gene_info <- gene_info |>
    dplyr::filter(X1 == 9606) |>
    dplyr::rename(
      entrezgene = X2,
      synonyms = X5,
      symbol = X3,
      name = X9,
      gene_biotype = X10,
      other_genename_designations = X14) |>
    dplyr::mutate(
      ensembl_gene_id = stringr::str_replace(
        stringr::str_match(X6,"Ensembl:ENSG[0-9]{1,}"), "Ensembl:", "")) |>
    dplyr::mutate(hgnc_id = stringr::str_replace(
      stringr::str_match(X6,"HGNC:HGNC:[0-9]{1,}"), "HGNC:HGNC:", "")) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::mutate(gene_biotype = dplyr::if_else(
      gene_biotype ==  "protin-coding",
      "protein_coding", as.character(gene_biotype))) |>
    dplyr::select(-c(X1,X6))

  ### for genes annotated with the same ensembl gene ids, ignore this annotation (set to NA)
  ensgene_id_count <- as.data.frame(
    dplyr::filter(gene_info, !is.na(ensembl_gene_id)) |>
      dplyr::group_by(ensembl_gene_id) |>
      dplyr::summarise(n = dplyr::n(),
                       .groups = "drop")
  )

  gene_info <- gene_info |>
    dplyr::left_join(
      ensgene_id_count,
      by=c("ensembl_gene_id")) |>
    dplyr::mutate(ensembl_gene_id = dplyr::if_else(
      !is.na(n) & n > 1,
      as.character(NA),
      as.character(ensembl_gene_id))) |>
    dplyr::select(-n)

  ## custom fixes
  gene_info[gene_info$symbol == "PTCSC3",]$ensembl_gene_id <-
    "ENSG00000259104"
  gene_info[gene_info$symbol == "RUNX1",]$ensembl_gene_id <-
    "ENSG00000159216"
  gene_info[gene_info$symbol == "TERF2IP",]$ensembl_gene_id <-
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

  lgr::lgr$info(paste0("Parsed n = ", nrow(gene_info),
                           " genes from NCBI's Gene resource (datestamp = ",
                           datestamp,")"))

  return(gene_info)

}

get_gene_aliases_ncbi <- function(gene_info,
                                  path_data_raw = NULL){

  primary_to_primary_all <-
    gene_info |>
    janitor::clean_names() |>
    dplyr::rename(alias = symbol) |>
    dplyr::select(entrezgene, alias) |>
    dplyr::distinct()

  gene_synonyms <-
    gene_info |>
    dplyr::select(symbol, synonyms, entrezgene) |>
    dplyr::rename(alias = synonyms) |>
    tidyr::separate_rows(alias, sep="\\|") |>
    dplyr::filter(!(symbol == "H3P10" & alias == "p16")) |>
    dplyr::filter(nchar(alias) > 2) |>
    dplyr::filter(!stringr::str_detect(alias,"^(-|[0-9])")) |>
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
      dplyr::inner_join(gene_synonyms, by = "alias") |>
      dplyr::select(alias, symbol, entrezgene,
                    n_primary_map) |>
      dplyr::mutate(ambiguous = T,
                    source = "NCBI") 
  )

  unambiguous_aliases <- gene_synonyms |>
    dplyr::inner_join(unique_aliases, by = "alias") |>
    dplyr::bind_rows(primary_to_primary_all) |>
    dplyr::filter(
      !(alias == "PD-1" & symbol != "PDCD1")) |>
    dplyr::filter(
      !(symbol == "HRAS" & 
          (alias == "KRAS" | 
             alias == "c-K-ras" | 
             alias == "c-Ki-ras"))) |>
    dplyr::arrange(symbol, alias) |>
    dplyr::mutate(symbol = dplyr::if_else(
      is.na(symbol) & !is.na(alias),
      as.character(alias),
      as.character(symbol)
    )) |>
    dplyr::mutate(ambiguous = F, 
                  source = "NCBI",
                  n_primary_map = 1) |>
    dplyr::select(alias, symbol, entrezgene,
                  n_primary_map, ambiguous, source)
  
  alias_custom <- as.data.frame(readr::read_tsv(
    file = file.path(path_data_raw,"custom_gene_aliases.tsv"),
    show_col_types = F)) |>
    dplyr::mutate(ambiguous = FALSE,
                  source = "custom",
                  n_primary_map = 1)

  gene_alias_info <- dplyr::bind_rows(
    ambiguous_aliases,
    alias_custom,
    unambiguous_aliases) |>
    dplyr::arrange(alias, symbol)

  return(gene_alias_info)
}



get_function_summary_ncbi <- function(
    gene_df = NULL,
    update = F){

  # assertable::assert_colnames(
  #   gene_df, c("entrezgene","gene_biotype"), only_colnames = F, quiet = T)
  #
  # pcg <- gene_df |>
  #   dplyr::filter(gene_biotype == "protein-coding") |>
  #   dplyr::select(entrezgene) |>
  #   dplyr::distinct()

  pcg <- gene_df

  i <- 1
  ncbi_gene_summary <- data.frame()
  while(i < NROW(pcg) - 300){

    queryset_stop <- i + 299
    queryset <- pcg[i:queryset_stop,"entrezgene"]

    summary_results <- suppressMessages(as.data.frame(
      mygene::queryMany(queryset,
                        scopes = "entrezgene",
                        fields= "summary",
                        species = "human")
    ))


    if("summary" %in% colnames(summary_results)){
      res <- summary_results |>
        dplyr::select(query, summary) |>
        dplyr::rename(entrezgene = query,
                      ncbi_function_summary = summary) |>
        dplyr::mutate(entrezgene = as.integer(entrezgene))

      ncbi_gene_summary <-
        dplyr::bind_rows(ncbi_gene_summary,
                         res)
    }
    i <- i + 300

  }

  queryset <- pcg[i:nrow(pcg),"entrezgene"]

  summary_results <- suppressMessages(
    as.data.frame(
      mygene::queryMany(queryset,
                        scopes = "entrezgene",
                        fields = "summary",
                        species = "human")
    ))

  if("summary" %in% colnames(summary_results)){
    res <- summary_results |>
      dplyr::select(query, summary) |>
      dplyr::rename(entrezgene = query,
                    ncbi_function_summary = summary) |>
      dplyr::mutate(entrezgene = as.integer(entrezgene))

    ncbi_gene_summary <-
      dplyr::bind_rows(ncbi_gene_summary,
                       res)
  }

  ncbi_gene_summary <- ncbi_gene_summary |>
    dplyr::mutate(
      ncbi_function_summary = stringr::str_replace(
        ncbi_function_summary,
        " \\[provided by RefSeq, [A-Za-z]{3} [0-9]{4}\\]\\.",
        "")
    ) |>
    dplyr::mutate(
      ncbi_function_summary = stringr::str_replace(
        ncbi_function_summary,
        "\\[supplied by OMIM, [A-Za-z]{3} [0-9]{4}\\]\\.",
        "")
    ) |>
    dplyr::mutate(
      ncbi_function_summary = stringr::str_replace(
        ncbi_function_summary,
        " \\[PubMed [0-9]{1,}\\]",
        "")
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

get_tso500 <- function(gene_info = NULL,
                       gene_alias = NULL){

  lgr::lgr$info("Retrieve genes covered by Illumina's TSO500 panel")
  tso500_all <- openxlsx::read.xlsx(
    "data-raw/tso500/journal.pone.0260089.s001.xlsx", sheet = "ST3")

  tso500_snv_indel <-
    data.frame('symbol' = stringr::str_split(tso500_all[1,], pattern = ", "),
               stringsAsFactors = F)
  colnames(tso500_snv_indel) <- c('symbol')
  tso500_snv_indel <- tso500_snv_indel |>
    dplyr::mutate(symbol = stringr::str_trim(symbol)) |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol), by = "symbol")
  tso500_snv_indel_complete <- tso500_snv_indel |>
    dplyr::filter(is.na(entrezgene)) |>
    dplyr::select(-c(entrezgene)) |>
    dplyr::rename(alias = symbol) |>
    dplyr::left_join(
      dplyr::filter(gene_alias$records, ambiguous == F),
      by = "alias") |>
    dplyr::select(-c(alias, n_primary_map, ambiguous)) |>
    dplyr::bind_rows(
      dplyr::filter(tso500_snv_indel, !is.na(entrezgene))) |>
    dplyr::mutate(
      TSO500 = "SNV_INDEL"
    )

  tso500_cna <-
    data.frame('symbol' = stringr::str_split(tso500_all[3,], pattern = ", "),
               stringsAsFactors = F)
  colnames(tso500_cna) <- c('symbol')
  tso500_cna <- tso500_cna |>
    dplyr::mutate(symbol = stringr::str_trim(symbol)) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "MYCL1","MYCL", as.character(symbol)
    )) |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol), by = "symbol") |>
    dplyr::mutate(
      TSO500 = "CNA_GAIN") |>
    dplyr::mutate(TSO500 = dplyr::if_else(
      symbol == "BRCA1" | symbol == "BRCA2" | symbol == "PTEN" | symbol == "ATM",
      "CNA_LOSS",
      as.character(TSO500)
    ))


  tso500_snv_indel_complete <- tso500_snv_indel |>
    dplyr::filter(is.na(entrezgene)) |>
    dplyr::select(-c(entrezgene)) |>
    dplyr::rename(alias = symbol) |>
    dplyr::left_join(
      dplyr::filter(gene_alias$records, ambiguous == F),
      by = "alias") |>
    dplyr::select(-c(alias, n_primary_map, ambiguous)) |>
    dplyr::bind_rows(
      dplyr::filter(tso500_snv_indel, !is.na(entrezgene))) |>
    dplyr::mutate(
      TSO500 = "SNV_INDEL"
    )

  tso500_rna <-
    data.frame('symbol' = stringr::str_split(tso500_all[5,], pattern = ", "),
               stringsAsFactors = F)
  colnames(tso500_rna) <- c('symbol')
  tso500_rna <- tso500_rna |>
    dplyr::mutate(symbol = stringr::str_trim(symbol)) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "KMT2A(MLL)","KMT2A", as.character(symbol))) |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol), by = "symbol") |>
    dplyr::mutate(
      TSO500 = "RNA_FUSION"
    )

  tso500_df <- tso500_snv_indel_complete |>
    dplyr::bind_rows(tso500_rna) |>
    dplyr::bind_rows(tso500_cna) |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(illumina_tso500 = paste(sort(unique(TSO500)), collapse=","),
                     .groups = "drop")


  return(tso500_df)
}

get_dna_repair_genes <- function(gene_info = NULL){

  lgr::lgr$info("Retrieve genes in DNA repair genes database")
  all_genes <- readr::read_tsv(
    file.path(
      "data-raw",
      "predisposition",
      "ge_panelapp",
      "dna_repair.tsv"
    ), show_col_types = F
  ) |>
    janitor::clean_names() |>
    dplyr::rename(symbol = gene_symbol) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "H2AFX","H2AX",as.character(symbol)
    )) |>
    dplyr::select(symbol, phenotypes) |>
    tidyr::separate_rows(phenotypes, sep = ";") |>
    dplyr::filter(stringr::str_detect(phenotypes, "^(Class|Activity)"))


  dnarepair_class <- all_genes |>
    dplyr::filter(stringr::str_detect(phenotypes, "^(Class)")) |>
    dplyr::rename(woods_dnarepair_class = phenotypes) |>
    dplyr::mutate(woods_dnarepair_class = stringr::str_replace(
      woods_dnarepair_class, "Class: ",""
    )) |>
    dplyr::bind_rows(
      data.frame(symbol = 'MNAT1', woods_dnarepair_class = 'Unknown')
    )

  dnarepair_activity <- all_genes |>
    dplyr::filter(stringr::str_detect(phenotypes, "^(Activity)")) |>
    dplyr::rename(woods_dnarepair_activity = phenotypes) |>
    dplyr::mutate(woods_dnarepair_activity = stringr::str_replace(
      woods_dnarepair_activity, "Activity: ",""
    ))

  dna_repair_all <-
    dnarepair_class |>
    dplyr::left_join(dnarepair_activity, by = "symbol") |>
    dplyr::left_join(
      dplyr::select(gene_info, entrezgene, symbol),
      by = "symbol"
    ) |>
    dplyr::select(-symbol)

  return(dna_repair_all)

}

get_dbnsfp_gene_annotations <- function(){

  lgr::lgr$info("Retrieving gene damage scores/OMIM annotation from dbNSFP_gene")
  dbnsfp_gene <- read.table(
    file=gzfile(file.path("data-raw","dbnsfp", "dbNSFP_gene.gz")), sep="\t",
    header = T, stringsAsFactors = F, na.strings = c(".",""), comment.char="",
    quote = NULL) |>
    janitor::clean_names() |>
    dplyr::rename(prob_gnomad_lof_intolerant = gnom_ad_p_li,
                  prob_gnomad_lof_intolerant_hom = gnom_ad_p_rec,
                  prob_gnomad_lof_tolerant_null = gnom_ad_p_null,
                  prob_exac_lof_intolerant = ex_ac_p_li,
                  prob_haploinsuffiency = p_hi,
                  prob_exac_lof_intolerant_hom = ex_ac_p_rec,
                  prob_exac_lof_tolerant_null = ex_ac_p_null,
                  prob_exac_nontcga_lof_intolerant = ex_ac_non_tcga_p_li,
                  prob_exac_nontcga_lof_intolerant_hom = ex_ac_non_tcga_p_rec,
                  prob_exac_nontcga_lof_tolerant_null = ex_ac_non_tcga_p_null) |>
    dplyr::select(entrez_gene_id,
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
                  prob_exac_nontcga_lof_tolerant_null)

  for(e in c('gene_indispensability_score',
             'prob_gnomad_lof_intolerant',
             'prob_gnomad_lof_intolerant_hom',
             'prob_gnomad_lof_tolerant_null',
             'prob_exac_lof_intolerant',
             'prob_exac_lof_intolerant_hom',
             'prob_exac_lof_tolerant_null',
             'prob_exac_nontcga_lof_intolerant',
             'prob_exac_nontcga_lof_intolerant_hom',
             'prob_exac_nontcga_lof_tolerant_null')){
    dbnsfp_gene[,e] <- suppressWarnings(
      signif(as.numeric(dbnsfp_gene[,e]),digits = 4))
  }

  # pmid_hits <- stringr::str_match_all(
  #   dbnsfp_gene$function_description,"PubMed:[0-9]{1,}")
  # i <- 1
  # pmid_all <- data.frame()
  # while(i <= length(pmid_hits)){
  #   pmid_df <- data.frame(
  #     'function_description_pmid' = paste(unique(sort(
  #       stringr::str_replace_all(pmid_hits[[i]][,1],"PubMed:",""))),
  #       collapse=","),
  #     stringsAsFactors = F) |>
  #     dplyr::mutate(function_description_pmid = dplyr::if_else(
  #       function_description_pmid == "",
  #       as.character(NA),
  #       as.character(function_description_pmid)))
  #   pmid_all <- dplyr::bind_rows(pmid_all, pmid_df)
  #   i <- i + 1
  #
  # }

  dbnsfp_gene <- dbnsfp_gene |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description,"HAMAP- Rule","HAMAP-Rule")) |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description,"FUNCTION: |\\{ECO:[0-9]{1,}(\\|(PubMed|HAMAP-Rule|UniProtKB):\\S+){0,1}(, ECO:[0-9]{1,}(\\|(PubMed|UniProtKB|HAMAP-Rule):\\S+){0,}){0,}\\}\\.","")) |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description,"( )?; $","")) |>
    dplyr::mutate(function_description = stringr::str_replace_all(
      function_description,"( )?\\(PubMed:[0-9]{1,}(, PubMed:[0-9]{1,}){0,}\\)","")) |>
    dplyr::mutate(dbnsfp_function_description = stringr::str_replace_all(
      function_description,"( )?\\(PubMed:[0-9]{1,}(, PubMed:[0-9]{1,}){0,}\\)","")) |>
    dplyr::rename(entrezgene = entrez_gene_id) |>
    dplyr::select(-function_description) |>
    #dplyr::bind_cols(pmid_all) |>
    dplyr::filter(!is.na(entrezgene) & stringr::str_detect(entrezgene,"^[0-9]{1,}$")) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::mutate(mim_phenotype_id = stringr::str_replace_all(
      stringr::str_replace(mim_phenotype_id,";$",""),";","&")) |>
    dplyr::mutate(mim_id = stringr::str_replace_all(
      stringr::str_replace(mim_id,";$",""),";","&")) |>
    dplyr::mutate(prob_exac_nontcga_lof_intolerant = dplyr::if_else(
      !is.na(prob_exac_nontcga_lof_intolerant),
      format(prob_exac_nontcga_lof_intolerant,
             scientific = T, trim = T, digits = 4),
      as.character(prob_exac_nontcga_lof_intolerant)
    )) |>
    dplyr::mutate(prob_exac_nontcga_lof_intolerant_hom = dplyr::if_else(
      !is.na(prob_exac_nontcga_lof_intolerant_hom),
      format(prob_exac_nontcga_lof_intolerant_hom,
             scientific = T, trim = T, digits = 4),
      as.character(prob_exac_nontcga_lof_intolerant_hom)
    )) |>
    dplyr::mutate(prob_exac_lof_tolerant_null = dplyr::if_else(
      !is.na(prob_exac_lof_tolerant_null),
      format(prob_exac_lof_tolerant_null,
             scientific = T, trim = T, digits = 4),
      as.character(prob_exac_lof_tolerant_null)
    )) |>
    dplyr::mutate(prob_exac_lof_intolerant = dplyr::if_else(
      !is.na(prob_exac_lof_intolerant),
      format(prob_exac_lof_intolerant,
             scientific = T, trim = T, digits = 4),
      as.character(prob_exac_lof_intolerant)
    )) |>
    dplyr::mutate(prob_exac_lof_intolerant_hom = dplyr::if_else(
      !is.na(prob_exac_lof_intolerant_hom),
      format(prob_exac_lof_intolerant_hom,
             scientific = T, trim = T, digits = 4),
      as.character(prob_exac_lof_intolerant_hom)
    )) |>
    dplyr::mutate(prob_exac_nontcga_lof_tolerant_null = dplyr::if_else(
      !is.na(prob_exac_nontcga_lof_tolerant_null),
      format(prob_exac_nontcga_lof_tolerant_null,
             scientific = T, trim = T, digits = 4),
      as.character(prob_exac_nontcga_lof_tolerant_null)
    )) |>
    dplyr::mutate(prob_gnomad_lof_intolerant = dplyr::if_else(
      !is.na(prob_gnomad_lof_intolerant),
      format(prob_gnomad_lof_intolerant,
             scientific = T, trim = T, digits = 4),
      as.character(prob_gnomad_lof_intolerant)
    )) |>
    dplyr::mutate(prob_gnomad_lof_intolerant_hom = dplyr::if_else(
      !is.na(prob_gnomad_lof_intolerant_hom),
      format(prob_gnomad_lof_intolerant_hom,
             scientific = T, trim = T, digits = 4),
      as.character(prob_gnomad_lof_intolerant_hom)
    )) |>
    dplyr::mutate(prob_gnomad_lof_tolerant_null = dplyr::if_else(
      !is.na(prob_gnomad_lof_tolerant_null),
      format(prob_gnomad_lof_tolerant_null,
             scientific = T, trim = T, digits = 4),
      as.character(prob_gnomad_lof_tolerant_null)
    ))


  return(dbnsfp_gene)
}
