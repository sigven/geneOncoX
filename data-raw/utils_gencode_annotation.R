

get_transcript_appris_annotation <-
  function(build = "grch37") {
    appris_base_url <-
      paste0(
        "http://apprisws.bioinfo.cnio.es/pub/",
        "current_release/datafiles/homo_sapiens/"
      )

    remote_appris <- paste0(
      appris_base_url,
      paste0(toupper(substr(build, 0, 3)), substr(build, 4, 6)),
      # "/appris_data.principal.txt")
      "/appris_data.principal.forENSEMBL.txt"
    )

    appris <- readr::read_tsv(
      file = remote_appris,
      show_col_types = F,
      col_names = F
    ) |>
      purrr::set_names(
        c(
          "ensembl_gene_id",
          "ensembl_transcript_id",
          "principal_isoform_flag"
        )
      )

    lgr::lgr$info(
      paste0(
        "A total of ",
        nrow(appris),
        " transcripts with APPRIS principal isoform ",
        "annotations were parsed"
      )
    )

    return(appris)
  }


gencode_get_transcripts <-
  function(build = "grch38",
           append_regulatory_region = TRUE,
           gencode_version = 46,
           ensembl_version = 112,
           uniprot_version = "2024_04",
           gene_info = NULL,
           gene_alias = NULL) {
    gencode_ftp_url <-
      "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    lgr::lgr$info(paste0(
      "Retrieving GENCODE transcripts - version ",
      gencode_version, ", build ", build
    ))

    remote_gtf_all <-
      paste0(
        gencode_ftp_url, "release_", gencode_version, "/",
        "gencode.v", gencode_version, 
        ".chr_patch_hapl_scaff",
        ".annotation.gtf.gz"
      )
    
    remote_gtf <-
      paste0(
        gencode_ftp_url, "release_", gencode_version, "/",
        "gencode.v", gencode_version, 
        ".annotation.gtf.gz"
      )
    
    #https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz

    options(timeout = 10000000)
    gencode_gtf_transcripts <- valr::read_gtf(remote_gtf_all) |>
      dplyr::filter(type == "transcript")
    
    if(gencode_version == 19){
      gencode_gtf_transcripts_chr <- valr::read_gtf(remote_gtf) |>
        dplyr::filter(type == "transcript") |>
        dplyr::select(transcript_id, protein_id)
      
      gencode_gtf_transcripts <- gencode_gtf_transcripts |>
        dplyr::left_join(
          gencode_gtf_transcripts_chr, by = "transcript_id"
        )
    }

    gencode_gtf <- gencode_gtf_transcripts |>
      dplyr::rename(
        transcript_start = start,
        transcript_end = end,
        ensembl_gene_id = gene_id,
        symbol = gene_name,
        gene_biotype = gene_type,
        transcript_biotype = transcript_type,
        ensembl_transcript_id_full = transcript_id,
        ensembl_protein_id = protein_id
      ) |>
      dplyr::select(
        chrom,
        transcript_start,
        transcript_end,
        strand,
        symbol,
        ensembl_gene_id,
        ensembl_transcript_id_full,
        gene_biotype,
        transcript_biotype,
        ensembl_protein_id
      ) |>
      dplyr::mutate(
        ensembl_transcript_id =
          stringr::str_replace(
            ensembl_transcript_id_full,
            "\\.[0-9]{1,}$", ""
          )
      ) |>
      dplyr::mutate(
        ensembl_gene_id =
          stringr::str_replace(
            ensembl_gene_id,
            "\\.[0-9]{1,}$", ""
          )
      ) |>
      dplyr::distinct()

    ## Temporary fix: Need additional processing to get multiple tag-value
    ## pairs in a proper fashion
    gencode_gtf_fix <- as.data.frame(
      data.table::fread(remote_gtf_all, skip = 5, verbose = F) |>
        dplyr::filter(V3 == "transcript")
    )
    tag_cols <- stringr::str_match_all(
      gencode_gtf_fix$V9, "tag \\\"(\\S+)\\\";"
    ) |>
      purrr::map_chr(~ stringr::str_c(.x[, ncol(.x)], collapse = "&"))

    gencode_gtf$tag <- tag_cols
    gencode <- gencode_gtf |>
      dplyr::mutate(tag = dplyr::if_else(
        nchar(tag) == 0,
        as.character(NA),
        as.character(tag)
      )) |>
      dplyr::filter(!is.na(ensembl_gene_id) &
        !is.na(ensembl_transcript_id)) |>
      dplyr::distinct()

    #
    lgr::lgr$info(paste0(
      "A total of ", nrow(gencode),
      " transcripts parsed"
    ))

    ## include regulatory region (for VEP annotation)
    if (append_regulatory_region == TRUE) {
      suppressPackageStartupMessages(
        library(BSgenome.Hsapiens.UCSC.hg19)
      )
      suppressPackageStartupMessages(
        library(BSgenome.Hsapiens.UCSC.hg38)
      )
      lgr::lgr$info(
        paste0(
          "Parameter 'append_regulatory_region' is TRUE: ",
          "expanding transcript start/end with 5kb ",
          "(for compliance with VEP upstream/downstream region)"
        )
      )
      chromosome_lengths <- data.frame(
        "chrom" = head(
          names(
            GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19)
          ),
          25
        ),
        "chrom_length" = head(
          GenomeInfoDb::seqlengths(
            BSgenome.Hsapiens.UCSC.hg19
          ), 25
        ),
        stringsAsFactors = FALSE, row.names = NULL
      )
      if (build == "grch38") {
        chromosome_lengths <- data.frame(
          "chrom" = head(
            names(
              GenomeInfoDb::seqlengths(
                BSgenome.Hsapiens.UCSC.hg38
              )
            ), 25
          ),
          "chrom_length" = head(
            GenomeInfoDb::seqlengths(
              BSgenome.Hsapiens.UCSC.hg38
            ), 25
          ),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }
      gencode <- as.data.frame(
        gencode |>
          dplyr::left_join(chromosome_lengths, 
                           by = c("chrom"), multiple = "all") |>
          dplyr::rowwise() |>
          dplyr::mutate(
            start = dplyr::if_else(
              startsWith(chrom,"chr") &
                chrom != "chrM",
              as.integer(
                max(1, transcript_start - 5001)
              ),
              as.integer(transcript_start)
            ),
            end = dplyr::if_else(
              startsWith(chrom,"chr") &
                chrom != "chrM",
              as.integer(
                min(chrom_length, transcript_end + 5001)
              ),
              as.integer(transcript_end)
            )
          ) |>
          dplyr::select(-chrom_length)
      )
    }

    attr(gencode, "groups") <- NULL

    transcript_df <- gencode |>
      dplyr::select(
        ensembl_transcript_id,
        ensembl_gene_id,
        symbol
      ) |>
      dplyr::distinct()
    
    

    ## get gene and transcript cross-references (biomart + gene_info)
    gencode_xrefs <- gencode_resolve_xrefs(
      transcript_df = transcript_df,
      ensembl_version = ensembl_version,
      build = build,
      gene_info = gene_info,
      gene_alias = gene_alias
    )
    
    #if(build == "grch37"){
    gene_info2 <- gene_info |>
      dplyr::select(entrezgene, symbol) |>
      dplyr::distinct()
    gencode_xrefs1 <- gencode_xrefs |>
      dplyr::filter(!is.na(entrezgene)) |>
      dplyr::select(-c("symbol")) |>
      dplyr::left_join(
        gene_info2,
        by = c("entrezgene")
      ) |>
      dplyr::distinct()
    
    gencode_xrefs2 <- gencode_xrefs |>
      dplyr::filter(is.na(entrezgene))
    gencode_xrefs3 <- dplyr::bind_rows(
      gencode_xrefs1, gencode_xrefs2)

    appris_annotation <- get_transcript_appris_annotation(
      build = build
    )

    uniprot_map <- get_uniprot_map(uniprot_version = uniprot_version) |>
      dplyr::select(-uniprot_version)

    gencode <- gencode |>
      dplyr::rename(symbol_gencode = symbol) |>
      dplyr::left_join(
        gencode_xrefs,
        by = c(
          #"symbol", "ensembl_gene_id",
          "ensembl_gene_id",
          "ensembl_transcript_id"
        ), multiple = "all"
      ) |>
      dplyr::mutate(symbol = dplyr::case_when(
        is.na(symbol) ~ symbol_gencode,
        TRUE ~ as.character(symbol)
      )) |>
      dplyr::left_join(
        dplyr::select(
          appris_annotation,
          principal_isoform_flag,
          ensembl_transcript_id
        ),
        by = "ensembl_transcript_id", 
        multiple = "all"
      ) |>
      gencode_expand_basic() |>
      dplyr::left_join(
        uniprot_map,
        by = c(
          "uniprot_acc",
          "ensembl_transcript_id"
        ), multiple = "all"
      ) |>
      sort_bed_regions() |>
      dplyr::rename(
        mane_select = transcript_mane_select,
        mane_plus_clinical = transcript_mane_plus_clinical,
      ) |>
      dplyr::mutate(
        chrom = as.character(chrom),
        entrezgene = as.integer(entrezgene),
        ensembl_version = ensembl_version,
        gencode_version = gencode_version,
        uniprot_version = uniprot_version
      ) |>
      dplyr::select(
        chrom,
        start,
        end,
        transcript_start,
        transcript_end,
        strand,
        ensembl_gene_id,
        ensembl_transcript_id,
        ensembl_transcript_id_full,
        ensembl_protein_id,
        symbol,
        symbol_gencode,
        hgnc_id,
        entrezgene,
        name,
        gene_biotype,
        transcript_biotype,
        tag,
        refseq_protein_id,
        refseq_transcript_id,
        mane_select,
        mane_plus_clinical,
        principal_isoform_flag,
        uniprot_acc,
        uniprot_id,
        ensembl_version,
        gencode_version,
        uniprot_version
      )
    
    if(build == "grch37"){
      gencode <- gencode |>
        dplyr::mutate(gene_biotype = dplyr::case_when(
          gene_biotype == "lincRNA" ~ "lncRNA",
          TRUE ~ as.character(gene_biotype)
        ))
    }
    
    if(build == "grch38"){
      gencode <- gencode |>
        dplyr::mutate(gene_biotype = dplyr::case_when(
          stringr::str_detect(gene_biotype, "processed") &
            stringr::str_detect(gene_biotype, "pseudogene") ~ "pseudogene",
          TRUE ~ as.character(gene_biotype)
        ))
    }
    

    lgr::lgr$info(paste0(
      "A total of ",
      nrow(gencode), " valid transcripts remaining"
    ))

    return(gencode)
  }

gencode_expand_basic <- function(gencode) {
  ## get all genes that has an annotated basic transcript

  ## perform anti_join against the set defined below

  gencode_basicset <- gencode |>
    dplyr::select(ensembl_gene_id, tag) |>
    dplyr::filter(!is.na(tag) &
      stringr::str_detect(tag, "basic"))


  single_transcripts_per_gene <- as.data.frame(
    gencode |>
      dplyr::filter(
        stringr::str_detect(
          gene_biotype,
          "^(IG_|TR_)|scaRNA|snRNA|snoRNA|sRNA|scRNA|pseudogene"
        )
      ) |>
      dplyr::select(
        ensembl_gene_id, symbol, tag,
        ensembl_transcript_id
      ) |>
      dplyr::anti_join(
        gencode_basicset, by = "ensembl_gene_id") |>
      dplyr::group_by(symbol, ensembl_gene_id) |>
      dplyr::summarise(
        n = dplyr::n(),
        tag = paste(tag, collapse = ","),
        .groups = "drop"
      ) |>
      dplyr::filter(n == 1) |>
      dplyr::filter(!stringr::str_detect(tag, "basic")) |>
      dplyr::select(ensembl_gene_id, symbol, tag) |>
      dplyr::mutate(basic_expanded = TRUE)
  )

  gencode <- gencode |>
    dplyr::left_join(
      dplyr::select(
        single_transcripts_per_gene,
        ensembl_gene_id, basic_expanded
      ),
      by = "ensembl_gene_id", multiple = "all"
    ) |>
    dplyr::mutate(tag = dplyr::case_when(
      is.na(tag) & basic_expanded == TRUE ~ "basic",
      !is.na(tag) &
        !stringr::str_detect(tag, "basic") &
        basic_expanded == TRUE ~ paste0("basic&", tag),
      TRUE ~ as.character(tag)
    )) |>
    dplyr::select(-basic_expanded)

  return(gencode)
}



gencode_resolve_xrefs <- function(transcript_df = NULL,
                                  build = "grch38",
                                  ensembl_version = 111,
                                  gene_info = NULL,
                                  gene_alias = NULL) {
  invisible(assertable::assert_colnames(
    transcript_df,
    colnames = c(
      "ensembl_transcript_id",
      "ensembl_gene_id",
      "symbol"
    ),
    quiet = TRUE
  ))

  options(timeout = 10000000)
  ensembl_mart <- list()
  ensembl_mart[["grch38"]] <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "hsapiens_gene_ensembl",
    #host = "https://www.ensembl.org",
    version = ensembl_version
  )

  ensembl_mart[["grch37"]] <- biomaRt::useEnsembl(
    biomart = "genes",
    GRCh = "37",
    #host = "https://www.ensembl.org",
    dataset = "hsapiens_gene_ensembl"
  )

  queryAttributes1 <- c(
    "ensembl_transcript_id",
    "refseq_mrna",
    "ensembl_gene_id",
    "hgnc_id",
    "entrezgene_id"
  )
  queryAttributes2 <- c(
    "ensembl_transcript_id",
    "uniprotswissprot",
    "refseq_peptide",
    "description"
  )

  queryAttributes3 <- c(
    "ensembl_transcript_id",
    "transcript_mane_select",
    "transcript_mane_plus_clinical"
  )

  queryAttributes4 <- c(
    "ensembl_transcript_id",
    "refseq_ncrna"
  )
  

  xref_biomart_1 <- biomaRt::getBM(
    attributes = queryAttributes1,
    mart = ensembl_mart[[build]]) |>
    dplyr::rename(entrezgene = entrezgene_id)
  
  xref_biomart_2 <- biomaRt::getBM(
    attributes = queryAttributes2,
    mart = ensembl_mart[[build]]) |>
    dplyr::rename(uniprot_acc = uniprotswissprot)

  ## map MANE xrefs with grch38 mart (not available for grch37)
  xref_biomart_3 <- biomaRt::getBM(
    attributes = queryAttributes3,
    mart = ensembl_mart[["grch38"]]
  )

  xref_biomart_4 <- biomaRt::getBM(
    attributes = queryAttributes4,
    mart = ensembl_mart[[build]])
  
  xref_biomart <- xref_biomart_1 |>
    dplyr::left_join(
      xref_biomart_2,
      by = "ensembl_transcript_id", 
      multiple = "all", 
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      xref_biomart_3,
      by = "ensembl_transcript_id", 
      multiple = "all",
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      xref_biomart_4,
      by = "ensembl_transcript_id", 
      multiple = "all",
      relationship = "many-to-many"
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(hgnc_id = as.integer(
      stringr::str_replace(hgnc_id, "HGNC:", "")
    )) |>
    dplyr::rename(
      refseq_transcript_id = refseq_mrna,
      refseq_protein_id = refseq_peptide
    ) |>
    dplyr::mutate(refseq_transcript_id = dplyr::if_else(
      refseq_transcript_id == "", 
      as.character(refseq_ncrna),
      as.character(refseq_transcript_id)
    ))
  
  valid_entrez_recs <- xref_biomart |> 
    dplyr::filter(!is.na(entrezgene) & 
                    !is.na(hgnc_id)) |> 
    #dplyr::select(ensembl_gene_id, 
    #              ensembl_transcript_id, 
    #              entrezgene, 
    #              hgnc_id) |> 
    dplyr::distinct()
  
  valid_other_recs <- xref_biomart |> 
    dplyr::filter(is.na(entrezgene)) |> 
    #dplyr::select(ensembl_gene_id, 
    #              ensembl_transcript_id, 
    #              entrezgene, 
    #              hgnc_id) |> 
    dplyr::distinct()
  
  invalid_entrez_recs <- xref_biomart |> 
    dplyr::filter(!is.na(entrezgene) & is.na(hgnc_id)) |> 
    #dplyr::select(
    #  ensembl_gene_id, ensembl_transcript_id, 
    #  entrezgene, hgnc_id) |> 
    dplyr::distinct() |> 
    dplyr::semi_join(
      valid_entrez_recs, by = "entrezgene") |> 
    dplyr::distinct() |>
    dplyr::mutate(entrezgene = as.integer(NA))
  
  xref_biomart <- rbind(
    valid_entrez_recs,
    valid_other_recs,
    invalid_entrez_recs
  )

  for (n in c(
    "refseq_protein_id",
    "refseq_transcript_id",
    "uniprot_acc",
    "transcript_mane_select",
    "transcript_mane_plus_clinical",
    "description"
  )) {
    xref_biomart[!is.na(xref_biomart[, n]) &
      (xref_biomart[, n] == "" |
        xref_biomart[, n] == "NA"), ][, n] <- NA
  }

  for (xref in c(
    "refseq_protein_id",
    "refseq_transcript_id",
    "uniprot_acc",
    "transcript_mane_select",
    "transcript_mane_plus_clinical",
    "description"
  )) {
    ensXref <- as.data.frame(
      xref_biomart |>
        dplyr::select(
          ensembl_gene_id,
          ensembl_transcript_id,
          !!rlang::sym(xref)
        ) |>
        dplyr::distinct() |>
        dplyr::group_by(
          ensembl_transcript_id,
          ensembl_gene_id
        ) |>
        dplyr::summarise(
          !!rlang::sym(xref) := paste(unique(!!rlang::sym(xref)),
            collapse = "&"
          ),
          .groups = "drop"
        ) |>
        dplyr::mutate(!!rlang::sym(xref) := dplyr::if_else(
          !!rlang::sym(xref) == "NA", as.character(NA),
          as.character(!!rlang::sym(xref))
        ))
    )
    transcript_df <- transcript_df |>
      dplyr::left_join(ensXref, by = c(
        "ensembl_gene_id",
        "ensembl_transcript_id"
      ), multiple = "all", 
      relationship = "many-to-many")
  }

  for (xref in c(
    "entrezgene",
    "hgnc_id"
  )) {
    ensXref <- as.data.frame(
      xref_biomart |>
        dplyr::select(
          ensembl_gene_id,
          !!rlang::sym(xref)
        ) |>
        dplyr::distinct() |>
        dplyr::group_by(
          ensembl_gene_id
        ) |>
        dplyr::summarise(
          !!rlang::sym(xref) := paste(unique(!!rlang::sym(xref)),
            collapse = "&"
          ),
          .groups = "drop"
        ) |>
        dplyr::mutate(!!rlang::sym(xref) := dplyr::if_else(
          !!rlang::sym(xref) == "NA", as.character(NA),
          as.character(!!rlang::sym(xref))
        ))
    )
    transcript_df <- transcript_df |>
      dplyr::left_join(
        ensXref, by = c("ensembl_gene_id"), 
        multiple = "all", 
        relationship = "many-to-many")
  }


  gene_xrefs_maps <- list()

  ## map gene cross-references (name) by hgnc identifier
  gene_xrefs_maps[["by_hgnc"]] <- as.data.frame(
    transcript_df |>
      dplyr::filter(!is.na(hgnc_id)) |>
      dplyr::select(-entrezgene) |>
      dplyr::left_join(
        dplyr::select(
          gene_info,
          hgnc_id,
          entrezgene,
          name
        ),
        by = "hgnc_id", 
        multiple = "all", 
        relationship = "many-to-many"
      )
  )

  ## map gene cross-references by Entrez gene identifier
  ## given that hgnc id is not provided
  gene_xrefs_maps[["by_entrezgene"]] <- as.data.frame(
    transcript_df |>
      dplyr::filter(is.na(hgnc_id) &
        !is.na(entrezgene) &
        !stringr::str_detect(entrezgene, "&")) |>
      dplyr::select(-hgnc_id) |>
      dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
      dplyr::left_join(
        dplyr::select(
          gene_info,
          hgnc_id,
          entrezgene,
          name
        ),
        by = "entrezgene", multiple = "all"
      )
  )

  ## map gene cross-references by symbol identifier
  ## given that hgnc id and entrezgene is not provided
  gene_xrefs_maps[["by_symbol"]] <- as.data.frame(
    transcript_df |>
      dplyr::filter(is.na(hgnc_id) &
        is.na(entrezgene)) |>
      dplyr::select(-c(entrezgene, hgnc_id)) |>
      dplyr::left_join(
        dplyr::select(
          gene_info,
          hgnc_id,
          entrezgene,
          symbol,
          name
        ),
        by = "symbol", multiple = "all", 
        relationship = "many-to-many"
      )
  )

  gene_xrefs_maps[["remain"]] <- as.data.frame(
    transcript_df |>
      dplyr::filter(is.na(hgnc_id) &
        !is.na(entrezgene) &
        stringr::str_detect(entrezgene, "&")) |>
      dplyr::mutate(name = NA)
  )


  gencode_transcripts_xref <-
    do.call(rbind, gene_xrefs_maps) |>
    dplyr::mutate(entrezgene = as.integer(entrezgene)) |>
    dplyr::mutate(name = dplyr::if_else(
      is.na(name) & !is.na(description),
      as.character(description),
      as.character(name)
    )) |>
    dplyr::select(-description)

  ## Map remaining gene cross-refs against unambiguous gene alias
  unambiguous_aliases <-
    gene_alias$records |>
    dplyr::filter(ambiguous == FALSE) |>
    dplyr::select(-c(symbol, n_primary_map, ambiguous)) |>
    dplyr::rename(symbol = alias) |>
    dplyr::inner_join(
      dplyr::select(gene_info, entrezgene, hgnc_id, name),
      by = "entrezgene", multiple = "all"
    )

  gene_xrefs_maps[["by_alias"]] <- gencode_transcripts_xref |>
    dplyr::filter(is.na(entrezgene) & is.na(name) & !is.na(symbol)) |>
    dplyr::select(-c(entrezgene, name, hgnc_id)) |>
    dplyr::left_join(unambiguous_aliases, by = "symbol", multiple = "all")

  gencode_transcripts_xref_final <- as.data.frame(
    dplyr::anti_join(
      gencode_transcripts_xref,
      gene_xrefs_maps[["by_alias"]],
      by = c("ensembl_gene_id", "ensembl_transcript_id", "symbol")
    ) |>
      dplyr::bind_rows(gene_xrefs_maps[["by_alias"]])
  )

  rownames(gencode_transcripts_xref_final) <- NULL
  attr(gencode_transcripts_xref_final, "groups") <- NULL

  return(gencode_transcripts_xref_final)
}

get_uniprot_map <- function(uniprot_version = "2024_04") {
  lgr::lgr$info("Retrieving UniProtKB annotation")
  withr::local_options(
    timeout = max(30000000, getOption("timeout")))

  remote_idmapping_dat_fname <-
    paste0(
      "https://ftp.ebi.ac.uk/pub/databases/uniprot/",
      "current_release/knowledgebase/idmapping/by_organism/",
      "HUMAN_9606_idmapping.dat.gz"
    )

  idmapping_up_kb <- readr::read_tsv(
    file = remote_idmapping_dat_fname,
    col_names = FALSE, quote = "", 
    show_col_types = FALSE
  ) |>
    dplyr::filter(X2 == "Ensembl_TRS" |
      X2 == "UniProtKB-ID" |
      X2 == "Ensembl" |
      X2 == "GeneID" |
      X2 == "RefSeq_NT") |>
    magrittr::set_colnames(c("acc", "type", "name")) |>
    dplyr::mutate(acc = stringr::str_replace(acc, "-[0-9]{1,}", "")) |>
    dplyr::filter(nchar(acc) == 6)


  ## UniProt accession to Ensembl transcript ID
  ensembl_up_acc <- dplyr::filter(
    idmapping_up_kb, type == "Ensembl_TRS") |>
    dplyr::select(acc, name) |>
    dplyr::rename(uniprot_acc = acc, ensembl_transcript_id = name) |>
    dplyr::distinct()
  ## remove transcripts mapped to more than one UniProt acc
  unique_ensembl_trans <- as.data.frame(
    dplyr::group_by(ensembl_up_acc, ensembl_transcript_id) |>
      dplyr::summarise(n_trans = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n_trans == 1)
  )
  ensembl_up_acc <- dplyr::semi_join(
    ensembl_up_acc,
    dplyr::select(
      unique_ensembl_trans, 
      ensembl_transcript_id),
    by = c("ensembl_transcript_id")
  )

  ## UniProt accession to UniProt ID
  up_id_acc <- as.data.frame(
    dplyr::filter(idmapping_up_kb, type == "UniProtKB-ID") |>
      dplyr::group_by(acc) |>
      dplyr::summarise(
        uniprot_id = paste(unique(name), collapse = "&"),
        .groups = "drop"
      ) |>
      dplyr::rename(uniprot_acc = acc)
  )
  lgr::lgr$info(
    paste0(
      "A total of ",
      nrow(up_id_acc), 
      " UniProt protein accessions were parsed"
    )
  )


  ## UniProt accession to refseq mRNA
  refseq_id_acc <- as.data.frame(
    dplyr::filter(idmapping_up_kb, type == "RefSeq_NT") |>
      dplyr::mutate(name = stringr::str_replace(
        name, "\\.[0-9]{1,}$", ""
      )) |>
      dplyr::group_by(acc) |>
      dplyr::summarise(
        refseq_transcript_id = paste(unique(name), collapse = "&"),
        .groups = "drop"
      ) |>
      dplyr::rename(uniprot_acc = acc)
  )


  uniprot_map <-
    dplyr::full_join(
      dplyr::left_join(
        refseq_id_acc, up_id_acc, 
        by = c("uniprot_acc"), multiple = "all"),
      dplyr::left_join(
        ensembl_up_acc, up_id_acc, 
        by = c("uniprot_acc"), multiple = "all"),
      by = c("uniprot_acc", "uniprot_id"), multiple = "all"
    ) |>
    dplyr::filter(!is.na(uniprot_id)) |>
    dplyr::mutate(ensembl_transcript_id_full = ensembl_transcript_id) |>
    dplyr::mutate(ensembl_transcript_id = stringr::str_replace(
      ensembl_transcript_id, "\\.[0-9]{1,}$", ""
    )) |>
    dplyr::mutate(
      uniprot_version = uniprot_version
    ) |>
    dplyr::select(
      uniprot_acc,
      uniprot_id,
      uniprot_version,
      ensembl_transcript_id
    ) |>
    dplyr::filter(!is.na(ensembl_transcript_id))


  return(uniprot_map)
}

sort_bed_regions <- function(unsorted_regions) {
  sorted_regions <- NULL
  all_regions <- NULL
  assertable::assert_colnames(
    unsorted_regions, c("start", "end"),
    only_colnames = FALSE, quiet = TRUE
  )
  if ("chrom" %in% colnames(unsorted_regions) &
    "start" %in% colnames(unsorted_regions) &
    "end" %in% colnames(unsorted_regions)) {
    chrOrder <- paste0(
      "chr", c(as.character(c(1:22)), "X", "Y", "M"))
   
    scaffold_regions <- unsorted_regions |>
      dplyr::filter(!(chrom %in% chrOrder))
    
    chrom_unsorted_regions <- unsorted_regions |>
      dplyr::filter(chrom %in% chrOrder)
    
    chrom_unsorted_regions$chrom <- factor(
      chrom_unsorted_regions$chrom,
      levels = chrOrder
    )
    chrom_unsorted_regions <- chrom_unsorted_regions[order(chrom_unsorted_regions$chrom), ]

    sorted_regions <- data.frame()
    for (chrom in chrOrder) {
      if (nrow(
        chrom_unsorted_regions[chrom_unsorted_regions$chrom == chrom, ]
      ) > 0) {
        chrom_regions <-
          chrom_unsorted_regions[chrom_unsorted_regions$chrom == chrom, ]
        chrom_regions_sorted <-
          chrom_regions[with(chrom_regions, order(start, end)), ]
        sorted_regions <- dplyr::bind_rows(
          sorted_regions, chrom_regions_sorted
        )
      }
    }

    sorted_regions$start <- as.integer(sorted_regions$start)
    sorted_regions$end <- as.integer(sorted_regions$end)
    
    all_regions <- sorted_regions |>
      dplyr::bind_rows(scaffold_regions)
  }
  
  
  
  return(all_regions)
}
