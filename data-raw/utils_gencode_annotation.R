
ensembl_cache_dir <- function() {
  ## Relative to the working directory, same convention as
  ## data-raw/gd_local/ elsewhere in this script - not here::here(),
  ## which caches the project root from whatever the working directory
  ## was when the 'here' package was first loaded in the session, and
  ## can silently drift from the actual working directory afterwards.
  cache_dir <- file.path("data-raw", "ensembl_cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  cache_dir
}

## Downloads 'remote_url' into data-raw/ensembl_cache/ (gitignored, kept
## separate from data-raw/gd_local/ which holds the final files released
## with geneOncoX) and returns the local path, skipping the download if
## a cached copy already exists. GENCODE/Ensembl URLs are release-pinned,
## so a cached file's content never goes stale for a given release.
download_cached <- function(remote_url, filename = basename(remote_url)) {
  local_path <- file.path(ensembl_cache_dir(), filename)
  if (!file.exists(local_path)) {
    utils::download.file(
      remote_url, destfile = local_path, mode = "wb", quiet = TRUE)
  }
  local_path
}

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
      file = download_cached(
        remote_appris,
        filename = paste0("appris_data.principal.forENSEMBL_", build, ".txt")
      ),
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


gencode_get_exons <- function(
    build = "grch38"){
  
  chrom_names <- NULL
  if(build == "grch38"){
    suppressMessages(
      library(BSgenome.Hsapiens.UCSC.hg38))
    chrom_names <- 
      head(
        BSgenome::seqnames(
          GenomeInfoDb::seqinfo(
            x = BSgenome.Hsapiens.UCSC.hg38)), 24)
  }
  if(build == "grch37"){
    suppressMessages(
      library(BSgenome.Hsapiens.UCSC.hg19))
    chrom_names <- 
      head(
        BSgenome::seqnames(
          GenomeInfoDb::seqinfo(
            x = BSgenome.Hsapiens.UCSC.hg19)), 24)
  }
  
  
  
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
  gencode_gtf_all <- read_gtf(remote_gtf_all)
  
  canonical_exons <- gencode_gtf_all |> 
    dplyr::filter(
      type == "exon" & 
        gene_type == "protein_coding" & 
        transcript_type == "protein_coding") |>
    dplyr::filter(
      transcript_type == "protein_coding",
      tag %in% c("GENCODE_Primary", "CCDS", "basic")
    ) |>
    dplyr::mutate(
      exon_number = as.integer(exon_number)
    )
  
  # Per-gene transcript count — needed to correctly label exons in
  # genes with only one canonical transcript (all exons => constitutive)
  transcripts_per_gene <- canonical_exons |>
    dplyr::mutate(
      gene_id = stringr::str_replace(gene_id, "\\.[0-9]+$", ""),
    ) |>
    dplyr::group_by(gene_name, gene_id, strand) |>
    dplyr::summarise(
      n_canonical_transcripts = dplyr::n_distinct(transcript_id),
      .groups = "drop"
    )

  all_canonical_protein_coding_exons <- canonical_exons |>
    dplyr::filter(chrom %in% chrom_names) |>
    dplyr::mutate(
      gene_id = stringr::str_replace(gene_id, "\\.[0-9]+$", ""),
    ) |>
    dplyr::group_by(chrom, gene_name, gene_id, exon_number, strand) |>
    dplyr::summarise(
      start              = min(start),
      end                = max(end),
      n_transcripts_with_exon = dplyr::n_distinct(transcript_id),
      exon_id                 = paste(unique(exon_id), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      transcripts_per_gene, 
      by = c("gene_name","gene_id","strand")) |>
    dplyr::arrange(gene_name, exon_number) |>
    dplyr::mutate(
      exon_type = dplyr::case_when(
        n_canonical_transcripts == 1     ~ "constitutive",
        n_transcripts_with_exon >= 2 ~ "constitutive",
        .default                     = "alternative"
      )
    ) |>
    dplyr::mutate(
      name = paste(
        gene_name, gene_id, 
        n_canonical_transcripts,
        n_transcripts_with_exon,
        paste("exon", exon_number, sep = "_"),
        exon_type,
        sep = "|"
      )) |>
    dplyr::mutate(
      score = 1000
    ) |>
    dplyr::select(
      chrom, 
      start, 
      end,
      name,
      score,
      strand
    )
}

gencode_get_transcripts <-
  function(build = "grch38",
           append_regulatory_region = TRUE,
           gencode_version = 50,
           ensembl_version = 116,
           uniprot_version = "2026_04",
           gene_info = NULL,
           gene_alias = NULL) {
    gencode_ftp_url <-
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
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
    gencode_gtf_all <- read_gtf(remote_gtf_all)
    
    cds_start_positions <- gencode_gtf_all |>
      dplyr::filter(.data$type == "start_codon" & phase == 0) |>
      dplyr::mutate(cds_start = dplyr::if_else(
        .data$strand == "-", 
        as.numeric(.data$end),
        as.numeric(.data$start)
      )) |>
      dplyr::select(c("transcript_id", "cds_start")) |>
      dplyr::rename(
        ensembl_transcript_id_full = "transcript_id") |>
      dplyr::distinct()
    
    gencode_gtf_transcripts <- gencode_gtf_all |>
      dplyr::filter(type == "transcript")
    
    if(gencode_version == 19){
      gencode_gtf_transcripts_chr <- read_gtf(remote_gtf) |>
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
        ensembl_gene_id_full = gene_id,
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
        ensembl_gene_id_full,
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
            ensembl_gene_id_full,
            "\\.[0-9]{1,}$", ""
          )
      ) |>
      dplyr::distinct()

    ## Temporary fix: Need additional processing to get multiple tag-value
    ## pairs in a proper fashion. data.table::fread() can read straight
    ## from a URL, but does its own download each call, bypassing the
    ## read_gtf() cache above - go through download_cached() so this
    ## reuses the same local copy instead of re-fetching the 100+MB file.
    gencode_gtf_fix <- as.data.frame(
      data.table::fread(download_cached(remote_gtf_all), skip = 5, verbose = F) |>
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
      dplyr::distinct() |>
      dplyr::left_join(
        cds_start_positions, 
        by = "ensembl_transcript_id_full"
      )

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
      ## pmax()/pmin() instead of rowwise()+max()/min() - vectorized,
      ## same result, ~50x faster at this row count (rowwise() evaluates
      ## one row at a time).
      gencode <- as.data.frame(
        gencode |>
          dplyr::left_join(chromosome_lengths,
                           by = c("chrom"), multiple = "all") |>
          dplyr::mutate(
            start = dplyr::if_else(
              startsWith(chrom,"chr") &
                chrom != "chrM",
              as.integer(
                pmax(1, transcript_start - 5001)
              ),
              as.integer(transcript_start)
            ),
            end = dplyr::if_else(
              startsWith(chrom,"chr") &
                chrom != "chrM",
              as.integer(
                pmin(chrom_length, transcript_end + 5001)
              ),
              as.integer(transcript_end)
            )
          ) |>
          dplyr::filter(!is.na(chrom_length)) |>
          dplyr::select(-chrom_length)
      )
    }

    attr(gencode, "groups") <- NULL

    transcript_df <- gencode |>
      dplyr::select(
        ensembl_transcript_id,
        ensembl_transcript_id_full,
        ensembl_gene_id_full,
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
          "ensembl_gene_id",
          "ensembl_transcript_id",
          "ensembl_gene_id_full",
          "ensembl_transcript_id_full"
        ), relationship = "one-to-one"
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
        cds_start,
        strand,
        ensembl_gene_id,
        ensembl_gene_id_full,
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
        refseq_select,
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



get_ensembl_gene_descriptions <- function(ensembl_version = 116,
                                          build = "grch38") {
  if (build == "grch37") {
    ## GRCh37 is frozen/legacy and lives on a separate archive site
    ## (ftp.ensembl.org/pub/grch37/...), not under /pub/release-<N>/ -
    ## 'ensembl_version' doesn't apply there, and the archive's own
    ## release number can change over time, so discover the current
    ## GRCh37 GFF3 filename from its directory listing rather than
    ## hardcoding a release.
    index_url <-
      "https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/"
    listing <- readLines(index_url, warn = FALSE)
    fname <- stringr::str_extract(
      listing,
      "Homo_sapiens\\.GRCh37\\.[0-9]+\\.chr_patch_hapl_scaff\\.gff3\\.gz"
    )
    fname <- fname[!is.na(fname)][1]
    remote_url <- paste0(index_url, fname)
  } else {
    remote_url <- paste0(
      "https://ftp.ensembl.org/pub/release-", ensembl_version,
      "/gff3/homo_sapiens/Homo_sapiens.GRCh38.",
      ensembl_version, ".chr_patch_hapl_scaff.gff3.gz"
    )
  }

  local_gz <- download_cached(remote_url)

  ## Pull only 'gene' feature rows (column 3) before loading into R -
  ## this file also carries every exon/CDS/UTR record, which we don't
  ## need here.
  gene_lines <- data.table::fread(
    cmd = paste0(
      "gzip -dc ", shQuote(local_gz), " | awk -F'\\t' '$3==\"gene\"'"
    ),
    header = FALSE, sep = "\t", quote = "", data.table = FALSE
  )

  gene_lines |>
    dplyr::transmute(
      ensembl_gene_id = stringr::str_match(
        V9, "gene_id=([^;]+)")[, 2],
      description = stringr::str_match(
        V9, "description=([^;]+)")[, 2]
    ) |>
    dplyr::filter(!is.na(description)) |>
    dplyr::mutate(
      description = vapply(description, utils::URLdecode, character(1))
    ) |>
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
}

read_ensembl_tsv_xref <- function(ensembl_version = 116,
                                  build = "grch38",
                                  xref_type = "refseq") {
  assembly <- if (build == "grch37") "GRCh37" else "GRCh38"
  remote_url <- paste0(
    "https://ftp.ensembl.org/pub/release-", ensembl_version,
    "/tsv/homo_sapiens/Homo_sapiens.", assembly, ".",
    ensembl_version, ".", xref_type, ".tsv.gz"
  )
  readr::read_tsv(
    download_cached(remote_url), show_col_types = FALSE, progress = FALSE,
    col_types = readr::cols(.default = "c")
  )
}

## GRCh37's transcript-level RefSeq/UniProt cross-references have no
## file-based equivalent to Ensembl's per-release TSV dumps used for
## GRCh38 - the GRCh37 archive (ftp.ensembl.org/pub/grch37/) only
## publishes 'ena', 'karyotype' and 'uniprot' TSVs, not 'refseq' or
## 'entrez'. grch37.ensembl.org's live BioMart service is a dedicated,
## standalone host though - it was never part of the www/asia/useast
## mirror pool that was retiring/failing, so it's used here instead,
## scoped only to build == "grch37".
get_grch37_biomart_xrefs <- function() {
  ensembl_mart_grch37 <- suppressWarnings(
    biomaRt::useEnsembl(
      biomart = "genes",
      GRCh = "37",
      dataset = "hsapiens_gene_ensembl"
    )
  )

  xref1 <- suppressWarnings(biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "refseq_mrna", "ensembl_gene_id"),
    mart = ensembl_mart_grch37
  )) |>
    dplyr::mutate(refseq_mrna = dplyr::if_else(
      nchar(refseq_mrna) == 0, NA_character_, refseq_mrna)) |>
    dplyr::filter(!is.na(refseq_mrna)) |>
    dplyr::distinct()

  refseq_mrna_map <- xref1 |>
    dplyr::select(ensembl_gene_id, ensembl_transcript_id, refseq_mrna) |>
    dplyr::distinct()

  xref2 <- suppressWarnings(biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "uniprotswissprot", "refseq_peptide"),
    mart = ensembl_mart_grch37
  )) |>
    dplyr::mutate(
      uniprotswissprot = dplyr::if_else(
        nchar(uniprotswissprot) == 0, NA_character_, uniprotswissprot),
      refseq_peptide = dplyr::if_else(
        nchar(refseq_peptide) == 0, NA_character_, refseq_peptide)
    )

  refseq_peptide_map <- xref2 |>
    dplyr::filter(!is.na(refseq_peptide)) |>
    dplyr::select(ensembl_transcript_id, refseq_peptide) |>
    dplyr::distinct()

  uniprot_swissprot_map <- xref2 |>
    dplyr::filter(!is.na(uniprotswissprot)) |>
    dplyr::transmute(ensembl_transcript_id, uniprot_acc = uniprotswissprot) |>
    dplyr::distinct()

  xref4 <- suppressWarnings(biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "refseq_ncrna"),
    mart = ensembl_mart_grch37
  )) |>
    dplyr::mutate(refseq_ncrna = dplyr::if_else(
      nchar(refseq_ncrna) == 0, NA_character_, refseq_ncrna))

  refseq_ncrna_map <- xref4 |>
    dplyr::filter(!is.na(refseq_ncrna)) |>
    dplyr::select(ensembl_transcript_id, refseq_ncrna) |>
    dplyr::distinct()

  list(
    refseq_mrna_map = refseq_mrna_map,
    refseq_peptide_map = refseq_peptide_map,
    refseq_ncrna_map = refseq_ncrna_map,
    uniprot_swissprot_map = uniprot_swissprot_map
  )
}

gencode_resolve_xrefs <- function(transcript_df = NULL,
                                  build = "grch38",
                                  ensembl_version = 114,
                                  gene_info = NULL,
                                  gene_alias = NULL) {
  invisible(assertable::assert_colnames(
    transcript_df,
    colnames = c(
      "ensembl_transcript_id",
      "ensembl_gene_id",
      "symbol"
    ),
    quiet = TRUE,
    only_colnames = F
  ))

  options(timeout = 10000000)

  refseq_data <- get_refseq_map(build = build)

  ## HGNC id and Entrez gene id come from NCBI's gene_info file (already
  ## supplied as 'gene_info'). Gene description text is taken primarily
  ## from Ensembl's own per-release GFF3 annotation (same underlying
  ## field BioMart's 'description' attribute used to surface, just read
  ## from a static file), falling back to gene_info's NCBI-sourced name
  ## for genes Ensembl's GFF3 doesn't carry a description for.
  ensembl_gene_descriptions <- get_ensembl_gene_descriptions(
    ensembl_version = ensembl_version, build = build)

  ensembl2entrez_hgnc <- gene_info |>
    dplyr::filter(!is.na(ensembl_gene_id)) |>
    dplyr::select(
      ensembl_gene_id, hgnc_id, entrezgene,
      description_ncbi = name
    ) |>
    dplyr::distinct() |>
    dplyr::full_join(
      ensembl_gene_descriptions, by = "ensembl_gene_id") |>
    dplyr::mutate(
      description = dplyr::if_else(
        is.na(description), description_ncbi, description
      )
    ) |>
    dplyr::select(-description_ncbi) |>
    dplyr::distinct()

  ## Transcript-level cross-references (RefSeq mRNA/ncRNA/peptide,
  ## UniProt SwissProt accession). For GRCh38 these come from Ensembl's
  ## static, per-release TSV xref dumps
  ## (ftp.ensembl.org/pub/release-<v>/tsv/homo_sapiens/) instead of the
  ## live BioMart web service - a plain HTTPS file download, unaffected
  ## by BioMart's mirror outages/retirement. GRCh37 has no such TSV
  ## dumps at all, so it uses grch37.ensembl.org's BioMart instead (see
  ## get_grch37_biomart_xrefs()).
  if (build == "grch37") {
    grch37_xrefs <- get_grch37_biomart_xrefs()
    refseq_mrna_map <- grch37_xrefs$refseq_mrna_map
    refseq_peptide_map <- grch37_xrefs$refseq_peptide_map
    refseq_ncrna_map <- grch37_xrefs$refseq_ncrna_map
    uniprot_swissprot_map <- grch37_xrefs$uniprot_swissprot_map
  } else {
    xref_refseq <- read_ensembl_tsv_xref(
      ensembl_version = ensembl_version, build = build, xref_type = "refseq")
    xref_uniprot <- read_ensembl_tsv_xref(
      ensembl_version = ensembl_version, build = build, xref_type = "uniprot")

    refseq_mrna_map <- xref_refseq |>
      dplyr::filter(db_name == "RefSeq_mRNA") |>
      dplyr::transmute(
        ensembl_transcript_id = transcript_stable_id,
        refseq_mrna = xref
      ) |>
      dplyr::distinct()

    refseq_peptide_map <- xref_refseq |>
      dplyr::filter(db_name == "RefSeq_peptide") |>
      dplyr::transmute(
        ensembl_transcript_id = transcript_stable_id,
        refseq_peptide = xref
      ) |>
      dplyr::distinct()

    refseq_ncrna_map <- xref_refseq |>
      dplyr::filter(db_name == "RefSeq_ncRNA") |>
      dplyr::transmute(
        ensembl_transcript_id = transcript_stable_id,
        refseq_ncrna = xref
      ) |>
      dplyr::distinct()

    uniprot_swissprot_map <- xref_uniprot |>
      dplyr::filter(db_name == "Uniprot/SWISSPROT") |>
      dplyr::transmute(
        ensembl_transcript_id = transcript_stable_id,
        uniprot_acc = xref
      ) |>
      dplyr::distinct()
  }

  ## 'transcript_df' already covers every GENCODE/Ensembl transcript for
  ## this assembly, so use its (gene, transcript) pairs as the backbone -
  ## matching what an unfiltered BioMart getBM() call used to return -
  ## rather than the (much smaller) set of transcripts that happen to
  ## have a RefSeq match.
  xref_biomart_1 <- transcript_df |>
    dplyr::select(ensembl_gene_id, ensembl_transcript_id) |>
    dplyr::distinct() |>
    dplyr::left_join(
      refseq_mrna_map, by = "ensembl_transcript_id") |>
    dplyr::left_join(
      ensembl2entrez_hgnc, by = "ensembl_gene_id",
      relationship = "many-to-many") |>
    dplyr::left_join(
       refseq_data,
       by = c("refseq_mrna" = "rs_transcript_id_stripped")) |>
    dplyr::mutate(
       refseq_transcript_id = dplyr::if_else(
         is.na(rs_transcript_id) & !is.na(refseq_mrna),
         refseq_mrna,
         rs_transcript_id
       )
    ) |>
    dplyr::select(-c("rs_transcript_id", "refseq_mrna")) |>
    dplyr::distinct()

  xref_biomart <- xref_biomart_1 |>
    dplyr::left_join(
      refseq_peptide_map,
      by = "ensembl_transcript_id",
      multiple = "all",
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      uniprot_swissprot_map,
      by = "ensembl_transcript_id",
      multiple = "all",
      relationship = "many-to-many"
    ) |>
    dplyr::left_join(
      refseq_ncrna_map,
      by = "ensembl_transcript_id",
      multiple = "all",
      relationship = "many-to-many"
    ) |>
    dplyr::distinct() |>
    dplyr::rename(
      refseq_protein_id = refseq_peptide
    ) |>
    dplyr::mutate(refseq_transcript_id = dplyr::if_else(
      is.na(refseq_transcript_id) & !is.na(refseq_ncrna),
      as.character(refseq_ncrna),
      as.character(refseq_transcript_id)
    )) |>
    dplyr::select(-refseq_ncrna) |>
    dplyr::distinct()

  valid_entrez_recs <- xref_biomart |>
    dplyr::filter(!is.na(entrezgene) & 
                    !is.na(hgnc_id)) |>
    dplyr::distinct()
  
  valid_other_recs <- xref_biomart |> 
    dplyr::filter(is.na(entrezgene)) |> 
    dplyr::distinct()
  
  ## A gene with an Entrez id but no HGNC id (common for very recently
  ## catalogued "uncharacterized LOC*" loci) has its entrezgene nulled
  ## out below as untrustworthy - but the row itself (and its transcript-
  ## level RefSeq/UniProt columns) is kept regardless of whether another
  ## gene happens to share that Entrez id. Dropping unmatched rows here
  ## used to silently delete their RefSeq/UniProt accessions too, which
  ## has nothing to do with Entrez/HGNC trustworthiness.
  invalid_entrez_recs <- xref_biomart |>
    dplyr::filter(!is.na(entrezgene) & is.na(hgnc_id)) |>
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
    "refseq_select",
    "hgnc_id",
    "transcript_mane_select",
    "transcript_mane_plus_clinical",
    "description")) {
    if(n %in% colnames(xref_biomart)){
      if(NROW(
        xref_biomart[!is.na(xref_biomart[, n]) &
                     (xref_biomart[, n] == "" |
                      xref_biomart[, n] == "NA"), ]) > 0){
        xref_biomart[!is.na(xref_biomart[, n]) &
                       (xref_biomart[, n] == "" |
                          xref_biomart[, n] == "NA"), ][, n] <- NA
      }
    }
  }
  
  for (xref in c(
    "refseq_protein_id",
    "refseq_transcript_id",
    "uniprot_acc",
    "transcript_mane_select",
    "transcript_mane_plus_clinical",
    "refseq_select"
    #"description"
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
          !!rlang::sym(xref) := paste(
            unique(sort(!!rlang::sym(xref), decreasing = T)),
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
    "hgnc_id",
    "description"
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
          !!rlang::sym(xref) := paste(
            unique(sort(!!rlang::sym(xref), decreasing = T)),
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

  transcript_df <- transcript_df |>
    dplyr::mutate(
      transcript_mane_select = stringr::str_replace(
        .data$transcript_mane_select,"^NA&|&NA$",""
      )
    ) |>
    dplyr::mutate(
      transcript_mane_plus_clinical = stringr::str_replace(
        .data$transcript_mane_plus_clinical,"^NA&|&NA$",""
      )
    ) |>
    dplyr::mutate(
      refseq_select = stringr::str_replace(
        .data$refseq_select,"^NA&|&NA$",""
      )
    )
  

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
  
  
  
  for (n in c(
    "refseq_protein_id",
    "refseq_transcript_id",
    "uniprot_acc",
    "refseq_select",
    "hgnc_id",
    "transcript_mane_select",
    "transcript_mane_plus_clinical",
    "name")) {
    if(n %in% colnames(gencode_transcripts_xref_final)){
      if(NROW(
        gencode_transcripts_xref_final[!is.na(gencode_transcripts_xref_final[, n]) &
                     (gencode_transcripts_xref_final[, n] == "" |
                      gencode_transcripts_xref_final[, n] == "NA"), ]) > 0){
        gencode_transcripts_xref_final[!is.na(gencode_transcripts_xref_final[, n]) &
                       (gencode_transcripts_xref_final[, n] == "" |
                          gencode_transcripts_xref_final[, n] == "NA"), ][, n] <- NA
      }
    }
  }
  
  gene_info_complete1 <- gene_info |>
    dplyr::filter(!is.na(entrezgene) & !is.na(ensembl_gene_id)) |>
    #dplyr::mutate(entrezgene = as.character(entrezgene)) |>
    dplyr::select(c("hgnc_id","entrezgene",
                     "symbol","ensembl_gene_id")) |>
    dplyr::rename(symbol2 = symbol)
  
  #gene_info_complete2 <- gene_info |>
  #dplyr::filter(!is.na(entrezgene) & is.na(ensembl_gene_id)) |>
  #   dplyr::mutate(entrezgene = as.character(entrezgene)) |>
  #   dplyr::select(c("hgnc_id","entrezgene",
  #                    "symbol"))
  
  missing_entrez1 <- gencode_transcripts_xref_final |>
    dplyr::filter(is.na(entrezgene)) |>
    dplyr::select(-c("entrezgene","hgnc_id")) |>
    dplyr::left_join(
      dplyr::select(
        gene_info_complete1,
        c("hgnc_id","entrezgene","symbol2","ensembl_gene_id")),
      by = c("ensembl_gene_id")
    ) |>
    dplyr::mutate(symbol = dplyr::if_else(
      !is.na(symbol2) &
        stringr::str_detect(symbol,"^ENSG"),
      symbol2,
      as.character(symbol)
    )) |>
    dplyr::select(-c("symbol2"))
  
  not_missing_entrez <- gencode_transcripts_xref_final |>
    dplyr::filter(!is.na(entrezgene)) |>
    dplyr::left_join(
      dplyr::select(
        gene_info_complete1,
        c("entrezgene","symbol2")),
      by = c("entrezgene")
    ) |>
    dplyr::mutate(symbol = dplyr::if_else(
      !is.na(symbol2) &
        stringr::str_detect(symbol,"^ENSG"),
      symbol2,
      as.character(symbol)
    )) |>
    dplyr::select(-c("symbol2"))
  
  gencode_transcripts_xref_final <- dplyr::bind_rows(
    missing_entrez1,not_missing_entrez) |>
    dplyr::arrange(entrezgene)
  
  
  set1 <- gencode_transcripts_xref_final |>
    dplyr::filter(!is.na(entrezgene) & is.na(name)) |>
    dplyr::select(-c("name")) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        c("entrezgene","name")),
      by = c("entrezgene")
    )
  
  set1 <- gencode_transcripts_xref_final |>
    dplyr::filter(!is.na(entrezgene))
  set2 <- gencode_transcripts_xref_final |>
    dplyr::filter(is.na(entrezgene) & 
                    !is.na(name) & !is.na(symbol)) |>
    dplyr::select(-c("name")) |>
    dplyr::left_join(
      dplyr::select(
        gene_info,
        c("symbol","name")),
      by = c("symbol" = "symbol")
    )
  
  set3 <- gencode_transcripts_xref_final |> dplyr::filter(
    is.na(entrezgene) & is.na(name))
  
  gencode_transcripts_xref_final <- 
    dplyr::bind_rows(set1,set2,set3) |>
    dplyr::distinct()
  
  gene_info_valid <- 
    gene_info |> 
    dplyr::select(
      c("entrezgene", "hgnc_id", 
        "symbol")) |> 
    dplyr::filter(!is.na(entrezgene))
  
  gencode_missing_entrez_hgnc <- 
    gencode_transcripts_xref_final |>
    dplyr::filter(is.na(entrezgene) &
                    is.na(hgnc_id)) |>
    dplyr::select(
      -c("entrezgene","hgnc_id")) |>
    dplyr::left_join(
      gene_info_valid, by = "symbol"
    )
  
  
  gencode_nonmissing_entrez_hgnc <- 
    gencode_transcripts_xref_final |>
    dplyr::filter(!(is.na(entrezgene) &
                      is.na(hgnc_id)))
  
  gencode_transcripts_xref_final <- 
    dplyr::bind_rows(
      gencode_missing_entrez_hgnc,
      gencode_nonmissing_entrez_hgnc
    ) |>
    dplyr::arrange(entrezgene)
  
  
  
  return(gencode_transcripts_xref_final)
}

get_uniprot_map <- function(uniprot_version = "2026_04") {
  lgr::lgr$info("Retrieving UniProtKB annotation")
  withr::local_options(
    timeout = max(30000000, getOption("timeout")))

  remote_idmapping_dat_fname <-
    paste0(
      "https://ftp.ebi.ac.uk/pub/databases/uniprot/",
      "current_release/knowledgebase/idmapping/by_organism/",
      "HUMAN_9606_idmapping.dat.gz"
    )

  ## readr::read_tsv() on a remote URL goes through vroom/curl, which
  ## has its own connection timeout independent of options(timeout=) -
  ## that option only applies to utils::download.file(). Downloading
  ## first (via download_cached(), which uses download.file()) avoids
  ## the ~10s curl connect-timeout this large file was hitting, and
  ## caches the file for re-reads.
  idmapping_up_kb <- readr::read_tsv(
    file = download_cached(remote_idmapping_dat_fname),
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


get_refseq_map <- function(build = "grch37"){
  
  refseq_base_url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation"
  remote_urls <- list()
  remote_urls[['grch37']] <- 
    file.path(refseq_base_url, 
              "GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz")
  remote_urls[['grch38']] <- 
    file.path(refseq_base_url, 
              "GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gtf.gz")
  
  
  mane_select_map <- data.frame()
  mane_plus_clinical_map <- data.frame()
  refseq_map <- data.frame()
  
  refseq_map_grch38 <- read_gtf(remote_urls[['grch38']]) |> 
    dplyr::filter(type == "transcript" & transcript_biotype == "mRNA") |>
    dplyr::rename(rs_transcript_id = transcript_id) |>
    dplyr::filter(
      !stringr::str_detect(rs_transcript_id, "\\.[0-9]{1,}_[0-9]$")
    ) |>
    dplyr::mutate(
      rs_transcript_id_stripped = stringr::str_replace(
        .data$rs_transcript_id, "\\.[0-9]{1,}$", ""
      )
    ) |>
    dplyr::filter(stringr::str_detect(
      rs_transcript_id_stripped, "^NM_"
    )) |>
    dplyr::mutate(
      transcript_mane_select = dplyr::if_else(
        tag == "MANE Select", 
        rs_transcript_id, 
        as.character(NA)
      )
    ) |>
    dplyr::mutate(
      transcript_mane_plus_clinical = dplyr::if_else(
        tag == "MANE Plus Clinical", 
        rs_transcript_id, 
        as.character(NA)
      )
    ) |>
    dplyr::select(
      c("rs_transcript_id", 
        "rs_transcript_id_stripped",
        "transcript_mane_select", 
        "transcript_mane_plus_clinical"
      )) |>
    dplyr::mutate(refseq_select = NA)
  
  mane_select_map <- refseq_map_grch38 |>
    dplyr::filter(!is.na(transcript_mane_select)) |>
    dplyr::select(transcript_mane_select) |>
    dplyr::mutate(mane_select = TRUE) |>
    dplyr::distinct()
  
  mane_plus_clinical_map <- refseq_map_grch38 |>
    dplyr::filter(!is.na(transcript_mane_plus_clinical)) |>
    dplyr::select(transcript_mane_plus_clinical) |>
    dplyr::mutate(mane_plus_clinical = TRUE) |>
    dplyr::distinct()
  
  refseq_map <- refseq_map_grch38
  
  
  if(build == "grch37"){
    refseq_map_grch37 <- read_gtf(remote_urls[['grch37']]) |> 
      dplyr::filter(.data$type == "transcript" & 
                      .data$transcript_biotype == "mRNA") |>
      dplyr::rename(rs_transcript_id = transcript_id,
                    rs_select = tag) |>
      dplyr::filter(
        !stringr::str_detect(
          .data$rs_transcript_id, "\\.[0-9]{1,}_[0-9]$")
      ) |>
      dplyr::mutate(refseq_select = dplyr::if_else(
        !is.na(.data$rs_select), TRUE, FALSE
      )) |>
      dplyr::mutate(
        rs_transcript_id_stripped = stringr::str_replace(
          .data$rs_transcript_id, "\\.[0-9]{1,}$", ""
        )
      ) |>
      dplyr::select(
        c("rs_transcript_id", 
          "rs_transcript_id_stripped",
          "refseq_select")) |>
      dplyr::left_join(
        mane_select_map,
        c("rs_transcript_id" = "transcript_mane_select")
      ) |>
      dplyr::left_join(
        mane_plus_clinical_map,
        c("rs_transcript_id" = "transcript_mane_plus_clinical")
      ) |>
      dplyr::mutate(
        transcript_mane_select = dplyr::if_else(
          !is.na(mane_select), 
          rs_transcript_id, 
          as.character(NA)
        )
      ) |>
      dplyr::mutate(refseq_select = dplyr::if_else(
        !is.na(refseq_select), rs_transcript_id, 
        as.character(NA)
      )) |>
      dplyr::mutate(
        transcript_mane_plus_clinical = dplyr::if_else(
          !is.na(mane_plus_clinical), 
          rs_transcript_id, 
          as.character(NA)
        )
      ) |>
      dplyr::distinct() |>
      dplyr::select(-c("mane_select", "mane_plus_clinical"))
    
    refseq_map <- refseq_map_grch37

  }
  
  return(refseq_map)
  
  
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

read_gtf <- function(path, zero_based = TRUE) {
  if (grepl("^https?://", path)) {
    path <- download_cached(path)
  }
  gtf <- rtracklayer::import(path)
  gtf <- as.data.frame(gtf)
  gtf <- dplyr::mutate_if(gtf, is.factor, as.character)
  res <- dplyr::rename(gtf, chrom = seqnames)
  
  if (zero_based) {
    res <- dplyr::mutate(res, start = start - 1L)
  }
  
  tibble::as_tibble(res)
}
  

