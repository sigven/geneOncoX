get_intogen_driver_genes <- function(gene_info = NULL){

  intogen_drivers <- as.data.frame(
    read.table(file = file.path(
      "data-raw", "intogen",
      "Compendium_Cancer_Genes.tsv"),sep="\t",
      header=T, quote=NULL, stringsAsFactors = F) |>
      janitor::clean_names() |>
      dplyr::group_by(symbol) |>
      dplyr::mutate(role = dplyr::if_else(
        role == "",as.character(NA),as.character(role))) |>
      dplyr::mutate(role = dplyr::if_else(
        role == "Act","Activating",as.character(role))) |>
      dplyr::mutate(role = dplyr::if_else(
        role == "LoF","Loss_of_Function",as.character(role))) |>
      dplyr::summarise(intogen_role = paste(unique(Hmisc::capitalize(role)),collapse="&"),
                       intogen_phenotype = paste(unique(cancer_type),collapse="&"),
                       .groups = "drop") |>
      dplyr::mutate(intogen_role = dplyr::if_else(
        is.na(intogen_role) | intogen_role == "NA",
        "Unknown",
        as.character(intogen_role)
      )) |>
      dplyr::mutate(symbol = dplyr::case_when(
        symbol == "CARS" ~ "CARS1",
        symbol == "FAM46C" ~ "TENT5C",
        symbol == "H3F3A" ~ "H3-3A",
        symbol == "HIST1H3B" ~ "H3C2",
        symbol == "HIST1H4I" ~ "H4C9",
        symbol == "SEPT9" ~ "SEPTIN9",
        TRUE ~ as.character(symbol)
      )) |>
      dplyr::left_join(
        dplyr::select(gene_info, symbol, entrezgene), 
        by = c("symbol")) |>
      dplyr::filter(!is.na(entrezgene)) |>
      dplyr::select(-symbol) |>
      dplyr::distinct()
  )
  lgr::lgr$info(paste0("Parsed n = ", nrow(intogen_drivers),
                           " predicted driver genes (Intogen)"))

  return(intogen_drivers)
}


get_curated_fp_cancer_genes <- function(
    gene_info = NULL){

  invisible(assertthat::assert_that(!is.null(gene_info), msg = "'gene_info' object is NULL"))
  assertable::assert_colnames(gene_info,c("symbol","entrezgene"), only_colnames = F, quiet = T)
  assertable::assert_coltypes(gene_info, list(symbol = character(), entrezgene = integer()), quiet = T)
  xlsx_fname <- file.path(
    "data-raw", "tcga_driver_genes", "bailey_2018_cell.xlsx")
  invisible(assertthat::assert_that(file.exists(xlsx_fname),msg = paste0("File ",xlsx_fname," does not exist")))
  tmp <- openxlsx::read.xlsx(xlsx_fname,sheet = 8,startRow = 4)
  fp_cancer_drivers <- data.frame('symbol' = tmp[,2], stringsAsFactors = F) |>
    dplyr::mutate(bailey2018_fp_driver = TRUE) |>
    dplyr::filter(!is.na(symbol)) |>
    dplyr::left_join(dplyr::select(gene_info, symbol, entrezgene), by = c("symbol")) |>
    dplyr::filter(!is.na(entrezgene)) |>
    dplyr::select(-symbol)
  return(fp_cancer_drivers)

  lgr::lgr$info(paste0("Parsed n = ", nrow(fp_cancer_drivers),
                           " false positive driver genes (TCGA, Bailey et al., Cell, 2018"))


}

get_tcga_driver_genes <- function(){

  tcga_projects <- readRDS(
    paste0("data-raw/tcga_projects.rds")) |>
    dplyr::bind_rows(
      data.frame('tumor' = 'PANCAN',
                 'name' = 'Pancancer', stringsAsFactors = F))

  tcga_drivers <- as.data.frame(
    openxlsx::read.xlsx(
      file.path(
        "data-raw",
        "tcga_driver_genes",
        "bailey_2018_cell.xlsx"),
      sheet = 2,startRow = 4,cols = c(1,2,4,6,7)) |>
      janitor::clean_names() |>
      dplyr::mutate(tissue_frequency =
                      dplyr::if_else(is.na(tissue_frequency) |
                                       tissue_frequency == "NA",
                                     as.character(pancan_frequency),
                                     as.character(tissue_frequency))) |>
      dplyr::rename(symbol = gene, tumor = cancer, frequency = tissue_frequency) |>
      dplyr::mutate(frequency = as.numeric(frequency)) |>
      dplyr::left_join(tcga_projects, by=c("tumor")) |>
      dplyr::select(symbol, tumor, name, frequency) |>
      dplyr::mutate(frequency = as.numeric(frequency)) |>
      dplyr::mutate(name_freq = paste0(tumor,":",round(frequency,3))) |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(tcga_driver = paste(name_freq, collapse = "&"), .groups = "drop") |>
      dplyr::filter(tcga_driver != 'PANCAN:0'))

  lgr::lgr$info(paste0("Parsed n = ", nrow(tcga_drivers),
                           " predicted driver genes (TCGA, Bailey et al., Cell, 2018"))


  return(tcga_drivers)
}

get_signaling_pathway_genes <- function(gene_info){

  signaling_genes_tcga_cell <-
    file.path("data-raw",
              "oncogenic_signaling",
              "tcga_pancan_2018_sanchez_vega.xlsx")

  cell_cycle_genes <-
    openxlsx::read.xlsx(signaling_genes_tcga_cell,
                        sheet = 1,startRow = 3) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "Cell Cycle", signaling_pathway ="CC")

  hippo_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell, sheet = 2, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "Hippo Signaling Cascade", signaling_pathway ="HIPPO")

  myc_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell, sheet = 3, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "MYC", signaling_pathway ="MYC")

  notch_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell, sheet = 4, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(symbol = dplyr::if_else(
      symbol == "NOV","CCN3",as.character(symbol)
    )) |>
    dplyr::mutate(pathway = "Notch Signaling Pathway", signaling_pathway ="NOTCH")

  nrf2_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
                                    sheet = 5, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "Oxidative Stress Response", signaling_pathway ="NRF2")

  pi3k_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
                                    sheet = 6, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "PI-3-Kinase signaling", signaling_pathway ="PI3K")

  rtkras_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
                                      sheet = 8, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "Receptor-Tyrosine Kinase (RTK)/RAS/MAP-Kinase Signaling", signaling_pathway ="RTKRAS")

  tgfbeta_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
                                       sheet = 7, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "TGFβ Signaling", signaling_pathway ="TFGBETA")

  tp53_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,
                                    sheet = 9, startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "TP53 Signaling", signaling_pathway ="TP53")

  wnt_genes <- openxlsx::read.xlsx(signaling_genes_tcga_cell,sheet = 10,startRow = 1) |>
    dplyr::select(Gene) |>
    dplyr::rename(symbol = Gene) |>
    dplyr::mutate(pathway = "β-catenin/Wnt Signaling", signaling_pathway ="WNT")

  signaling_genes <- as.data.frame(
    dplyr::bind_rows(cell_cycle_genes, hippo_genes,
                     nrf2_genes, myc_genes,
                     wnt_genes, tp53_genes,
                     rtkras_genes, tgfbeta_genes,
                     pi3k_genes, notch_genes) |>
      dplyr::group_by(symbol) |>
      dplyr::summarise(sanchezvega2018_signaling_pathway_short = paste(unique(signaling_pathway),collapse="&"),
                       sanchezvega2018_signaling_pathway = stringi::stri_enc_toascii(paste(unique(pathway),collapse="&")), .groups = "drop")) |>
    dplyr::left_join(
      dplyr::select(gene_info, symbol, entrezgene), by = "symbol"
    ) |>
    dplyr::select(-c(symbol, sanchezvega2018_signaling_pathway_short))

  lgr::lgr$info(paste0("Parsed n = ", nrow(signaling_genes),
                           " genes with pathway signalling annotation (TCGA, Sanchez-Vega et al., Cell, 2018"))

  return(signaling_genes)

}

get_cancer_gene_census <- function(
    origin = "somatic",
    cgc_version = "96"){

  cosmic_cgc <- readr::read_csv(
    file = file.path(
      "data-raw", "cancer_gene_census",
      "cancer_gene_census_96.csv"),
    show_col_types = F) |>
    janitor::clean_names()

  if(origin == "somatic"){
    cosmic_cgc <- cosmic_cgc |>
      dplyr::filter(somatic == "yes")
  }else{
    cosmic_cgc <- cosmic_cgc |>
      dplyr::filter(germline == "yes")
  }

  cosmic_cgc <- cosmic_cgc |>
    dplyr::select(gene_symbol,
                  entrez_gene_id,
                  germline,
                  somatic,
                  tier,
                  role_in_cancer,
                  molecular_genetics,
                  cancer_syndrome,
                  tumour_types_somatic,
                  tumour_types_germline) |>
    dplyr::rename(symbol = gene_symbol,
                  entrezgene = entrez_gene_id,
                  cgc_tier = tier,
                  moi = molecular_genetics) |>
    dplyr::mutate(cgc_oncogene = dplyr::if_else(
      !is.na(role_in_cancer) & stringr::str_detect(role_in_cancer,"oncogene"),
      TRUE,FALSE
    )) |>
    dplyr::mutate(cgc_tsg = dplyr::if_else(
      !is.na(role_in_cancer) & stringr::str_detect(role_in_cancer,"TSG"),
      TRUE,FALSE
    )) |>
    dplyr::mutate(moi = stringr::str_trim(moi)) |>
    dplyr::mutate(moi = dplyr::if_else(stringr::str_detect(moi,"^Rec$"),"AR",as.character(moi))) |>
    dplyr::mutate(moi = dplyr::if_else(stringr::str_detect(moi,"^Dom$"),"AD",as.character(moi))) |>
    dplyr::mutate(moi = dplyr::if_else(stringr::str_detect(moi,"^Dom/Rec$"),"AD/AR",as.character(moi))) |>
    dplyr::mutate(cgc_moi = dplyr::if_else(stringr::str_detect(moi,"Rec/X"),as.character(NA),as.character(moi))) |>
    dplyr::select(symbol, cgc_moi, entrezgene, cgc_tsg, cgc_oncogene,
                  tumour_types_somatic, cancer_syndrome, tumour_types_germline) |>
    ## bug in cancer gene census for ERCC5
    dplyr::mutate(entrezgene = dplyr::if_else(
      symbol == 'ERCC5',as.integer(2073),as.integer(entrezgene))
    ) |>
    dplyr::select(-symbol)

  if(origin == "somatic"){
    cosmic_cgc <- cosmic_cgc |>
      dplyr::mutate(cgc_somatic = T) |>
      dplyr::mutate(cgc_phenotype_somatic = tumour_types_somatic) |>
      dplyr::select(-c(cancer_syndrome,
                       tumour_types_germline,
                       tumour_types_somatic))

  }else{
    cosmic_cgc <- cosmic_cgc |>
      dplyr::mutate(cgc_germline = T) |>
      dplyr::mutate(cgc_phenotype_germline = paste(tumour_types_germline,
                                                   cancer_syndrome, sep="|")) |>
      dplyr::select(-c(cancer_syndrome,
                       tumour_types_germline,
                       tumour_types_somatic))
  }

  lgr::lgr$info(paste0("Parsed n = ", nrow(cosmic_cgc),
                           " genes in COSMIC's Cancer Gene Census (version ",
                           cgc_version,") - ", origin, " context"))

  return(cosmic_cgc)

}

get_network_of_cancer_genes <- function(
    ncg_version = "7.0"){

  ncg <- read.table(file = file.path(
    "data-raw","ncg","ncg.tsv"),header = T,
    stringsAsFactors = F, sep = "\t", 
    quote = "", comment.char = "") |>
    janitor::clean_names() |>
    dplyr::mutate(ncg_driver = dplyr::if_else(
      !is.na(type) & type == "\"Canonical Cancer Driver\"",
      as.logical(TRUE), as.logical(FALSE)
    )) |>
    dplyr::mutate(cancer_type = stringr::str_replace_all(
      cancer_type,"\\\"","")) |>
    dplyr::mutate(
      ncg_tsg = dplyr::if_else(
        stringr::str_detect(cgc_annotation,"TSG") | 
          ncg_tsg == 1,
        as.logical(TRUE),as.logical(FALSE))) |>
    dplyr::mutate(
      ncg_oncogene = dplyr::if_else(
        stringr::str_detect(cgc_annotation,"oncogene") | 
          ncg_oncogene == 1,
        as.logical(TRUE),as.logical(FALSE))) |>
    dplyr::rename(entrezgene = entrez) |>
    dplyr::filter(!(ncg_tsg == F &
                      ncg_oncogene == F &
                      ncg_driver == F)) |>
    dplyr::distinct() |>
    dplyr::group_by(entrezgene) |>
    dplyr::summarise(
      ncg_phenotype = paste(
        sort(unique(cancer_type)), collapse=","),
      ncg_pmid = paste(
        sort(unique(pubmed_id)), collapse=";"),
      ncg_tsg = paste(
        unique(ncg_tsg), collapse=","),
      ncg_oncogene = paste(
        unique(ncg_oncogene), collapse=","),
      ncg_driver = paste(
        unique(ncg_driver), collapse=",")) |>
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
      ncg_oncogene, ncg_pmid, ncg_phenotype) |>
    dplyr::mutate(ncg_phenotype = stringr::str_replace_all(
      ncg_phenotype, ",,",","
    )) |>
    dplyr::mutate(ncg_phenotype = stringr::str_replace_all(
      ncg_phenotype, "^,",""
    )) |>
    dplyr::distinct()

  lgr::lgr$info(paste0("Found n = ", nrow(ncg),
                           " cancer-relvant genes (drivers, proto-oncogenes, tumor suppressors) in Network of Cancer Genes (version ",
                           ncg_version,")"))
  return(ncg)

}

get_cancermine_genes <- function(cancermine_version = "48"){

  cancermine_sentences_fname <-
    paste0('data-raw/cancermine/data-raw/cancermine_sentences.v',
           cancermine_version,'.tsv.gz')
  cancermine_citations_fname <-
    paste0("data-raw/cancermine/output/cancermine_citations.v",
           cancermine_version,".tsv.gz")
  cancermine_collated_fname <-
    paste0("data-raw/cancermine/data-raw/cancermine_collated.v",
           cancermine_version,".tsv.gz")

  pmids <- as.data.frame(
    read.table(file = gzfile(cancermine_sentences_fname),
               header = T, comment.char = "", quote = "",
               sep="\t", stringsAsFactors = F) |>
      dplyr::filter(predictprob >= 0.8) |>
      
      ## some entries wrongly captured, are in
      ## fact mentions of anti-sense non-coding genes
      dplyr::filter(
        !stringr::str_detect(
          formatted_sentence, "-<b>AS1|b>-AS1"
        )
      ) |>
      dplyr::group_by(role, gene_entrez_id, pmid) |>
      dplyr::summarise(doid = paste(
        unique(cancer_id), collapse=","), .groups = "drop") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::distinct()
  )

  all_citations <- data.frame()
  if(file.exists(cancermine_citations_fname)){
    all_citations <- read.table(
      file = gzfile(cancermine_citations_fname),
      sep = "\t", header = F, quote = "",
      comment.char = "", stringsAsFactors = F) |>
      magrittr::set_colnames(c('pmid','citation','citation_link')) |>
      dplyr::mutate(
        citation = stringi::stri_enc_toascii(citation)) |>
      dplyr::mutate(
        citation_link = stringi::stri_enc_toascii(citation_link)) |>
      dplyr::distinct() |>
      dplyr::arrange(desc(pmid))

  }

  pmids <- pmids |>
    dplyr::inner_join(all_citations, by=c("pmid"))

  pmids_oncogene <- as.data.frame(
    pmids |>
      dplyr::filter(role == "Oncogene") |>
      dplyr::arrange(entrezgene, desc(pmid)) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        pmids_oncogene = paste(pmid,collapse=";"),
        citations_oncogene = paste(
          citation,collapse="; "),
        citation_links_oncogene = paste(
          head(citation_link,50),collapse=", "),
        .groups = "drop") |>
      dplyr::mutate(n_citations_oncogene = as.integer(
        stringr::str_count(pmids_oncogene,";")) + 1) |>
      dplyr::filter(n_citations_oncogene > 1) |>
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(
        n_citations_oncogene >= 6,
        "MC",
        as.character("LC"))) |>
      dplyr::mutate(oncogene_cancermine = dplyr::if_else(
        n_citations_oncogene >= 15,
        "HC",
        as.character(oncogene_cancermine)))
  )

  pmids_tsgene <- as.data.frame(
    pmids |>
      dplyr::filter(role == "Tumor_Suppressor") |>
      dplyr::arrange(entrezgene, desc(pmid)) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        pmids_tsgene = paste(unique(pmid),collapse=";"),
        citations_tsgene = paste(
          citation,collapse="; "),
        citation_links_tsgene = paste(head(
          citation_link,50),collapse=", "),
        .groups = "drop") |>
      dplyr::ungroup() |>
      dplyr::mutate(n_citations_tsgene = as.integer(
        stringr::str_count(pmids_tsgene,";")) + 1) |>
      dplyr::filter(n_citations_tsgene > 1) |>
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(
        n_citations_tsgene >= 6,
        "MC",
        as.character("LC"))) |>
      dplyr::mutate(tumor_suppressor_cancermine = dplyr::if_else(
        n_citations_tsgene >= 15,
        "HC",
        as.character(tumor_suppressor_cancermine)))
  )

  pmids_cdriver <- as.data.frame(
    pmids |>
      dplyr::filter(role == "Driver") |>
      dplyr::arrange(entrezgene, desc(pmid)) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(
        pmids_cdriver = paste(pmid,collapse=";"),
        citations_cdriver = paste(citation,collapse="; "),
        citation_links_cdriver = paste(
          head(citation_link,50),collapse=", "),
        .groups = "drop") |>
      dplyr::mutate(n_citations_cdriver = as.integer(
        stringr::str_count(pmids_cdriver,";")) + 1) |>
      dplyr::filter(n_citations_cdriver > 1) |>
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(
        n_citations_cdriver >= 6,
        "MC",
        as.character("LC"))) |>
      dplyr::mutate(cancer_driver_cancermine = dplyr::if_else(
        n_citations_cdriver >= 10,
        "HC",
        as.character(cancer_driver_cancermine)))
  )


  lgr::lgr$info("Retrieving known proto-oncogenes/tumor suppressor genes from CancerMine")
  oncogene <- as.data.frame(
    readr::read_tsv(cancermine_collated_fname,
                    col_names = T,na ="-",
                    comment = "#", quote = "",
                    show_col_types = F) |>
      dplyr::filter(role == "Oncogene") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(doid_oncogene = paste(unique(cancer_id), collapse=","), .groups = "drop") |>
      dplyr::inner_join(pmids_oncogene, by = "entrezgene") |>
      dplyr::distinct())

  n_hc_oncogene <- oncogene |>
    dplyr::filter(oncogene_cancermine == "HC") |>
    nrow()


  tsgene <- as.data.frame(
    readr::read_tsv(cancermine_collated_fname,
                    col_names = T,na ="-",
                    comment = "#", quote = "",
                    show_col_types = F) |>
      dplyr::filter(role == "Tumor_Suppressor") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(doid_tsgene = paste(
        unique(cancer_id),collapse=","), .groups = "drop") |>
      dplyr::inner_join(pmids_tsgene, by = "entrezgene") |>
      dplyr::distinct()
  )
  n_hc_tsgene <- tsgene |>
    dplyr::filter(tumor_suppressor_cancermine == "HC") |>
    nrow()


  cdriver <- as.data.frame(
    readr::read_tsv(cancermine_collated_fname,
                    col_names = T,na ="-",
                    comment = "#", quote = "",
                    show_col_types = F) |>
      dplyr::filter(role == "Driver") |>
      dplyr::rename(entrezgene = gene_entrez_id) |>
      dplyr::group_by(entrezgene) |>
      dplyr::summarise(doid_cdriver = paste(
        unique(cancer_id),collapse=","), .groups = "drop") |>
      dplyr::inner_join(pmids_cdriver, by = "entrezgene") |>
      dplyr::distinct()
  )
  n_hc_cdriver <- cdriver |>
    dplyr::filter(cancer_driver_cancermine == "HC") |>
    nrow()

  cancermine_full <- cdriver |>
    dplyr::full_join(tsgene, by = "entrezgene") |>
    dplyr::full_join(oncogene, by = "entrezgene") |>
    dplyr::rename(cancermine_pmid_driver = pmids_cdriver,
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
                  cancermine_cit_links_driver = citation_links_cdriver) |>
    dplyr::select(entrezgene,
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
                  cancermine_cit_links_tsg) |>
    dplyr::mutate(cancermine_n_cit_oncogene = dplyr::if_else(
      is.na(cancermine_n_cit_oncogene),
      as.integer(0),
      as.integer(cancermine_n_cit_oncogene))) |>
    dplyr::mutate(cancermine_n_cit_tsg = dplyr::if_else(
      is.na(cancermine_n_cit_tsg),
      as.integer(0),
      as.integer(cancermine_n_cit_tsg))) |>
    dplyr::mutate(cancermine_n_cit_driver = dplyr::if_else(
      is.na(cancermine_n_cit_driver),
      as.integer(0),
      as.integer(cancermine_n_cit_driver)))


  return(cancermine_full)


}
