#' Function to write virtual panel BED files (CPSR)
#'
#' @keywords internal
#'
create_virtual_panel_bed <- function(
    gwas_bed_fpath = NA,
    gox_gencode = NULL,
    gox_panels = NULL,
    gox_predisposition = NULL,
    super_panel = F,
    #include_secondary_findings = TRUE,
    build = "grch37",
    dest_dir = NA){
  
  
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))
  
  if(is.na(dest_dir)){
    lgr::lgr$error(paste0("Argument dest_dir = '",
                          dest_dir, "' is not defined"))
  }
  
  if(!dir.exists(dest_dir)){
    lgr::lgr$error(paste0("Argument dest_dir = '",
                          dest_dir, "' does not exist"))
  }
  
  if(!(file.exists(gwas_bed_fpath) & endsWith(gwas_bed_fpath,".bed.gz"))){
    lgr::lgr$error(paste0("File gwas_bed_file = '",
                          gwas_bed_fpath, "' does not exist"))
  }
  
  ## read build-specific GWAS BED file (argument)
  gwas_bed_fpath_uncomp <- file.path(
    dest_dir, stringr::str_replace(basename(gwas_bed_fpath),".gz",""))
  system(paste0('bgzip -dc ',gwas_bed_fpath," > ",gwas_bed_fpath_uncomp))
  gwas_bed_data <- readr::read_tsv(gwas_bed_fpath_uncomp, col_names = F,
                                   show_col_types = F) 
  colnames(gwas_bed_data) <- c('chrom','start','end','name')
  system(paste0('rm -f ', gwas_bed_fpath_uncomp))
  
  ## Get all build-specific Genomics England PanelApp records
  all_panel_records <- gox_panels$records |>
    dplyr::filter(genome_build == build)
  
  ## Get build-specific GENCODE transcripts
  gencode_transcripts <- gox_gencode$records[[build]] |>
    dplyr::mutate(chrom = as.character(stringr::str_replace(chrom,"chr","")))
  
  ## Get ACMG secondary-findings genes and intersect
  ## them with GENCODE tracks -> BED data (chrom, start, end, name)
  acmg_sf_bed_data <- gox_predisposition$records |>
    dplyr::filter(stringr::str_detect(predisp_source,"ACMG_SF")) |>
    dplyr::select(entrezgene, moi) |>
    dplyr::inner_join(
      dplyr::select(gencode_transcripts, chrom, start, end,
                    ensembl_gene_id, ensembl_transcript_id,
                    entrezgene, symbol),
      by = c("entrezgene")) |>
    dplyr::mutate(name = paste(
      symbol, ensembl_gene_id, entrezgene, ensembl_transcript_id,
      "ACMG_SF", moi, sep="|")) |>
    dplyr::select(chrom, start, end, name) |>
    dplyr::mutate(chrom = as.character(stringr::str_replace(chrom,"chr","")))
    
  
  unique_panel_ids <- unique(all_panel_records$id)
  i <- 1
  while(i <= length(unique_panel_ids)){
    panel_id <- unique_panel_ids[i]
    panel_name <- unique(
      all_panel_records[all_panel_records$id == i,"gepa_panel_name"])
    
    ## Get genes for the specific panel, and intersect them
    ## with GENCODE transcripts  --> BED format (chrom, start, end, name)
    panel_bed_data <- dplyr::filter(all_panel_records, id == i) |>
      dplyr::select(id, entrezgene, ensembl_gene_id,
                    gepa_panel_name, gepa_panel_id, 
                    gepa_confidence_level, gepa_moi,
                    gepa_panel_version) |>
      dplyr::inner_join(
        dplyr::select(gencode_transcripts, chrom, start, end,
                      ensembl_gene_id, ensembl_transcript_id,
                      entrezgene, symbol),
        by = c("ensembl_gene_id","entrezgene")) |>
      dplyr::mutate(name = paste(
        symbol, ensembl_gene_id, entrezgene, ensembl_transcript_id,
        stringr::str_replace_all(gepa_panel_name," ","_"),
        paste0("GEPA_PANEL_", id), gepa_confidence_level,
        gepa_panel_id, gepa_panel_version, gepa_moi, sep="|")) |>
      dplyr::select(chrom, start, end, name, gepa_confidence_level) |>
      dplyr::mutate(chrom = as.character(stringr::str_replace(chrom,"chr","")))
    
    
    ## Add BED records with ACMG genes and GWAS variants
    ## - For all and high-confidence genes only (green)
    
    virtual_panels <- list()
    virtual_panels[['green']] <- panel_bed_data |>
      dplyr::filter(gepa_confidence_level == 3) |>
      dplyr::select(-gepa_confidence_level)
    virtual_panels[['all']]  <- panel_bed_data |>
      dplyr::select(-gepa_confidence_level)
    
    cat(paste0("Panel ", i, ": level green: ", nrow(virtual_panels[['green']]), "\n"))
    
    if(nrow(virtual_panels[['green']]) == 0){
      lgr::lgr$warn(paste0("Panel - ", panel_name, " has zero GREEN entries - using all"))
      virtual_panels[['green']] <- virtual_panels[['all']]
    }
    
    ## Green
    acmg_sf_bed_data_missing <- acmg_sf_bed_data |>
      dplyr::anti_join(virtual_panels[['green']], 
                       by = c("chrom", "start", "end"))
    
    virtual_panels[['green']] <- virtual_panels[['green']] |>
      dplyr::bind_rows(acmg_sf_bed_data_missing) |>
      dplyr::bind_rows(gwas_bed_data)
    
    geneOncoX:::write_bed_file(bed_data = virtual_panels[['green']],
                               bed_fname = paste0(panel_id,'.',build,'.GREEN.bed'),
                               dest_dir = "~/Downloads/test_panels")
    
    ## All
    acmg_sf_bed_data_missing <- acmg_sf_bed_data |>
      dplyr::anti_join(virtual_panels[['all']], 
                       by = c("chrom", "start", "end"))
    
    virtual_panels[['all']] <- virtual_panels[['all']] |>
      dplyr::bind_rows(acmg_sf_bed_data_missing) |>
      dplyr::bind_rows(gwas_bed_data)
    
    geneOncoX:::write_bed_file(bed_data = virtual_panels[['all']],
                               bed_fname = paste0(panel_id,'.',build, '.bed'),
                               dest_dir = "~/Downloads/test_panels")
    
    i <- i + 1
    
      
  }
  
  
  
}

#' Function to sort BED data
#' 
#' @param unsorted_regions data frame with unsorted regions
#' @param chr_names logical indicating if chromosome names contain 'chr' or not
#'
#' @keywords internal
#' 
sort_bed_regions <- function(unsorted_regions, chr_names = F){
  sorted_regions <- NULL
  assertable::assert_colnames(
    unsorted_regions,c('start','end'),only_colnames = F, quiet = T)
  if("chrom" %in% colnames(unsorted_regions) & 
     "start" %in% colnames(unsorted_regions) & 
     "end" %in% colnames(unsorted_regions)){
    
    chrOrder <- c(as.character(c(1:22)),"X","Y","M")
    if(chr_names == T){
      chrOrder <- paste0('chr',c(as.character(c(1:22)),"X","Y","M"))
    }
    unsorted_regions$chrom <- factor(unsorted_regions$chrom, levels=chrOrder)
    unsorted_regions <- unsorted_regions[order(unsorted_regions$chrom),]
    
    sorted_regions <- data.frame()
    for(chrom in chrOrder){
      if(nrow(unsorted_regions[unsorted_regions$chrom == chrom,]) > 0){
        chrom_regions <- unsorted_regions[unsorted_regions$chrom == chrom,]
        chrom_regions_sorted <- chrom_regions[with(chrom_regions, order(start, end)),]
        sorted_regions <- dplyr::bind_rows(sorted_regions, chrom_regions_sorted)
      }
    }
    
    sorted_regions$start <- as.integer(sorted_regions$start)
    sorted_regions$end <- as.integer(sorted_regions$end)
  }
  return(sorted_regions)
  
}

#' Function write BED data to file
#' 
#' @param bed_data data frame with BED data
#' @param bed_fname file name for BED file
#' @param dest_dir directory to write BED file
#'
#' @keywords internal
#' 
write_bed_file <- function(bed_data, 
                           bed_fname, 
                           dest_dir = NA){
  
  if(!dir.exists(dest_dir)){
    lgr::lgr$error(paste0("Argument dest_dir = '",
                          dest_dir, "' does not exist"))
  }
  
  bed_fname_full <- file.path(dest_dir, bed_fname)
  bed_data_sorted <- geneOncoX:::sort_bed_regions(bed_data)
  write.table(bed_data_sorted, file = bed_fname_full, sep = "\t",
              row.names = F, col.names = F, quote = F)
  system(paste0("bgzip -f ",bed_fname_full))
  system(paste0("tabix -p bed ",bed_fname_full,'.gz'))
  
}

#' Function to classify genes as tumor suppressors/proto-oncogenes,
#' driver etc. based on multiple lines of evidence.
#'
#' @param gox_basic output from geneOncoX::get_basic()
#' @param min_citations_tsg minimum citations (CancerMine) for tsg support
#' @param min_citations_oncogene minimum citations (CancerMine) for oncogene support
#' @param min_citations_driver minimum citations (CancerMine) for driver support
#'
#' @keywords internal
#'
assign_cancer_gene_evidence <- function(gox_basic = NULL,
                                       min_citations_tsg = 15,
                                       min_citations_oncogene = 15,
                                       min_citations_driver = 10){
  
  
  tsg_oncogene_evidence <- gox_basic$records |>
    dplyr::filter(ncg_oncogene == T |
                    ncg_tsg == T |
                    cgc_oncogene == T |
                    cgc_tsg == T) |>
    dplyr::mutate(cancermine_onco_ts_citation_ratio = dplyr::if_else(
      !is.na(cancermine_n_cit_tsg),
      as.numeric(cancermine_n_cit_oncogene)/cancermine_n_cit_tsg,
      as.numeric(NA))) |>
    
    dplyr::mutate(oncogene = dplyr::if_else(
      ncg_oncogene == T | cgc_oncogene == T |
        (!is.na(cancermine_n_cit_oncogene) &
           cancermine_n_cit_oncogene >= min_citations_oncogene),
      TRUE, FALSE)) |>
    dplyr::mutate(oncogene = dplyr::if_else(
      is.na(oncogene), 
      FALSE, 
      as.logical(oncogene))) |>
    dplyr::mutate(tsg = dplyr::if_else(
      ncg_tsg == T | cgc_tsg == T |
        (!is.na(cancermine_n_cit_tsg) &
           cancermine_n_cit_tsg >= min_citations_tsg),
      TRUE, FALSE)) |>
    dplyr::mutate(tsg = dplyr::if_else(
      is.na(tsg),
      FALSE,
      as.logical(tsg))) |>
    dplyr::mutate(cancer_driver = dplyr::if_else(
      (!is.na(cancermine_n_cit_driver) &
         cancermine_n_cit_driver >= min_citations_driver),
      TRUE,FALSE)) |>
    dplyr::mutate(cancer_driver = dplyr::if_else(
      is.na(cancer_driver), 
      FALSE, 
      as.logical(cancer_driver))) |>

    ## Ignore tsg classification if considerable more support
    ## for oncogene
    dplyr::mutate(tsg = dplyr::if_else(
      tsg == T &
        !is.na(cancermine_onco_ts_citation_ratio) & 
        cancermine_onco_ts_citation_ratio > 3,
      FALSE,
      as.logical(tsg))) |>
    
    ## Ignore oncogene classification if considerable more support
    ## for tumor suppressors
    dplyr::mutate(oncogene = dplyr::if_else(
      oncogene == T &
        !is.na(cancermine_onco_ts_citation_ratio) & 
        cancermine_onco_ts_citation_ratio <= 0.33,
      FALSE,
      as.logical(oncogene))) |>
    
    
    dplyr::mutate(ncg_links_tsg = dplyr::if_else(
      ncg_tsg == T,
      paste0("NCG-TumorSuppressor: <a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",
             entrezgene,"\" target=\"_blank\">YES</a>"),
      as.character("NCG-TumorSuppressor: Not defined"))) |>
    dplyr::mutate(ncg_links_tsg = dplyr::if_else(
      is.na(ncg_links_tsg),
      as.character("NCG-TumorSuppressor: Not defined"),
      as.character(ncg_links_tsg))) |>
    dplyr::mutate(ncg_links_oncogene = dplyr::if_else(
      ncg_oncogene == T,
      paste0("NCG-OncoGene: <a href=\"http://ncg.kcl.ac.uk/query.php?gene_name=",
             entrezgene,"\" target=\"_blank\">Yes</a>"),
      as.character("NCG-OncoGene: Not defined"))) |>
    dplyr::mutate(ncg_links_oncogene = dplyr::if_else(
      is.na(ncg_links_oncogene),
      as.character("NCG-OncoGene: Not defined"),
      as.character(ncg_links_oncogene))) |>
    dplyr::mutate(ncg_link = paste(ncg_links_tsg,
                                   ncg_links_oncogene,
                                   sep=", ")) |>
    
    dplyr::mutate(cgc_links_tsg = dplyr::if_else(
      cgc_tsg == T,
      paste0("CGC-TumorSuppressor: <a href='https://cancer.sanger.ac.uk/census'",
             " target='_blank'>YES</a>"),
      as.character("CGC-TumorSuppressor: Not defined"))) |>
    dplyr::mutate(cgc_links_tsg = dplyr::if_else(
      is.na(cgc_links_tsg),
      as.character("CGC-TumorSuppressor: Not defined"),
      as.character(cgc_links_tsg))) |>
    dplyr::mutate(cgc_links_oncogene = dplyr::if_else(
      cgc_oncogene == T,
      paste0("CGC-OncoGene: <a href='https://cancer.sanger.ac.uk/census'",
             " target='_blank'>YES</a>"),
      as.character("CGC-OncoGene: Not defined"))) |>
    dplyr::mutate(cgc_links_oncogene = dplyr::if_else(
      is.na(cgc_links_oncogene),
      as.character("CGC-OncoGene: Not defined"),
      as.character(cgc_links_oncogene))) |>
    dplyr::mutate(cgc_link = paste(cgc_links_tsg,
                                   cgc_links_oncogene,
                                   sep=", ")) |>
    
    
    dplyr::mutate(
      citation_links_cdriver = dplyr::if_else(
        is.na(cancer_driver),
        as.character(NA),
        as.character(cancermine_cit_links_driver)),
      citations_cdriver = dplyr::if_else(
        is.na(cancer_driver), 
        as.character(NA), 
        as.character(cancermine_cit_driver)),
      citation_links_oncogene = dplyr::if_else(
        cancermine_n_cit_oncogene == 0,
        "Oncogenic role (CancerMine): No records",
        as.character(paste0("Oncogenic role (CancerMine): ",cancermine_cit_links_oncogene))),
      citations_oncogene = dplyr::if_else(
        cancermine_n_cit_oncogene == 0,
        as.character(NA),
        as.character(cancermine_cit_oncogene)),
      citation_links_tsg = dplyr::if_else(
        cancermine_n_cit_tsg == 0,
        "Tumor suppressive role (CancerMine): No records",
        as.character(paste0("Tumor suppressor role (CancerMine): ",cancermine_cit_links_tsg))),
      citations_tsg = dplyr::if_else(
        cancermine_n_cit_tsg == 0,
        as.character(NA),
        as.character(cancermine_cit_tsg))) |>
    dplyr::mutate(cancergene_support = stringr::str_replace(
      paste(citation_links_oncogene,
            citation_links_tsg,
            cgc_link,
            ncg_link, sep=", "),"^, ","")) |>
    dplyr::mutate(tsg_evidence = paste0(
      "NCG:",ncg_tsg,
      "&CGC:", cgc_tsg,
      "&CancerMine:",cancermine_n_cit_tsg)) |>
    dplyr::mutate(oncogene_evidence = paste0(
      "NCG:",ncg_oncogene,
      "&CGC:",cgc_oncogene,
      "&CancerMine:",cancermine_n_cit_oncogene)) |>
    dplyr::mutate(cancer_driver_evidence = paste0(
      "CancerMine:",
      cancermine_n_cit_driver)) |>
    dplyr::mutate(
      citation_links_oncogene =
        stringr::str_replace(
          citation_links_oncogene,"(NA, ){1,}","")) |>
    dplyr::mutate(
      citation_links_tsg =
        stringr::str_replace(
          citation_links_tsg,"(NA, ){1,}","")) |>
    dplyr::filter(bailey2018_fp_driver == F) |>
    dplyr::select(entrezgene, cancergene_support,
                  cancer_driver, cancer_driver_evidence,
                  oncogene, oncogene_evidence,
                  tsg, tsg_evidence) |>
    dplyr::rename(tumor_suppressor = tsg,
                  tumor_suppressor_evidence = tsg_evidence)
  
  
  n_onc_ts <- dplyr::filter(tsg_oncogene_evidence, oncogene == T & tumor_suppressor == T)
  n_ts <- dplyr::filter(tsg_oncogene_evidence, tumor_suppressor == T & oncogene == F)
  n_onc <- dplyr::filter(tsg_oncogene_evidence, oncogene == T & tumor_suppressor == F)
  
  lgr::lgr$info(paste0("A total of n = ",nrow(n_onc)," classified proto-oncogenes were retrieved from CancerMine/NCG/CGC"))
  lgr::lgr$info(paste0("A total of n = ",nrow(n_ts)," classified tumor suppressor genes were retrieved from CancerMine/NCG/CGC"))
  lgr::lgr$info(paste0("A total of n = ",nrow(n_onc_ts)," genes were annotated with dual roles as tumor suppressor genes and oncogenes from CancerMine/NCG/CGC"))
  
  
}