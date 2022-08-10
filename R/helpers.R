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

