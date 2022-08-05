#' Function that retrieves geneOncoXREF data from Google Drive
#'
#' @param cache_dir Local directory for data download
#' @param overwrite Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' @param db type of dataset to be retrieved
#'
#' @keywords internal
#'
#'
get_gox_data <- function(cache_dir = NULL,
                         overwrite = F,
                         db = "alias"){

  if(is.null(cache_dir) | !dir.exists(cache_dir)){
    lgr::lgr$error(paste0("Argument cache_dir = '",
                          cache_dir, "' is NULL or does not exist"))
  }

  fname_local <- file.path(
    cache_dir,
    paste0("gene_",db,"_v",
           geneOncoXREF:::db_id_ref[geneOncoXREF:::db_id_ref$name == db,]$pVersion,
           '.rds')
  )

  fname_gd <- googledrive::as_id(
    geneOncoXREF:::db_id_ref[geneOncoXREF:::db_id_ref$name == db,]$gid)

  md5checksum_package <-
    geneOncoXREF:::db_id_ref[geneOncoXREF:::db_id_ref$name == db,]$md5Checksum

  dat <- NULL
  if(file.exists(fname_local) & overwrite == F){
    dat <- readRDS(fname_local)
    if(!is.null(dat[['records']]) & !is.null(dat[['metadata']])){
      lgr::lgr$info(paste0(
        "Reading from cache_dir = '", cache_dir, "', argument overwrite = F"))
      lgr::lgr$info(paste0("Object 'gene_",db,"' sucessfully loaded"))
      if(db == 'gencode'){
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[['records']][['grch37']]), " records for",
          " genome build grch37"))
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[['records']][['grch38']]), " records for",
          " genome build grch38"))
      }else{
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[['records']]), " records"))
      }
    }

  }else{

    googledrive::drive_deauth()

    lgr::lgr$info("Downloading remote dataset from Google Drive to cache_dir")
    dl <- googledrive::with_drive_quiet(
      googledrive::drive_download(
        fname_gd,
        path = fname_local, overwrite = TRUE)
    )

    md5checksum_remote <- dl$drive_resource[[1]]$md5Checksum
    md5checksum_local <- tools::md5sum(fname_local)
    names(md5checksum_local) <- NULL

    if(md5checksum_remote == md5checksum_local){
      dat <- readRDS(fname_local)
      if(!is.null(dat[['records']]) & !is.null(dat[['metadata']])){

        lgr::lgr$info(paste0(
          "Reading from cache_dir = ' (", cache_dir, "'), argument overwrite = F"))
        lgr::lgr$info(paste0("Object 'gene_",db,"' sucessfully loaded"))
        lgr::lgr$info(paste0("md5 checksum is valid: ", md5checksum_remote))

        if(db == 'gencode'){
          lgr::lgr$info(paste0(
            "Retrieved n = ", nrow(dat[['records']][['grch37']]), " records for",
            " genome build grch37"))
          lgr::lgr$info(paste0(
            "Retrieved n = ", nrow(dat[['records']][['grch38']]), " records for",
            " genome build grch38"))
        }else{
          lgr::lgr$info(paste0(
            "Retrieved ", nrow(dat[['records']]), " records"))
        }
      }
    }else{
      lgr::lgr$error(paste0("md5 checksum of local file (", md5checksum_local,
                            ") is inconsistent with remote file (",
                            md5checksum_remote,")"))
    }

  }
  return(dat)
}
