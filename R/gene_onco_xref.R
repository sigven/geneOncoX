#' Function that retrieves geneOncoX data from Google Drive
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache
#'   should be overwritten (set to TRUE to re-download if file
#'   exists in cache)
#' @param db type of dataset to be retrieved
#'
#' @return pre-processed gene/transcript annotation records and metadata
#'
#' @keywords internal
#'
#'
get_gox_data <- function(cache_dir = NA,
                         force_download = FALSE,
                         db = "alias") {
  dat <- list()
  dat[["metadata"]] <- data.frame()
  dat[["records"]] <- data.frame()

  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T")
  )

  if (is.na(cache_dir)) {
    lgr::lgr$error(paste0(
      "Argument cache_dir = '",
      cache_dir, "' is not defined"
    ))
    return(dat)
  }
  if (!dir.exists(cache_dir)) {
    lgr::lgr$error(paste0(
      "Argument cache_dir = '",
      cache_dir, "' does not exist"
    ))
    return(dat)
  }
  fname_local <- file.path(
    cache_dir,
    paste0(
      "gene_", db, "_v",
      db_id_ref[db_id_ref$name == db, ]$pVersion,
      ".rds"
    )
  )

  fname_gd <- googledrive::as_id(
    db_id_ref[db_id_ref$name == db, ]$gid
  )

  md5checksum_package <-
    db_id_ref[db_id_ref$name == db, ]$md5Checksum

  dat <- NULL
  if (file.exists(fname_local) && force_download == FALSE) {
    dat <- readRDS(fname_local)
    dat$fpath <- fname_local
    if (!is.null(dat[["records"]]) && !is.null(dat[["metadata"]])) {
      lgr::lgr$info(paste0(
        "Reading from cache_dir = '",
        cache_dir,
        "', argument force_download = FALSE"
      ))
      lgr::lgr$info(paste0("Object 'gene_", db, "' sucessfully loaded"))
      if (startsWith(db,"gencode")) {
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[["records"]][["grch37"]]),
          " records for",
          " genome build grch37"
        ))
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[["records"]][["grch38"]]),
          " records for",
          " genome build grch38"
        ))
      } else {
        lgr::lgr$info(paste0(
          "Retrieved n = ", nrow(dat[["records"]]), " records"
        ))
      }
    }
  } else {
    googledrive::drive_deauth()

    lgr::lgr$info(
      "Downloading remote dataset from Google Drive to cache_dir"
    )
    dl <- googledrive::with_drive_quiet(
      googledrive::drive_download(
        fname_gd,
        path = fname_local,
        overwrite = TRUE
      )
    )

    md5checksum_remote <- dl$drive_resource[[1]]$md5Checksum
    md5checksum_local <- tools::md5sum(fname_local)
    names(md5checksum_local) <- NULL

    if (md5checksum_remote == md5checksum_local) {
      dat <- readRDS(fname_local)
      dat$fpath <- fname_local
      if (!is.null(dat[["records"]]) &&
        !is.null(dat[["metadata"]])) {
        lgr::lgr$info(paste0(
          "Reading from cache_dir = '", cache_dir,
          "', argument force_download = ", force_download
        ))
        lgr::lgr$info(paste0(
          "Object 'gene_", db, "' sucessfully loaded"
        ))
        lgr::lgr$info(paste0(
          "md5 checksum is valid: ", md5checksum_remote
        ))

        if (db == "gencode") {
          lgr::lgr$info(paste0(
            "Retrieved n = ",
            nrow(dat[["records"]][["grch37"]]), " records for",
            " genome build grch37"
          ))
          lgr::lgr$info(paste0(
            "Retrieved n = ",
            nrow(dat[["records"]][["grch38"]]), " records for",
            " genome build grch38"
          ))
        } else {
          lgr::lgr$info(paste0(
            "Retrieved ", nrow(dat[["records"]]), " records"
          ))
        }
      }
    } else {
      lgr::lgr$error(
        paste0(
          "md5 checksum of local file (", md5checksum_local,
          ") is inconsistent with remote file (",
          md5checksum_remote, ")"
        )
      )
    }
  }
  return(dat)
}
