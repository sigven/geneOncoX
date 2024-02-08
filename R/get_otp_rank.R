#' @title
#' Get rank of genes according to association to cancer (Open Targets Platform)
#'
#' @description
#' Downloads and returns a dataset that ranks genes according to aggregated 
#' association scores between genes and cancer phenotype terms 
#' from the Open Targets Platform
#'
#' * `metadata` - a data frame with metadata regarding annotation
#'                resources used
#' * `records` - a data frame with gene scores/ranks (one record per gene and primary site)
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#'
#' @returns
#' \strong{metadata} - A data frame with 10 rows and 6 columns:
#' \itemize{
#'   \item \emph{source} - gene annotation source
#'   \item \emph{annotation_data} - type of annotations used
#'   \item \emph{url} - URL of annotation resource
#'   \item \emph{citation} - publication to cite for annotation source
#'   (citation; PMID)
#'   \item \emph{version} - version used
#'   \item \emph{abbreviation} - abbreviation used in column names of records
#'  }
#'
#' \strong{records} - A data frame with 20,954 rows and 6 columns:
#' \itemize{
#'   \item \emph{entrezgene} - NCBI Entrez gene identifier
#'   \item \emph{primary_site} - Primary tumor site
#'   \item \emph{tissue_assoc_score} - tissue-specific aggregated association score
#'   \item \emph{tissue_assoc_rank} - tissue-specific cancer rank
#'   \item \emph{global_assoc_score} - pan-cancer aggregated association score
#'   \item \emph{global_assoc_rank} - cancer gene rank (pan-cancer, across tissues)
#' }
#'
#' @examples
#' \dontrun{
#' library(geneOncoX)
#' gene_otp_rank <- get_otp_rank(cache_dir = tempdir())
#' }
#'
#' @export
#'

get_otp_rank <- function(cache_dir = NA, force_download = FALSE) {
  dat <- get_gox_data(
    cache_dir = cache_dir,
    force_download = force_download,
    db = "otp_rank"
  )
  return(dat)
}
