#' @title
#' Get human gene aliases from NCBI
#'
#' @description
#' A dataset that indicate ambiguous and unambiguous gene aliases/synonyms
#' for human genes. The dataset comes as a `list` object,
#' with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation resources used
#' * `records` - a list with gene aliases indicating ambiguous/non-ambiguous state
#'
#' @param cache_dir Local directory for data download
#' @param overwrite Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#'
#' @returns
#' \bold{metadata} - A data frame with 1 row and 6 columns:
#' \itemize{
#'   \item \emph{source} - gene annotation source
#'   \item \emph{annotation_data} - type of annotations used
#'   \item \emph{url} - URL of annotation resource
#'   \item \emph{citation} - publication to cite for annotation source (citation; PMID)
#'   \item \emph{version} - version used/datestamp
#'   \item \emph{abbreviation} - abbreviation used in column names of records
#'  }
#'
#' \bold{records} - A data frame with 144,316 rows and 5 columns:
#' \itemize{
#'   \item \emph{alias} - gene alias/synonym
#'   \item \emph{symbol} - primary symbol
#'   \item \emph{entrezgene} - Entrez gene identifier
#'   \item \emph{n_primary_map} - number of primary symbols linked to the alias
#'   \item \emph{ambiguous} - logical indicating if alias is ambiguous or not
#
#' }
#'
#' @export
#'

get_alias <- function(cache_dir = NA, overwrite = F){

  dat <- get_gox_data(cache_dir = cache_dir,
                      overwrite = overwrite,
                      db = "alias")
  return(dat)
}
