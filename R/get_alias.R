#' @title
#' Get human gene aliases from NCBI
#'
#' @description
#' Downloads and returns a dataset that indicate ambiguous and unambiguous 
#' gene aliases/synonyms for human genes. Gene aliases with less than three 
#' characters have been ignored, and a few custom aliases have been added 
#' (source = `custom`). The dataset comes as a `list` object, with two 
#' elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation 
#'                resources used
#' * `records` - a list with gene aliases indicating ambiguous/ 
#'               non-ambiguous state
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be 
#'  overwritten (set to TRUE to re-download if file exists in cache)
#'
#' @returns
#' \strong{metadata} - A data frame with 1 row and 6 columns:
#' \itemize{
#'   \item \emph{source} - gene annotation source
#'   \item \emph{annotation_data} - type of annotations used
#'   \item \emph{url} - URL of annotation resource
#'   \item \emph{citation} - publication to cite for annotation source
#'   (citation; PMID)
#'   \item \emph{version} - version used/datestamp
#'   \item \emph{abbreviation} - abbreviation used in column names of records
#'  }
#'
#' \strong{records} - A data frame with 177,916 rows and 5 columns:
#' \itemize{
#'   \item \emph{alias} - gene alias/synonym
#'   \item \emph{symbol} - primary symbol
#'   \item \emph{entrezgene} - Entrez gene identifier
#'   \item \emph{n_primary_map} - number of primary symbols linked to the alias
#'   \item \emph{ambiguous} - logical indicating if alias is ambiguous or not
#'   \item \emph{source} - source for gene synonyms (NCBI, custom)
#
#' }
#' 
#' @examples
#' 
#' \dontrun{
#' library(geneOncoX)
#' gene_alias <- get_alias(cache_dir = tempdir())
#' }
#'
#' @export
#'

get_alias <- function(cache_dir = NA,
                      force_download = FALSE) {
    dat <- get_gox_data(
        cache_dir = cache_dir,
        force_download = force_download,
        db = "alias"
    )
    return(dat)
}
