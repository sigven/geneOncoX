#' @title
#' Get cancer gene panel collections from PanelApp
#'
#' @description
#' Downloads and returns a collection of >40 cancer gene panels from
#' \href{https://panelapp.genomicsengland.co.uk/}{Genomics England PanelApp}.
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation
#'                resources used
#' * `records` - a data frame with genes found in each panel
#'
#' @details \strong{NOTE:} Gene panel records are provided per genome build
#' (filter on column `genome_build` to get a build-specific set of panels)
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#'
#' @returns
#' \strong{metadata} - A data frame with 1 row and 6 columns:
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
#' \strong{records} - A data frame with 2,566 rows and 13 columns:
#' \itemize{
#'   \item \emph{genome_build} - human assembly build (grch37/grch38)
#'   \item \emph{id} - panel identifier (local)
#'   \item \emph{entrezgene} - Entrez gene identifier
#'   \item \emph{gene_biotype} - gene biotype
#'   \item \emph{genename} - Gene name
#'   \item \emph{ensembl_gene_id} - Ensembl gene identifier
#'   \item \emph{gepa_moi} - mechanism of inheritance
#'   (Genomics England PanelApp)
#'   \item \emph{gepa_penetrance} - penetrance (Genomics England PanelApp)
#'   \item \emph{gepa_confidence_level} - confidence level
#'   (Genomics England PanelApp)
#'   \item \emph{gepa_panel_name} - panel name (Genomics England PanelApp)
#'   \item \emph{gepa_panel_id} - panel identifier (Genomics England PanelApp)
#'   \item \emph{gepa_panel_version} - panel version (Genomics England PanelApp)
#'   \item \emph{gepa_phenotype} - associated cancer phenotypes
#'   (Genomics England PanelApp)
#'   \item \emph{gepa_panel_url} - panel URL (Genomics England PanelApp)
#' }
#'
#' @examples
#' \dontrun{
#' library(geneOncoX)
#' gene_panels <- get_panels(cache_dir = tempdir())
#' }
#'
#' @export
#'


get_panels <- function(cache_dir = NA, force_download = FALSE) {
  dat <- get_gox_data(
    cache_dir = cache_dir,
    force_download = force_download,
    db = "panels"
  )
  return(dat)
}
