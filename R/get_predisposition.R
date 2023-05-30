#' @title
#' Get human cancer predisposition genes
#'
#' @description
#' Downloads and retrieves a pre-processed dataset of human cancer
#' predisposition genes (CPGs). Genes deemed relevant for cancer predisposition
#' have been collected from multiple resources:
#'
#' * \href{https://panelapp.genomicsengland.co.uk/}{Genomics England PanelApp}
#' * \href{https://pubmed.ncbi.nlm.nih.gov/29625052/}{TCGA's PanCancer study}
#' * \href{https://cancer.sanger.ac.uk/census}{Cancer Gene Census}
#' * Other/curated
#'
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation
#'                resources used
#' * `records` - a list with records of cancer predisposing gene,
#'               one record per gene
#'
#' @details
#' \strong{NOTE}: The dataset also contains genes recommended for reporting of
#' incidental findings (ACMG_SF).
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
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
#' \strong{records} - A data frame with 594 rows and 8 columns:
#' \itemize{
#'   \item \emph{symbol} - official gene symbol
#'   \item \emph{entrezgene} - Entrez gene identifier
#'   \item \emph{cpg_moi} - mechanism of inheritance (AD, AR)
#'   \item \emph{cpg_syndrome_cui} - Concept unique identifiers
#'   (CUI, UMLS) - inherited cancer syndromes
#'   \item \emph{cpg_cancer_cui} - Concept unique identifiers (CUI, UMLS) -
#'   inherited cancer conditions
#'   \item \emph{cpg_source} - Sources supporting predisposition gene
#'   PANEL_APP, CGC, TCGA_PANCAN_2018, OTHER, ACMG_SF
#'   \item \emph{cpg_phenotypes} - associated cancer phenotypes
#'   \item \emph{cpg_mod} - mechanism of disease)
#' }
#'
#' @examples
#' \dontrun{
#' library(geneOncoX)
#' gene_predisp <- get_predisposition(cache_dir = tempdir())
#' }
#'
#' @export
#'

get_predisposition <- function(cache_dir = NA, force_download = FALSE) {
  dat <- get_gox_data(
    cache_dir = cache_dir,
    force_download = force_download,
    db = "predisposition"
  )
  return(dat)
}
