#' @title
#' Get human cancer predisposition genes
#'
#' @description
#' Downloads and retrieves a pre-processed dataset of human cancer predisposition
#' genes. Genes deemed relevant for cancer predisposition
#' have been collected from multiple resources:
#'
#' * \href{https://panelapp.genomicsengland.co.uk/}{Genomics England PanelApp}
#' * \href{https://pubmed.ncbi.nlm.nih.gov/29625052/}{TCGA's PanCancer study}
#' * \href{https://cancer.sanger.ac.uk/census}{Cancer Gene Census}
#' * \href{https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html}{Human DNA repair genes}
#' * Other/curated
#'
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation resources used
#' * `records` - a list with records of cancer predisposing gene, one record per gene
#'
#' @details
#' \bold{NOTE}: The dataset also contains genes recommended for reporting of incidental
#' findings (ACMG_SF).
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#'
#' @returns
#' \bold{metadata} - A data frame with 1 row and 6 columns:
#' \itemize{
#'   \item \emph{source} - gene annotation source
#'   \item \emph{annotation_data} - type of annotations used
#'   \item \emph{url} - URL of annotation resource
#'   \item \emph{citation} - publication to cite for annotation source (citation; PMID)
#'   \item \emph{version} - version used
#'   \item \emph{abbreviation} - abbreviation used in column names of records
#'  }
#'
#' \bold{records} - A data frame with 591 rows and 8 columns:
#' \itemize{
#'   \item \emph{symbol} - official gene symbol
#'   \item \emph{entrezgene} - Entrez gene identifier
#'   \item \emph{moi} - mechanism of inheritance (AD, AR)
#'   \item \emph{predisp_syndrome_cui} - Concept unique identifiers (CUI, UMLS) -
#'   inherited cancer syndromes
#'   \item \emph{predisp_cancer_cui} - Concept unique identifiers (CUI, UMLS) -
#'   inherited cancer conditions
#'   \item \emph{predisp_source} - Sources for cancer-predisposition gene
#'   PANEL_APP, CGC, TCGA_PANCAN_2018, OTHER, ACMG_SF
#'   \item \emph{phenotypes} - associated cancer phenotypes
#'   \item \emph{mechanism_of_disease} - mechanism of disease)
#' }
#'
#' @export
#'

get_predisposition <- function(cache_dir = NA, force_download = F){

  dat <- get_gox_data(cache_dir = cache_dir,
                      force_download = force_download,
                      db = "predisposition")
  return(dat)
}
