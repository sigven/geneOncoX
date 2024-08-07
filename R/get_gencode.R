#' @title
#' Get human GENCODE transcript dataset
#'
#' @description
#' Downloads and returns a pre-processed dataset of GENCODE transcripts.
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation
#'                resources used
#' * `records` - a list with two data frames of transcripts, one
#' for `grch37`, and one for `grch38` (see specific column format below)
#'
#' The GENCODE datasets provided were established by parsing the GENCODE
#' \href{https://www.gencodegenes.org/human/}{comprehensive
#' gene annotation GTF file} with the
#' \href{https://rnabioco.github.io/valr/}{valr} package, and
#' also providing cross-references to RefSeq (through the use of the
#' \href{https://bioconductor.org/packages/biomaRt/}{biomaRt}
#' package), and UniProt accessions/identifiers.
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should 
#' be overwritten (set to TRUE to re-download if file exists in 
#' cache)
#' @param ensembl_release version of Ensembl to use - this will
#' dictate the version of GENCODE used for grch38
#' @param chromosomes_only Logical indicating if transcripts returned
#' should belong to ordinary chromosomes only (scaffolds etc. 
#' should be ignored)
#'
#' @returns
#' \strong{metadata} - A data frame with 3 rows and 6 columns:
#' \itemize{
#'   \item \emph{source} - gene annotation source
#'   \item \emph{annotation_data} - type of annotations used
#'   \item \emph{url} - URL of annotation resource
#'   \item \emph{citation} - publication to cite for annotation source
#'   (citation; PMID)
#'   \item \emph{version} - version used
#'   \item \emph{abbreviation} - abbreviation (e.g. used in column
#'   names of records)
#'  }
#'
#' \strong{records} - A list with one data frame per assembly
#' \itemize{
#'   \item \emph{chrom} - chromosome
#'   \item \emph{start} - transcript start with 5kb padding (upstream)
#'   \item \emph{end} - transcript end with 5kb padding (downstream)
#'   \item \emph{transcript_start} - transcript start
#'   \item \emph{transcript_end} - transcript end
#'   \item \emph{strand} - strand
#'   \item \emph{ensembl_gene_id} - Ensembl gene identifier
#'   \item \emph{ensembl_transcript_id} - Ensembl transcript identifier
#'   \item \emph{ensembl_transcript_id_full} - Ensembl transcript
#'   identifier (with version)
#'   \item \emph{ensembl_protein_id} - Ensembl protein identifier
#'   \item \emph{symbol} - official gene symbol
#'   \item \emph{hgnc_id} - HGNC gene identifier
#'   \item \emph{entrezgene} - Entrez gene identifier
#'   \item \emph{name} - genename
#'   \item \emph{gene_biotype} -
#'   \href{https://www.gencodegenes.org/pages/biotypes.html}{gene biotype}
#'   (GENCODE)
#'   \item \emph{transcript_biotype} -
#'   \href{https://www.gencodegenes.org/pages/biotypes.html}{transcript biotype}
#'   (GENCODE))
#'   \item \emph{tag} -
#'   \href{https://www.gencodegenes.org/pages/biotypes.html}{tag} (GENCODE))
#'   \item \emph{refseq_protein_id} - RefSeq peptide identifier
#'   \item \emph{refseq_transcript_id} - RefSeq mRNA identifier
#'   \item \emph{mane_select} -
#'   \href{https://www.ncbi.nlm.nih.gov/refseq/MANE/}{MANE} Select
#'   \item \emph{mane_plus_clinical} -
#'   \href{https://www.ncbi.nlm.nih.gov/refseq/MANE/}{MANE} Plus Clinical
#'   \item \emph{principal_isoform_flag} - Principal
#'   \href{https://appris.bioinfo.cnio.es/#/help/scores}{isoform tag} (APPRIS)
#'   \item \emph{uniprot_acc} - UniProtKB protein accession
#'   \item \emph{uniprot_id} - UniProtKB protein identifier
#'   \item \emph{ensembl_version} - Ensembl version
#'   \item \emph{gencode_version} - GENCODE version
#'   \item \emph{uniprot_version} - UniProtKB version
#' }
#'
#' @source <https://www.gencodegenes.org/>
#' @source <https://uniprot.org>
#' @source <https://www.ncbi.nlm.nih.gov/gene>
#' @source <https://apprisws.bioinfo.cnio.es/landing_page/>
#'
#' @return list object with GENCODE transcripts (records) and associated
#'   metadata (metadata)
#'
#' @examples
#' \dontrun{
#' library(geneOncoX)
#' transcripts_gencode <- get_gencode(cache_dir = tempdir())
#' }
#'
#' @export
#'

get_gencode <- function(cache_dir = NA,
                        force_download = FALSE,
                        chromosomes_only = TRUE,
                        ensembl_release = 112) {
  
  if(ensembl_release > 112 | ensembl_release < 112){
    lgr::lgr$fatal(
      paste0("ERROR: Ensembl release must be equal to 112",
             " - exiting"))
      return(0)
  }
  
  gencode_release <- 46
  # if(ensembl_release == 112){
  #   gencode_release <- 46
  # }
  
  dat <- get_gox_data(
    cache_dir = cache_dir,
    force_download = force_download,
    db = paste0(
      "gencode_", gencode_release, 
      "_", ensembl_release)
  )
  
  if(chromosomes_only == TRUE){
    dat$records$grch37 <- dat$records$grch37 |> 
      dplyr::filter(
        stringr::str_detect(chrom, "^chr"))
    
    dat$records$grch38 <- dat$records$grch38 |> 
      dplyr::filter(
        stringr::str_detect(chrom, "^chr"))
  }
  return(dat)
}
