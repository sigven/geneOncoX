#' @title
#' Get basic cancer-relevant gene annotations
#'
#' @description
#' Downloads and returns a dataset that combines multiple human cancer gene
#' annotations,i.e. from IntOGen, CancerMine, Network of Cancer Genes,
#' Cancer Gene Census, NCBI, dbNSFP etc. The dataset comes as a `list` object,
#' with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation
#'                resources used
#' * `records` - a data frame with gene annotations (one record per gene)
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
#' \strong{records} - A data frame with 108,648 rows and 65 columns:
#' \itemize{
#'   \item \emph{entrezgene} - NCBI Entrez identifier
#'   \item \emph{symbol} - primary gene symbol
#'   \item \emph{gene_biotype} - type of gene (ncRNA, protein-coding, or pseudo)
#'   \item \emph{name} - gene name
#'   \item \emph{other_genename_designations} - other gene name designations
#'   \item \emph{hgnc_id} - HGNC gene identifier
#'   \item \emph{ncbi_function_summary} - gene function summary (NCBI Gene)
#'   \item \emph{cgc_hallmark} - Annotated with cancer gene hallmarks(s) by
#'   Cancer Gene Census (CGC)
#'   \item \emph{cgc_tier} - Cancer Gene Census tier (TIER1/TIER2)
#'   \item \emph{cgc_driver_tier1} - Logical indicating if gene is part of
#'   Cancer Gene Census - TIER1
#'   \item \emph{cgc_driver_tier2} - Logical indicating if gene is part of
#'   Cancer Gene Census - TIER2
#'   \item \emph{cgc_tsg} - tumor suppressor gene (Cancer Gene Census)
#'   \item \emph{cgc_oncogene} - proto-oncogene (Cancer Gene Census)
#'   \item \emph{cgc_somatic} - logical indicating whether the cancer relevance
#'   of this gene relates to the soma (Cancer Gene Census)
#'   \item \emph{cgc_phenotype_somatic} - cancer phenotypes relevant for
#'   somatic mutations of this gene (Cancer Gene Census)
#'   \item \emph{cgc_germline} - logical indicating whether the cancer
#'   relevance of this gene relates to the germline (Cancer Gene Census)
#'   \item \emph{cgc_phenotype_germline} - cancer phenotypes relevant for
#'   germline mutations of this gene (Cancer Gene Census)
#'   \item \emph{ncg_driver} - canonical cancer driver gene according to
#'   Network of Cancer Genes (NCG)
#'   \item \emph{ncg_tsg} - tumor suppressor gene (NCG)
#'   \item \emph{ncg_oncogene} - proto-oncogene (NCG)
#'   \item \emph{ncg_phenotype} - cancer phenotypes relevant for this gene (NCG)
#'   \item \emph{ncg_pmid} - supporting literature identifiers (PMIDs, NCG)
#'   \item \emph{intogen_role} - cancer driver role (IntOGen)
#'   \item \emph{intogen_phenotype} - cancer phenotypes relevant for
#'   this gene (IntOGen)
#'   \item \emph{intogen_driver} - logical indicatin if gene is predicted
#'   as cancer driver by IntOGen
#'   \item \emph{bailey2018_fp_driver} - logical indicating whether this
#'   gene is likely a false positive driver gene (Bailey et al., Cell, 2018)
#'   \item \emph{woods_dnarepair_class} - class of DNA repair
#'   (DNA repair database, Woods et al.)
#'   \item \emph{woods_dnarepair_activity} - type of DNA repair activity
#'   involved (DNA repair database, Woods et al.)
#'   \item \emph{illumina_tso500} - gene is part of Illumina's TSO500 panel
#'   (SNV_INDEL, CNA_GAIN, CNA_LOSS, RNA_FUSION)
#'   \item \emph{foundation_one_f1cdx} - gene is part of Foundation
#'   One's F1CDx panel (SNV_INDEL, CNA, FUSION, PROMOTER)
#'   \item \emph{cpic_pgx} - gene related to pharmacogenomics (CPIC, antineoplastic drugs)
#'   \item \emph{sanchezvega2018_signaling_pathway} - curated signalling
#'   pathways (Sanchez-Vega et al., Cell, 2018)
#'   \item \emph{cancermine_pmid_driver} - PMIDs that support
#'   (from text mining) a role for this gene as a driver (CancerMine)
#'   \item \emph{cancermine_pmid_oncogene} - PMIDs that support (from
#'   text mining) a role for this gene as a proto-oncogene (CancerMine)
#'   \item \emph{cancermine_pmid_tsg} - PMIDs that support (from text mining)
#'   a role for this gene as a tumor suppressor gene (CancerMine)
#'   \item \emph{cancermine_doid_driver} - cancer phenotypes relevant for
#'   the given role (Disease Ontology identifiers, CancerMine)
#'   \item \emph{cancermine_doid_oncogene} - cancer phenotypes relevant for
#'   the given role (Disease Ontology identifiers, CancerMine)
#'   \item \emph{cancermine_doid_tsg} - cancer phenotypes relevant for the
#'   given role (Disease Ontology identifiers, CancerMine)
#'   \item \emph{cancermine_n_cit_driver} - number of citations (PMIDs)
#'   that support a role for this as a driver (CancerMine)
#'   \item \emph{cancermine_n_cit_oncogene} - number of citations (PMIDs)
#'   that support a role for this as a proto-oncogene (CancerMine)
#'   \item \emph{cancermine_n_cit_tsg} - number of citations (PMIDs)
#'   that support a role for this as a tumor suppressor (CancerMine)
#'   \item \emph{cancermine_cit_tsg} - citations for cancer driver
#'   support (all with prob > 0.8, CancerMine)
#'   \item \emph{cancermine_cit_oncogene} - citations for proto-oncogene
#'   support (all with prob > 0.8, max 50, CancerMine)
#'   \item \emph{cancermine_cit_driver} - citations for tumor suppressor
#'   gene support (all with prob > 0.8, CancerMine)
#'   \item \emph{cancermine_cit_links_driver} - citation links for cancer
#'   driver support (50 most recent, CancerMine)
#'   \item \emph{cancermine_cit_links_oncogene} - citation links for
#'   proto-oncogene support (50 most recent, CancerMine)
#'   \item \emph{cancermine_cit_links_tsg} - citation links for tumor
#'   suppressor gene support (50 most recent, CancerMine)
#'   \item \emph{mim_id} - MIM gene id (from HGNC)
#'   \item \emph{mim_phenotype_id} - MIM (ids) of the phenotype the gene
#'   caused or associated with (dbNSFP, from Uniprot)
#'   \item \emph{prob_haploinsuffiency} - Estimated probability of
#'   haploinsufficiency of the gene
#'   (from doi:10.1371/journal.pgen.1001154) (dbNSFP)
#'   \item \emph{gene_indispensability_score} - A probability prediction of
#'   the gene being essential. From doi:10.1371/journal.pcbi.1002886 (dbNSFP)
#'   \item \emph{gene_indispensability_pred} - Essential ("E") or
#'   loss-of-function tolerant ("N") based on
#'   gene_indispensability_score (dbNSFP)
#'   \item \emph{essential_gene_crispr} - Essential ("E") or Non-essential
#'   phenotype-changing ("N") based on
#'   large scale CRISPR experiments. from doi: 10.1126/science.aac7041 (dbNSFP)
#'   \item \emph{essential_gene_crispr2} - Essential ("E"),
#'   context-Specific essential ("S"), or Non-essential phenotype-changing ("N")
#'   based on large scale CRISPR experiments.
#'   from http://dx.doi.org/10.1016/j.cell.2015.11.015 (dbNSFP)
#'   \item \emph{prob_gnomad_lof_intolerant} - the probability of being
#'   loss-of-function intolerant (intolerant of both heterozygous and
#'   homozygous lof variants) based on gnomAD 2.1 data (dbNSFP)
#'   \item \emph{prob_gnomad_lof_intolerant_hom} - the probability of being
#'   intolerant of homozygous, but not heterozygous lof variants based on
#'   gnomAD 2.1 data (dbNSFP)
#'   \item \emph{prob_gnomad_lof_tolerant_null} - the probability of being
#'   tolerant of both heterozygous and homozygous lof variants based on gnomAD
#'   2.1 data (dbNSFP)
#'   \item \emph{prob_exac_lof_intolerant} - the probability of being
#'   loss-of-function intolerant (intolerant of both heterozygous and
#'   homozygous lof variants) based on ExAC r0.3 data (dbNSFP)
#'   \item \emph{prob_exac_lof_intolerant_hom} - the probability of being
#'   intolerant of homozygous, but not heterozygous lof variants based on
#'   ExAC r0.3 data (dbNSFP)
#'   \item \emph{prob_exac_lof_tolerant_null} - the probability of being
#'   tolerant of both heterozygous and homozygous lof variants based on
#'   ExAC r0.3 data (dbNSFP)
#'   \item \emph{prob_exac_nontcga_lof_intolerant} - the probability of
#'   being loss-of-function intolerant (intolerant of both heterozygous and
#'   homozygous lof variants) based on ExAC r0.3 nonTCGA subset (dbNSFP)
#'   \item \emph{prob_exac_nontcga_lof_intolerant_hom} - the probability of
#'   being intolerant of homozygous, but not heterozygous lof variants based
#'   on ExAC 0.3 nonTCGA subset (dbNSFP)
#'   \item \emph{prob_exac_nontcga_lof_tolerant_null} - the probability of
#'   being tolerant of both heterozygous and homozygous lof variants based on
#'   ExAC r0.3 nonTCGA subset (dbNSFP)
#'   \item \emph{dbnsfp_function_description} - gene function description
#'   (dbNSFP/UniProtKB)
#' }
#'
#' @examples
#' \dontrun{
#' library(geneOncoX)
#' gene_basic <- get_basic(cache_dir = tempdir())
#' }
#'
#' @export
#'

get_basic <- function(cache_dir = NA, force_download = FALSE) {
  dat <- get_gox_data(
    cache_dir = cache_dir,
    force_download = force_download,
    db = "basic"
  )
  return(dat)
}
